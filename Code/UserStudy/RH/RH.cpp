#include "RH.h"



/**
 * @brief Algorithm RH
 * @param pset      The dataset
 * @param realSet   The real dataset
 * @param e         The expected point (for test)
 * @param TID       The id of the recommended point
 * @param fp        The file for recording results
 */
void RH(point_set *pset, point_set *realSet, point_t* e, int &TID, std::ofstream &fp)
{
    int dim = pset->points[0]->d, dimu = pset->points[0]->d_unorder, Qcount = 0;
    double M = pset->points.size();
    pset->random(0.5);
    hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    std::vector<point_t*> top_current;
    std::vector<point_t*> currentPoint, currentUse;
    point_t* point_result = NULL;
    currentUse.push_back(pset->points[0]);
    point_t *p1 = NULL, *p2 = NULL;

    //interaction
    for(int i = 1; i < M; ++i)
    {
        bool same_exist=false;
        for(int j = 0; j < currentUse.size(); ++j)
        {
            if(pset->points[i]->is_same(currentUse[j]))
            {
                same_exist=true;
                break;
            }
        }
        if(!same_exist)
        {
            for(int j = 0; j < currentUse.size(); j++)
                currentPoint.push_back(currentUse[j]);
            currentUse.push_back(pset->points[i]);
            p1 = pset->points[i];

            while(currentPoint.size() > 0)
            {
                int numPoint = currentPoint.size(), scan_index = 0;//the index where the points we have scanned
                bool need_ask = false; double distance = INF;
                int p_index = -1;// the index where the point we select to ask question
                //find the question asked user
                for(int j = 0; j < numPoint; j++)
                {
                    hyperplane *h = new hyperplane(p1, currentPoint[scan_index], R->expDim);
                    if(R->check_relation(h) == 0)
                    {
                        need_ask = true;
                        double d_h = h->distance(R->average_point());
                        if(d_h < distance)
                        {
                            distance = d_h;
                            p2 = currentPoint[scan_index];
                            p_index = scan_index;
                        }
                        scan_index++;
                    }
                    else
                    {
                        currentPoint.erase(currentPoint.begin() + scan_index);
                    }
                }
                if(need_ask)
                {
                    Qcount++;
                    if(Qcount >= 80)
                    {
                        point_t* Rpoint = R->approxNearest(pset);
                        realSet->points[Rpoint->id]->final_result("RH", Qcount, fp);
                        TID = Rpoint->id;
                        return;
                    }
                    double dist1 = p1->distance(e);
                    double dist2 = p2->distance(e);
                    cout << "\nThe " << Qcount <<"th question\n";
                    int option = realSet->show_to_user(p1->id, p2->id);
                    if(dimu > 1)
                    {
                        hyperplane *h;
                        if (option == 2)
                        {
                            h = new hyperplane(p1, p2, R->expDim);
                        }
                        else
                            h = new hyperplane(p2, p1, R->expDim);
                        R->hyperplanes.push_back(h);
                        bool update = R->set_ext_pts();
                        if(!update)
                            R->hyperplanes.pop_back();
                    }
                    else
                    {
                        double bd = p1->bound(p2, R->expDim->attr[0]);
                        if (option == 2)
                        {
                            if (p1->attr[1] < p2->attr[1])
                            {
                                R->ext_pts[1]->attr[0] = bd;
                            }
                            else
                            {
                                R->ext_pts[0]->attr[0] = bd;
                            }
                        }
                        else
                        {
                            if (p1->attr[1] < p2->attr[1])
                            {
                                R->ext_pts[0]->attr[0] = bd;
                            }
                            else
                            {
                                R->ext_pts[1]->attr[0] = bd;
                            }
                        }
                    }
                    currentPoint.erase(currentPoint.begin() + p_index);
                }

                point_t* Rpoint = R->findNearest(pset);
                if(Rpoint != NULL)
                {
                    realSet->points[Rpoint->id]->final_result("RH", Qcount, fp);
                    TID = Rpoint->id;
                    return;
                }
            }
        }
    }


    point_t* Rpoint = R->approxNearest(pset);
    realSet->points[Rpoint->id]->final_result("RH", Qcount, fp);
    TID = Rpoint->id;
    return;



}