#include "RH.h"



/**
 * @brief Algorithm RH (Random-Half)
 * @param pset 		The point set
 * @param u 		The expected point
 */
void RH(point_set *pset, point_t* e, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
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
    preTime = timeCost(t1);

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
                    if(Qcount == 50)
                        std::cout << "";
                    double dist1 = p1->distance(e);
                    double dist2 = p2->distance(e);
                    if(dimu > 1)
                    {
                        hyperplane *h;
                        if (dist1 > dist2)
                        {
                            h = new hyperplane(p1, p2, R->expDim);
                        }
                        else
                            h = new hyperplane(p2, p1, R->expDim);
                        R->hyperplanes.push_back(h);
                        R->set_ext_pts();
                    }
                    else
                    {
                        int d = p1->d;
                        double bd = p1->bound(p2, R->expDim);
                        if(dist1 > dist2)
                        {
                            if (p1->attr[d - 1] < p2->attr[d - 1])
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
                            if (p1->attr[d - 1] < p2->attr[d - 1])
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
                    printMiddleResult(out_cp, t1, preTime, Qcount, 100, 100, mem_baseline, type);
                }

                point_t* Rpoint = R->findNearest(pset);
                if(Rpoint != NULL)
                {
                    Rpoint->printResult(out_cp, "RH", Qcount, t1, preTime, mem_baseline, M, type);
                    out_cp.close();
                    return;
                }
            }
        }
    }
}