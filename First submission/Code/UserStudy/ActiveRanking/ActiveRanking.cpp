#include "ActiveRanking.h"

/**
 * @brief Algorithm Active-Ranking
 * @param pset      The dataset
 * @param realSet   The real dataset
 * @param e         The expected point (for test)
 * @param TID       The recommended point ID
 * @param fp        The file for recording results
 * @return
 */
int ActiveRanking(point_set* pset, point_set* realSet, point_t* e, int &TID, std::ofstream &fp)
{
    hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    int Qcount = 0, M = pset->points.size(), dimu = pset->points[0]->d_unorder;
    pset->random(0.5);
    std::vector<point_t*> current_use;
    current_use.push_back(pset->points[0]); //store all the points in order
    point_t *p1 = NULL, *p2 = NULL;

    //interaction
    for (int i = 1; i < M; i++) //compare: p_set contains all the points
    {
        bool same_exist = false;
        for (int j = 0; j < current_use.size(); j++)
        {
            if (pset->points[i]->is_same(current_use[j]))
            {
                same_exist = true;
                break;
            }
        }
        if (!same_exist)
        {
            int num_point = current_use.size(), place = 0; //the place of the point inserted into the current_use
            p1 = pset->points[i];
            //find the question asked user
            for (int j = 0; j < num_point; j++)
            {
                p2 = current_use[j];
                hyperplane *h = new hyperplane(p1, p2, R->expDim);
                int relation = R->check_relation(h);
                if (relation == 0)
                {
                    Qcount++;
                    if(Qcount >= 80)
                    {
                        realSet->points[current_use[0]->id]->final_result("ActiveRanking", Qcount, fp);
                        TID = current_use[0]->id;
                        return 1;
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
                            place = j + 1;
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
                            place = j + 1;
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
                }
                else if (relation == -1)
                {
                    place = j + 1;
                }
            }
            current_use.insert(current_use.begin() + place, pset->points[i]);
        }
    }
    realSet->points[current_use[0]->id]->final_result("ActiveRanking", Qcount, fp);
    TID = current_use[0]->id;
}

