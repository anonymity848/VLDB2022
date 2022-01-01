#include "ActiveRanking.h"

/**
 * @brief Ask user questions and give a ranking
 * @param original_set 		The original dataset
 * @param u 				The linear function
 * @param k 				The threshold top-k
 */
int ActiveRanking(point_set* pset, point_t* e, long mem_baseline)
{
    //std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
    hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    int Qcount = 0, M = pset->points.size(), dimu = pset->points[0]->d_unorder;
    pset->random(0.5);
    std::vector<point_t*> current_use;
    current_use.push_back(pset->points[0]); //store all the points in order
    point_t *p1 = NULL, *p2 = NULL;
    preTime = timeCost(t1);

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
                    double dist1 = p1->distance(e);
                    double dist2 = p2->distance(e);
                    if(dimu > 1)
                    {
                        hyperplane *h;
                        if (dist1 > dist2)
                        {
                            h = new hyperplane(p1, p2, R->expDim);
                            place = j + 1;
                        }
                        else
                            h = new hyperplane(p2, p1, R->expDim);
                        R->hyperplanes.push_back(h);
                        R->set_ext_pts();
                    }
                    else
                    {
                        double bd = p1->bound(p2, R->expDim->attr[0]);
                        if(dist1 > dist2)
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
                    //printMiddleResult(out_cp, t1, preTime, Qcount, 100, mem_baseline);
                }
                else if (relation == -1)
                {
                    place = j + 1;
                }
            }
            current_use.insert(current_use.begin() + place, pset->points[i]);
        }
    }
    current_use[0]->printResult("ActiveRanking", Qcount, t1, preTime, mem_baseline);
    //out_cp.close();
}
