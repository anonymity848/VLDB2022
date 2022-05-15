#include "BS.h"

/**
 * @brief Algorithm RS
 * @param pset  The point set
 * @param e     The user's expected point
 */
void BS(point_set *pset, point_t* e, double Beta, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
    int Qcount = 0, dimu = pset->points[0]->d_unorder, first = 0;
    double volume;
    point_t *p1 = NULL, *p2 =NULL;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size(), LastCsize = pointSize;
    std::vector<hyperplane*> candHyper;
    point_set *cPoint = new point_set();
    std::vector<hyperplane_set*> partitionSet;
    std::queue<int> candPoint;
    int tt = pset->findClosest(R->expDim);
    candPoint.push(tt);
    pset->points[tt]->value = 1;
    while (candPoint.size() > 0)
    {
        int index= candPoint.front(); candPoint.pop();
        if(R->find_boundary(pset, index, candPoint, candHyper, partitionSet))
        {
            //pset->points[index]->print();
            cPoint->points.push_back(pset->points[index]);
        }
    }

    double dis = (pow(volume/cPoint->points.size(), 1.0/dimu))/2;

    //cPoint->points[0]->printResult(out_cp, "BS", Qcount, t1, preTime, mem_baseline);

    while (cPoint->points.size() > 1)
    {
        int hIndex, heven = -1;
        for(int i = 0; i < candHyper.size(); ++i)
        {
            int pri = candHyper[i]->priority(cPoint, dis, Beta);
            /*
            std::cout <<"\n"<< i <<" "<<pri <<"  " << candHyper[i]->p_1->bound(candHyper[i]->p_2, R->expDim->attr[0]) <<"\n";
            candHyper[i]->p_1->internal->print();
            candHyper[i]->p_2->internal->print();
            */
            if(pri == -1)
            {
                candHyper.erase(candHyper.begin() + i);
                --i;
            }
            else if(pri >= heven)
            {
                heven = pri;
                hIndex = i;
            }
        }

        hyperplane *h;
        if(candHyper.size() > 0)
            h = candHyper[hIndex];
        else
            h = new hyperplane(cPoint->points[0], cPoint->points[1]);
        p1 = h->p_1, p2 = h->p_2;
        int relation = R->check_relation(h);// Check whether it intersects with R
        if (relation == 0) //if intersect, calculate the distance
        {
            if(first == 0)
            {
                preTime = timeCost(t1);
                first = 1;
            }

            Qcount++;
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);
            if(dimu > 1)
            {
                if (dist1 > dist2)
                {
                    h = new hyperplane(p1, p2, R->expDim);
                    p1->value = 2;
                    cPoint->prunePt(p1);
                }
                else
                {
                    h = new hyperplane(p2, p1, R->expDim);
                    p2->value = 2;
                    cPoint->prunePt(p2);
                }
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
                    } else
                    {
                        R->ext_pts[0]->attr[0] = bd;
                    }
                    p1->value = 2;
                    cPoint->prunePt(p1);
                }
                else
                {
                    if (p1->attr[d - 1] < p2->attr[d - 1])
                    {
                        R->ext_pts[0]->attr[0] = bd;
                    } else
                    {
                        R->ext_pts[1]->attr[0] = bd;
                    }
                    p2->value = 2;
                    cPoint->prunePt(p2);
                }
            }
            //R->nearestSkyline(cPoint);
            if(dimu > 1)
            {
                for (int a = 0; a < cPoint->points.size(); ++a)
                {
                    if (h->necessaryPrune(cPoint->points[a]))
                    {
                        if (R->is_prune(cPoint, a, partitionSet))
                        {
                            cPoint->points[a]->value = 2;
                            cPoint->points.erase(cPoint->points.begin() + a);
                            --a;
                        }
                    }
                }
            }
            else
            {
                for (int a = 0; a < cPoint->points.size(); ++a)
                {
                    if(cPoint->points[a]->ext[0]->attr[0] > R->ext_pts[0]->attr[0])
                        cPoint->points[a]->ext[0]->attr[0] = R->ext_pts[0]->attr[0];
                    if(cPoint->points[a]->ext[1]->attr[0] < R->ext_pts[1]->attr[0])
                        cPoint->points[a]->ext[1]->attr[0] = R->ext_pts[1]->attr[0];
                    cPoint->points[a]->internal->attr[0] = (cPoint->points[a]->ext[0]->attr[0] + cPoint->points[a]->ext[1]->attr[0])/2;
                    if(cPoint->points[a]->ext[0]->attr[0] <= cPoint->points[a]->ext[1]->attr[0])
                    {
                        cPoint->points[a]->value = 2;
                        cPoint->points.erase(cPoint->points.begin() + a);
                        --a;
                    }
                }
            }
            printMiddleResult(out_cp, t1, preTime, Qcount, (cPoint->points.size()/pointSize) * 100,
                              (1 - cPoint->points.size()/LastCsize) * 100, mem_baseline, type);
            LastCsize = cPoint->points.size();
        }
        else
        {
            if (relation == -1)
            {
                p1->value = 2;
                cPoint->prunePt(p1);
            }
            else if (relation == 1)
            {
                p2->value = 2;
                cPoint->prunePt(p2);
            }
        }
        if(candHyper.size() > 0)
            candHyper.erase(candHyper.begin() + hIndex);

        if (cPoint->points.size() == 1)
        {
            cPoint->points[0]->printResult(out_cp, "BS", Qcount, t1, preTime, mem_baseline, pointSize, type);
            out_cp.close();
            return;
        }
    }

}


/**
 * @brief Algorithm RS
 * @param pset  The point set
 * @param e     The user's expected point
 */
void BSTopk(point_set *origset, point_t* e, double Beta, int k, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1;
    gettimeofday(&t1, 0);
    double preTime;
    int Qcount = 0, dimu = origset->points[0]->d_unorder, first = 0;
    double volume;
    point_t *p1 = NULL, *p2 = NULL;
    hyperplane_set *R = new hyperplane_set(origset, volume);
    R->nearestSkyband(origset, k);
    double pointSize = origset->points.size(), LastCsize = pointSize;
    point_set *maintainSet = new point_set(origset);
    point_set *resultSet = new point_set();

    while (resultSet->points.size() < k)
    {

        point_set *pset = new point_set(maintainSet);
        R->nearestSkyline(pset);
        std::vector<hyperplane *> candHyper;
        point_set *cPoint = new point_set();
        std::vector<hyperplane_set *> partitionSet;
        std::queue<int> candPoint;
        R->update_expDim();
        //R->expDim->print();
        int tt = pset->findClosest(R->expDim);
        candPoint.push(tt);
        pset->points[tt]->value = 1;
        while (candPoint.size() > 0)
        {
            int index = candPoint.front();
            candPoint.pop();
            if (R->find_boundary(pset, index, candPoint, candHyper, partitionSet))
            {
                //pset->points[index]->print();
                cPoint->points.push_back(pset->points[index]);
            }
        }
        double dis = (pow(volume / cPoint->points.size(), 1.0 / dimu)) / 2;


        //interaction
        while (cPoint->points.size() > 1)
        {
            int hIndex, heven = -1;
            for (int i = 0; i < candHyper.size(); ++i)
            {
                int pri = candHyper[i]->priority(cPoint, dis, Beta);
                /*
                std::cout <<"\n"<< i <<" "<<pri <<"  " << candHyper[i]->p_1->bound(candHyper[i]->p_2, R->expDim->attr[0]) <<"\n";
                candHyper[i]->p_1->internal->print();
                candHyper[i]->p_2->internal->print();
                */
                if (pri == -1)
                {
                    candHyper.erase(candHyper.begin() + i);
                    --i;
                } else if (pri >= heven)
                {
                    heven = pri;
                    hIndex = i;
                }
            }

            hyperplane *h;
            if (candHyper.size() > 0)
                h = candHyper[hIndex];
            else
                h = new hyperplane(cPoint->points[0], cPoint->points[1]);
            p1 = h->p_1, p2 = h->p_2;
            int relation = R->check_relation(h);// Check whether it intersects with R
            if (relation == 0) //if intersect, calculate the distance
            {
                if (first == 0)
                {
                    preTime = timeCost(t1);
                    first = 1;
                }

                Qcount++;
                double dist1 = p1->distance(e);
                double dist2 = p2->distance(e);
                if (dimu > 1)
                {
                    if (dist1 > dist2)
                    {
                        h = new hyperplane(p1, p2, R->expDim);
                        p1->value = 2;
                        cPoint->prunePt(p1);
                    } else
                    {
                        h = new hyperplane(p2, p1, R->expDim);
                        p2->value = 2;
                        cPoint->prunePt(p2);
                    }
                    R->hyperplanes.push_back(h);
                    R->set_ext_pts();
                } else
                {
                    int d = p1->d;
                    double bd = p1->bound(p2, R->expDim);
                    if (dist1 > dist2)
                    {

                        if (p1->attr[d - 1] < p2->attr[d - 1])
                        {
                            R->ext_pts[1]->attr[0] = bd;
                        } else
                        {
                            R->ext_pts[0]->attr[0] = bd;
                        }
                        p1->value = 2;
                        cPoint->prunePt(p1);
                    } else
                    {
                        if (p1->attr[d - 1] < p2->attr[d - 1])
                        {
                            R->ext_pts[0]->attr[0] = bd;
                        } else
                        {
                            R->ext_pts[1]->attr[0] = bd;
                        }
                        p2->value = 2;
                        cPoint->prunePt(p2);
                    }
                }
                //R->nearestSkyline(cPoint);
                if (dimu > 1)
                {
                    for (int a = 0; a < cPoint->points.size(); ++a)
                    {
                        if (h->necessaryPrune(cPoint->points[a]))
                        {
                            if (R->is_prune(cPoint, a, partitionSet))
                            {
                                cPoint->points[a]->value = 2;
                                cPoint->points.erase(cPoint->points.begin() + a);
                                --a;
                            }
                        }
                    }
                } else
                {
                    for (int a = 0; a < cPoint->points.size(); ++a)
                    {
                        if (cPoint->points[a]->ext[0]->attr[0] > R->ext_pts[0]->attr[0])
                            cPoint->points[a]->ext[0]->attr[0] = R->ext_pts[0]->attr[0];
                        if (cPoint->points[a]->ext[1]->attr[0] < R->ext_pts[1]->attr[0])
                            cPoint->points[a]->ext[1]->attr[0] = R->ext_pts[1]->attr[0];
                        cPoint->points[a]->internal->attr[0] =
                                (cPoint->points[a]->ext[0]->attr[0] + cPoint->points[a]->ext[1]->attr[0]) / 2;
                        if (cPoint->points[a]->ext[0]->attr[0] <= cPoint->points[a]->ext[1]->attr[0])
                        {
                            cPoint->points[a]->value = 2;
                            cPoint->points.erase(cPoint->points.begin() + a);
                            --a;
                        }
                    }
                }
                printMiddleResult(out_cp, t1, preTime, Qcount, (cPoint->points.size() / pointSize) * 100,
                                  (1 - cPoint->points.size() / LastCsize) * 100, mem_baseline, type);
                LastCsize = cPoint->points.size();
            } else
            {
                if (relation == -1)
                {
                    p1->value = 2;
                    cPoint->prunePt(p1);
                } else if (relation == 1)
                {
                    p2->value = 2;
                    cPoint->prunePt(p2);
                }
            }
            if (candHyper.size() > 0)
                candHyper.erase(candHyper.begin() + hIndex);

        }

        resultSet = R->findPossibleNearestK(origset, k);
        int PID = cPoint->points[0]->id;
        for(int i = 0; i < maintainSet->points.size(); ++i)
        {
            if(maintainSet->points[i]->id == PID)
            {
                maintainSet->points.erase(maintainSet->points.begin() + i);
                break;
            }
        }
    }
    resultSet->printResult(out_cp, "BSTopk", Qcount, t1, preTime, mem_baseline, pointSize, type);
    out_cp.close();
    return;
}