#include "Baseline.h"


void baseline(point_set *orgSet, point_t *e, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;

    point_set *pset = new point_set(orgSet);
    int Qcount = 0, dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder, batch = 0, first = 0;
    double volume; point_t *p1 = NULL, *p2 =NULL;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size(), LastCsize = pointSize;
    point_set *cPoint = new point_set();
    std::vector<hyperplane_set*> candPartition;
    R->findRegion(orgSet, pset, cPoint, candPartition);
    std::vector<hyperplane*> candHyper;
    for(int i = 0; i < candPartition.size(); ++i)
    {
        for(int j = 0; j < candPartition[i]->hyperplanes.size(); ++j)
        {
            if(candPartition[i]->hyperplanes[j]->p_1 != NULL && candPartition[i]->hyperplanes[j]->p_2 != NULL)
            {
                bool exist = false;
                for (int k = 0; k < candHyper.size(); ++k)
                {
                    if (candHyper[k]->is_same(candPartition[i]->hyperplanes[j]))
                    {
                        exist = true;
                        break;
                    }
                }
                if (!exist)
                    candHyper.push_back(candPartition[i]->hyperplanes[j]);
            }
        }
    }


    double dis = (pow(volume/cPoint->points.size(), 1.0/dimu))/2;
    double Beta = 0.1;
    //interaction
    while (cPoint->points.size() > 1)
    {
        int hIndex, heven = -1;
        if(dimu == 1)
        {
            for (int i = 0; i < candHyper.size(); ++i)
            {
                int pri = candHyper[i]->priority(cPoint, dis, Beta);
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
        }
        else
        {
            hIndex = rand() % candHyper.size();
            while(candHyper[hIndex]->p_1->value == 2 || candHyper[hIndex]->p_2->value == 2)
                hIndex = rand() % candHyper.size();
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

            if(dimu > 1 && candHyper.size() == 0)
                std::cout << "wrong\n";

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
                //R->print();
            }
            else
            {
                double bd = p1->bound(p2, R->expDim->attr[0]);
                if(dist1 > dist2)
                {
                    if (p1->attr[1] < p2->attr[1])
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
                    if (p1->attr[1] < p2->attr[1])
                    {
                        R->ext_pts[0]->attr[0] = bd;
                    } else
                    {
                        R->ext_pts[1]->attr[0] = bd;
                    }
                    p2->value = 2;
                    cPoint->prunePt(p2);
                }
                //R->print();
            }

            //R->nearestSkyline(cPoint);
            if(dimu > 1)
            {
                for (int a = 0; a < cPoint->points.size(); ++a)
                {
                    if (h->necessaryPrune(cPoint->points[a]))
                    {
                        if (R->is_prune(cPoint, a, candPartition))
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
            cPoint->points[0]->printResult(out_cp, "Baseline", Qcount, t1, preTime, mem_baseline, pointSize, type);
            out_cp.close();
            return;
        }
    }

    /*
    std::vector<hyperplane_set*> sortPartition;

    for(int i = 0; i < candPartition.size(); ++i)
    {
        int postition = 0;
        for(int j = 0; j < sortPartition.size(); ++j)
        {
            if(sortPartition[j]->ext_pts[0]->attr[0] < candPartition[i]->ext_pts[0]->attr[0])
                postition++;
            else
                break;
        }
        sortPartition.insert(sortPartition.begin() + postition, candPartition[i]);
    }

    for(int i = 0; i < sortPartition.size(); ++i)
    {
        std::cout << "Partition: " << i << "\n";
        sortPartition[i]->print();
    }

    for(int i = 0; i < sortPartition.size() - 1; ++i)
    {
        if(sortPartition[i]->ext_pts[1]->attr[0] != sortPartition[i+1]->ext_pts[0]->attr[0])
            std::cout<<"wrong\n";

    }
    */
}
