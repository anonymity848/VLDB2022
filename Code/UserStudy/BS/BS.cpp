#include "BS.h"

/**
 * @brief Algorithm RS
 * @param pset  The point set
 * @param e     The user's expected point
 */
void BS(point_set* pset, point_set* realSet, point_t* e, double Beta, int &TID, std::ofstream &fp)
{
    int Qcount = 0, dimu = pset->points[0]->d_unorder;
    double volume;
    point_t *p1 = NULL, *p2 =NULL;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size();
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

    //interaction
    while (cPoint->points.size() > 1)
    {
        int hIndex, heven = -1;
        for(int i = 0; i < candHyper.size(); ++i)
        {
            int pri = candHyper[i]->priority(cPoint, dis, Beta);
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
            Qcount++;
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);
            cout << "\nThe " << Qcount <<"th question\n";
            int option = realSet->show_to_user(p1->id, p2->id);
            if(dimu > 1)
            {
                if (option == 2)
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
                double bd = p1->bound(p2, R->expDim->attr[0]);
                if(option == 2)
                {
                    if (p1->attr[1] < p2->attr[1])
                    {
                        R->ext_pts[1]->attr[0] = bd;
                    } else
                    {
                        R->ext_pts[0]->attr[0] = bd;
                    }
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
            TID = cPoint->points[0]->id;
            realSet->points[TID]->final_result("BS", Qcount, fp);
            return;
        }
    }

}




/**
 * @brief Algorithm RS
 * @param pset  The point set
 * @param e     The user's expected point
 */
void BS_unanswer(point_set* pset, point_set* realSet, point_t* e, double Beta, std::vector<int> &TID, std::ofstream &fp)
{
    int Qcount = 0, dimu = pset->points[0]->d_unorder;
    double volume;
    point_t *p1 = NULL, *p2 =NULL;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size();
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

    int Psize = cPoint->points.size();
    int **PIndex = new int*[Psize];
    for(int i = 0; i < Psize; ++i)
    {
        cPoint->points[i]->index = i;
        PIndex[i] = new int[Psize];
        for(int j = 0; j < Psize; ++j)
            if(i == j)
                PIndex[i][j] = 1;
            else
                PIndex[i][j] = 0;
    }


    //interaction
    while (cPoint->points.size() > 1)
    {
        int hIndex, heven = -1;
        for(int i = 0; i < candHyper.size(); ++i)
        {
            int pri = candHyper[i]->priority(cPoint, dis, Beta);
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
            Qcount++;
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);
            cout << "\nThe " << Qcount <<"th question\n";
            int option = realSet->show_to_user_unanswer(p1->id, p2->id);
            if(option != 3)
            {
                if (dimu > 1)
                {
                    if (option == 2)
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
                }
                else
                {
                    double bd = p1->bound(p2, R->expDim->attr[0]);
                    if (option == 2)
                    {
                        if (p1->attr[1] < p2->attr[1])
                        {
                            R->ext_pts[1]->attr[0] = bd;
                        } else
                        {
                            R->ext_pts[0]->attr[0] = bd;
                        }
                    } else
                    {
                        if (p1->attr[1] < p2->attr[1])
                        {
                            R->ext_pts[0]->attr[0] = bd;
                        } else
                        {
                            R->ext_pts[1]->attr[0] = bd;
                        }
                    }
                }


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
            }
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


        if (cPoint->points.size() == 1 || candHyper.size() < 1)
        {
            for(int k = 0; k < cPoint->points.size(); ++k)
                TID.push_back(cPoint->points[k]->id);
            realSet->realPrint("BS", Qcount, TID, fp);
            return;
        }
    }

}

