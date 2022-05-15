#include "EDI.h"

struct Cmp
{
    bool operator()(const point_t* left, const point_t* right) const
    {
        return left->traAttr[0] < right->traAttr[0];
    }
};


void EDI(point_set* pset, point_set* realSet, point_t* e, double Beta, int &TID, std::ofstream &fp)
{
    int Qcount = 0, dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder;
    double volume;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size();
    std::vector<hyperplane*> candHyper;
    point_set *cPoint = new point_set();
    std::vector<hyperplane_set*> partitionSet;
    hyperplane *h;
    std::queue<int> candP;
    int tt = pset->findClosest(R->expDim);
    candP.push(tt);
    pset->points[tt]->value = 1;
    while (candP.size() > 0)
    {
        int index= candP.front(); candP.pop();
        if(R->find_boundary(pset, index, candP, candHyper, partitionSet))
        {
            //pset->points[index]->print();
            cPoint->points.push_back(pset->points[index]);
        }
    }

    double dis = (pow(volume/cPoint->points.size(), 1.0/dimu))/2;

    while (cPoint->points.size() > 1)
    {
        double vLength = 0, bd;
        int index1, index2, numPt = -INF, eteration = 0, maxPriority = -INF;

        while(numPt < 2 || eteration < 0.5 * R->ext_pts.size())
        {
            ++eteration;
            int ind1 = rand() % R->ext_pts.size(), ind2 = ind1;
            while (ind1 == ind2)
                ind2 = rand() % R->ext_pts.size();

            point_t *Expt = new point_t(dimo, dimu);
            for (int j = 0; j < dimo; ++j)
                Expt->attr[j] = R->expDim->attr[j];
            for (int j = dimo; j < dimo + dimu; ++j)
                Expt->attr[j] = R->ext_pts[ind1]->attr[j - dimo];
            int extBest = cPoint->findClosest(Expt);
            for (int j = dimo; j < dimo + dimu; ++j)
                Expt->attr[j] = R->ext_pts[ind2]->attr[j - dimo];
            if (extBest != cPoint->findClosest(Expt))
            {
                //transform the coordinates of points based on a vector
                point_t *startExpt = new point_t(dimo, dimu);
                for (int i = 0; i < dimo; ++i)
                    startExpt->attr[i] = R->expDim->attr[i];
                for (int i = dimo; i < dimo + dimu; ++i)
                    startExpt->attr[i] = R->ext_pts[ind1]->attr[i - dimo];

                point_t *vec = new point_t(dimu), *vecNorm = new point_t(dimu);
                for (int i = 0; i < dimu; ++i)
                {
                    vec->attr[i] = R->ext_pts[ind2]->attr[i] - R->ext_pts[ind1]->attr[i];
                    vLength += vec->attr[i] * vec->attr[i];
                }
                vLength = sqrt(vLength);
                for (int i = 0; i < dimu; ++i)
                {
                    vecNorm->attr[i] = vec->attr[i] / vLength;
                }
                cPoint->transformAttr(startExpt, vec, vecNorm);

                //partition the expected space
                sort(cPoint->points.begin(), cPoint->points.end(), Cmp());
                std::vector<double> boundary;
                boundary.push_back(0);
                std::vector<point_t *> candPoint;
                candPoint.push_back(cPoint->points[0]);
                for (int i = 1; i < cPoint->points.size(); ++i)
                {
                    bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
                    while (bd < boundary[boundary.size() - 1])
                    {
                        candPoint.pop_back();
                        boundary.pop_back();
                        if (boundary.size() == 0)
                            break;
                        bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
                    }
                    if (bd < vLength)
                    {
                        if (bd < 0)
                            boundary.push_back(0);
                        else
                            boundary.push_back(bd);
                        candPoint.push_back(cPoint->points[i]);
                    }
                }

                int candPointSize = candPoint.size();
                if(candPointSize > 1)
                {
                    int candMiddle = candPointSize / 2 - 1;
                    hyperplane *hhh = new hyperplane(candPoint[candMiddle], candPoint[candMiddle + 1]);
                    int pri = hhh->priority(cPoint, dis, Beta);
                    if (pri > maxPriority) //(candPointSize > numPt || (candPointSize == numPt && pri > maxPriority))
                    {
                        index1 = ind1;
                        index2 = ind2;
                        numPt = candPointSize;
                        maxPriority = pri;
                    }
                }
            }
            else
            {
                if(R->ext_pts.size() == 2)
                {
                    TID = cPoint->points[0]->id;
                    realSet->points[TID]->final_result("EDI", Qcount, fp);
                    return;
                }
            }
        }


        //transform the coordinates of points based on a vector
        point_t *startExpt = new point_t(dimo, dimu);
        for(int i = 0; i < dimo; ++i)
            startExpt->attr[i] = R->expDim->attr[i];
        for(int i = dimo; i < dimo + dimu; ++i)
            startExpt->attr[i] = R->ext_pts[index1]->attr[i - dimo];

        point_t *vec = new point_t(dimu), *vecNorm = new point_t(dimu);
        for(int i = 0; i < dimu; ++i)
        {
            vec->attr[i] = R->ext_pts[index2]->attr[i] - R->ext_pts[index1]->attr[i];
            vLength += vec->attr[i] * vec->attr[i];
        }
        vLength = sqrt(vLength);
        for(int i = 0; i < dimu; ++i)
        {
            vecNorm->attr[i] = vec->attr[i] / vLength;
        }
        cPoint->transformAttr(startExpt, vec, vecNorm);

        //partition the expected space
        sort(cPoint->points.begin(), cPoint->points.end(), Cmp());
        std::vector<double> boundary;
        boundary.push_back(0);
        std::vector<point_t *> candPoint;
        candPoint.push_back(cPoint->points[0]);
        for (int i = 1; i < cPoint->points.size(); ++i)
        {
            bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
            while (bd < boundary[boundary.size() - 1])
            {
                candPoint.pop_back();
                boundary.pop_back();
                if (boundary.size() == 0)
                    break;
                bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
            }
            if (bd < vLength)
            {
                if (bd < 0)
                    boundary.push_back(0);
                else
                    boundary.push_back(bd);
                candPoint.push_back(cPoint->points[i]);
            }
        }

        //interaction
        double beforesize = candPoint.size();
        double num_interaction = 0;
        while (candPoint.size() > 4 || num_interaction == 0)
        {
            ++Qcount; ++num_interaction;
            int middle = candPoint.size() / 2 - 1;
            point_t *p1 = candPoint[middle], *p2 = candPoint[middle + 1];
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);
            cout << "\nThe " << Qcount <<"th question\n";
            int option = realSet->show_to_user(p1->id, p2->id);

            if (option == 1)
            {
                point_set *pps = new point_set();
                for (int i = 0; i <= middle; ++i)
                {
                    pps->points.push_back(candPoint[i]);
                }
                candPoint = pps->points;
            }
            else
            {
                point_set *pps = new point_set();
                for (int i = middle + 1; i < candPoint.size(); ++i)
                {
                    pps->points.push_back(candPoint[i]);
                }
                candPoint = pps->points;
            }

            if(dimu > 1)
            {

                if (option == 1)
                {
                    h = new hyperplane(p2, p1, R->expDim);
                    cPoint->prunePt(p2);
                }
                else
                {
                    h = new hyperplane(p1, p2, R->expDim);
                    cPoint->prunePt(p1);
                }
                R->hyperplanes.push_back(h);
                R->set_ext_pts();

            }
            else
            {
                double BD = p1->bound(p2, R->expDim->attr[0]);
                if(dist1 > dist2) //(option == 2)
                {
                    if (p1->attr[1] < p2->attr[1])
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    }
                }
                else
                {
                    if (p1->attr[1] < p2->attr[1])
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    }
                }
            }

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
    }
    TID = cPoint->points[0]->id;
    realSet->points[TID]->final_result("EDI", Qcount, fp);
    return;
}




void EDI_unasnwer(point_set* pset, point_set* realSet, point_t* e, double Beta, std::vector<int> &TID, std::ofstream &fp)
{
    int Qcount = 0, dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder;
    double volume;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size();
    std::vector<hyperplane*> candHyper;
    point_set *cPoint = new point_set();
    std::vector<hyperplane_set*> partitionSet;
    hyperplane *h;
    std::queue<int> candP;
    int tt = pset->findClosest(R->expDim);
    candP.push(tt);
    pset->points[tt]->value = 1;
    while (candP.size() > 0)
    {
        int index= candP.front(); candP.pop();
        if(R->find_boundary(pset, index, candP, candHyper, partitionSet))
        {
            //pset->points[index]->print();
            cPoint->points.push_back(pset->points[index]);
        }
    }

    double dis = (pow(volume/cPoint->points.size(), 1.0/dimu))/2;

    int **expIndex = new int*[R->ext_pts.size()];
    for(int i = 0; i < R->ext_pts.size(); ++i)
    {
        expIndex[i] = new int[R->ext_pts.size()];
        for(int j = 0; j < R->ext_pts.size(); ++j)
            expIndex[i][j] = 0;
    }


    while (cPoint->points.size() > 1)
    {
        double vLength = 0, bd;
        int index1, index2, numPt = -INF, eteration = 0, maxPriority = -INF;

        while(numPt < 2 || eteration < 0.5 * R->ext_pts.size())
        {
            ++eteration;
            int ind1 = rand() % R->ext_pts.size(), ind2 = ind1;
            while (ind1 == ind2)
                ind2 = rand() % R->ext_pts.size();

            point_t *Expt = new point_t(dimo, dimu);
            for (int j = 0; j < dimo; ++j)
                Expt->attr[j] = R->expDim->attr[j];
            for (int j = dimo; j < dimo + dimu; ++j)
                Expt->attr[j] = R->ext_pts[ind1]->attr[j - dimo];
            int extBest = cPoint->findClosest(Expt);
            for (int j = dimo; j < dimo + dimu; ++j)
                Expt->attr[j] = R->ext_pts[ind2]->attr[j - dimo];
            if (extBest != cPoint->findClosest(Expt) && expIndex[ind1][ind2] == 0)
            {
                //transform the coordinates of points based on a vector
                point_t *startExpt = new point_t(dimo, dimu);
                for (int i = 0; i < dimo; ++i)
                    startExpt->attr[i] = R->expDim->attr[i];
                for (int i = dimo; i < dimo + dimu; ++i)
                    startExpt->attr[i] = R->ext_pts[ind1]->attr[i - dimo];

                point_t *vec = new point_t(dimu), *vecNorm = new point_t(dimu);
                for (int i = 0; i < dimu; ++i)
                {
                    vec->attr[i] = R->ext_pts[ind2]->attr[i] - R->ext_pts[ind1]->attr[i];
                    vLength += vec->attr[i] * vec->attr[i];
                }
                vLength = sqrt(vLength);
                for (int i = 0; i < dimu; ++i)
                {
                    vecNorm->attr[i] = vec->attr[i] / vLength;
                }
                cPoint->transformAttr(startExpt, vec, vecNorm);

                //partition the expected space
                sort(cPoint->points.begin(), cPoint->points.end(), Cmp());
                std::vector<double> boundary;
                boundary.push_back(0);
                std::vector<point_t *> candPoint;
                candPoint.push_back(cPoint->points[0]);
                for (int i = 1; i < cPoint->points.size(); ++i)
                {
                    bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
                    while (bd < boundary[boundary.size() - 1])
                    {
                        candPoint.pop_back();
                        boundary.pop_back();
                        if (boundary.size() == 0)
                            break;
                        bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
                    }
                    if (bd < vLength)
                    {
                        if (bd < 0)
                            boundary.push_back(0);
                        else
                            boundary.push_back(bd);
                        candPoint.push_back(cPoint->points[i]);
                    }
                }

                int candPointSize = candPoint.size();
                if(candPointSize > 1)
                {
                    int candMiddle = candPointSize / 2 - 1;
                    hyperplane *hhh = new hyperplane(candPoint[candMiddle], candPoint[candMiddle + 1]);
                    int pri = hhh->priority(cPoint, dis, Beta);
                    if (pri > maxPriority) //(candPointSize > numPt || (candPointSize == numPt && pri > maxPriority))
                    {
                        index1 = ind1;
                        index2 = ind2;
                        numPt = candPointSize;
                        maxPriority = pri;
                    }
                }
            }

        }

        //transform the coordinates of points based on a vector
        point_t *startExpt = new point_t(dimo, dimu);
        for(int i = 0; i < dimo; ++i)
            startExpt->attr[i] = R->expDim->attr[i];
        for(int i = dimo; i < dimo + dimu; ++i)
            startExpt->attr[i] = R->ext_pts[index1]->attr[i - dimo];

        point_t *vec = new point_t(dimu), *vecNorm = new point_t(dimu);
        for(int i = 0; i < dimu; ++i)
        {
            vec->attr[i] = R->ext_pts[index2]->attr[i] - R->ext_pts[index1]->attr[i];
            vLength += vec->attr[i] * vec->attr[i];
        }
        vLength = sqrt(vLength);
        for(int i = 0; i < dimu; ++i)
        {
            vecNorm->attr[i] = vec->attr[i] / vLength;
        }
        cPoint->transformAttr(startExpt, vec, vecNorm);

        //partition the expected space
        sort(cPoint->points.begin(), cPoint->points.end(), Cmp());
        std::vector<double> boundary;
        boundary.push_back(0);
        std::vector<point_t *> candPoint;
        candPoint.push_back(cPoint->points[0]);
        for (int i = 1; i < cPoint->points.size(); ++i)
        {
            bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
            while (bd < boundary[boundary.size() - 1])
            {
                candPoint.pop_back();
                boundary.pop_back();
                if (boundary.size() == 0)
                    break;
                bd = candPoint[candPoint.size() - 1]->tranBound(cPoint->points[i]);
            }
            if (bd < vLength)
            {
                if (bd < 0)
                    boundary.push_back(0);
                else
                    boundary.push_back(bd);
                candPoint.push_back(cPoint->points[i]);
            }
        }

        //interaction
        double beforesize = candPoint.size();
        double num_interaction = 0;
        while (candPoint.size() > 4 || num_interaction == 0)
        {
            ++Qcount;
            ++num_interaction;
            int middle = candPoint.size() / 2 - 1;
            point_t *p1 = candPoint[middle], *p2 = candPoint[middle + 1];
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);
            cout << "\nThe " << Qcount << "th question\n";
            int option = realSet->show_to_user_unanswer(p1->id, p2->id);
            if (option == 3)
            {
                expIndex[index1][index2] = 1;
                expIndex[index2][index1] = 1;

                int dataPoint1, dataPoint2;
                point_t *Expt1 = new point_t(dimo, dimu);
                point_t *Expt2 = new point_t(dimo, dimu);
                for (int j = 0; j < dimo; ++j)
                {
                    Expt1->attr[j] = R->expDim->attr[j];
                    Expt2->attr[j] = R->expDim->attr[j];
                }
                for (int j = dimo; j < dimo + dimu; ++j)
                {
                    Expt1->attr[j] = R->ext_pts[index1]->attr[j - dimo];
                    Expt2->attr[j] = R->ext_pts[index2]->attr[j - dimo];
                }
                dataPoint1 = cPoint->findClosest(Expt1);
                dataPoint2 = cPoint->findClosest(Expt2);

                for(int i = 0; i < R->ext_pts.size(); ++i)
                {
                    for(int j = 0; j < R->ext_pts.size(); ++j)
                    {
                        if(i != j)
                        {
                            int dataP1, dataP2;
                            for (int k = dimo; k < dimo + dimu; ++k)
                            {
                                Expt1->attr[k] = R->ext_pts[i]->attr[k - dimo];
                                Expt2->attr[k] = R->ext_pts[j]->attr[k - dimo];
                            }
                            dataP1 = cPoint->findClosest(Expt1);
                            dataP2 = cPoint->findClosest(Expt2);
                            if(dataP1 == dataP2)
                            {
                                expIndex[i][j] = 1;
                                expIndex[j][i] = 1;
                            }
                            else if(dataPoint1 == dataP1 && dataPoint2 == dataP2 || dataPoint1 == dataP2 && dataPoint2 == dataP1)
                            {
                                expIndex[i][j] = 1;
                                expIndex[j][i] = 1;
                            }
                        }
                    }
                }

                /**
                for(int i = 0; i < R->ext_pts.size(); ++i)
                {
                    for(int j = 0; j < R->ext_pts.size(); ++j)
                    {
                        cout << expIndex[i][j] << "  ";
                    }
                    cout <<"\n";
                }
                */
                break;
            }
            else
            {
                if (option == 1)
                {
                    point_set *pps = new point_set();
                    for (int i = 0; i <= middle; ++i)
                    {
                        pps->points.push_back(candPoint[i]);
                    }
                    candPoint = pps->points;
                } else
                {
                    point_set *pps = new point_set();
                    for (int i = middle + 1; i < candPoint.size(); ++i)
                    {
                        pps->points.push_back(candPoint[i]);
                    }
                    candPoint = pps->points;
                }

                if (dimu > 1)
                {

                    if (option == 1)
                    {
                        h = new hyperplane(p2, p1, R->expDim);
                        cPoint->prunePt(p2);
                    } else
                    {
                        h = new hyperplane(p1, p2, R->expDim);
                        cPoint->prunePt(p1);
                    }
                    R->hyperplanes.push_back(h);
                    R->set_ext_pts();

                    expIndex = new int*[R->ext_pts.size()];
                    for(int i = 0; i < R->ext_pts.size(); ++i)
                    {
                        expIndex[i] = new int[R->ext_pts.size()];
                        for(int j = 0; j < R->ext_pts.size(); ++j)
                            expIndex[i][j] = 0;
                    }

                }
                else
                {
                    double BD = p1->bound(p2, R->expDim->attr[0]);
                    if (dist1 > dist2) //(option == 2)
                    {
                        if (p1->attr[1] < p2->attr[1])
                        {
                            R->ext_pts[1]->attr[0] = BD;
                        } else
                        {
                            R->ext_pts[0]->attr[0] = BD;
                        }
                    } else
                    {
                        if (p1->attr[1] < p2->attr[1])
                        {
                            R->ext_pts[0]->attr[0] = BD;
                        } else
                        {
                            R->ext_pts[1]->attr[0] = BD;
                        }
                    }

                    expIndex = new int*[R->ext_pts.size()];
                    for(int i = 0; i < R->ext_pts.size(); ++i)
                    {
                        expIndex[i] = new int[R->ext_pts.size()];
                        for(int j = 0; j < R->ext_pts.size(); ++j)
                            expIndex[i][j] = 0;
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

        bool is_exist = false;
        for(int i = 0; i < R->ext_pts.size(); ++i)
        {
            for(int j = 0; j < R->ext_pts.size(); ++j)
            {
                if(i != j && expIndex[i][j] == 0)
                {
                    is_exist = true;
                    break;
                }
            }
            if(is_exist)
                break;
        }
        if(!is_exist)
            break;

    }
    for(int k = 0; k <cPoint->points.size(); ++k)
        TID.push_back(cPoint->points[k]->id);
    realSet->realPrint("EDI", Qcount, TID, fp);
    return;
}

