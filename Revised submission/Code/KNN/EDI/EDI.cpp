#include "EDI.h"

struct Cmp
{
    bool operator()(const point_t* left, const point_t* right) const
    {
        return left->traAttr[0] < right->traAttr[0];
    }
};


void EDI(point_set *pset, point_t *e, double Beta, double gamma, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
    int Qcount = 0, dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder, batch = 0, first = 0;
    double volume;
    hyperplane_set *R = new hyperplane_set(pset, volume);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size(), LastCsize = pointSize;
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
        ++batch;
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
                    //cout << "candPointSize: " << candPointSize << "  " << "priority: " << pri << "\n";
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
                    cPoint->points[0]->printResult(out_cp, "EDI", Qcount, t1, preTime, mem_baseline, pointSize, type);
                    cout << batch << "\n";
                    out_cp.close();
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
        if(first == 0) //finish the preprocessing
        {
            preTime = timeCost(t1);
            first = 1;
        }

        //interaction
        double beforesize = candPoint.size();
        double num_interaction = 0;
        while (candPoint.size() > gamma || num_interaction == 0)//candPoint.size()/beforesize > 0.5 &&
        {
            ++Qcount; ++num_interaction;
            int middle = candPoint.size() / 2 - 1;
            point_t *p1 = candPoint[middle], *p2 = candPoint[middle + 1];
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);

            if (dist1 < dist2)
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

                if (dist1 < dist2)
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
                int d = p1->d;
                double BD = p1->bound(p2, R->expDim);
                if(dist1 > dist2)
                {
                    if (p1->attr[d - 1] < p2->attr[d - 1])
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    }
                    cPoint->prunePt(p1);
                }
                else
                {
                    if (p1->attr[d - 1] < p2->attr[d - 1])
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    }
                    cPoint->prunePt(p2);
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
            //cout << candPoint.size() << "\n";
            printMiddleResult(out_cp, t1, preTime, Qcount, (cPoint->points.size()/pointSize) * 100,
                              (1 - cPoint->points.size()/LastCsize) * 100, mem_baseline, type);
            LastCsize = cPoint->points.size();
        }


    }
    cPoint->points[0]->printResult(out_cp, "EDI", Qcount, t1, preTime, mem_baseline, pointSize, type);
    out_cp.close();
}


void EDITopk(point_set *origset, point_t *e, double Beta, double gamma, int k, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1;
    gettimeofday(&t1, 0);
    double preTime;
    int Qcount = 0, dimo = origset->points[0]->d_order, dimu = origset->points[0]->d_unorder, batch = 0, first = 0;
    double volume;
    hyperplane_set *R = new hyperplane_set(origset, volume);
    R->nearestSkyband(origset, k);
    double pointSize = origset->points.size(), LastCsize = pointSize;
    point_set *maintainSet = new point_set(origset);
    point_set *resultSet = new point_set();

    while (resultSet->points.size() < k)
    {
        point_set* pset = new point_set(maintainSet);
        R->nearestSkyline(pset);
        std::vector<hyperplane *> candHyper;
        point_set *cPoint = new point_set();
        std::vector<hyperplane_set *> partitionSet;
        hyperplane *h;
        std::queue<int> candP;
        R->update_expDim();
        int tt = pset->findClosest(R->expDim);
        candP.push(tt);
        pset->points[tt]->value = 1;
        while (candP.size() > 0)
        {
            int index = candP.front();
            candP.pop();
            if (R->find_boundary(pset, index, candP, candHyper, partitionSet))
            {
                //pset->points[index]->print();
                cPoint->points.push_back(pset->points[index]);
            }
        }

        double dis = (pow(volume / cPoint->points.size(), 1.0 / dimu)) / 2;

       while (cPoint->points.size() > 1)
    {
        ++batch;
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
                    //cout << "candPointSize: " << candPointSize << "  " << "priority: " << pri << "\n";
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
                    cPoint->points[0]->printResult(out_cp, "EDI", Qcount, t1, preTime, mem_baseline, pointSize, type);
                    cout << batch << "\n";
                    out_cp.close();
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
        if(first == 0) //finish the preprocessing
        {
            preTime = timeCost(t1);
            first = 1;
        }

        //interaction
        double beforesize = candPoint.size();
        double num_interaction = 0;
        while (candPoint.size() > gamma || num_interaction == 0)//candPoint.size()/beforesize > 0.5 &&
        {
            ++Qcount; ++num_interaction;
            int middle = candPoint.size() / 2 - 1;
            point_t *p1 = candPoint[middle], *p2 = candPoint[middle + 1];
            double dist1 = p1->distance(e);
            double dist2 = p2->distance(e);

            if (dist1 < dist2)
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

                if (dist1 < dist2)
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
                int d = p1->d;
                double BD = p1->bound(p2, R->expDim);
                if(dist1 > dist2)
                {
                    if (p1->attr[d - 1] < p2->attr[d - 1])
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    }
                    cPoint->prunePt(p1);
                }
                else
                {
                    if (p1->attr[d - 1] < p2->attr[d - 1])
                    {
                        R->ext_pts[0]->attr[0] = BD;
                    } else
                    {
                        R->ext_pts[1]->attr[0] = BD;
                    }
                    cPoint->prunePt(p2);
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
            //cout << candPoint.size() << "\n";
            printMiddleResult(out_cp, t1, preTime, Qcount, (cPoint->points.size()/pointSize) * 100,
                              (1 - cPoint->points.size()/LastCsize) * 100, mem_baseline, type);
            LastCsize = cPoint->points.size();
        }


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
    resultSet->printResult(out_cp, "EDITopk", Qcount, t1, preTime, mem_baseline, pointSize, type);
    out_cp.close();
}

