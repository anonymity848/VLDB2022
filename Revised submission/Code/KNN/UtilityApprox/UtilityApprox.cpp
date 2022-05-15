#include "UtilityApprox.h"


/**
 * @brief Calculated the regret ratio of utility range
 * @param L the lower bound of each dimension
 * @param U the upper bound of each dimension
 * @param D the dimension
 * @return the regret ratio
 */
double rr_bound(double *L, double *U, int D)
{
    double bound = 0;

    for (int i = 0; i < D; i++)
    {
        bound += U[i] - L[i];
    }
    return bound;
}


void UtilityApprox(point_set *pset, point_t *e, long mem_baseline, int type)
{
    std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
    int dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder;
    hyperplane_set *R = new hyperplane_set(pset);
    int M = pset->points.size();
    //double bestUtility, currUtility;
    double Qcount = 0;
    point_t *p1 = new point_t(dimo, dimu), *p2 = new point_t(dimo, dimu), *ap;
    int indexDimension = -1;
    point_t* Rpoint = R->findNearest(pset);
    preTime = timeCost(t1);
    double L1dst = R->findL1Dis(dimu);

    //Interaction
    while (L1dst >= dimu && Qcount < 10000)
    {
        ++Qcount;
        indexDimension = (indexDimension + 1) % dimu;
        double max, min, gap;
        R->findMinMax(indexDimension, min, max);
        gap = (max - min)/4;


        point_t *v = new point_t(dimo, dimu);
        for(int i = 0; i < dimo; ++i)
            v->attr[i] = R->expDim->attr[i];
        ap = R->average_point();
        for(int i = dimo; i < dimo + dimu; ++i)
            v->attr[i] = ap->attr[i - dimo];
        //v->print();
        int index = pset->findClosest(v);
        for(int i = 0; i < dimo + dimu; ++i)
        {
            p1->attr[i] = pset->points[index]->attr[i];
            p2->attr[i] = pset->points[index]->attr[i];
        }
        p1->attr[dimo + indexDimension] = min + gap;
        p2->attr[dimo + indexDimension] = min + gap * 3;


        double dist1 = p1->distance(e);
        double dist2 = p2->distance(e);
        if(dimu > 1)
        {
            hyperplane *hh;
            if (dist1 > dist2)
                hh = new hyperplane(p1, p2, R->expDim);
            else
                hh = new hyperplane(p2, p1, R->expDim);
            R->hyperplanes.push_back(hh);
            R->set_ext_pts();
            //R->print();
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
            }
        }
        //Rpoint = R->findNearest(pset);
        L1dst = R->findL1Dis(dimu);
        printMiddleResult(out_cp, t1, preTime, Qcount, 100, 100, mem_baseline, type);
    }

    Rpoint = R->findNearest(pset);
    Rpoint->printResult(out_cp, "UtilityApprox", Qcount, t1, preTime, mem_baseline, M, type);
    out_cp.close();
    return;
}