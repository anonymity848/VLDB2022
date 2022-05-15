#include "DI.h"


struct aCmp
{
    bool operator()(const point_t* left, const point_t* right) const
    {
        return left->attr[1] < right->attr[1];
    }
};


/**
 * @brief Algorithm DI
 * @param pset The point set
 * @param e    The user's expected point
 */
void DI(point_set *pset, point_t *e, long mem_baseline)
{
    //std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0); double preTime;
    int Qcount = 0;
    double maxx = -1, maxy = -1, bd;
    hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    double pointSize = pset->points.size();

    for(int i = 0; i < pset->points.size(); ++i)
    {
        if(pset->points[i]->attr[0] > maxx)
            maxx = pset->points[i]->attr[0];
        if(pset->points[i]->attr[1] > maxy)
            maxy = pset->points[i]->attr[1];
    }

    //divide space
    sort(pset->points.begin(), pset->points.end(), aCmp());
    std::vector<double> boundary; boundary.push_back(0);
    std::vector<point_t*> candPoint; candPoint.push_back(pset->points[0]);
    for(int i = 1; i < pset->points.size(); ++i)
    {
        bd = candPoint[candPoint.size() - 1]->bound(pset->points[i], maxx);
        while(bd < boundary[boundary.size() - 1])
        {
            candPoint.pop_back();
            boundary.pop_back();
            if(boundary.size() == 0)
                break;
            bd = candPoint[candPoint.size() - 1]->bound(pset->points[i], maxx);
        }
        if(bd < maxy)
        {
            if(bd < 0)
                boundary.push_back(0);
            else
                boundary.push_back(bd);
            candPoint.push_back(pset->points[i]);
        }

    }
    preTime = timeCost(t1);

    //interaction
    while (candPoint.size() > 1)
    {
        ++Qcount;
        int middle = candPoint.size() / 2 - 1;
        point_t *p1 = candPoint[middle], *p2 = candPoint[middle + 1];
        double dist1 = p1->distance(e), dist2 = p2->distance(e);
        if(dist1 < dist2)
        {
            point_set *pps = new point_set();
            for(int i = 0; i <= middle; ++i)
            {
                pps->points.push_back(candPoint[i]);
            }
            candPoint = pps->points;
        }
        else
        {
            point_set *pps = new point_set();
            for(int i = middle + 1; i < candPoint.size(); ++i)
            {
                pps->points.push_back(candPoint[i]);
            }
            candPoint = pps->points;
        }
        //printMiddleResult(out_cp, t1, preTime, Qcount, (candPoint.size()/pointSize) * 100, mem_baseline);
    }
    candPoint[0]->printResult("DI", Qcount, t1, preTime, mem_baseline);
    //out_cp.close();
}