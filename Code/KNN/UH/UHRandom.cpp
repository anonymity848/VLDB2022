#include "UHRandom.h"

/**
 * @brief Algorithm UH-Random
 * @param pset  The point set
 * @param e     The user's expected point
 */
void UH_Random(point_set* pset, point_t* e, long mem_baseline)
{
    //std::ofstream out_cp("../../result.txt");
    timeval t1; gettimeofday(&t1, 0);
    double preTime;
	vector<int> C_idx; //index of points (position of point)
	int dim = pset->points[0]->d, dimu = pset->points[0]->d_unorder;
	hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    for(int i = 0; i < pset->points.size(); i++)
        C_idx.push_back(i);
    int index1, index2, Qcount = 0;;
	point_t *point1 = NULL, *point2 = NULL;
    preTime = timeCost(t1);
    double pointSize = pset->points.size();


	// interactively reduce the candidate set and shrink R
	while (C_idx.size() > 1)
	{
		Qcount++;
		index1 = rand() % C_idx.size();
		point1 = pset->points[C_idx[index1]];
		index2 = rand() % C_idx.size();
        point2 = pset->points[C_idx[index2]];
        while (point1->id == point2->id)
        {
            index2 = rand() % C_idx.size();
            point2 = pset->points[C_idx[index2]];
        }

        double dist1 = point1->distance(e);
        double dist2 = point2->distance(e);
        if(dimu > 1)
        {
            if (dist1 > dist2)
            {
                R->hyperplanes.push_back(new hyperplane(point1, point2, R->expDim));
                C_idx.erase(C_idx.begin() + index1);
            } else
            {
                R->hyperplanes.push_back(new hyperplane(point2, point1, R->expDim));
                C_idx.erase(C_idx.begin() + index2);
            }
            R->set_ext_pts();
        }
        else
        {
            double bd = point1->bound(point2, R->expDim->attr[0]);
            if(dist1 > dist2)
            {
                if (point1->attr[1] < point2->attr[1])
                {
                    R->ext_pts[1]->attr[0] = bd;
                }
                else
                {
                    R->ext_pts[0]->attr[0] = bd;
                }
                C_idx.erase(C_idx.begin() + index1);
            }
            else
            {
                if (point1->attr[1] < point2->attr[1])
                {
                    R->ext_pts[0]->attr[0] = bd;
                } else
                {
                    R->ext_pts[1]->attr[0] = bd;
                }
                C_idx.erase(C_idx.begin() + index2);
            }
        }

		//update candidate set
		rtree_pruning(pset, C_idx, R);
		//printMiddleResult(out_cp, t1, preTime, Qcount, (C_idx.size()/pointSize) * 100, mem_baseline);
    }

	pset->points[C_idx[0]]->printResult("UH-Random", Qcount, t1, preTime, mem_baseline);
    //out_cp.close();
}






























