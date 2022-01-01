#include "UHRandom.h"

/**
 * @brief Algorithm UH-Random
 * @param pset      The dataset
 * @param realSet   The real dataset
 * @param e         The expected point (for test)
 * @param TID       The recommended point ID
 * @param fp        The file for recording results
 * @return
 */
void UH_Random(point_set* pset, point_set* realSet, point_t* e, int &TID, std::ofstream &fp)
{
	vector<int> C_idx; //index of points (position of point)
	int dim = pset->points[0]->d, dimu = pset->points[0]->d_unorder;
	hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    for(int i = 0; i < pset->points.size(); i++)
        C_idx.push_back(i);
    int index1, index2, Qcount = 0;;
	point_t *point1 = NULL, *point2 = NULL;
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
        cout << "\nThe " << Qcount <<"th question\n";
        int option = realSet->show_to_user(point1->id, point2->id);
        if(dimu > 1)
        {
            if (option == 2)
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
            if (option == 2)
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
    }
    TID = pset->points[C_idx[0]]->id;
	realSet->points[TID]->final_result("UH-Random", Qcount, fp);
}




/**
 * @brief Algorithm UH-Random
 * @param pset      The dataset
 * @param realSet   The real dataset
 * @param e         The expected point (for test)
 * @param TID       The recommended point ID
 * @param fp        The file for recording results
 * @return
 */
void UH_Random_unanswer(point_set* pset, point_set* realSet, point_t* e, std::vector<int> &TID, std::ofstream &fp)
{
    vector<int> C_idx; //index of points (position of point)
    int dim = pset->points[0]->d, dimu = pset->points[0]->d_unorder;
    hyperplane_set *R = new hyperplane_set(pset);
    R->nearestSkyline(pset);
    for(int i = 0; i < pset->points.size(); i++)
        C_idx.push_back(i);
    int index1, index2, Qcount = 0;;
    point_t *point1 = NULL, *point2 = NULL;

    int Psize = pset->points.size();
    int **PIndex = new int*[Psize];
    for(int i = 0; i < Psize; ++i)
    {
        PIndex[i] = new int[Psize];
        for(int j = 0; j < Psize; ++j)
            if(i == j)
                PIndex[i][j] = 1;
            else
                PIndex[i][j] = 0;
    }

    // interactively reduce the candidate set and shrink R
    while (C_idx.size() > 1)
    {
        bool is_exist = false;
        for(int i = 0; i < C_idx.size(); ++i)
        {
            for (int j = 0; j < C_idx.size(); ++j)
            {
                if (i != j && PIndex[pset->points[C_idx[i]]->index][pset->points[C_idx[j]]->index] == 0)
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

        Qcount++;
        index1 = rand() % C_idx.size();
        point1 = pset->points[C_idx[index1]];
        index2 = rand() % C_idx.size();
        point2 = pset->points[C_idx[index2]];
        while (point1->id == point2->id || PIndex[point1->index][point2->index] == 1)
        {
            index2 = rand() % C_idx.size();
            point2 = pset->points[C_idx[index2]];
        }

        cout << "\nThe " << Qcount << "th question\n";
        int option = realSet->show_to_user_unanswer(point1->id, point2->id);
        if (option == 3)
        {
            PIndex[point1->index][point2->index] = 1;
            PIndex[point2->index][point1->index] = 1;
        }
        else
        {
            if (dimu > 1)
            {
                if (option == 2)
                {
                    R->hyperplanes.push_back(new hyperplane(point1, point2, R->expDim));
                    C_idx.erase(C_idx.begin() + index1);
                } else
                {
                    R->hyperplanes.push_back(new hyperplane(point2, point1, R->expDim));
                    C_idx.erase(C_idx.begin() + index2);
                }
                R->set_ext_pts();
            } else
            {
                double bd = point1->bound(point2, R->expDim->attr[0]);
                if (option == 2)
                {
                    if (point1->attr[1] < point2->attr[1])
                    {
                        R->ext_pts[1]->attr[0] = bd;
                    } else
                    {
                        R->ext_pts[0]->attr[0] = bd;
                    }
                    C_idx.erase(C_idx.begin() + index1);
                } else
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

        }
    }
    for(int i = 0; i < C_idx.size(); ++i)
    {
        TID.push_back(pset->points[C_idx[i]]->id);
    }
    realSet->realPrint("UH-Random", Qcount, TID, fp);
}



























