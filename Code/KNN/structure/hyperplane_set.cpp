#include "hyperplane_set.h"
#include "../Others/lp.h"
#include "../Others/pruning.h"
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <algorithm>

/**
 * @brief Constructor
 */
hyperplane_set::hyperplane_set(){}

/**
 * @brief       Constructor: initial some hyperplanes so that min <= u[i] <= max
 * @param       The point set
 */
hyperplane_set::hyperplane_set(point_set* p_set)
{
    int dimo = p_set->points[0]->d_order, dimu = p_set->points[0]->d_unorder, dim = dimo + dimu;

    double *max =new double[dim], *min = new double[dim];
    //find min, max of each dimension
    //**********************************************************
    for(int i = 0; i < dim; ++i)
    {
        max[i] = 0; min[i] = INF;
    }
    for(int i = 0; i < p_set->points.size(); ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            if(p_set->points[i]->attr[j] < min[j])
                min[j] = p_set->points[i]->attr[j];
            if(p_set->points[i]->attr[j] > max[j])
                max[j] = p_set->points[i]->attr[j];
        }
    }
    //**********************************************************

    //initialize expDim
    //**********************************************************
    expDim = new point_t(dimo, dimu);
    for(int j = 0; j < dimo; ++j)
        expDim->attr[j] = max[j];
    for(int j = dimo; j < dim; ++j)
        expDim->attr[j] = (min[j] + max[j])/2;
    //**********************************************************


    if(dimu > 1)
    {
        //initialize hyperplanes
        hyperplane *h;
        for(int i = 0; i < dimu; ++i)
        {
            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = -1;
                else
                    h->norm[j] = 0;
            }
            h->offset = min[dimo + i];
            hyperplanes.push_back(h);

            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = 1;
                else
                    h->norm[j] = 0;
            }
            h->offset = -max[dimo + i];
            hyperplanes.push_back(h);
        }

        //set extreme points
        this->set_ext_pts();
    }
    else if (dimu == 1)
    {
        //upper
        point_t *p = new point_t(dimu);
        p->attr[0] = max[dimo];
        ext_pts.push_back(p);
        //lower
        p = new point_t(dimu);
        p->attr[0] = min[dimo];
        ext_pts.push_back(p);
    }
}

/**
 * @brief       Constructor: initial some hyperplanes so that min <= u[i] <= max
 * @param       The point set
 */
hyperplane_set::hyperplane_set(point_set* p_set, double &volume)
{
    int dimo = p_set->points[0]->d_order, dimu = p_set->points[0]->d_unorder, dim = dimo + dimu;

    double *max =new double[dim], *min = new double[dim];
    //find min, max of each dimension
    //**********************************************************
    for(int i = 0; i < dim; ++i)
    {
        max[i] = 0; min[i] = INF;
    }
    for(int i = 0; i < p_set->points.size(); ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            if(p_set->points[i]->attr[j] < min[j])
                min[j] = p_set->points[i]->attr[j];
            if(p_set->points[i]->attr[j] > max[j])
                max[j] = p_set->points[i]->attr[j];
        }
    }
    //**********************************************************

    //initialize expDim
    //**********************************************************
    expDim = new point_t(dimo, dimu);
    for(int j = 0; j < dimo; ++j)
        expDim->attr[j] = max[j];
    for(int j = dimo; j < dim; ++j)
        expDim->attr[j] = (min[j] + max[j])/2;
    //**********************************************************

    //calculate the volume
    volume = 1;
    for(int j = dimo; j < dim; ++j)
    {
        volume *= (max[j] - min[j]);
    }

    if(dimu > 1)
    {
        //initialize hyperplanes
        hyperplane *h;
        for(int i = 0; i < dimu; ++i)
        {
            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = -1;
                else
                    h->norm[j] = 0;
            }
            h->offset = min[dimo + i];
            hyperplanes.push_back(h);

            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = 1;
                else
                    h->norm[j] = 0;
            }
            h->offset = -max[dimo + i];
            hyperplanes.push_back(h);
        }

        //set extreme points
        this->set_ext_pts();
    }
    else if (dimu == 1)
    {
        //upper
        point_t *p = new point_t(dimu);
        p->attr[0] = max[dimo];
        ext_pts.push_back(p);
        //lower
        p = new point_t(dimu);
        p->attr[0] = min[dimo];
        ext_pts.push_back(p);
    }
}


/**
 * @brief       Constructor: initial some hyperplanes so that min <= u[i] <= max
 * @param       The point set
 */
hyperplane_set::hyperplane_set(point_set* p_set, double *max, double *min)
{
    int dimo = p_set->points[0]->d_order, dimu = p_set->points[0]->d_unorder, dim = dimo + dimu;

    //find min, max of each dimension
    //**********************************************************
    for(int i = 0; i < dim; ++i)
    {
        max[i] = 0; min[i] = INF;
    }
    for(int i = 0; i < p_set->points.size(); ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            if(p_set->points[i]->attr[j] < min[j])
                min[j] = p_set->points[i]->attr[j];
            if(p_set->points[i]->attr[j] > max[j])
                max[j] = p_set->points[i]->attr[j];
        }
    }
    //**********************************************************

    //initialize expDim
    //**********************************************************
    expDim = new point_t(dimo, dimu);
    for(int j = 0; j < dimo; ++j)
        expDim->attr[j] = max[j];
    for(int j = dimo; j < dim; ++j)
        expDim->attr[j] = (min[j] + max[j])/2;
    //**********************************************************


    if(dimu > 1)
    {
        //initialize hyperplanes
        hyperplane *h;
        for(int i = 0; i < dimu; ++i)
        {
            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = -1;
                else
                    h->norm[j] = 0;
            }
            h->offset = min[dimo + i];
            hyperplanes.push_back(h);

            h = new hyperplane(dimu);
            for (int j = 0; j < dimu; ++j)
            {
                if (i == j)
                    h->norm[j] = 1;
                else
                    h->norm[j] = 0;
            }
            h->offset = -max[dimo + i];
            hyperplanes.push_back(h);
        }

        //set extreme points
        this->set_ext_pts();
    }
    else if (dimu == 1)
    {
        //upper
        point_t *p = new point_t(dimu);
        p->attr[0] = max[dimo];
        ext_pts.push_back(p);
        //lower
        p = new point_t(dimu);
        p->attr[0] = min[dimo];
        ext_pts.push_back(p);
    }
}


hyperplane_set::hyperplane_set(point_set* pset, int index, int index1) {
    int dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder, dim = dimo + dimu;
    if (pset->points[0]->d_unorder > 0) {
        //find min, max of each dimension
        double *max = new double[dim];
        double *min = new double[dim];
        for (int i = 0; i < dim; ++i) {
            max[i] = 0;
            min[i] = INF;
        }
        for (int i = 0; i < pset->points.size(); ++i) {
            for (int j = 0; j < dim; ++j) {
                if (pset->points[i]->attr[j] < min[j])
                    min[j] = pset->points[i]->attr[j];
                if (pset->points[i]->attr[j] > max[j])
                    max[j] = pset->points[i]->attr[j];
            }
        }


        //initialize expDim
        expDim = new point_t(dimo, dimu);
        for(int j = 0; j < dimo; ++j)
            expDim->attr[j] = max[j];


        //initialize hyperplanes
        hyperplane *h;
        for (int i = 0; i < dim; ++i) {
            h = new hyperplane(dim);
            for (int j = 0; j < dim; ++j) {
                if (j == i)
                    h->norm[j] = -1;
                else
                    h->norm[j] = 0;
            }
            h->offset = min[i];
            hyperplanes.push_back(h);

            h = new hyperplane(dim);
            for (int j = 0; j < dim; ++j) {
                if (j == i)
                    h->norm[j] = 1;
                else
                    h->norm[j] = 0;
            }
            h->offset = -max[i];
            hyperplanes.push_back(h);
        }

        for (int i = 0; i < pset->points.size(); ++i) {
            if (i != index && i != index1)
                hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[index], expDim));
        }

        //set extreme points
        this->set_ext_pts();
    }
}


/**
 * @brief Constructor: Initialize all the possible hyperplanes as questions
 * @param p_set
 */
hyperplane_set::hyperplane_set(point_set *p_set, point_t *exp, double RandRate)
{
    int M = p_set->points.size();
    for(int i = 0; i < M; ++i)
    {
        for(int j = i + 1; j < M; ++j)
        {
            hyperplane *h = new hyperplane(p_set->points[i], p_set->points[j], exp);
            hyperplanes.push_back(h);
        }
    }

    //reinsert
    for (int i = 0; i < M * RandRate; i++)
    {
        int n = ((int) rand()) % M;
        hyperplane *h = hyperplanes[n];
        hyperplanes.erase(hyperplanes.begin() + n);
        hyperplanes.push_back(h);
    }
}


/**
 * @brief Destructor
 *        Delete the hyperplanes and extreme points
 */
hyperplane_set::~hyperplane_set()
{
    //delete the hyperplanes
    int i = hyperplanes.size();
    while(i > 0)
    {
        delete hyperplanes[i-1];
        i--;
    }
    hyperplanes.clear();

    //delete the extreme points
    i = ext_pts.size();
    while(i > 0)
    {
        delete ext_pts[i-1];
        i--;
    }
    ext_pts.clear();

}


/**
    * @brief Prepare the file for computing the convex hull (the utility range R)
    *        via halfspace interaction
    * @param feasible_pt   A points inside R
    * @param filename      The name of the file written
    */
void hyperplane_set::write(point_t *feasible_pt, char *filename)
{
    //char filename[MAX_FILENAME_LENG];
    int dim = feasible_pt->d, size = hyperplanes.size();
    ofstream wPtr;
    wPtr.open(filename, std::ios::out);
    wPtr.setf(ios::fixed, ios::floatfield);  // set as fixed model
    wPtr.precision(6);  // set precision to 6

    // write the feasible point
    wPtr << dim <<"\n" << 1 <<"\n";
    for(int i = 0; i < dim; i++)
        wPtr << feasible_pt->attr[i] << " ";
    wPtr << "\n";

    // write the hyperplane
    wPtr << dim + 1 <<"\n" << size << "\n";//record the offset as one dimension
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            wPtr << hyperplanes[i]->norm[j] <<" ";
        }
        wPtr << hyperplanes[i]->offset <<"\n";
    }
    wPtr.close();
}

/**
 * @brief   Print the information of the hyperplane set
 *          Including hyperplanes and extreme points of the intersection(convex hull)
 */
void hyperplane_set::print()
{
    //print hyperplanes
    for (int i = 0; i < hyperplanes.size(); i++)
    {
        printf("hyperplane: ");
        hyperplanes[i]->print();
    }
    //print extreme points
    for (int i = 0; i < ext_pts.size(); i++)
    {
        printf("extreme point: ");
        ext_pts[i]->print();
    }
}

/**
 * @brief   Check whether the intersection of the hyperplane exists
 *          Set the extreme points of the hyperplane set.
 *          Refine the bounded hyperplanes
 * @return  Whether R is updated
 *          true  R is the different
 *          false R is the same as before
 */
bool hyperplane_set::set_ext_pts() {
    int M = hyperplanes.size(), size = 0;
    if (M <= 0) {
        printf("%s\n", "Error: No hyperplane in the set.");
        return false;
    }
    int dim = hyperplanes[0]->dim;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");
    point_t *feasible_pt = find_feasible();//the feasible point

    if (feasible_pt == NULL) {
        cout << "The intersection is infeasible.\n";
        return false;
    }
    write(feasible_pt, file1);//write hyperplanes and the feasible point to file1,

    //conduct half space intersection and write results to file2
    FILE *rPtr, *wPtr;
    if ((rPtr = fopen(file1, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file.\n");
        exit(0);
    }
    wPtr = (FILE *) fopen(file2, "w");
    halfspace(rPtr, wPtr);
    fclose(rPtr);
    fclose(wPtr);

    //read extreme points in file2
    if ((rPtr = fopen(file2, "r")) == NULL) {
        fprintf(stderr, "Cannot open the data file %s.\n", file2);
        exit(0);
    }

    //update the set of extreme points
    fscanf(rPtr, "%i%i", &dim, &size);
    if (size > 0) {
        //delete all the original extreme points
        while (ext_pts.size() > 0) {
            point_t *pt = ext_pts[ext_pts.size() - 1];
            ext_pts.pop_back();
            delete pt;
        }
        ext_pts.clear();

        //input new extreme points
        for (int i = 0; i < size; ++i) {
            point_t *p = new point_t(dim);
            p->d_unorder = dim;
            for (int j = 0; j < dim; ++j)
                fscanf(rPtr, "%lf", &p->attr[j]);
            ext_pts.push_back(p);
        }

        // update the set of hyperplanes
        fscanf(rPtr, "%i", &size);
        int *hid = new int[size + 1];
        for (int i = 1; i <= size; ++i)
            fscanf(rPtr, "%i", &hid[i]);
        sort(hid, hid + size + 1);
        for (int i = M - 1, count = size; i > -1; --i) {
            if (hid[count] < i)
                hyperplanes.erase(hyperplanes.begin() + i);
            else
                --count;
        }
        delete[]hid;
    }
    fclose(rPtr);
    return true;
}



bool hyperplane_set::find_boundary(point_set *pset, int index, std::queue<int> &candPoint, std::vector<hyperplane*> &candHyper, std::vector<hyperplane_set*> &partitionSet)
{
    int dim = pset->points[0]->d, dimu = pset->points[0]->d_unorder;
    //Initial a test
    if(dimu > 1) {
        int M = hyperplanes.size(), size = 0;
        char file1[MAX_FILENAME_LENG];
        sprintf(file1, "../output/hyperplane_data.txt");
        char file2[MAX_FILENAME_LENG];
        sprintf(file2, "../output/ext_pt.txt");

        hyperplane_set *cHull = new hyperplane_set();
        cHull->expDim = new point_t(expDim);
        for (int i = 0; i < M; ++i)
            cHull->hyperplanes.push_back(new hyperplane(hyperplanes[i]));
        for (int i = 0; i < pset->points.size() * 0.05; ++i) {
            if (i != index)
                cHull->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[index], expDim));
        }
        if (cHull->set_ext_pts()) {
            //cHull->print();
            for (int i = pset->points.size() * 0.05; i < pset->points.size(); ++i) {
                if (i != index) {
                    hyperplane *h = new hyperplane(pset->points[i], pset->points[index], expDim);
                    if (cHull->check_relation(h) == 0)
                        cHull->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[index], expDim));
                }
            }
            point_t *feasible_pt = cHull->find_feasible();//the feasible point

            if (feasible_pt == NULL) {
                //std::cout << "Intersection is infeasible. \n";
                return false;
            }
            cHull->write(feasible_pt, file1);//write hyperplanes and the feasible point to file1,

            //conduct half space intersection and write results to file2
            FILE *rPtr, *wPtr;
            if ((rPtr = fopen(file1, "r")) == NULL) {
                fprintf(stderr, "Cannot open the data file.\n");
                exit(0);
            }
            wPtr = (FILE *) fopen(file2, "w");
            halfspace(rPtr, wPtr);
            fclose(rPtr);
            fclose(wPtr);

            //read extreme points in file2
            if ((rPtr = fopen(file2, "r")) == NULL) {
                fprintf(stderr, "Cannot open the data file %s.\n", file2);
                exit(0);
            }

            //update the set of extreme points
            fscanf(rPtr, "%i%i", &dim, &size);
            if (size > 0)
            {
                //extreme points
                point_t *p = new point_t(dim);
                for (int i = 0; i < size; ++i)
                {
                    for (int j = 0; j < dim; ++j)
                        fscanf(rPtr, "%lf", &p->attr[j]);
                    pset->points[index]->ext.push_back(p);
                }


                // update the set of hyperplanes
                fscanf(rPtr, "%i", &size);
                int *hid = new int[size + 1];
                for (int i = 1; i <= size; ++i)
                    fscanf(rPtr, "%i", &hid[i]);
                sort(hid, hid + size + 1);
                M = cHull->hyperplanes.size();
                for (int i = M - 1, count = size; i > -1; --i)
                {
                    if (hid[count] < i)
                        cHull->hyperplanes.erase(cHull->hyperplanes.begin() + i);
                    else
                        --count;
                }
                delete[]hid;

                for (int i = 0; i < size; ++i)
                {
                    if (cHull->hyperplanes[i]->p_1 != NULL && cHull->hyperplanes[i]->p_2 != NULL)
                    {
                        if (cHull->hyperplanes[i]->p_1->hyper == 0 && cHull->hyperplanes[i]->p_2->hyper == 0)
                        {
                            candHyper.push_back(new hyperplane(cHull->hyperplanes[i]->p_1, cHull->hyperplanes[i]->p_2, expDim));
                        }
                        if (cHull->hyperplanes[i]->p_1->value == 0)
                        {
                            cHull->hyperplanes[i]->p_1->value = 1;
                            candPoint.push(cHull->hyperplanes[i]->p_1->index);
                        }
                        if (cHull->hyperplanes[i]->p_2->value == 0)
                        {
                            cHull->hyperplanes[i]->p_2->value = 1;
                            candPoint.push(cHull->hyperplanes[i]->p_2->index);
                        }
                    }
                }
                pset->points[index]->hyper = 1;
                pset->points[index]->internal = feasible_pt;
                partitionSet.push_back(cHull);
                pset->points[index]->boundaryHyperID = partitionSet.size() - 1;
            }
            fclose(rPtr);
            return true;
        } else {
            return false;
        }
    }
    else
    {
        double min = -INF, max = INF;
        int minPoint = -1, maxPoint = -1;
        for(int i = 0; i < pset->points.size(); ++i)
        {
            if(i != index)
            {
                double bd = pset->points[index]->bound(pset->points[i], expDim->attr[0]);
                if(pset->points[index]->attr[1] > pset->points[i]->attr[1] && min < bd)
                {
                    min = bd;
                    minPoint = i;
                }
                else if (pset->points[index]->attr[1] < pset->points[i]->attr[1] && max > bd)
                {
                    max = bd;
                    maxPoint = i;
                }

                if(min >= max)
                    return false;
            }
        }

        //update the set of extreme points
        point_t *p = new point_t(dimu);
        p->attr[0] = max;
        pset->points[index]->ext.push_back(p);
        p = new point_t(dimu);
        p->attr[0] = min;
        pset->points[index]->ext.push_back(p);

        // update the set of hyperplanes
        if(minPoint < pset->points.size() && minPoint >= 0)
        {
            if (pset->points[index]->hyper == 0 && pset->points[minPoint]->hyper == 0)
                candHyper.push_back(new hyperplane(pset->points[index], pset->points[minPoint], expDim));
            if (pset->points[minPoint]->value == 0) {
                pset->points[minPoint]->value = 1;
                candPoint.push(pset->points[minPoint]->index);
            }
        }

        if(maxPoint < pset->points.size() && maxPoint >= 0)
        {
            if (pset->points[index]->hyper == 0 && pset->points[maxPoint]->hyper == 0)
                candHyper.push_back(new hyperplane(pset->points[index], pset->points[maxPoint], expDim));
            if (pset->points[maxPoint]->value == 0) {
                pset->points[maxPoint]->value = 1;
                candPoint.push(pset->points[maxPoint]->index);
            }
        }

        pset->points[index]->hyper = 1;
        p = new point_t(dimu);
        p->attr[0] = (min + max)/2;
        pset->points[index]->internal = p;
        return true;
    }


}


bool hyperplane_set::find_boundary(point_set *pset, int index)
{
    int M = hyperplanes.size(), size = 0, dim = hyperplanes[0]->dim;
    char file1[MAX_FILENAME_LENG];
    sprintf(file1, "../output/hyperplane_data.txt");
    char file2[MAX_FILENAME_LENG];
    sprintf(file2, "../output/ext_pt.txt");

    //list all the necessary hyperplanes
    hyperplane_set *cHull = new hyperplane_set();
    cHull->expDim = new point_t(expDim);
    for(int i = 0; i < M; ++i)
        cHull->hyperplanes.push_back(new hyperplane(hyperplanes[i]));

    //cHull->print();
    for (int i = 0; i < pset->points.size(); ++i)
    {
        if (i != index)
            cHull->hyperplanes.push_back(new hyperplane(pset->points[i], pset->points[index], expDim));
    }
    point_t *feasible_pt = cHull->find_feasible();//the feasible point

    if (feasible_pt == NULL)
    {
        std::cout << "Intersection is infeasible. \n";
        return false;
    }
    return true;
}

/**
 * @brief   Find a point inside the convex hull
 * @return  The point inside the convex hull
 */
point_t *hyperplane_set::find_feasible() {
    int M = hyperplanes.size();
    int D = hyperplanes[0]->dim;

    // D + 2variables: D for dim, 2 for additional var for feasible
    int* ia = new int[1 + (D + 2) * M];  //TODO: delete
    int* ja = new int[1 + (D + 2) * M];  //TODO: delete
    double* ar = new double[1 + (D + 2) * M];   //TODO: delete
    int i, j;

    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "find_feasible");
    glp_set_obj_dir(lp, GLP_MAX);


    glp_add_rows(lp, M);  // add D rows: q_1...q_D
    // Add rows q_1 ... q_D
    for (i = 1; i <= M; i++) {
        char buf[10];
        sprintf(buf, "q%d", i);
        glp_set_row_name(lp, i, buf);
        glp_set_row_bnds(lp, i, GLP_UP, 0, 0); // qi = 0
    }


    glp_add_cols(lp, D + 2);    // add D columns: v[1] ... v[D]
    // Add col v[1] ... v[D]
    for (i = 1; i <= D + 2; i++) {
        char buf[10];
        sprintf(buf, "v%d", i);

        glp_set_col_name(lp, i, buf);

        if(i <= D)
            glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0); // -infty <= v[i] < infty
        else if (i == D + 1)
            glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0); // 0 <= v[i] < infty
        else
            glp_set_col_bnds(lp, i, GLP_UP, 0.0, D+1);

        if(i == D + 2)
            glp_set_obj_coef(lp, i, 1);  // objective: 0
        else
            glp_set_obj_coef(lp, i, 0.0);  // objective: 0
    }


    int counter = 1;
    // set value on row q1 ... qD
    for (i = 1; i <= M; i++) {
        for (j = 1; j <= D + 2; j++) {

            ia[counter] = i; ja[counter] = j;

            if(j <= D)
            {
                ar[counter++] = hyperplanes[i-1]->norm[j-1];
                //printf("%lf ", hyperplane[i-1]->normal->coord[j-1]);
            }
            else if (j == D+1)
            {
                ar[counter++] = hyperplanes[i-1]->offset;
                //printf("%lf ", hyperplane[i-1]->offset);
            }
            else if (j == D+2)
            {
                ar[counter++] = 1;
                //printf("1.00000\n");
            }
        }
    }

    // loading data
    glp_load_matrix(lp, counter - 1, ia, ja, ar);

    // running simplex
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF; // turn off all message by glp_simplex

    glp_simplex(lp, &parm);


    point_t* feasible_pt = new point_t(D);
    double w1, w2;
    w1 = glp_get_col_prim(lp, D+1);
    w2 = glp_get_col_prim(lp, D+2);
    double status = glp_get_prim_stat(lp);

    if(w1 < EQN3 || w2 < EQN3 || isZero(w1) || isZero(w2))
    {
        //printf("LP feasible error.\n");
        return NULL;
    }
    for (i = 0; i < D; i++)
    {
        double v = glp_get_col_prim(lp, i + 1);
        //printf("w%d = %lf\n", i + 1, v);
        feasible_pt->attr[i] = v / w1;
    }
    feasible_pt->d = D;

    glp_delete_prob(lp); // clean up
    delete[]ia;
    delete[]ja;
    delete[]ar;

    return feasible_pt;
}


/**
 * @brief Check the relation between the hyperplane and the intersection of the hyperplane set
 *        Since the extreme points of the half_set can not be accurate enough, we set "Precision" to solve the error
 * @param h     The hyperplane
 * @return      1: half_set on the positive side of the hyperplane
 *              -1: half_set on the negative side of the hyperplane
 *              0: half_set intersects with the hyperplane
 */
int hyperplane_set::check_relation(hyperplane *h)
{
    int M = ext_pts.size();
    if (M < 1)
    {
        printf("%s\n", "None of the ext_pts in the set.");
        return -2;
    }
    int positive = 0, negative = 0;
    for (int i = 0; i < M; i++)
    {
        int relation = h->check_position(ext_pts[i]);
        if (relation == -1)
            ++negative;
        else if (relation == 1)
            ++positive;
        if (positive > 0 && negative > 0)
        {
            return 0;
        }
    }
    if (negative > 0)
        return -1;
    else
        return 1;
}


/**
 * @brief Check whether p1 R-dominates p2
 * @param p1 The first point
 * @param p2 The second point
 * @return   1 R-dominates
 *          -1 Does not R-dominate
 */
bool hyperplane_set::R_dominate(point_t *p1, point_t *p2)
{
    int size = ext_pts.size(), d = ext_pts[0]->d;
    hyperplane *h = new hyperplane(p1, p2, expDim);
    for (int i = 0; i < size; i++)
    {
        if (h->check_positive(ext_pts[i]) != 1)
        {
            delete h;
            return false;
        }
    }
    delete h;
    return true;
}

/**
 * @brief Prune all the points which are R-dominated by at least one point
 *        There is a point which the closer to any expected point in R
 * @param The point set
 */
void hyperplane_set::nearestSkyline(point_set *pset)
{
    int i = 0;
    while (i < pset->points.size())
    {
        //cout << "ID: " <<pset->points[i]->id << "\n";
        //pset->points[i]->print();
        for(int j = 0; j < pset->points.size(); ++j)
        {
            //pset->points[j]->print();
            if(i != j && R_dominate(pset->points[j], pset->points[i]))
            {
                pset->points.erase(pset->points.begin() + i);
                --i;
                break;
            }
        }
        ++i;
    }
    //record the index of the point in the list
    for(int i = 0; i < pset->points.size(); ++i)
        pset->points[i]->index = i;
}


/**
 * @brief Calculate the average point of the extreme points
 * @param ap The average point
 */
point_t* hyperplane_set::average_point()
{
    int size = ext_pts.size(), dim = ext_pts[0]->d;
    point_t* avgPoint = new point_t(dim);
    for(int j = 0; j < dim; ++j)
        avgPoint->attr[j] = 0;
    //calculate the sum
    for(int i = 0; i < size; ++i)
    {
        for (int j = 0; j < dim; ++j)
            avgPoint->attr[j] += ext_pts[i]->attr[j];
    }
    for(int j = 0; j < dim; ++j)
        avgPoint->attr[j] /= size;
    return avgPoint;
}


/**
 * @brief Find the nearest point of any expected point in R
 * @param p_set     The point set
 * @return          The nearest point
 *                  If return NULL, there does not exist the nearest point
 */
point_t* hyperplane_set::findNearest(point_set *p_set)
{
    int M = ext_pts.size(), dim = p_set->points[0]->d, dimo = p_set->points[0]->d_order, dimu = p_set->points[0]->d_unorder;
    point_t *p = new point_t(dimo, dimu);
    for(int i = 0; i < dimo; ++i)
        p->attr[i] = expDim->attr[i];

    point_t *pp = average_point();
    for(int i = dimo; i < dim; ++i)
        p->attr[i] = pp->attr[i-dimo];
    int nearestPoint = p_set->findClosest(p);
    //p_set->points[nearestPoint]->print();

    for(int i = 0; i < p_set->points.size(); ++i)
    {
        if(nearestPoint != i)
        {
            hyperplane *h = new hyperplane(p_set->points[nearestPoint], p_set->points[i], expDim);
            for(int j = 0; j < M; ++j)
            {
                if(h->check_positive(ext_pts[j]) != 1)
                    return NULL;
            }
        }
    }
    return p_set->points[nearestPoint];
}


bool hyperplane_set::is_prune(point_set *pset, int index, std::vector<hyperplane_set*> &partitionSet)
{
    int M = hyperplanes.size(), size = 0, dim = hyperplanes[0]->dim;
    //Initial a test
    hyperplane_set *cHull = new hyperplane_set();
    cHull->expDim = new point_t(expDim);
    for(int i = 0; i < M; ++i)
        cHull->hyperplanes.push_back(new hyperplane(hyperplanes[i]));

    //cHull->print();
    hyperplane_set *hs = partitionSet[pset->points[index]->boundaryHyperID];
    for (int i = 0; i < hs->hyperplanes.size(); ++i)
    {
        cHull->hyperplanes.push_back(new hyperplane(hs->hyperplanes[i]));
        //cHull->hyperplanes.push_back(new hyperplane(hs->hyperplanes[i]->p_1, hs->hyperplanes[i]->p_2, expDim, 0.05));
    }
    point_t *feasible_pt = cHull->find_feasible();//the feasible point
    pset->points[index]->internal = feasible_pt;
    if(feasible_pt == NULL)
        return true;
    else
        return false;
}

/**
 * @brief Find two extreme points in the expected space as vertices of the vector
 * @param pset      The point set
 * @param index1    The index of the first extreme points
 * @param index2
 */
void hyperplane_set::findExpforVec(point_set *pset, int &index1, int &index2)
{
    int dimo = pset->points[0]->d_order, dimu = pset->points[0]->d_unorder;
    int *extBest = new int[ext_pts.size()];
    point_t *Expt = new point_t(dimo, dimu);
    for(int j = 0; j < dimo; ++j)
        Expt->attr[j] = expDim->attr[j];
    for(int i = 0; i < ext_pts.size(); ++i)
    {
        for(int j = dimo; j < dimo + dimu; ++j)
            Expt->attr[j] = ext_pts[i]->attr[j - dimo];
        extBest[i] = pset->findClosest(Expt);
    }

    double dist = -1;
    for(int i = 0; i < ext_pts.size() - 1; ++i)
    {
        for(int j = i + 1; j < ext_pts.size(); ++j)
        {
            if(extBest[i] != extBest[j])
            {
                double ds = ext_pts[i]->distance(ext_pts[j]);
                if(ds > dist)
                {
                    dist = ds;
                    index1 = i;
                    index2 = j;
                }
            }
        }
    }
}

/**
 * @brief Find the upper and lower bound of the indexDimension-th dimension among the extreme points
 * @param indexDimenion     The index of the dimension
 * @param min               The maximum value
 * @param max               The minimum value
 */
void hyperplane_set::findMinMax(int indexDimenion, double &min, double &max)
{
    int size = ext_pts.size();
    min = INF; max = -INF;
    for(int i = 0; i < ext_pts.size(); ++i)
    {
        if(min > ext_pts[i]->attr[indexDimenion])
            min = ext_pts[i]->attr[indexDimenion];
        if(max < ext_pts[i]->attr[indexDimenion])
            max = ext_pts[i]->attr[indexDimenion];
    }
}


double hyperplane_set::findL1Dis(int dim)
{
    double L1 = 0, min, max;
    for(int i = 0; i < dim; ++i)
    {
        findMinMax(i, min, max);
        L1 = L1 + max - min;
    }
    return L1;
}





/**
 * @brief draw a line crossing the average points and calculate the intersecting points of the line and R
 * @param leftPts   One intersecting point
 * @param rightPts  One intersecting point
 * @param dimInex
 */
void hyperplane_set::findendPts(point_t *lowerPts, point_t *upperPts, int dimInex)
{
    int dim = hyperplanes[0]->dim;
    point_t* avgPts = average_point();
    for(int i = 0; i < dim; ++i)
    {
        lowerPts->attr[i] = avgPts->attr[i];
        upperPts->attr[i] = avgPts->attr[i];
    }


    lowerPts->attr[dimInex] = -INF;
    upperPts->attr[dimInex] = INF;

    for(int i = 0; i < hyperplanes.size(); ++i)
    {
        double sum = hyperplanes[i]->offset;
        for(int j = 0; j < dim; ++j)
        {
            if(j != dimInex)
                sum += hyperplanes[i]->norm[j] * avgPts->attr[j];
        }
        if(hyperplanes[i]->norm[dimInex] > 0)
        {
            sum = sum / (-hyperplanes[i]->norm[dimInex]);
            if(upperPts->attr[dimInex] > sum)
                upperPts->attr[dimInex] = sum;
        }
        else if(hyperplanes[i]->norm[dimInex] < 0)
        {
            sum = sum / (-hyperplanes[i]->norm[dimInex]);
            if(lowerPts->attr[dimInex] < sum)
                lowerPts->attr[dimInex] = sum;
        }
    }
}

/**
 * @brief Find a extreme point which makes the point (index) is the cloest point
 * @param e         The extreme points
 * @param CPoint    The candidate set
 * @param index     The point
 */
point_t* hyperplane_set::findExePt(point_set *CPoint, int index)
{
    int dimo = expDim->d_order, dimu = expDim->d_unorder;
    point_t *e = new point_t(dimo, dimu);
    for(int i = 0; i < dimo; ++i)
        e->attr[i] = expDim->attr[i];

    for(int i = 0; i < ext_pts.size(); ++i)
    {
        for(int j = dimo; j < dimu + dimo; ++j)
            e->attr[j] = ext_pts[i]->attr[j - dimo];
        if (index == CPoint->findClosest(e))
            return ext_pts[i];
    }
}























