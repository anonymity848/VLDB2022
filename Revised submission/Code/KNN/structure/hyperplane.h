#ifndef HYPERPLANE_H
#define HYPERPLANE_H

#include "point_t.h"
#include "point_set.h"

class hyperplane
{
public:
    int dim;          //dimension
    double* norm;   //parameter of the norm vector
    double offset;
    point_t* p_1;
    point_t* p_2;

    hyperplane();
    explicit hyperplane(int dim);
    explicit hyperplane(hyperplane* h);
    hyperplane(int dim, double *n, double offset);
    hyperplane(point_t *p1, point_t *p2);
    hyperplane(point_t* p1, point_t* p2, point_t* exp);
    hyperplane(point_t *p1, point_t *p2, point_t* exp, double error);
    bool is_same(hyperplane *h);
    hyperplane(int d, double* p1, double* p2);
    int check_position(point_t *p);
    int check_positive(point_t *p);
    double check_distance(point_t *p);
    ~hyperplane();
    void print();
    int priority(point_set *pset);
    int priority2(point_set *pset);
    int priority(point_set *pset, double dNN, double Beta);
    bool is_boundary(point_set *pset);
    bool necessaryPrune(point_t *p);
    double distance(point_t *p);
    void transform(point_t* p);
    double bound(point_t *p1, point_t *p2);
};


#endif //U_2_HYPERPLANE_H
