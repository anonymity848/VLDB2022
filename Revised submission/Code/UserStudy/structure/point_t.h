#ifndef POINT_T_H
#define POINT_T_H
#include "define.h"

class point_t
{
public:
    int	id, index, d, d_order, d_unorder;
    double *attr, *traAttr;
    double value, hyper;
    point_t *internal;
    std::vector<point_t*> ext;
    int boundaryHyperID;

    explicit point_t(int dim);
    explicit point_t(int d_order, int d_unorder);
    point_t(int d_order, int d_unorder, int id);
    explicit point_t(point_t *p);
    point_t(point_t *p1, point_t *p2);
    ~point_t();
    void print(); //print the point

    bool is_same(point_t *p);//check whether the point is the same as p
    bool is_zero();
    double dot_product(point_t* p);
    double dot_product(double *v);
    point_t* sub(point_t* p);
    point_t* add(point_t* p);
    point_t* scale(double c);
    double distance(point_t* p);
    double distance(double *p);
    double bound(point_t* p, double y);
    double tranBound(point_t *p);
    void printResult(char *name, int Qcount, timeval t1);
    void printResult(std::ofstream &out_cp, char *name, int Qcount, timeval t1, double preTime, long mem_baseline);
    void transformAttr(point_t *start, point_t *vec, point_t *vecNorm);
    point_t* generateEXP(point_t *p2);
    void final_result(std::string s, int numOfQuestion, std::ofstream &fp);
};










#endif //U_2_POINT_T_H
