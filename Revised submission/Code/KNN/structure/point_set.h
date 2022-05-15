#ifndef POINT_SET_H
#define POINT_SET_H

#include "point_t.h"
#include <string>


class point_set
{
public:
    std::vector<point_t*> points;

    point_set();
    explicit point_set(point_set *p_set);
    explicit point_set(const char* input);
    ~point_set();

    void print();
    void random(double RandRate);
    point_set* sort(point_t *u);
    int findClosest(point_t *e);
    point_set* findNearestK(point_t* e, int k);
    void possibleNearest(double *max, double *min, point_set *top, point_t *e, int level);
    point_t* initializeExpectedPoint();
    void transformAttr(point_t *start, point_t *vec, point_t *vecNorm);
    void write(std::string fileName);
    void prunePt(point_t* p);
    void printResult(std::ofstream &out_cp, char *name, int Qcount, timeval t1, double preTime, long mem_baseline,
                                double Pcize, int type);

};


#endif //U_2_POINT_SET_H
