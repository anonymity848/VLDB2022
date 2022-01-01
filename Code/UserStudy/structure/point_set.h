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
    void possibleNearest(double *max, double *min, point_set *top, point_t *e, int level);
    point_t* initializeExpectedPoint();
    void transformAttr(point_t *start, point_t *vec, point_t *vecNorm);
    void write(std::string fileName);
    void prunePt(point_t* p);
    int show_to_user(int idx_1, int idx_2);
    int show_to_user_unanswer(int idx_1, int idx_2);
    void realPrint(std::string s, int numOfQuestion, std::vector<int> pIndex, std::ofstream &fp);

};


#endif //U_2_POINT_SET_H
