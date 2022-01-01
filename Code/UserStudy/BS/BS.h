#ifndef BS_H
#define BS_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

void BS(point_set* pset, point_set* realSet, point_t* e, double Beta, int &TID, std::ofstream &fp);


void BS_unanswer(point_set* pset, point_set* realSet, point_t* e, double Beta, std::vector<int> &TID, std::ofstream &fp);
#endif //RUN_BS_H
