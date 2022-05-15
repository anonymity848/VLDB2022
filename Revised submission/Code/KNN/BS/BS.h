#ifndef BS_H
#define BS_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

void BS(point_set *pset, point_t* e, double Beta, long mem_baseline, int type);

void BSTopk(point_set *origset, point_t* e, double Beta, int k, long mem_baseline, int type);

#endif //RUN_BS_H
