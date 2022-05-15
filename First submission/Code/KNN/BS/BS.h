#ifndef BS_H
#define BS_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

void BS(point_set *pset, point_t* e, double Beta, long mem_baseline);

void BS_sampling(point_set *pset, point_t* e);


#endif //RUN_BS_H
