#ifndef ACTIVERANKING_H
#define ACTIVERANKING_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

int ActiveRanking(point_set* pset, point_t* e, long mem_baseline);

#endif //RUN_ACTIVERANKING_H
