#ifndef ACTIVERANKING_H
#define ACTIVERANKING_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

int ActiveRanking(point_set* pset, point_t* e, long mem_baseline, int type);

int ActiveRankingTopk(point_set* pset, point_t* e, int k, long mem_baseline, int type);

#endif //RUN_ACTIVERANKING_H
