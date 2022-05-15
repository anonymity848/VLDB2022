#ifndef EDI_H
#define EDI_H
#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "DI/DI.h"
#include "Others/operation.h"


void EDI(point_set *pset, point_t *e, long mem_baseline);

void EDI2(point_set *pset, point_t *e, long mem_baseline);

void EDI3(point_set *pset, point_t *e, long mem_baseline);

void EDI4(point_set *pset, point_t *e, long mem_baseline);

void EDI(point_set *pset, point_t *e, double Beta, long mem_baseline);

#endif //RUN_EDI_H
