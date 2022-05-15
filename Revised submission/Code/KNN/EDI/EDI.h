#ifndef EDI_H
#define EDI_H
#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "DI/DI.h"
#include "Others/operation.h"


void EDI1(point_set *pset, point_t *e, double Beta, double gamma, long mem_baseline, int type);
void EDI(point_set *pset, point_t *e, double Beta, double gamma, long mem_baseline, int type);

void EDITopk(point_set *origset, point_t *e, double Beta, double gamma, int k, long mem_baseline, int type);

#endif //RUN_EDI_H
