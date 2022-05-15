#ifndef RH_H
#define RH_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"


void RH(point_set *pset, point_set *realSet, point_t* e, int &TID, std::ofstream &fp);

#endif //RUN_RH_H
