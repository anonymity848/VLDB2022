#ifndef EDI_H
#define EDI_H
#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "DI/DI.h"
#include "Others/operation.h"


void EDI(point_set* pset, point_set* realSet, point_t* e, double Beta, int &TID, std::ofstream &fp);


void EDI_unasnwer(point_set* pset, point_set* realSet, point_t* e, double Beta, std::vector<int> &TID, std::ofstream &fp);
#endif //RUN_EDI_H
