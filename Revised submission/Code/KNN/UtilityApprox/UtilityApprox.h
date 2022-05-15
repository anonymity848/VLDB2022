#ifndef RUN_UTILITYAPPROX_H
#define RUN_UTILITYAPPROX_H

#include "data_utility.h"
#include "point_set.h"
#include "hyperplane_set.h"
#include "Others/operation.h"

/**
 * @brief Danupon's Fake-Points algorithm UtilityApprox
 * @param P the point dataset
 * @param u the expected point
 */
void UtilityApprox(point_set *pset, point_t *e, long mem_baseline, int type);


#endif //RUN_UTILITYAPPROX_H
