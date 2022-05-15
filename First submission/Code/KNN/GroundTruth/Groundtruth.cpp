#include "Groundtruth.h"

/**
 * @brief Find the nearest point of e
 * @param p_skyline The point set
 * @param e         The expected point
 */
void ground_truth(point_set *pset, point_t *e)
{
    point_t *p = pset->points[pset->findClosest(e)];
    std::cout << "-----------------------------------------------------------------------------------\n";
    printf("|%15s |%15s |%15s |%15s |%10d |\n", "Ground Truth", "-", "-", "-", p->id);
    std::cout << "-----------------------------------------------------------------------------------\n";
}
