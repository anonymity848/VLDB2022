#ifndef MAXUTILITY_H
#define MAXUTILITY_H

#include "structure/define.h"
#include "structure/data_struct.h"
#include "structure/data_utility.h"

#include <vector>
#include <algorithm>
#include "structure/rtree.h"
#include "Others/lp.h"
#include "Others/pruning.h"
#include "Others/operation.h"
#include <queue>

using namespace std;

void UH_Random(point_set* pset, point_set* realSet, point_t* e, int &TID, std::ofstream &fp);

void UH_Random_unanswer(point_set* pset, point_set* realSet, point_t* e, std::vector<int> &TID, std::ofstream &fp);

#endif