#include "structure/data_utility.h"
#include "structure/data_struct.h"
#include "structure/point_set.h"
#include "structure/define.h"
#include <vector>
#include <ctime>
#include <sys/time.h>
#include "UH/UHRandom.h"
#include "GroundTruth/Groundtruth.h"
#include "BS/BS.h"
#include "DI/DI.h"
#include "ActiveRanking/ActiveRanking.h"
#include "RH/RH.h"
#include "EDI/EDI.h"
#include "UtilityApprox/UtilityApprox.h"
#include "Baseline/Baseline.h"

int main(int argc, char *argv[])
{
    long mem_baseline = get_mem_usage();
    ifstream config("../config.txt");
    char data_name[MAX_FILENAME_LENG]; double Beta, gamma; int type, k;
    config >>  data_name >> type >> Beta >> gamma >> k;
    cout  << data_name << "  Type:" << type << "   Beta:" << Beta << "   gamma:" << gamma << "   k:" << k << "\n";
    //initialization
    //point set
    point_set *p_set = new point_set(data_name);
    point_t *e = new point_t(p_set->points[0]->d_order, p_set->points[0]->d_unorder);
    for(int i = 0; i < p_set->points[0]->d; ++i)
        config >> e->attr[i];
    //e->print();

    std::cout << "--------------------------------------------------------------------------------------\n";
    printf("|%20s |%15s |%15s |%15s |%10s |\n", "Algorithm", "# of Questions", "Preprocessing",
           "Interaction", "Point #ID");
    std::cout << "--------------------------------------------------------------------------------------\n";
    ground_truth(p_set, e, k); //look for the ground truth maximum utility point

    point_set *pps = new point_set(p_set);

    if(k == 1)
    {
        baseline(pps, e, mem_baseline, type);

        pps = new point_set(p_set);
        DI(pps, e, mem_baseline, type);

        pps = new point_set(p_set);
        EDI(pps, e, Beta, gamma, mem_baseline, type);

        pps = new point_set(p_set);
        BS(pps, e, Beta, mem_baseline, type);

        pps = new point_set(p_set);
        RH(pps, e, mem_baseline, type);

        pps = new point_set(p_set);
        ActiveRanking(pps, e, mem_baseline, type);

        pps = new point_set(p_set);
        UtilityApprox(pps, e, mem_baseline, type);

        pps = new point_set(p_set);
        UH_Random(pps, e, mem_baseline, type);
    }

    pps = new point_set(p_set);
    BSTopk(pps, e, Beta, k, mem_baseline, type);

    pps = new point_set(p_set);
    EDITopk(pps, e, Beta, gamma, k, mem_baseline, type);

    pps = new point_set(p_set);
    ActiveRankingTopk(pps, e, k, mem_baseline, type);

    //cout << get_mem_usage() - mem_baseline << "\n";
    delete p_set;
    return 0;
}
