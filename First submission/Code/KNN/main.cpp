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


int main(int argc, char *argv[])
{
    long mem_baseline = get_mem_usage();
    ifstream config("../config.txt");
    char data_name[MAX_FILENAME_LENG]; double Beta;
    config >> data_name >> Beta;
    //initialization
    //point set
    point_set *p_set = new point_set(data_name);
    point_t *e = new point_t(p_set->points[0]->d_order, p_set->points[0]->d_unorder);
    for(int i = 0; i < p_set->points[0]->d; ++i)
        config >> e->attr[i];

    printf("-----------------------------------------------------------------------------------\n");
    printf("|%15s |%15s |%15s |%15s |%10s |\n", "Algorithm", "# of Questions", "Preprocessing", "Interaction", "Point #ID");
    printf("-----------------------------------------------------------------------------------\n");
    ground_truth(p_set, e); //look for the ground truth maximum utility point

    point_set *pps = new point_set(p_set);
    DI(pps, e, mem_baseline);
    
    pps = new point_set(p_set);
    EDI(pps, e, Beta, mem_baseline);
    
    pps = new point_set(p_set);
    BS(pps, e, Beta, mem_baseline);
    
    pps = new point_set(p_set);
    RH(pps, e, mem_baseline);
    
    pps = new point_set(p_set);
    ActiveRanking(pps, e, mem_baseline);
    
    pps = new point_set(p_set);
    UtilityApprox(pps, e, mem_baseline);
    
    pps = new point_set(p_set);
    UH_Random(pps, e, mem_baseline);

    delete p_set;
    return 0;
}
