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
    double Beta = 4;
    //point set
    point_set *pSet = new point_set("car.txt");
    point_set *RealSet = new point_set("carReal.txt");
    point_t *e = new point_t(pSet->points[0]->d_order, pSet->points[0]->d_unorder);
    for(int i = 0; i < pSet->points[0]->d_order; ++i)
        e->attr[i] = 1000;
    for(int i = 0; i < pSet->points[0]->d_unorder; ++i)
        e->attr[i + pSet->points[0]->d_order] = 500;

    //look for the ground truth maximum utility point
    //ground_truth(pSet, e);


    string userName;
    // the welcome message
    cout <<"Please input your name: "; cin >> userName;
    cout << "-------------------------Welcome to the recommending car system------------------------------\n";
    cout <<"1. In our research project, we want to ask as few questions as possible so that we could help\n"
           "you to find your favorite used car in our car database. \n"
           "2. The final car returned by our system should be your favorite car.\n"
           "3. There are 8 rounds in the system. We will ask you a list of consecutive questions for each\n"
           "round. Each round involving consecutive questions corresponds to a method in our system.\n"
           "4. You will be presented two cars each time and you need to choose the one you favor more.\n"
           "5. (1). Price(¥) 27k-5215k                  "
           "   (2). Year(Year of production) 1997-2020  "
           "   (3). Used KM 0-1000,000\n"
           "   (4). Length(mm) 2499-6297                "
           "   (5). Width(mm) 1275-2360                 "
           "   (6). Seats 1-9                     ";

    point_set *pps;
    int TID; vector<int> TID_list;
    ofstream fp("../Result/" + userName + ".txt");


    //Active-Ranking Algorithm
    cout << "\n\n\n\n===============================================Round 1======================================================\n";
    pps = new point_set(pSet);
    ActiveRanking(pps, RealSet, e, TID, fp);
    insert_intlist(TID_list, TID);

    //Algorithm EDI
    cout << "\n\n===============================================Round 2======================================================\n";
    pps = new point_set(pSet);
    EDI(pps, RealSet, e, Beta, TID, fp);
    insert_intlist(TID_list, TID);

    //Algorithm RH
    cout << "\n\n===============================================Round 3======================================================\n";
    pps = new point_set(pSet);
    RH(pps, RealSet, e, TID, fp);
    insert_intlist(TID_list, TID);

    //Algorithm BS
    cout << "\n\n===============================================Round 4======================================================\n";
    pps = new point_set(pSet);
    BS(pps, RealSet, e, Beta, TID, fp);
    insert_intlist(TID_list, TID);

    //Algorithm UH-Random
    cout << "\n\n===============================================Round 5======================================================\n";
    pps = new point_set(pSet);
    UH_Random(pps, RealSet, e, TID, fp);
    insert_intlist(TID_list, TID);


    cout<< "\n\n\nThe recommended tuples: \n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::setw(10) << " " << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
              << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
              << std::setw(16) << "Seats\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < TID_list.size(); ++i)
    {
        std::cout << std::setw(10) << i + 1;
        std::cout << std::setw(14) << RealSet->points[TID_list[i]]->attr[0]/1000 << "k";
        for (int j = 1; j < pSet->points[0]->d_order + pSet->points[0]->d_unorder; ++j)
            std::cout << std::setw(15) << RealSet->points[TID_list[i]]->attr[j];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";
    }
    cout <<"Please give an order for the shown used car(s) (e.g., ";
    for(int i = 0; i < TID_list.size() - 1; ++i)
        cout << i + 1 << " ";
    cout << TID_list.size();
    cout <<"), ";
    cout << "where the first one \nis the most preferred tuple and the last one is the least preferred tuple: ";
    int *order = new int[TID_list.size()];
    for(int i = 0; i < TID_list.size(); ++i)
        cin >> order[i];
    fp << "order: \n";
    fp << "------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < TID_list.size(); ++i)
    {
        fp << std::setw(10) << i + 1;
        fp << std::setw(14) << RealSet->points[TID_list[order[i] - 1]]->attr[0]/1000 << "k";
        for(int j = 1; j < pSet->points[0]->d_order + pSet->points[0]->d_unorder; ++j)
            fp << std::setw(15) << RealSet->points[TID_list[order[i] - 1]]->attr[j];
        fp << "\n------------------------------------------------------------------------------------------------------------\n";
    }
    fp << "\n\n\n";



    //unanswered questions

    //Algorithm EDI
    cout << "\n\n===============================================Round 6======================================================\n";
    std::vector<int> TIDset1;
    pps = new point_set(pSet);
    EDI_unasnwer(pps, RealSet, e, Beta, TIDset1, fp);

    //Algorithm BS
    cout << "\n\n===============================================Round 7======================================================\n";
    std::vector<int> TIDset2;
    pps = new point_set(pSet);
    BS_unanswer(pps, RealSet, e, Beta, TIDset2, fp);


    //Algorithm UH-Random
    cout << "\n\n===============================================Round 8======================================================\n";
    std::vector<int> TIDset3;
    pps = new point_set(pSet);
    UH_Random_unanswer(pps, RealSet, e, TIDset3, fp);

    cout<< "\n\n\nThe recommended tuple sets: \n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::setw(10) << "Set 1" << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
              << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
              << std::setw(16) << "Seats\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < TIDset1.size(); ++i)
    {
        std::cout << std::setw(10) << i + 1;
        std::cout << std::setw(14) << RealSet->points[TIDset1[i]]->attr[0]/1000 << "k";
        for (int j = 1; j < pSet->points[0]->d_order + pSet->points[0]->d_unorder; ++j)
            std::cout << std::setw(15) << RealSet->points[TIDset1[i]]->attr[j];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";
    }
    std::cout << "\n\n------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::setw(10) << "Set 2" << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
              << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
              << std::setw(16) << "Seats\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < TIDset2.size(); ++i)
    {
        std::cout << std::setw(10) << i + 1;
        std::cout << std::setw(14) << RealSet->points[TIDset2[i]]->attr[0]/1000 << "k";
        for (int j = 1; j < pSet->points[0]->d_order + pSet->points[0]->d_unorder; ++j)
            std::cout << std::setw(15) << RealSet->points[TIDset2[i]]->attr[j];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";
    }
    std::cout << "\n\n------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::setw(10) << "Set 3" << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
              << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
              << std::setw(16) << "Seats\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    for(int i = 0; i < TIDset3.size(); ++i)
    {
        std::cout << std::setw(10) << i + 1;
        std::cout << std::setw(14) << RealSet->points[TIDset3[i]]->attr[0]/1000 << "k";
        for (int j = 1; j < pSet->points[0]->d_order + pSet->points[0]->d_unorder; ++j)
            std::cout << std::setw(15) << RealSet->points[TIDset3[i]]->attr[j];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";
    }

    cout <<"Please give an order for the shown sets (e.g., ";
    for(int i = 0; i < 2; ++i)
        cout << i + 1 << " ";
    cout <<"3), ";
    cout << "where the first one \nis the most preferred set and the last one is the least preferred set: ";
    int order1[3];
    for(int i = 0; i < 3; ++i)
        cin >> order1[i];
    fp << "order: ";
    for(int i = 0; i < 3; ++i)
    {
        fp << order1[i] << "  ";
    }
    fp << "\n\n";





    fp.close();
    return 0;
}
