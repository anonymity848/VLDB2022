#include "point_set.h"
#include "../qhull/qhull_build.h"


/**
 * @brief Constructor
 */
point_set::point_set(){}

/**
 * @brief Constructor
 *        Create a point set the same as p_set, all the points are re-created
 * @param p_set     Point set
 */
point_set::point_set(point_set *p_set)
{
    int dimo = p_set->points[0]->d_order, dimu = p_set->points[0]->d_unorder,
        num = p_set->points.size(), dim = dimo + dimu;
    point_t *p;
    for (int i = 0; i < num; i++)
    {
        p = new point_t(dimo, dimu, p_set->points[i]->id);
        for (int j = 0; j < dim; j++)
        {
            p->attr[j] = p_set->points[i]->attr[j];
        }
        p->index = p_set->points[i]->index;
        p->value = p_set->points[i]->value;
        for(int a = 0; a < p_set->points[i]->ext.size(); ++a)
        {
            p->ext.push_back(p_set->points[i]->ext[a]);
        }
        points.push_back(p);
    }

}

/**
 * @brief Constructor
 *        Record all the points in the input file to the point set
 * @param input     Name of the data file.
 */
point_set::point_set(const char* input)
{
    FILE* c_fp;
    char filename[MAX_FILENAME_LENG];
    sprintf(filename, "../input/%s", input);
    //printf("%s\n", filename);
    if ((c_fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Cannot open the data file %s.\n", filename);
        exit(0);
    }

    int num, dimo, dimu;
    point_t* p;
    fscanf(c_fp, "%i%i%i", &num, &dimo, &dimu);

    // read points line by line
    for (int i = 0; i < num; i++)
    {
        p = new point_t(dimo, dimu, i);
        for (int j = 0; j < dimo + dimu; j++)
            fscanf(c_fp, "%lf", &p->attr[j]);
        points.push_back(p);
    }
    fclose(c_fp);
}

/**
 *@brief  Destructor
 *        Delete the points in the array
 */
point_set::~point_set()
{
    int i = points.size();
    point_t *p;
    while(i>0)
    {
        p = points[i-1];
        points.pop_back();
        delete p;
        i--;
    }
    points.clear();
}

/*
 *	For debug purpose, print all the points in the set
 */
void point_set::print()
{
    for (int i = 0; i < points.size(); i++)
        points[i]->print();
    printf("\n");
}

/**
 * @brief           Reload the points Randomly
 *                  Define "RandRate" to control how many point reinserted
 * @param p_set     The returned dataset where the points are inserted randomly
 */
void point_set::random(double RandRate)
{
    int size = points.size();
    //reinsert
    for (int i = 0; i < size * RandRate; i++)
    {
        int n = ((int) rand()) % size;
        point_t *p = points[n];
        points.erase(points.begin() + n);
        points.push_back(p);
    }
}

/**
 * @brief       Sort points based on their distance to e
 * @param e     The expected point
 * @return      The point set which contains all the point in order
 */
point_set* point_set::sort(point_t *e)
{
    int size = points.size();
    if(size <=0)
        return NULL;
    point_set *return_set = new point_set();
    return_set->points.push_back(points[0]);
    for (int i = 1; i < size; i++)
    {
        double v0 = points[i]->distance(e);
        int left = 0, right = return_set->points.size() - 1;
        //find the place for p_set[i] in return_point and record the place index in "left"
        while (left <= right)
        {
            int middle = (left + right) / 2;
            double v = return_set->points[middle]->distance(e);
            if (v0 < v)
            {
                right = middle - 1;
            }
            else
            {
                left = middle + 1;
            }
        }
        return_set->points.insert(return_set->points.begin() + left, points[i]);
    }
    /*
    for(int i=0; i<return_set->points.size();i++)
    {
        return_set->points[i]->print();
    }
    */
    return return_set;
}

/**
 * @brief Find the nearest point of e in the dataset
 * @param e     The expected point
 * @return      The nearest point of e
 */
int point_set::findClosest(point_t *e) {

    int M = points.size();
    double mindist = INF;
    int position = -1;
    for(int i = 0; i < M; ++i)
    {
        double dist = points[i]->distance(e);
        if(dist < mindist)
        {
            mindist = dist;
            position = i;
        }
    }
    return position;
}

/**
 * @brief Initialize the user's expected point
 * @return  The expected point
 */
point_t *point_set::initializeExpectedPoint()
{
    int dimo = points[0]->d_order, dimu = points[0]->d_unorder, dim = dimo + dimu;
    point_t* expectedPoint = new point_t(dimo, dimu);
    double *max = new double[dim];
    double *min = new double[dim];
    for(int i = 0; i < dim; ++i)
    {
        max[i] = 0; min[i] = INF;
    }
    for(int i = 0; i < points.size(); ++i)
    {
        for(int j = 0; j < dim; ++j)
        {
            if(points[i]->attr[j] < min[j])
                min[j] = points[i]->attr[j];
            if(points[i]->attr[j] > max[j])
                max[j] = points[i]->attr[j];
        }
    }

    for(int i = 0; i < dimo; ++i)
        expectedPoint->attr[i] = max[i];
    for(int i = dimo; i < dim; ++i)
    {
        double difference = (max[i] - min[i])/100;
        expectedPoint->attr[i] = min[i] + (rand() % 100) * difference;
    }
    return expectedPoint;
}


/**
 * @brief Find all the possible nearest points
 * @param max   The upper bound of each dimension
 * @param min   The lower bound of each dimension
 * @param top   The point set which contains all the possible nearest points
 * @param level Record the dimension searched (initial level = dimo)
 */
void point_set::possibleNearest(double *max, double *min, point_set *top, point_t *e, int level)
{
    int dim = points[0]->d, M = points.size();
    double segment = 10;
    if (level >= dim - 1)
    {
        for (int i = 0; i <= segment; ++i)
        {
            e->attr[level] = min[level] + i * (max[level] - min[level]) / segment;
            //Find the point with max utility w.r.t u
            e->print();
            int maxId = findClosest(e);
            if (points[maxId]->value == 0)
            {
                points[maxId]->value = 1;
                top->points.push_back(points[maxId]);
            }
        }
    }
    else
    {
        for (int i = 0; i <= segment; ++i)
        {
            e->attr[level] = min[level] + i * (max[level] - min[level]) / segment;
            possibleNearest(max, min, top, e, level + 1);
        }
    }
}

/**
 * @brief Transform the coordinate of the points to parallel and vertical to vec
 * @param vec The direction vector
 */
void point_set::transformAttr(point_t *start, point_t *vec, point_t *vecNorm)
{
    int M = points.size();
    for(int i = 0; i < M; ++i)
    {
        points[i]->transformAttr(start, vec, vecNorm);
    }
}


/**
 * @brief   Write the dataset to the txt file
 * @param   fileName  The name of the txt file
 */
void point_set::write(std::string fileName)
{
    ofstream wPtr;
    wPtr.open(fileName, std::ios::out);
    wPtr.setf(ios::fixed, ios::floatfield);  // set as fixed model
    wPtr.precision(6);  // set precision to 6


    int size = points.size(), dimo = points[0]->d_order, dimu = points[0]->d_unorder;
    // write the points
    wPtr << size << "   " << dimo << "   " << dimu << " \n";//record the offset as one dimension
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < dimo + dimu; j++)
        {
            wPtr << points[i]->attr[j] <<" ";
        }
        wPtr <<"\n";
    }
    wPtr.close();
}

/**
 * @param Delete point p in the point set
 * @param p The point
 */
void point_set::prunePt(point_t *p)
{
    for(int i = 0; i < points.size(); ++i)
    {
        if(points[i]->id == p->id)
        {
            points.erase(points.begin() + i);
            return;
        }
    }
}



/**
 * Show a user two points and ask which one is more preferred by the user
 * @param idx_1 The index of the first point
 * @param idx_2 The index of the second point
 * @return
 */
int point_set::show_to_user(int idx_1, int idx_2)
{

    int option = 0;
    // ask the user for the better car among two given options
    while (option != 1 && option != 2)
    {
        std::cout << "Please choose the car you favor more:\n";
        std::cout << "------------------------------------------------------------------------------------------------------------\n";
        std::cout << std::setw(10) << " " << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
                  << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
                  << std::setw(16) << "Seats\n";
        std::cout << "------------------------------------------------------------------------------------------------------------\n";
        std::cout << std::setw(10) << "Option 1";
        cout << std::setw(14) << points[idx_1]->attr[0]/1000 << "k";
        for(int i = 1; i < points[idx_1]->d; ++i)
            cout << std::setw(15) << points[idx_1]->attr[i];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";

        std::cout << std::setw(10) << "Option 2";
        cout << std::setw(14) << points[idx_2]->attr[0]/1000 << "k";
        for(int i = 1; i < points[idx_2]->d; ++i)
            cout << std::setw(15) << points[idx_2]->attr[i];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";

        printf("Your choice: ");
        std::cin >> option;
    }
    return option;
}


/**
 * Show a user two points and ask which one is more preferred by the user
 * @param idx_1 The index of the first point
 * @param idx_2 The index of the second point
 * @return
 */
int point_set::show_to_user_unanswer(int idx_1, int idx_2)
{

    int option = -1;
    // ask the user for the better car among two given options
    while (option != 1 && option != 2 && option != 3)
    {
        std::cout << "Please choose the car you favor more:\n";
        std::cout << "------------------------------------------------------------------------------------------------------------\n";
        std::cout << std::setw(10) << " " << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
                  << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
                  << std::setw(16) << "Seats\n";
        std::cout << "------------------------------------------------------------------------------------------------------------\n";
        std::cout << std::setw(10) << "Option 1";
        cout << std::setw(14) << points[idx_1]->attr[0]/1000 << "k";
        for(int i = 1; i < points[idx_1]->d; ++i)
            cout << std::setw(15) << points[idx_1]->attr[i];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";

        std::cout << std::setw(10) << "Option 2";
        cout << std::setw(14) << points[idx_2]->attr[0]/1000 << "k";
        for(int i = 1; i < points[idx_2]->d; ++i)
            cout << std::setw(15) << points[idx_2]->attr[i];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";

        std::cout << std::setw(10) << "Option 3";
        cout << "           Hard to decide which one is more preferred.";
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";

        printf("Your choice: ");
        std::cin >> option;
    }
    return option;
}





/**
 * @brief   The user evaluation step
 * @return  The boredness given by user
 */
void point_set::realPrint(std::string s, int numOfQuestion, std::vector<int> pIndex, std::ofstream &fp)
{
    std::cout << "\nRound Finish.        The number of questions asked: " << numOfQuestion << "\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";
    std::cout << std::setw(10) << "ID" << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
              << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
              << std::setw(16) << "Seats\n";
    std::cout << "------------------------------------------------------------------------------------------------------------\n";

    for(int k = 0; k < pIndex.size(); ++k)
    {
        std::cout << std::setw(10) << points[pIndex[k]]->id;
        cout << std::setw(14) << points[pIndex[k]]->attr[0] / 1000 << "k";
        for (int i = 1; i < points[pIndex[k]]->d; ++i)
            cout << std::setw(15) << points[pIndex[k]]->attr[i];
        std::cout << "\n------------------------------------------------------------------------------------------------------------\n";
    }

    printf("\nPlease give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how \n"
           "satisfy you are when seeing the recommended cars. (Note: 10 denotes that \n"
           "you are very satisfied with the recommended cars and 1 denotes that you are\n"
           "unsatisfied with the recommended cars): ");
    int satisfy = 0;
    while (satisfy > 10 || satisfy < 1)
        std::cin >> satisfy;

    printf("\nPlease give a number from 1 to 10 (i.e., 1, 2, .., 10) to indicate how \n"
           "bored you feel when you are asked with %d questions in this round in \n"
           "order to obtain the recommended cars (Note: 10 denotes that you feel the \n"
           "most bored and 1 denotes that you feel the least bored.): ", numOfQuestion);
    int bor = 0;
    while (bor > 10 || bor < 1)
        std::cin >> bor;

    fp << "Algorithm: " << s << "   Question: " << numOfQuestion << "\n";
    fp << "------------------------------------------------------------------------------------------------------------\n";
    fp << std::setw(10) << "ID" << std::setw(18) << "Price(¥)" << std::setw(13) << "Year"
       << std::setw(16) << "Used KM" << std::setw(16) << "Length(mm)" << std::setw(15) << "Width(mm)"
       << std::setw(16) << "Seats\n";
    fp << "------------------------------------------------------------------------------------------------------------\n";
    for(int k = 0; k < pIndex.size(); ++k)
    {
        fp << std::setw(10) << points[pIndex[k]]->id;
        fp << std::setw(14) << points[pIndex[k]]->attr[0] / 1000 << "k";
        for (int i = 1; i < points[pIndex[k]]->d; ++i)
            fp << std::setw(15) << points[pIndex[k]]->attr[i];
        fp << "\n------------------------------------------------------------------------------------------------------------\n";
    }

    fp << "Satisfaction: " << satisfy << "\n";
    fp << "Boredness: " << bor << "\n\n\n";
}
























