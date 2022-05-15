#include <operation.h>
#include "point_t.h"

/**
 * @brief Constructor
 * @param dim   The number of dimensions
 */
point_t::point_t(int dim)
{
    d = dim;
    attr = new double[d];
    d_order = 0;
    d_unorder = 0;
    id = -1;
    this->value = 0;
    this->hyper = 0;
}

/**
 * @brief Constructor
 * @param d_order   The number of ordered dimensions
 * @param d_unorder The number of unordered dimensions
 */
point_t::point_t(int d_order, int d_unorder)
{
    this->d_order = d_order;
    this->d_unorder = d_unorder;
    d = d_order + d_unorder;
    attr = new double[d_order + d_unorder];
    id = -1;
    this->value = 0;
    this->hyper = 0;
}

/**
 * @brief Constructor
 * @param d_order   The number of ordered dimensions
 * @param d_unorder The number of unordered dimensions
 * @param id        ID
 */
point_t::point_t(int d_order, int d_unorder, int id)
{
    this->d_order = d_order;
    this->d_unorder = d_unorder;
    d = d_order + d_unorder;
    attr = new double[d_order + d_unorder];
    this->id = id;
    this->value = 0;
    this->hyper = 0;
}

/**
 * @brief Constructor
 * @param p   The point
 */
point_t::point_t(point_t *p)
{
    d = p->d;
    d_order = p->d_order;
    d_unorder = p->d_unorder;
    id = p->id;
    attr = new double[d];
    for(int i = 0;i < d; ++i)
        attr[i] = p->attr[i];
    value = p->value;
    hyper = p->hyper;
}

/*
 * @brief Constructor
 * @param p1   point
 */
point_t::point_t(point_t *p1, point_t *p2)
{
    d = p1->d;
    this->d_order = p1->d_order;
    this->d_unorder = p1->d_unorder;
    attr = new double[d];
    for(int i = 0; i < d; ++i)
        attr[i] = p1->attr[i] - p2->attr[i];
    id = -1;
    this->value = 0;
    this->hyper = 0;
}

/*
 * @brief Destructor
 *        delete the memory of the array
 */
point_t::~point_t()
{
    delete []attr;
}

/*
 * @brief For debug purpose, print the coordinates for a given point
 */
void point_t::print()
{
    std::cout<< id <<"  ";
    for (int i = 0; i < d; i++)
        std::cout<< attr[i] << " ";
    std::cout << "\n";
}

/**
 * @brief       Check whether the point is the same as p
 * @param p     Point
 * @return      1 same
 *              0 different
 */
bool point_t::is_same(point_t *p)
{
    if(d != p->d)
        return false;
    for (int i = 0; i < d; ++i)
    {
        if (attr[i] - p->attr[i] < -EQN2 || attr[i] - p->attr[i] > EQN2)
            return false;
    }
    return true;
}

/**
 * @brief   Check whether all the attribute values are 0
 * @return  -true   all attributes values are 0
 *          -false  there exists attribute value which is not 0
 */
bool point_t::is_zero()
{
    for(int i = 0; i < d; ++i)
    {
        if(attr[i] < -EQN2 || attr[i] > EQN2)
            return false;
    }
    return true;
}

/**
 * @brief	    Calculate the dot product between two points
 * @param p     One point
 */
double point_t::dot_product(point_t *p)
{
    double result = 0;
    for(int i = 0; i < d; i++)
    {
        result += attr[i] * p->attr[i];
    }
    return result;
}

/**
 * @brief	    Calculate the dot product between two points
 * @param v     One array
 */
double point_t::dot_product(double *v)
{
    double result = 0;
    for(int i = 0; i < d; i++)
    {
        result += attr[i] * v[i];
    }
    return result;
}

/**
 * @brief	Calculate the subtraction between two points.
 *          Remember to release the returned point to save memory.
 * @param p The subtractor
 * @return  The subtraction(new points)
 */
point_t *point_t::sub(point_t *p)
{
    int dim = d_order + d_unorder;
    point_t* result = new point_t(d_order, d_unorder);
    for(int i = 0; i < dim; i++)
    {
        result->attr[i] = attr[i] - p->attr[i];
    }
    return result;
}

/**
 * @brief	Calculate the addition between two points.
 *          Remember to release the returned point to save memory.
 * @param p The point
 * @return  The addition(new points)
 */
point_t *point_t::add(point_t *p)
{
    point_t* result = new point_t(d_order, d_unorder);
    for(int i = 0; i < d; i++)
    {
        result->attr[i] = attr[i] + p->attr[i];
    }
    return result;
}

/**
 * @brief	Scale the point
 *          Remember to release the returned point to save memory.
 * @param c The scaled coefficient
 * @return  The scaled point
 */
point_t *point_t::scale(double c)
{
    point_t* result = new point_t(d_order, d_unorder);
    for(int i = 0; i < d; i++)
    {
        result->attr[i] = attr[i] * c;
    }
    return result;
}


/**
 * @brief       Calculate the distance between two points
 * @param p     The points
 * @return      The distance
 */
double point_t::distance(point_t *p)
{
    double diff = 0;
    for(int i = 0; i < d; i++)
    {
        diff += (double) pow(attr[i]/100 - p->attr[i]/100, 2);
    }
    return 100 * sqrt(diff);
}

/**
 * @brief       Calculate the distance between a points and a vector
 * @param p     The vector
 * @return      The distance
 */
double point_t::distance(double *p)
{
    double diff = 0;
    for(int i = 0; i < d; i++)
    {
        diff += (double) pow(attr[i] - p[i], 2);
    }
    return sqrt(diff);
}

/**
 * @brief Print the result of the algorithm
 * @param out_cp    The name of the output file
 * @param name      The name of the algorithm
 * @param Qcount    The number of question asked
 * @param t1        The start time
 */
void point_t::printResult(char *name, int Qcount, timeval t1)
{
    timeval t2;
    std::ofstream out_cp("../../result.txt");
    gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << "-----------------------------------------------------------------\n";
    printf("|%15s |%15d |%15lf |%10d |\n", name, Qcount, time_cost, id);
    std::cout << "-----------------------------------------------------------------\n";
    out_cp << Qcount << "       " << time_cost << "\n";
    out_cp.close();
}


/**
 * @brief Print the result of the algorithm
 * @param out_cp    The name of the output file
 * @param name      The name of the algorithm
 * @param Qcount    The number of question asked
 * @param t1        The start time
 * @param preTime   The preprocessing time cost
 */
void point_t::printResult(std::ofstream &out_cp, char *name, int Qcount, timeval t1, double preTime, long mem_baseline,
                          double Pcize, int type)
{
    timeval t2; gettimeofday(&t2, 0);
    double time_cost = (double) t2.tv_sec + (double) t2.tv_usec / 1000000 - (double) t1.tv_sec - (double) t1.tv_usec / 1000000;
    std::cout << "--------------------------------------------------------------------------------------\n";
    printf("|%20s |%15d |%15lf |%15lf |%10d |\n", name, Qcount, preTime, time_cost - preTime, id);// get_mem_usage() - mem_baseline);
    std::cout << "--------------------------------------------------------------------------------------\n";
    if(type == 2)
        out_cp << Qcount << "       " << preTime << "       " << time_cost - preTime << "\n";
    //"      " << Pcize << "\n" << "    " << get_mem_usage() - mem_baseline << "\n";
}


/**
 * @brief calculate the place where the order of the two points change
 * @param p  The second point
 * @return
 */
double point_t::bound(point_t *p, double x)
{
    double sum = 0;
    sum =  (this->attr[0] * this->attr[0]) - (p->attr[0] * p->attr[0])
             + (this->attr[1] * this->attr[1]) - (p->attr[1] * p->attr[1]);
    sum = sum / 2 - x * (this->attr[0] - p->attr[0]);
    sum = sum /  (this->attr[1] - p->attr[1]);
    return sum;
}

/**
 * @brief calculate the place where the order of the two points change
 * @param p  The second point
 * @return
 */
double point_t::bound(point_t *p, point_t* exp)
{
    double sum = 0;
    for(int i = 0; i < p->d; ++i)
        sum +=  (this->attr[i] * this->attr[i]) - (p->attr[i] * p->attr[i]);
    sum = sum / 2;
    for(int i = 0; i < d_order; ++i)
        sum -= exp->attr[i] * (this->attr[i] - p->attr[i]);
    sum = sum /  (this->attr[d - 1] - p->attr[d - 1]);
    return sum;
}

/**
 * @brief calculate the place where the order of the two points change
 * @param p  The second point
 * @return
 */
double point_t::tranBound(point_t *p)
{
    double sum = 0;
    sum =  (this->traAttr[0] * this->traAttr[0]) - (p->traAttr[0] * p->traAttr[0])
           + (this->traAttr[1] * this->traAttr[1]) - (p->traAttr[1] * p->traAttr[1]);
    sum = sum / (2 * (this->traAttr[0] - p->traAttr[0]));
    return sum;
}


/**
 * @brief Transform the attributes of the point based on a vector
 * @param vec       The vector
 * @param vecNorm   The vector which is normalized
 */
void point_t::transformAttr(point_t *start, point_t *vec, point_t *vecNorm)
{
    traAttr = new double[2];
    traAttr[0] = 0; traAttr[1] = 0;
    double ttLength = 0;

    double *p = new double[d];
    for(int i = 0; i < d; ++i)
    {
        p[i] = attr[i] - start->attr[i];
        ttLength += p[i] * p[i];
    }

    for(int i = 0; i < d_unorder; ++i)
        traAttr[0] += p[i + d_order] * vecNorm->attr[i];
    traAttr[1] = sqrt(ttLength - traAttr[0] * traAttr[0]);
}

point_t *point_t::generateEXP(point_t* p2)
{
    int dimo = d_order, dimu = d_unorder;
    point_t *e = new point_t(dimo, dimu);
    for(int i = 0; i < dimo; ++i)
    {
        e->attr[i] = attr[i];
    }
    for(int i = dimo; i < dimo + dimu; ++i)
    {
        e->attr[i] = p2->attr[i - dimo];
    }
    return e;
}












