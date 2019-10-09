//
// Created by lurker on 9/26/19.
//

#ifndef TREECODE_RTE_H
#define TREECODE_RTE_H

#include "utility/Profiler.h"
#include "treecode.h"

class treecode_rte {
public:
    index_t time_steps;
    scalar_t T;
    scalar_t h;
    index_t  N;
    scalar_t MAC;
    index_t nChebshev;
    scalar_t kappa;



    TreeCode tc;
    index_t nSource;
    index_t rank;
    index_t maxLevel;

    vector<scalar_t> mu_t; // predefined.
    vector<scalar_t> mu_s; // predefined.

    std::function<scalar_t (scalar_t, point&)> sourceFunc; // only mid-point samples are needed

    vector<TreeCode> rhs;
    vector<Vector> solution;

    treecode_rte(index_t _time_steps, scalar_t _T, index_t _N, index_t _nChebyshev, vector<point> &_source,
                 index_t _nSource, index_t _rank, index_t _maxLevel, scalar_t _MAC);
    ~treecode_rte();

    void compute();

    scalar_t interaction(index_t step, index_t pointId, index_t  rootId);
    index_t getRow(scalar_t y);
    index_t getCol(scalar_t x);
    scalar_t getIntegral(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1);
    scalar_t getLocalIntegral(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1);
    scalar_t getAttribute(scalar_t x, scalar_t y);
    scalar_t integral_block(scalar_t a, scalar_t b);
    scalar_t kernel(point& r0, point& r1);

    void output(std::string& filename);

};


#endif //TREECODE_RTE_H
