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
    scalar_t MAC{};

    TreeCode tc;

    std::function<scalar_t (point&, point&)> eval;
    std::function<scalar_t (scalar_t, point&)> f; // only mid-point samples are needed

    vector<TreeCode> rhs;

    treecode_rte(index_t _time_steps, scalar_t _T, index_t _nChebyshev, vector<point> &_source,
                 index_t _nSource, index_t _rank, index_t _maxLevel);
    ~treecode_rte();
};


#endif //TREECODE_RTE_H
