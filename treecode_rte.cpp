//
// Created by lurker on 9/26/19.
//

#include "treecode_rte.h"

treecode_rte::treecode_rte(index_t _time_steps, scalar_t _T, index_t _nChebyshev, vector<point> &_source,
        index_t _nSource, index_t _rank, index_t _maxLevel) {

    /*
     * at each time step, solve a forward stepping based on the treecode algorithm.
     * the previous treecode structures are stored.
     */
    time_steps = _time_steps;
    T = _T;
    h = T/time_steps; // h is the time step, should be less than the spatial resolution.

    Vector charge(_nSource);
    setValue(charge, 1.0);

    tc.initalize(_nChebyshev, _source, charge, _nSource, _rank, _maxLevel);
    tc.upPass(0);

    rhs.resize(time_steps);

    for (index_t step = 0; step < time_steps; ++step) {
        // each step calculates the solution on the time interval [step * dt, (step + 1) * dt]
        for (index_t point_id = 0; point_id < _nSource; ++point_id) {
            // for each point, compute the interaction between this point and other point clusters.





        }

    }



}


treecode_rte::~treecode_rte() = default;

