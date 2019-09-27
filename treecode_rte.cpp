//
// Created by lurker on 9/26/19.
//

#include "treecode_rte.h"

treecode_rte::treecode_rte(index_t _time_steps, scalar_t _T, index_t _nChebyshev, vector<point> &_source,
        index_t _nSource, index_t _rank, index_t _maxLevel, scalar_t _MAC) {

    /*
     * at each time step, solve a forward stepping based on the treecode algorithm.
     * the previous treecode structures are stored.
     */
    time_steps = _time_steps;
    T = _T;
    h = T/time_steps; // h is the time step, should be less than the spatial resolution.
    MAC = _MAC;

    Vector charge(_nSource);
    setValue(charge, 1.0);

    tc.initalize(_nChebyshev, _source, charge, _nSource, _rank, _maxLevel);
    tc.upPass(0);

    rhs.resize(time_steps);
    solution.resize(time_steps);

    for (index_t step = 0; step < time_steps; ++step) {
        // each step calculates the solution on the time interval [step * dt, (step + 1) * dt]

        solution[step].resize(_nSource);

        for (index_t point_id = 0; point_id < _nSource; ++point_id) {
            // for each point, compute the interaction between this point and other point clusters.
            solution[step](point_id) = interaction(point_id, 0);
            // update rhs for the next step.
        }
    }
}

scalar_t treecode_rte::interaction(index_t pointId, index_t rootId) {
    // compare the distance and diameter to match MAC.
    auto x0 = tc.t.sourceTree[pointId].x;
    auto y0 = tc.t.sourceTree[pointId].y;
    auto x1 = tc.t.dict[rootId].center.x;
    auto y1 = tc.t.dict[rootId].center.y;
    auto dist = sqrt(SQR(x0-x1) + SQR(y0-y1));
    auto diam = sqrt(SQR(tc.t.dict[rootId].radius.x) + SQR(tc.t.dict[rootId].radius.y));

    auto potential = 0.;

    if (tc.t.dict[rootId].isLeaf) {
        // direct summation.
        return 0.; // todo
    }
    else if (diam/dist < MAC) {
        // far field interaction.
        return 0.; // todo
    }
    else {
        // compute the interaction on the children and take summation.
        for (int i = 0; i < 4; ++i) {
            potential += interaction(pointId, tc.t.dict[rootId].child[i]);
        }
        return potential;
    }
}


treecode_rte::~treecode_rte() = default;

