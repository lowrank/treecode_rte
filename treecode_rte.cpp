//
// Created by lurker on 9/26/19.
//

#include "treecode_rte.h"

treecode_rte::treecode_rte(index_t _time_steps, scalar_t _T,  index_t _N, index_t _nChebyshev, vector<point> &_source,
        index_t _nSource, index_t _rank, index_t _maxLevel, scalar_t _MAC) {

    /*
     * at each time step, solve a forward stepping based on the treecode algorithm.
     * the previous treecode structures are stored.
     */
    time_steps = _time_steps + 1;
    T = _T;
    h = T/(time_steps - 1); // h is the time step, should be less than the spatial resolution.

    assert(h <= 1.0/N);

    N = _N;
    MAC = _MAC;

    nChebshev = _nChebyshev;
    nSource = _nSource;
    rank = _rank;
    maxLevel = _maxLevel;

    Vector charge(nSource);
    setValue(charge, 0.0);

    tc.initalize(_nChebyshev, _source, charge, nSource, rank, maxLevel);
    tc.upPass(0);

    rhs.resize(time_steps);
    solution.resize(time_steps);

    // set mu_t
    mu_t.resize(nSource);
    mu_s.resize(nSource);


    for (int i = 0; i < nSource; ++i) {
        mu_t[i] = 5.2;
        mu_s[i] = 5.0;
//        mu_t[i] = 2.2 + 20 * SQR(tc.t.sourceTree[i].y - 0.5);
//        mu_s[i] = 2. + 20 * SQR(tc.t.sourceTree[i].y - 0.5);; // or set according to point_id
    }

    sourceFunc = [&](scalar_t t, point& p) {
        if (t == 0.) {
            return 0.;
        }

        scalar_t cx = 0.5 + 0.2 * cos(4 * M_PI * t);
        scalar_t cy = 0.5 + 0.2 * sin(4 * M_PI * t);
        scalar_t r2 = (SQR(p.x-cx) + SQR(p.y-cy));

        return 4 * t * exp(-r2 * 40);

    };


//    sourceFunc = [&](scalar_t t, point& p) {
//        scalar_t cx = t;
//        scalar_t cy = 0.5;
//        scalar_t r2 = (SQR(p.x-cx) + SQR(p.y-cy));
//
//        return  4 * SQR(t) * exp(-r2 * 40);
//
//    };



    kappa = kernel(tc.t.sourceTree[0], tc.t.sourceTree[0]);

}

void treecode_rte::compute() {
    for (index_t step = 0; step < time_steps; ++step) {
        // each step calculates the solution on the time interval {step * dt}
        std::cout << "time step: " << step << std::endl;

        Vector charge(nSource);

        for (index_t point_id = 0; point_id < nSource;++point_id) {
            charge(point_id) = sourceFunc(step * h, tc.t.sourceTree[point_id]);
        }

        rhs[step].initalize(nChebshev, tc.t.sourceTree, charge, nSource, rank, maxLevel);
        rhs[step].upPass(0);
        solution[step].resize(nSource);

#ifdef RUN_OMP
#pragma omp parallel for shared(step) default (none)
#endif
        for (index_t point_id = 0; point_id < nSource; ++point_id) {
            // for each point, compute the interaction between this point and other point clusters.
            solution[step](point_id) = interaction(step, point_id, 0) / (1 - kappa * mu_s[point_id]); // multiply a factor
            // update rhs for the next step.
        }

        rhs[step].reset(0);
        for (index_t point_id = 0; point_id < nSource;++point_id) {
            rhs[step].chargeTree(point_id) = charge(point_id) +  mu_s[point_id] * solution[step](point_id);
        }
        rhs[step].upPass(0);
    }
}

scalar_t treecode_rte::interaction(index_t step, index_t pointId, index_t rootId) {
    // compare the distance and diameter to match MAC.
    auto x0 = tc.t.sourceTree[pointId].x;
    auto y0 = tc.t.sourceTree[pointId].y;
    auto x1 = tc.t.dict[rootId].center.x;
    auto y1 = tc.t.dict[rootId].center.y;
    auto dist = sqrt(SQR(x0-x1) + SQR(y0-y1));
    auto diam = sqrt(SQR(tc.t.dict[rootId].radius.x) + SQR(tc.t.dict[rootId].radius.y));


    if (tc.t.dict[rootId].isLeaf) {
        // direct summation or use vector inner product from ddot.
        Vector l_weights(tc.t.dict[rootId].nSource);
        Vector r_weights(tc.t.dict[rootId].nSource);
        Vector l_charge(tc.t.dict[rootId].nSource);
        Vector r_charge(tc.t.dict[rootId].nSource);

        Vector m_weights(tc.t.dict[rootId].nSource);
        Vector m_charge(tc.t.dict[rootId].nSource);

        for (int i = 0; i < tc.t.dict[rootId].nSource;++i) {
            auto id = tc.t.dict[rootId].sourceIndex[i];
            scalar_t _dist = sqrt(SQR(x0 - tc.t.sourceTree[id].x) + SQR(y0 - tc.t.sourceTree[id].y));

            index_t l_step = step - index_t(floor(_dist / h));
            scalar_t alpha = _dist / h - floor(_dist / h);

            scalar_t theta = kernel(tc.t.sourceTree[pointId], tc.t.sourceTree[id]);
//            l_weights(i) = theta * (1-alpha);
//            r_weights(i) = theta * (alpha);
            l_weights(i) = theta * (alpha - 1) * (alpha - 2) / 2.0;
            r_weights(i) = -theta * alpha * (alpha - 2);
            m_weights(i) = theta * alpha * (alpha - 1) / 2.0;

            if (l_step >= 0) {
                l_charge(i) = rhs[l_step].t.dict[rootId].charge(i);
            }
            else{
                l_charge(i) = 0.;
            }

            if (l_step - 1 >= 0) {
                r_charge(i) = rhs[l_step - 1].t.dict[rootId].charge(i);
            }
            else {
                r_charge(i) = 0.;
            }

            if (l_step - 2 >= 0) {
                m_charge(i) = rhs[l_step - 2].t.dict[rootId].charge(i);
            } else {
                m_charge(i) = 0.;
            }

        }
        return ddot(l_weights, l_charge) + ddot(r_weights, r_charge) \
 + ddot(m_weights, m_charge);

    }
    else if (diam/dist < MAC) {
        // far field interaction.
        vector<point> clusterCnode;
        for (int i = 0; i < tc.nChebyshev; ++i) {
            for (int j = 0; j < tc.nChebyshev; ++j) {
                point chebyPoint; // order is important
                chebyPoint.x = tc.t.dict[rootId].scaledCnode[i].x;
                chebyPoint.y = tc.t.dict[rootId].scaledCnode[j].y;
                clusterCnode.push_back(chebyPoint);
            }
        }

        Vector l_weights(tc.rank);
        Vector r_weights(tc.rank);
        Vector l_charge(tc.rank);
        Vector r_charge(tc.rank);


        for (int i = 0; i < tc.rank;++i) {
            scalar_t _dist = sqrt(SQR(x0 - clusterCnode[i].x) + SQR(y0 - clusterCnode[i].y));

            index_t l_step = step - index_t(floor(_dist / h));
            scalar_t alpha = _dist / h - floor(_dist / h);

            scalar_t theta = kernel(tc.t.sourceTree[pointId], clusterCnode[i]);
            l_weights(i) = theta * (1-alpha);
            r_weights(i) = theta * (alpha);

            if (l_step >= 0) {
                l_charge(i) = rhs[l_step].t.dict[rootId].nodeCharge(i);
            }
            else{
                l_charge(i) = 0.;
            }

            if (l_step - 1 >= 0) {
                r_charge(i) = rhs[l_step - 1].t.dict[rootId].nodeCharge(i);
            }
            else {
                r_charge(i) = 0.;
            }

        }
        return ddot(l_weights, l_charge) + ddot(r_weights, r_charge);
    }
    else {
        // compute the interaction on the children and take summation.
        scalar_t potential = 0.;
        // can use **reduce**
        for (int i = 0; i < 4; ++i) {
            potential += interaction(step, pointId, tc.t.dict[rootId].child[i]);
        }
        return potential;
    }
}

index_t treecode_rte::getRow(scalar_t y) {
    return index_t (floor(y * N));
}

index_t treecode_rte::getCol(scalar_t x) {
    return index_t (floor(x * N));
}

scalar_t treecode_rte::getIntegral(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) {
    return getLocalIntegral(x0, y0, x1, y1);
    auto col0 = getCol(x0);
    auto col1 = getCol(x1);
    auto row0 = getRow(y0);
    auto row1 = getRow(y1);

    if ((row0 == row1) && (col0 == col1)) {
        return getLocalIntegral(x0, y0, x1, y1);
    }
    else if ((row0 == row1 + 1) && (col0 == col1)) {
        auto row0 = getCol(y0);
        scalar_t ybar = scalar_t(row0)/N;
        scalar_t xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return getLocalIntegral(x0, y0, xbar, ybar) + getLocalIntegral(xbar, ybar, x1, y1);
    }
    else if ((row0 == row1 - 1) && (col0 == col1)) {
        auto row1 = getCol(y1);
        scalar_t ybar = scalar_t(row1)/N;
        scalar_t xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return getLocalIntegral(x0, y0, xbar, ybar) + getLocalIntegral(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 + 1) && (row0 == row1)) {
        auto col0 = getCol(x0);
        scalar_t xbar = scalar_t(col0)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        return getLocalIntegral(x0, y0, xbar, ybar) + getLocalIntegral(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 - 1) && (row0 == row1)) {
        auto col1 = getCol(x1);
        scalar_t xbar = scalar_t(col1)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        return getLocalIntegral(x0, y0, xbar, ybar) + getLocalIntegral(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 + 1) && (row0 == row1 + 1)) {
        scalar_t xbar = scalar_t(col0)/N;
        scalar_t ybar2 = scalar_t(row0)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar < xbar2) {
            return getLocalIntegral(x1, y1, xbar, ybar) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar2, ybar2, x0, y0);
        }
        else {
            return getLocalIntegral(x1, y1, xbar2, ybar2) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar, ybar, x0, y0);

        }
    }
    else if ((col0 == col1 + 1) && (row0 == row1 - 1)) {
        scalar_t xbar = scalar_t(col0)/N;
        scalar_t ybar2 = scalar_t(row1)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar < xbar2) {
            return getLocalIntegral(x1, y1, xbar, ybar) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar2, ybar2, x0, y0);
        }
        else {
            return getLocalIntegral(x1, y1, xbar2, ybar2) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar, ybar, x0, y0);

        }

    }
    else if ((col0 == col1 - 1) && (row0 == row1 + 1)) {
        scalar_t xbar = scalar_t(col1)/N;
        scalar_t ybar2 = scalar_t(row0)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar > xbar2) {
            return getLocalIntegral(x1, y1, xbar, ybar) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar2, ybar2, x0, y0);
        }
        else {
            return getLocalIntegral(x1, y1, xbar2, ybar2) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar, ybar, x0, y0);
        }
    }
    else if ((col0 == col1 - 1) && (row0 == row1 - 1)) {
        scalar_t xbar = scalar_t(col1)/N;
        scalar_t ybar2 = scalar_t(row1)/N;
        scalar_t ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        scalar_t xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar > xbar2) {
            return getLocalIntegral(x1, y1, xbar, ybar) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar2, ybar2, x0, y0);
        }
        else {
            return getLocalIntegral(x1, y1, xbar2, ybar2) + \
                        getLocalIntegral(xbar, ybar, xbar2, ybar2) + getLocalIntegral(xbar, ybar, x0, y0);

        }
    }
    else {
        scalar_t xm = (x0 + x1) * 0.5;
        scalar_t ym = (y0 + y1) * 0.5;
        return getIntegral(x0, y0, xm, ym) + getIntegral(xm, ym, x1, y1);
    }

}

scalar_t treecode_rte::getLocalIntegral(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) {
    // mid point rule for local
    scalar_t cx = (x0 + x1) * 0.5;
    scalar_t cy = (y0 + y1) * 0.5;
    scalar_t l = sqrt(SQR(x0-x1) + SQR(y0-y1));
    scalar_t mu = getAttribute(cx, cy);
    return l * mu;
}

scalar_t treecode_rte::getAttribute(scalar_t x, scalar_t y) {
    // closet point.
    auto col = index_t (N * x);
    auto row = index_t (N * y);
    return mu_t[row * N + col]; // row major.
}

scalar_t treecode_rte::integral_block(scalar_t a, scalar_t b) {
        int sgnA = a > 0?1:-1;
        int sgnB = b > 0?1:-1;
        scalar_t absA = fabs(a);
        scalar_t absB = fabs(b);
        return sgnA * sgnB * (-absA * log(absA) - absB * log(absB) + \
                absB * log(absA + sqrt(absA * absA + absB * absB)) + \
                absA * log(absB + sqrt(absA * absA + absB * absB)));
    }

scalar_t treecode_rte::kernel(point &r0, point &r1) {
    scalar_t dx = r1.x - r0.x;
    scalar_t dy = r1.y - r0.y;
    scalar_t dist = sqrt(SQR(dx) + SQR(dy));
    scalar_t l = 1.0/N;
    scalar_t f    = fabs(integral_block(dx + l / 2, dy + l / 2) + integral_block(dx - l / 2, dy - l / 2) - \
    integral_block(dx - l / 2, dy + l / 2) - integral_block(dx + l / 2, dy - l / 2));
    return exp(-getIntegral(r0.x, r0.y, r1.x, r1.y)) * f / (2 * M_PI);
}

void treecode_rte::output(const std::string filename) {
    std::ofstream File;
    File.open(filename);

    for (index_t step = 0; step < time_steps; ++step) {
        for (index_t i = 0; i < tc.t.nSource; ++i) {
            File << solution[step](i) << " ";
        }
        File << "\n";
    }
    File.close();
}




treecode_rte::~treecode_rte() = default;

