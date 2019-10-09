#include "treecode_rte.h"

/*
 * Computational domain is set as [0, 1]^2. light speed is normalized as 1.
 * The time interval is h = 1/512. Total traveling time is set as 1 as well, hence there are 512 steps.
 * The spatial resolution is set as dx = 1/64, 1/128, 1/256, 1/512.
 * The Chebyshev nodes are using nc = 4, 6, 8.
 *
 * The parameter can be changed through config file eventually. h < dx is required to keep maximum principle.
 */

int main() {
    /*
     * time setting
     */
    index_t time_steps = 150;
    scalar_t T = 1.0;

    /*
     * domain setting
     */
    int N = 128;
    scalar_t dx = 1.0/N;
    index_t nSource = N * N;
    vector<point> source(nSource);
    for (index_t i = 0; i < N; ++i) {
        for (index_t j = 0; j < N; ++j) {
            source[i*N+j].x = (i+0.5)*dx;
            source[i*N+j].y = (j+0.5)*dx;
        }
    }

    /*
     * algorithm setting
     */
    index_t nChebyshev = 4;
    index_t rank = nChebyshev * nChebyshev;
    index_t maxLevel = 10;
    scalar_t MAC = 0.5;


    auto tc_rte = treecode_rte(time_steps, T, N, nChebyshev, source, nSource, rank, maxLevel, MAC);

    RUN("COMPUTE", tc_rte.compute());

    tc_rte.output("solution.txt");


}