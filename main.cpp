#include "treecode_rte.h"
#include "utility/config.h"

/*
 * Computational domain is set as [0, 1]^2. light speed is normalized as 1.
 * The time interval is h = 1/512. Total traveling time is set as 1 as well, hence there are 512 steps.
 * The spatial resolution is set as dx = 1/64, 1/128, 1/256, 1/512.
 * The Chebyshev nodes are using nc = 4, 6, 8.
 *
 * The parameter can be changed through config file eventually. h <= dx is required to keep maximum principle.
 */

int main(int argc, char *argv[]) {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());

#endif

    if (argc <= 1) {
        std::cout << "USE " << argv[0] << " PATH OF CONFIG FILE " << std::endl;
        exit(0);
    }


    config cfg;
    std::ifstream cfgFile;

    cfgFile.open(argv[1], std::ifstream::in);
    cfg.parse(cfgFile);
    cfgFile.close();

    cfg.print();


    /*
     * time setting
     */
    index_t time_steps = atoi(cfg.options["time_step"].c_str());
    scalar_t T = atof(cfg.options["T"].c_str());

    /*
     * domain setting
     */
    int N = atoi(cfg.options["N"].c_str());
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
    index_t nChebyshev = atoi(cfg.options["Cheb"].c_str());
    index_t rank = nChebyshev * nChebyshev;
    index_t maxLevel = atoi(cfg.options["maxLevel"].c_str());
    scalar_t MAC = atof(cfg.options["MAC"].c_str());


    auto tc_rte = treecode_rte(time_steps, T, N, nChebyshev, source, nSource, rank, maxLevel, MAC);

    RUN("COMPUTE", tc_rte.compute());


    std::string outputFile = cfg.options["output_dir"] + "solution-T-" \
                + std::to_string(T) + "-step-"+std::to_string(time_steps)+"-N-"+std::to_string(N) + "-nC-" + std::to_string(nChebyshev) \
                + "-MAC-"+std::to_string(MAC) + ".txt";

    tc_rte.output(outputFile);


}