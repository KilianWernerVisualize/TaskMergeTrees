#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/align/aligned_alloc.hpp>
#include <boost/program_options.hpp>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Log.h"
#include "TreeConstructor.h"

/**
 * @brief HPX main entry.
 * @param vm
 * @return
 */
int hpx_main(boost::program_options::variables_map& vm)
{
    /*
     * Parse command line
     */
    std::string input;

    try {
        // Input (positional)
        if (vm.count("hpx:positional")) {
            std::vector<std::string> positionals = vm["hpx:positional"].as<std::vector<std::string>>();

            if (positionals.size() >= 1)
                input = positionals[0];
        }
        if (input.empty())
            throw std::runtime_error("No input specified");

    } catch (const std::exception& e) {
        std::cout << "Parsing error: " << e.what() << std::endl
                  << std::endl;
        std::cout << "Usage: hpxct [options] input" << std::endl
                  << std::endl;
        std::cout << "Check --help (-h) for details." << std::endl;
        return hpx::finalize();
    }

    Log() << "Input: " << input;

    // Print HPX info
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    int numThreads = hpx::get_os_thread_count();
    Log().tag("HPX") << "Localities: " << localities.size();
    Log().tag("HPX") << "Threads: " << numThreads;

    /*
     * Initialize components
     */
    hpx::util::high_resolution_timer timer;

    // Create tree constructor on each locality
    std::vector<hpx::id_type> treeConstructors;
    for (hpx::id_type locality : localities)
        treeConstructors.push_back(hpx::new_<TreeConstructor>(locality).get());

    // Initialize
    timer.restart();
    std::vector<hpx::shared_future<void>> initFutures;
    for (hpx::id_type treeConstructor : treeConstructors)
        initFutures.push_back(hpx::async<TreeConstructor::init_action>(treeConstructor, treeConstructors, input));
    hpx::lcos::wait_all(initFutures);
    Log() << "Initialization: " << timer.elapsed() << " s";

    // Construction
    timer.restart();
    std::vector<hpx::shared_future<void>> constructFutures;
    for (hpx::id_type treeConstructor : treeConstructors)
        constructFutures.push_back(hpx::async<TreeConstructor::construct_action>(treeConstructor));
    hpx::lcos::wait_all(constructFutures);
    Log() << "Construction: " << timer.elapsed() << " s";

    // Destroy
    std::vector<hpx::shared_future<void>> destroyFutures;
    for (hpx::id_type treeConstructor : treeConstructors)
        destroyFutures.push_back(hpx::async<TreeConstructor::destroy_action>(treeConstructor));
    hpx::lcos::wait_all(destroyFutures);

    return hpx::finalize();
}

/**
 * @brief Application entry point.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[])
{
    // Command line description
    // ...

    boost::program_options::options_description allOptionDescriptions("hpxct [options] input");
    //    allOptionDescriptions.add(inputOutputOptionsDescription).add(parallelizationOptionsDescription).add(imageOptionsDescription).add(cameraOptionsDescription).add(debuggingOptionsDescription);

    // HPX config
    std::vector<std::string> const cfg = {
        //        "hpx.stacks.small_size=0x80000"
        "hpx.stacks.use_guard_pages=0"
    };

    // Run HPX
    return hpx::init(allOptionDescriptions, argc, argv, cfg);
}
