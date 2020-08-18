#define HPX_HPX_MAIN_IMPL_HPP

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>

#include <boost/program_options.hpp>

#include <cstring>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include <glm/glm.hpp>

#include "DataManager.h"
#include "SweepQueue.h"
#include "TreeConstructor.h"

#include "Arc.h"

vtkStandardNewMacro(MouseInteractorStyle2);
Datamanager* DataValueComparator::data = nullptr;
Datamanager* ArcDataValueComparator::data = nullptr;

int main(int argc, char* argv[])
{
    boost::program_options::options_description desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

    // clang-format off
    desc_commandline.add_options()
            ("resample", boost::program_options::value<std::string>()->default_value("100"), "Resamples the input file with the given axis scale percentages (format: X,Y,Z or XxYxZ or N).")
            ("input", boost::program_options::value<std::string>()->default_value("ctBones.vti"), "Path to input data that should be used.")
            ("visualize", "Resulting tree embedded in input data will be rendered with vtk")
            ("no-trunkskip", "Perform explicit trunk computation and instead of collecting dangling saddles.");

    return hpx::init(desc_commandline, argc, argv);
}

int hpx_main(boost::program_options::variables_map& vm)
{
    std::cout << "Threads: " << hpx::get_num_worker_threads() << std::endl;

    // Read Command Line Config

    Options opt;
    opt.visualize = false;
    opt.trunkSkipping = true;

    if (vm.count("visualize"))
        opt.visualize = true;

    if (vm.count("no-trunkskip"))
        opt.trunkSkipping = false;

    opt.inputFile = vm["input"].as<std::string>();

    glm::ivec3 resample;

    std::vector<std::string> resampleComponents;
    boost::split(resampleComponents, vm["resample"].as<std::string>(), boost::is_any_of(",x"));

    if (resampleComponents.size() == 3)
        resample = glm::ivec3(std::stoi(resampleComponents[0]), std::stoi(resampleComponents[1]), std::stoi(resampleComponents[2]));
    else if (resampleComponents.size() == 1)
        resample = glm::ivec3(std::stoi(resampleComponents[0]));
    else
        throw std::runtime_error("Invalid resample size");

    opt.resampleX = resample.x;
    opt.resampleY = resample.y;
    opt.resampleZ = resample.z;

    // Create and Initialize local Tree Constructor
    hpx::id_type joinTree = hpx::new_<::server::TreeConstructor>(hpx::find_here()).get();
    hpx::async<::server::TreeConstructor::init_action>(joinTree, opt).get();

    // Global Barrier
    // or set event in init_action and wait for event in every remote action

    // Start Tree Construction

    hpx::async<::server::TreeConstructor::construct_action>(joinTree).get();

    return hpx::finalize();
}
