#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"
#include "core/Spectrum.hh"
#include "core/Coordinates.hh"
#include "core/MultiSource.hh"

/**
 * @file edge_complete.cpp
 * @brief 完整的EDGE C++程序 - 与Python版本功能对等
 */

class EdgeApplication {
public:
    EdgeApplication() : config_() {
        // calculator_ will be created later
    }

    void run(int argc, char* argv[]) {
        setupDefaultConfig();
        parseCommandLine(argc, argv);
        parseCommandLine(argc, argv);

        if (config_.name.empty() || show_help_) {
            printUsage();
            return;
        }

        if (!config_.validate()) {
            throw std::runtime_error("Invalid configuration parameters");
        }

        executeCalculation();
    }

private:
    edge::EdgeConfig config_;
    std::unique_ptr<edge::EdgeCalculator> calculator_;
    std::unique_ptr<edge::PulsarCatalog> pulsar_catalog_;
    bool show_help_ = false;

    void setupDefaultConfig() {
        // 设置与Python版本相同的默认参数
        config_.name = "Source";
        config_.profile_file = "Data/GemingaProfile.dat";
        config_.alpha = 2.2;
        config_.distance = 0.25;
        config_.delta = 0.33;
        config_.age = 3.42e5;
        config_.emax = 500.0;
        config_.emin = 0.001;
        config_.d0 = 4.0e27;
        config_.size = 5.0;
        config_.tot_e_dens = 1.06;
        config_.bfield = 3.0;
        config_.edot = 3.2e34;
        config_.brind = 3.0;
        config_.tau0 = 1.2e4;
        config_.period = 237.0;
        config_.p0 = 40.5;
        config_.ebins = 100;
        config_.rbins = 400;
    }

    void parseCommandLine(int argc, char* argv[]) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            if (arg == "-h" || arg == "--help") {
                show_help_ = true;
            } else if (arg == "-n" || arg == "--name") {
                if (++i < argc) config_.name = argv[i];
            } else if (arg == "-f" || arg == "--file") {
                if (++i < argc) config_.profile_file = argv[i];
            } else if (arg == "-al" || arg == "--alpha") {
                if (++i < argc) config_.alpha = std::stod(argv[i]);
            } else if (arg == "-d" || arg == "--distance") {
                if (++i < argc) config_.distance = std::stod(argv[i]);
            } else if (arg == "-del" || arg == "--delta") {
                if (++i < argc) config_.delta = std::stod(argv[i]);
            } else if (arg == "-a" || arg == "--age") {
                if (++i < argc) config_.age = std::stod(argv[i]);
            } else if (arg == "-emax" || arg == "--emax") {
                if (++i < argc) config_.emax = std::stod(argv[i]);
            } else if (arg == "-emin" || arg == "--emin") {
                if (++i < argc) config_.emin = std::stod(argv[i]);
            } else if (arg == "-kn" || arg == "--kn") {
                config_.use_kn = true;
            } else if (arg == "-d0" || arg == "--d0") {
                if (++i < argc) config_.d0 = std::stod(argv[i]);
            } else if (arg == "-s" || arg == "--size") {
                if (++i < argc) config_.size = std::stod(argv[i]);
            } else if (arg == "-mu" || arg == "--mu") {
                if (++i < argc) config_.mu = std::stod(argv[i]);
            } else if (arg == "-edens" || arg == "--edens") {
                if (++i < argc) config_.tot_e_dens = std::stod(argv[i]);
            } else if (arg == "-bfield" || arg == "--bfield") {
                if (++i < argc) config_.bfield = std::stod(argv[i]);
            } else if (arg == "-edot" || arg == "--edot") {
                if (++i < argc) config_.edot = std::stod(argv[i]);
            } else if (arg == "-brind" || arg == "--brind") {
                if (++i < argc) config_.brind = std::stod(argv[i]);
            } else if (arg == "-tsupr" || arg == "--tsupr") {
                if (++i < argc) config_.timesupr = std::stod(argv[i]);
            } else if (arg == "-eb" || arg == "--ebins") {
                if (++i < argc) config_.ebins = std::stoi(argv[i]);
            } else if (arg == "-rb" || arg == "--rbins") {
                if (++i < argc) config_.rbins = std::stoi(argv[i]);
            } else if (arg == "-birth_period" || arg == "--birth_period") {
                config_.birth_period = true;
            } else if (arg == "-all_pulsar" || arg == "--all_pulsar") {
                config_.all_pulsar = true;
            } else if (arg == "-only_flux_earth" || arg == "--only_flux_earth") {
                config_.only_flux_earth = true;
            } else if (arg == "-eps" || arg == "--eps") {
                config_.fig_eps = true;
            } else if (arg == "-o" || arg == "--output") {
                if (++i < argc) output_dir_ = argv[i];
            }
        }
    }

    void printUsage() {
        std::cout << "\n";
        std::cout << "#######################################################\n";
        std::cout << "#                                                     #\n";
        std::cout << "#         Calculation of electron spectra,            #\n";
        std::cout << "#         gamma-ray spectra and electrons             #\n";
        std::cout << "#         flux at the Earth for different             #\n";
        std::cout << "#              initial parameters                     #\n";
        std::cout << "#                                                     #\n";
        std::cout << "#######################################################\n";
        std::cout << "#                                                     #\n";
        std::cout << "#     C++ Version - Migrated from Python EDGE.py      #\n";
        std::cout << "#     Using GAMERA C++ Library for Physics            #\n";
        std::cout << "#                                                     #\n";
        std::cout << "#######################################################\n";
        std::cout << "\n";

        std::cout << "Usage: edge_complete [options]\n\n";
        std::cout << "Source Parameters:\n";
        std::cout << "  -n, --name <string>              Source name (default: Source)\n";
        std::cout << "  -f, --file <string>              Angular profile file (default: Data/GemingaProfile.dat)\n";
        std::cout << "  -al, --alpha <float>             Spectral index (default: 2.2)\n";
        std::cout << "  -d, --distance <float>           Distance [kpc] (default: 0.25)\n";
        std::cout << "  -del, --delta <float>            Diffusion index (default: 0.33)\n";
        std::cout << "  -a, --age <float>                Age [years] (default: 3.42e5)\n";
        std::cout << "  -emax, --emax <float>            EMAX [TeV] (default: 500.0)\n";
        std::cout << "  -emin, --emin <float>            EMIN [TeV] (default: 0.001)\n";
        std::cout << "  -d0, --d0 <float>                Diffusion coefficient [cm^2/s] (default: 4e27)\n";
        std::cout << "  -s, --size <float>              Size of source [kpc] (default: 5.0)\n";
        std::cout << "  -mu, --mu <float>                Fraction of energy in electrons (default: 0.5)\n";

        std::cout << "\nPhysics Parameters:\n";
        std::cout << "  -edens, --edens <float>          Total energy density [eV/cm^3] (default: 1.06)\n";
        std::cout << "  -bfield, --bfield <float>       Magnetic field [G] (default: 3.0)\n";
        std::cout << "  -kn, --kn                       Use Klein-Nishina energy losses\n";

        std::cout << "\nPulsar Parameters:\n";
        std::cout << "  -edot, --edot <float>           Spin-down power [erg/s] (default: 3.2e34)\n";
        std::cout << "  -brind, --brind <float>         Breaking index (default: 3.0)\n";
        std::cout << "  -p, --p <float>                 Pulsar period [ms] (default: 237.0)\n";
        std::cout << "  -p0, --p0 <float>                Initial period [ms] (default: 40.5)\n";
        std::cout << "  -tau0, --tau0 <float>           Initial spin-down timescale [years] (default: 1.2e4)\n";

        std::cout << "\nComputational Parameters:\n";
        std::cout << "  -eb, --ebins <int>              Energy bins (default: 100)\n";
        std::cout << "  -rb, --rbins <int>              Radial bins (default: 400)\n";

        std::cout << "\nRuntime Options:\n";
        std::cout << "  -birth_period, --birth_period  Calculate from birth period\n";
        std::cout << "  -all_pulsar, --all_pulsar      Calculate all pulsar contribution\n";
        std::cout << "  -only_flux_earth, --only_flux_earth  Only calculate flux at Earth\n";
        std::cout << "  -eps, --eps                    Save figures in EPS format\n";
        std::cout << "  -o, --output <dir>             Output directory (default: Results)\n";

        std::cout << "\nExamples:\n";
        std::cout << "  edge_complete --name Geminga --distance 0.25 --alpha 2.2\n";
        std::cout << "  edge_complete --all_pulsar --output Results/Geminga\n";
        std::cout << "  edge_complete --name Test --kn --emax 100\n";

        std::cout << "\nFor more information, see the documentation.\n";
    }

    void executeCalculation() {
        std::cout << "\n=== Starting EDGE C++ Calculation ===" << std::endl;
        std::cout << "Configuration:\n";
        config_.print();

        // 创建输出目录
        if (output_dir_.empty()) {
            output_dir_ = "Results/" + config_.name;
        }
        std::string mkdir_cmd = "mkdir -p " + output_dir_;
        system(mkdir_cmd.c_str());

        // 创建计算器
        calculator_ = std::make_unique<edge::EdgeCalculator>(config_);

        // 加载脉冲星目录（如果需要）
        if (config_.all_pulsar) {
            std::cout << "\n=== Generating Pulsar Catalog ===" << std::endl;
            edge::EdgeUtils utils;
            pulsar_catalog_ = std::make_unique<edge::PulsarCatalog>(&utils);

            // 生成均匀分布的脉冲星
            double max_age = 1.0e7; // years
            double sn_rate = 2.0 / 100.0; // per year
            pulsar_catalog_->generateHomogeneousDistribution(max_age, sn_rate, true);

            calculator_->setPulsarCatalog(*pulsar_catalog_);
        }

        // 创建能量数组
        std::vector<double> energies = edge::InterpolationUtils::logspace(
            std::log10(config_.emin),
            std::log10(config_.emax),
            config_.ebins
        );

        // 运行完整计算
        auto result = calculator_->runCompleteCalculation(energies);

        // 保存结果
        calculator_->saveResults(result, output_dir_);

        // 打印结果摘要
        printResultsSummary(result);
    }

    void printResultsSummary(const edge::EdgeCalculator::CalculationResult& result) {
        std::cout << "\n=== Results Summary ===" << std::endl;

        std::cout << "Energy Range: " << config_.emin << " - " << config_.emax << " TeV" << std::endl;
        std::cout << "Computation Time: " << result.calculation_time << " seconds" << std::endl;

        if (!result.earth_flux.empty()) {
            std::cout << "\nEarth Flux [1/(TeV*cm^2*s)]:\n";
            std::cout << "Energy [TeV] | Flux" << std::endl;
            std::cout << "-------------|----------------------" << std::endl;

            for (size_t i = 0; i < std::min(size_t(10), result.earth_flux.size()); ++i) {
                std::cout << std::fixed << std::setprecision(3)
                          << std::setw(11) << result.energies[i] << " | "
                          << std::scientific << std::setprecision(2)
                          << result.earth_flux[i] << std::fixed << std::endl;
            }

            if (result.earth_flux.size() > 10) {
                std::cout << "... (showing first 10 of " << result.earth_flux.size() << " energy bins)" << std::endl;
            }
        }

        if (!result.pulsars.empty()) {
            std::cout << "\nPulsar Statistics:" << std::endl;
            std::cout << "Total pulsars: " << result.pulsars.size() << std::endl;

            // 找出最近的脉冲星
            auto nearby = findNearestPulsars(result.pulsars, 5);
            if (!nearby.empty()) {
                std::cout << "Nearest pulsars:" << std::endl;
                std::cout << "Distance [kpc] | Name" << std::endl;
                std::cout << "---------------|------" << std::endl;

                for (const auto& pulsar : nearby) {
                    std::cout << std::fixed << std::setprecision(3)
                              << std::setw(13) << pulsar.distanceToEarth() << " | "
                              << pulsar.name << std::endl;
                }
            }
        }

        std::cout << "\nResults saved to: " << output_dir_ << std::endl;
        std::cout << "=== Calculation Complete ===" << std::endl;
    }

    std::vector<edge::PulsarData> findNearestPulsars(
        const std::vector<edge::PulsarData>& pulsars, int n) {

        std::vector<edge::PulsarData> sorted = pulsars;
        std::sort(sorted.begin(), sorted.end(),
            [](const edge::PulsarData& a, const edge::PulsarData& b) {
                return a.distanceToEarth() < b.distanceToEarth();
            });

        n = std::min(n, static_cast<int>(sorted.size()));
        return std::vector<edge::PulsarData>(sorted.begin(), sorted.begin() + n);
    }

    std::string output_dir_ = "Results";
};

int main(int argc, char* argv[]) {
    try {
        EdgeApplication app;
        app.run(argc, argv);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}