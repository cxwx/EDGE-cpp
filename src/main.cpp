#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"
#include "core/GammaRay.hh"

/**
 * @brief EDGE C++ 主程序
 *
 * 计算电子扩散、伽马射线光谱和地球位置的电子流量
 * 这是原 Python EDGE.py 的 C++ 版本
 */

void printHeader() {
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
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help                 Show this help message\n";
    std::cout << "  -n, --name <string>        Source name (default: Source)\n";
    std::cout << "  -al, --alpha <float>       Spectral index (default: 2.2)\n";
    std::cout << "  -d, --distance <float>     Distance [kpc] (default: 0.25)\n";
    std::cout << "  -del, --delta <float>      Diffusion index (default: 0.33)\n";
    std::cout << "  -a, --age <float>          Age [years] (default: 3.42e5)\n";
    std::cout << "  -emax, --emax <float>      EMAX [TeV] (default: 500.0)\n";
    std::cout << "  -emin, --emin <float>      EMIN [TeV] (default: 0.001)\n";
    std::cout << "  -kn, --kn                  Use Klein-Nishina energy losses\n";
    std::cout << "  -d0, --d0 <float>          Diffusion coefficient [cm^2/s] (default: 4e27)\n";
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    printHeader();

    // 默认配置
    edge::EdgeConfig config;

    // 简单的命令行解析（实际应用中可以使用更复杂的库）
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-n" || arg == "--name") {
            if (i + 1 < argc) config.name = argv[++i];
        } else if (arg == "-al" || arg == "--alpha") {
            if (i + 1 < argc) config.alpha = std::stod(argv[++i]);
        } else if (arg == "-d" || arg == "--distance") {
            if (i + 1 < argc) config.distance = std::stod(argv[++i]);
        } else if (arg == "-del" || arg == "--delta") {
            if (i + 1 < argc) config.delta = std::stod(argv[++i]);
        } else if (arg == "-a" || arg == "--age") {
            if (i + 1 < argc) config.age = std::stod(argv[++i]);
        } else if (arg == "-emax" || arg == "--emax") {
            if (i + 1 < argc) config.emax = std::stod(argv[++i]);
        } else if (arg == "-emin" || arg == "--emin") {
            if (i + 1 < argc) config.emin = std::stod(argv[++i]);
        } else if (arg == "-kn" || arg == "--kn") {
            config.use_kn = true;
        } else if (arg == "-d0" || arg == "--d0") {
            if (i + 1 < argc) config.d0 = std::stod(argv[++i]);
        }
    }

    // 验证配置
    if (!config.validate()) {
        std::cerr << "Invalid configuration parameters" << std::endl;
        return 1;
    }

    // 打印配置
    config.print();

    try {
        // 初始化 EDGE 工具
        std::cout << "\n=== Initializing EDGE C++ ===" << std::endl;
        edge::EdgeUtils utils;

        // 创建计算器
        std::cout << "\n=== Creating Calculator Objects ===" << std::endl;
        edge::DiffusionCalculator diffusion_calc(config, &utils);
        edge::ElectronSpectrum electron_spectrum(config, &utils);
        edge::PulsarEvolution pulsar_evolution(config);

        // 计算能量轨迹
        std::cout << "\n=== Calculating Energy Trajectory ===" << std::endl;
        auto trajectory = diffusion_calc.calculateEnergyTrajectory();

        std::cout << "\nTrajectory Summary:" << std::endl;
        std::cout << "  Number of steps: " << trajectory.time.size() << std::endl;
        std::cout << "  Time range: 0 - " << trajectory.time.back() << " years" << std::endl;
        std::cout << "  Energy range: "
                  << edge::EdgeUtils::ergToTeV(trajectory.energy.front()) << " - "
                  << edge::EdgeUtils::ergToTeV(trajectory.energy.back()) << " TeV" << std::endl;

        // 创建网格
        std::cout << "\n=== Creating Energy-Radius Grid ===" << std::endl;
        auto [energies, radii] = diffusion_calc.createEnergyRadiusGrid();

        // 计算电子密度分布（示例：只计算几个点以节省时间）
        std::cout << "\n=== Calculating Electron Density (sample) ===" << std::endl;
        std::vector<double> sample_energies = {
            edge::EdgeUtils::teVToErg(1.0),   // 1 TeV
            edge::EdgeUtils::teVToErg(10.0),  // 10 TeV
            edge::EdgeUtils::teVToErg(100.0)  // 100 TeV
        };
        std::vector<double> sample_radii = {1.0, 10.0, 50.0}; // pc

        for (double energy : sample_energies) {
            for (double radius : sample_radii) {
                double density = diffusion_calc.calculateElectronDensity(
                    energy, radius, config.age);

                std::cout << "E = " << std::setw(6) << edge::EdgeUtils::ergToTeV(energy) << " TeV, "
                          << "r = " << std::setw(4) << radius << " pc: "
                          << "n = " << std::scientific << density << std::fixed
                          << " 1/(erg*cm^3)" << std::endl;
            }
        }

        // 计算脉冲星光度演化
        std::cout << "\n=== Calculating Pulsar Evolution ===" << std::endl;
        std::vector<double> time_points = {0, 100, 1000, 10000, config.age};
        std::cout << "Time [years] | Luminosity [erg/s]" << std::endl;
        std::cout << "-------------|---------------------" << std::endl;

        for (double t : time_points) {
            double luminosity = pulsar_evolution.calculateLuminosity(t);
            std::cout << std::setw(12) << static_cast<int>(t) << " | "
                      << std::scientific << std::setprecision(3)
                      << luminosity << std::fixed << std::endl;
        }

        double total_energy = pulsar_evolution.calculateTotalEnergy(config.age);
        std::cout << "\nTotal energy injected over " << config.age << " years: "
                  << std::scientific << total_energy << std::fixed
                  << " erg" << std::endl;

        // 计算注入光谱
        std::cout << "\n=== Calculating Injection Spectrum ===" << std::endl;
        std::cout << "Energy [TeV] | Injection Spectrum [arbitrary units]" << std::endl;
        std::cout << "-------------|--------------------------------------" << std::endl;

        for (double energy_TeV : {0.1, 1.0, 10.0, 100.0, 500.0}) {
            double energy = edge::EdgeUtils::teVToErg(energy_TeV);
            double spectrum = electron_spectrum.calculateInjectionSpectrum(energy);

            std::cout << std::setw(11) << energy_TeV << " | "
                      << std::scientific << spectrum << std::fixed << std::endl;
        }

        std::cout << "\n=== EDGE C++ Calculation Complete ===" << std::endl;
        std::cout << "\nAll calculations completed successfully!" << std::endl;
        std::cout << "Results can be compared with Python EDGE.py for validation." << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "\nUNKNOWN ERROR occurred" << std::endl;
        return 1;
    }

    return 0;
}