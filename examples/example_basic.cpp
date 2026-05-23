/**
 * @file example_basic.cpp
 * @brief EDGE C++ 基础使用示例
 *
 * 这个示例展示了如何使用 EDGE C++ 进行基本的电子扩散计算
 */

#include <iostream>
#include <iomanip>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"

int main() {
    std::cout << "EDGE C++ Basic Usage Example" << std::endl;
    std::cout << "=============================" << std::endl;

    // 1. 创建配置
    edge::EdgeConfig config;
    config.name = "Geminga";
    config.distance = 0.25;     // kpc
    config.alpha = 2.2;         // 光谱指数
    config.age = 3.42e5;        // years
    config.d0 = 4.0e27;         // cm^2/s
    config.delta = 0.33;        // 扩散指数

    std::cout << "\nConfiguration:" << std::endl;
    std::cout << "Source: " << config.name << std::endl;
    std::cout << "Distance: " << config.distance << " kpc" << std::endl;
    std::cout << "Age: " << config.age << " years" << std::endl;

    try {
        // 2. 初始化 EDGE 工具
        edge::EdgeUtils utils;

        // 3. 创建计算器
        edge::DiffusionCalculator diff_calc(config, &utils);
        edge::PulsarEvolution pulsar(config);

        // 4. 计算能量轨迹
        std::cout << "\nCalculating energy trajectory..." << std::endl;
        auto trajectory = diff_calc.calculateEnergyTrajectory();

        std::cout << "Trajectory calculated:" << std::endl;
        std::cout << "  Steps: " << trajectory.time.size() << std::endl;
        std::cout << "  Max time: " << trajectory.time.back() << " years" << std::endl;

        // 5. 计算一些样本点的电子密度
        std::cout << "\nSample electron densities:" << std::endl;
        std::cout << "Energy[TeV] | Radius[pc] | Density[1/(erg*cm^3)]" << std::endl;
        std::cout << "-------------|------------|------------------------" << std::endl;

        std::vector<double> test_energies = {
            edge::EdgeUtils::teVToErg(1.0),
            edge::EdgeUtils::teVToErg(10.0),
            edge::EdgeUtils::teVToErg(100.0)
        };

        std::vector<double> test_radii = {1.0, 10.0, 50.0};  // pc

        for (double energy : test_energies) {
            for (double radius : test_radii) {
                double density = diff_calc.calculateElectronDensity(
                    energy, radius, config.age);

                std::cout << std::setw(11) << edge::EdgeUtils::ergToTeV(energy) << " | "
                          << std::setw(10) << radius << " | "
                          << std::scientific << std::setprecision(2)
                          << density << std::fixed << std::endl;
            }
        }

        // 6. 计算脉冲星光度演化
        std::cout << "\nPulsar luminosity evolution:" << std::endl;
        std::cout << "Time[years] | Luminosity[erg/s]" << std::endl;
        std::cout << "------------|------------------" << std::endl;

        std::vector<double> time_points = {
            0, 1000, 10000, 100000, config.age
        };

        for (double t : time_points) {
            double lum = pulsar.calculateLuminosity(t);
            std::cout << std::setw(11) << static_cast<int>(t) << " | "
                      << std::scientific << std::setprecision(2)
                      << lum << std::fixed << std::endl;
        }

        std::cout << "\n=== Example completed successfully! ===" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}