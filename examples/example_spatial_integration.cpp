#include <iostream>
#include <iomanip>
#include <vector>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"
#include "core/Spectrum.hh"
#include "core/Coordinates.hh"

/**
 * @file example_spatial_integration.cpp
 * @brief 演示空间积分和地球流量计算的完整示例
 */

int main() {
    std::cout << "\n";
    std::cout << "========================================================\n";
    std::cout << "  EDGE C++ Spatial Integration Example              \n";
    std::cout << "  Computing Electron Flux at Earth Position         \n";
    std::cout << "========================================================\n";

    try {
        // 1. 设置源参数（例如：Geminga脉冲星）
        edge::EdgeConfig config;
        config.name = "Geminga";
        config.distance = 0.25;     // kpc
        config.alpha = 2.2;         // 光谱指数
        config.age = 3.42e5;        // years
        config.d0 = 4.0e27;         // cm^2/s
        config.delta = 0.33;        // 扩散指数
        config.emin = 0.001;        // TeV
        config.emax = 100.0;        // TeV (为了演示速度，缩小范围)
        config.ebins = 20;          // 减少bins以加快计算
        config.rbins = 100;

        std::cout << "\nSource Configuration:" << std::endl;
        std::cout << "Name: " << config.name << std::endl;
        std::cout << "Distance: " << config.distance << " kpc" << std::endl;
        std::cout << "Age: " << config.age << " years" << std::endl;
        std::cout << "Spectral index: " << config.alpha << std::endl;

        // 2. 初始化工具和计算器
        edge::EdgeUtils utils;

        // 3. 计算源位置（银道坐标）
        edge::AstroCoordinates::GalacticCoords source_coords;
        source_coords.l = 195.1;    // Geminga的银经（度）
        source_coords.b = -4.2;     // Geminga的银纬（度）
        source_coords.distance = config.distance;

        std::cout << "Galactic coordinates: (l=" << source_coords.l << "°, b="
                  << source_coords.b << "°, d=" << source_coords.distance << " kpc)" << std::endl;

        // 4. 设置地球位置
        edge::AstroCoordinates::CartesianCoords earth_pos(8.5, 0.0, 0.0); // 银河系中心坐标系

        // 5. 初始化计算器
        edge::DiffusionCalculator diff_calc(config, &utils);
        edge::SpectrumCalculator spectrum_calc(config, &utils);
        edge::AstroCoordinates coords(&utils);

        // 6. 计算能量轨迹
        std::cout << "\n=== Computing Energy Trajectory ===" << std::endl;
        auto trajectory = diff_calc.calculateEnergyTrajectory();

        // 7. 准备插值数据
        edge::EnergyTrajectoryData traj_data;
        traj_data.time = trajectory.time;
        traj_data.energy = trajectory.energy;
        traj_data.lambda = trajectory.lambda;
        traj_data.prepareInterpolationArrays();

        spectrum_calc.setEnergyTrajectory(traj_data);

        // 8. 创建计算网格
        std::cout << "=== Creating Computational Grid ===" << std::endl;
        auto [energies, radii] = diff_calc.createEnergyRadiusGrid();

        std::cout << "Energy bins: " << energies.size() << std::endl;
        std::cout << "Radial bins: " << radii.size() << std::endl;

        // 9. 计算电子密度分布
        std::cout << "\n=== Computing Electron Density Distribution ===" << std::endl;
        std::vector<std::vector<double>> density_2d(energies.size(), std::vector<double>(radii.size()));

        // 对有限的能量点计算（演示目的）
        size_t max_energies = std::min(size_t(10), energies.size());
        for (size_t i = 0; i < max_energies; ++i) {
            for (size_t j = 0; j < radii.size(); ++j) {
                density_2d[i][j] = diff_calc.calculateElectronDensity(energies[i], radii[j], config.age);
            }

            if ((i + 1) % 5 == 0) {
                std::cout << "Progress: " << (i + 1) << "/" << max_energies << " energy bins" << std::endl;
            }
        }

        // 10. 计算源到地球的方向
        std::cout << "\n=== Computing Source Direction ===" << std::endl;
        auto source_cartesian = coords.galacticToCartesian(source_coords, earth_pos);

        double dx = source_cartesian.x - earth_pos.x;
        double dy = source_cartesian.y - earth_pos.y;
        double dz = source_cartesian.z - earth_pos.z;

        double l_deg = std::atan2(dy, dx) * 180.0 / M_PI;
        double dist_xy = std::sqrt(dx*dx + dy*dy);
        double b_deg = std::atan2(dz, dist_xy) * 180.0 / M_PI;

        std::cout << "Source direction from Earth: (l=" << l_deg << "°, b=" << b_deg << "°)" << std::endl;

        // 11. 进行视线积分计算
        std::cout << "\n=== Computing Line of Sight Integration ===" << std::endl;
        edge::LineOfSightCalculator los_calc(config, &utils);

        std::vector<double> fluxes = los_calc.lineOfSightIntegration(
            l_deg, b_deg, density_2d, energies, radii, earth_pos);

        // 12. 显示结果
        std::cout << "\n=== Results: Electron Flux at Earth ===" << std::endl;
        std::cout << "Energy [TeV] | Electron Flux [1/(erg*cm^2)]" << std::endl;
        std::cout << "-------------|-------------------------------" << std::endl;

        for (size_t i = 0; i < std::min(max_energies, fluxes.size()); ++i) {
            double energy_TeV = edge::EdgeUtils::ergToTeV(energies[i]);

            std::cout << std::fixed << std::setprecision(3)
                      << std::setw(11) << energy_TeV << " | "
                      << std::scientific << std::setprecision(2)
                      << fluxes[i] << std::fixed << std::endl;
        }

        // 13. 计算总流量
        double total_flux = 0.0;
        for (size_t i = 0; i < std::min(max_energies, fluxes.size()); ++i) {
            total_flux += fluxes[i];
        }

        std::cout << "\n=== Summary ===" << std::endl;
        std::cout << "Source: " << config.name << std::endl;
        std::cout << "Distance to Earth: " << source_coords.distance << " kpc" << std::endl;
        std::cout << "Angular position: (l=" << l_deg << "°, b=" << b_deg << "°)" << std::endl;
        std::cout << "Total integrated electron flux: "
                  << std::scientific << total_flux << std::fixed
                  << " 1/(erg*cm^2)" << std::endl;

        std::cout << "\n=== Example Completed Successfully ===" << std::endl;
        std::cout << "\nThis demonstrates the complete workflow from source parameters" << std::endl;
        std::cout << "to electron flux calculation at Earth position using:" << std::endl;
        std::cout << "✓ GAMERA physics engine" << std::endl;
        std::cout << "✓ Coordinate transformations" << std::endl;
        std::cout << "✓ Energy trajectory calculations" << std::endl;
        std::cout << "✓ Electron density modeling" << std::endl;
        std::cout << "✓ Line-of-sight integration" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << std::endl;
        return 1;
    }
}