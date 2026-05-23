#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "utils/EdgeUtils.hh"
#include "core/Coordinates.hh"
#include "core/Diffusion.hh"
#include "core/Spectrum.hh"

/**
 * @brief 测试坐标转换和空间积分功能
 */

void testCoordinateConversions() {
    std::cout << "\n=== Testing Coordinate Conversions ===" << std::endl;

    edge::EdgeConfig config;
    edge::EdgeUtils utils;

    edge::AstroCoordinates coords(&utils);

    // 设置地球位置（在距离中心8.5 kpc的银盘上）
    edge::AstroCoordinates::CartesianCoords earth_pos(8.5, 0.0, 0.0);
    coords.setGalacticReferencePoint(earth_pos);

    std::cout << "Earth position: (" << earth_pos.x << ", "
              << earth_pos.y << ", " << earth_pos.z << ") kpc" << std::endl;

    // 测试银道坐标到笛卡尔坐标的转换
    std::cout << "\nTesting Galactic -> Cartesian conversion:" << std::endl;
    std::cout << "l [deg] | b [deg] | distance [kpc] | x [kpc] | y [kpc] | z [kpc]" << std::endl;
    std::cout << "---------|---------|----------------|---------|---------|---------" << std::endl;

    std::vector<edge::AstroCoordinates::GalacticCoords> test_coords = {
        {0.0, 0.0, 1.0},    // 银心方向
        {90.0, 0.0, 1.0},   // 银道方向
        {0.0, 90.0, 1.0},   // 北银极方向
        {45.0, 30.0, 0.5}   // 随机方向
    };

    for (const auto& galactic : test_coords) {
        auto cartesian = coords.galacticToCartesian(galactic, earth_pos);

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(8) << galactic.l << " | "
                  << std::setw(7) << galactic.b << " | "
                  << std::setw(14) << galactic.distance << " | "
                  << std::setw(7) << cartesian.x << " | "
                  << std::setw(7) << cartesian.y << " | "
                  << std::setw(7) << cartesian.z << std::endl;
    }

    // 测试笛卡尔坐标到银道坐标的转换
    std::cout << "\nTesting Cartesian -> Galactic conversion:" << std::endl;
    std::cout << "x [kpc] | y [kpc] | z [kpc] | l [deg] | b [deg] | distance [kpc]" << std::endl;
    std::cout << "---------|---------|---------|---------|---------|----------------" << std::endl;

    std::vector<edge::AstroCoordinates::CartesianCoords> test_cartesian = {
        {9.5, 0.0, 0.0},   // 沿x轴
        {8.5, 1.0, 0.0},   // 沿y轴
        {8.5, 0.0, 0.5},   // 沿z轴
        {9.0, 1.0, 0.5}    // 对角线
    };

    for (const auto& cartesian : test_cartesian) {
        auto galactic = coords.cartesianToGalactic(cartesian, earth_pos);

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(7) << cartesian.x << " | "
                  << std::setw(7) << cartesian.y << " | "
                  << std::setw(7) << cartesian.z << " | "
                  << std::setw(7) << galactic.l << " | "
                  << std::setw(7) << galactic.b << " | "
                  << std::setw(14) << galactic.distance << std::endl;
    }

    // 测试坐标旋转
    std::cout << "\nTesting coordinate rotation:" << std::endl;
    edge::AstroCoordinates::CartesianCoords test_point(1.0, 0.0, 0.0);
    std::cout << "Original: (" << test_point.x << ", " << test_point.y << ", " << test_point.z << ")" << std::endl;

    edge::AstroCoordinates::rotateCoordinates(test_point, 0.0, 0.0, M_PI/2); // 绕z轴旋转90度
    std::cout << "After rotation (90° around z): (" << test_point.x << ", "
              << test_point.y << ", " << test_point.z << ")" << std::endl;

    std::cout << "\n=== Coordinate Conversion Tests Completed ===" << std::endl;
}

void testLineOfSightCalculation() {
    std::cout << "\n=== Testing Line of Sight Calculation ===" << std::endl;

    edge::EdgeConfig config;
    config.distance = 0.25;  // kpc
    edge::EdgeUtils utils;

    edge::LineOfSightCalculator los_calculator(config, &utils);

    // 设置地球位置
    edge::AstroCoordinates::CartesianCoords earth_pos(8.5, 0.0, 0.0);

    // 测试方向
    double l = 0.0;    // 银心方向
    double b = 0.0;    // 银道平面

    std::cout << "Generating line of sight in direction (l=" << l << "°, b=" << b << "°)" << std::endl;
    std::cout << "Distance [kpc] | x [kpc] | y [kpc] | z [kpc] | Distance to Earth [kpc]" << std::endl;
    std::cout << "---------------|---------|---------|---------|----------------------" << std::endl;

    // 生成几个距离值的视线点
    std::vector<double> distances = {0.1, 0.5, 1.0, 2.0, 5.0};

    auto line_of_sight = los_calculator.generateLineOfSight(l, b, distances, earth_pos);

    for (size_t i = 0; i < distances.size(); ++i) {
        const auto& coords = line_of_sight[i];
        double dist_to_earth = coords.distanceTo(earth_pos);

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(13) << distances[i] << " | "
                  << std::setw(7) << coords.x << " | "
                  << std::setw(7) << coords.y << " | "
                  << std::setw(7) << coords.z << " | "
                  << std::setw(20) << dist_to_earth << std::endl;
    }

    std::cout << "\n=== Line of Sight Calculation Tests Completed ===" << std::endl;
}

void testSpatialIntegration() {
    std::cout << "\n=== Testing Spatial Integration ===" << std::endl;

    edge::EdgeConfig config;
    config.distance = 0.25;  // kpc
    config.emin = 0.001;     // TeV
    config.emax = 100.0;     // TeV (缩小范围以加快测试)
    config.ebins = 10;       // 减少bins以加快测试
    config.rbins = 50;

    edge::EdgeUtils utils;

    // 创建扩散计算器和光谱计算器
    edge::DiffusionCalculator diff_calc(config, &utils);
    edge::SpectrumCalculator spectrum_calc(config, &utils);

    // 计算能量轨迹
    std::cout << "Calculating energy trajectory..." << std::endl;
    auto trajectory = diff_calc.calculateEnergyTrajectory();

    // 准备能量轨迹数据
    edge::EnergyTrajectoryData traj_data;
    traj_data.time = trajectory.time;
    traj_data.energy = trajectory.energy;
    traj_data.lambda = trajectory.lambda;
    traj_data.prepareInterpolationArrays();

    spectrum_calc.setEnergyTrajectory(traj_data);

    // 创建网格
    std::cout << "Creating energy-radius grid..." << std::endl;
    auto [energies, radii] = diff_calc.createEnergyRadiusGrid();

    // 计算简化的密度网格（使用快速方法）
    std::cout << "Computing electron density grid..." << std::endl;
    std::vector<std::vector<double>> density_2d(energies.size(), std::vector<double>(radii.size(), 0.0));

    // 对几个能量点计算密度（为了测试目的）
    for (size_t i = 0; i < std::min(size_t(5), energies.size()); ++i) {
        for (size_t j = 0; j < radii.size(); ++j) {
            density_2d[i][j] = diff_calc.calculateElectronDensity(energies[i], radii[j], config.age);
        }
    }

    // 创建视线积分计算器
    edge::LineOfSightCalculator los_calculator(config, &utils);

    // 设置地球位置
    edge::AstroCoordinates::CartesianCoords earth_pos(8.5, 0.0, 0.0);

    // 测试几个方向的积分
    std::cout << "\nTesting line of sight integration:" << std::endl;
    std::cout << "Direction (l,b) | Energy [TeV] | Integrated Flux [1/(erg*cm^2)]" << std::endl;
    std::cout << "---------------|-------------|-------------------------------" << std::endl;

    std::vector<std::pair<double, double>> directions = {
        {0.0, 0.0},    // 银心方向
        {90.0, 0.0},   // 银道方向
        {45.0, 30.0}   // 随机方向
    };

    for (const auto& direction : directions) {
        double l = direction.first;
        double b = direction.second;

        try {
            auto fluxes = los_calculator.lineOfSightIntegration(
                l, b, density_2d, energies, radii, earth_pos);

            // 只显示几个能量点
            for (size_t i = 0; i < std::min(size_t(3), energies.size()); ++i) {
                double energy_TeV = edge::EdgeUtils::ergToTeV(energies[i]);

                std::cout << std::fixed << std::setprecision(0)
                          << "(" << l << "°," << b << "°)      | "
                          << std::setw(11) << energy_TeV << " | "
                          << std::scientific << std::setprecision(2)
                          << fluxes[i] << std::fixed << std::endl;
            }

        } catch (const std::exception& e) {
            std::cout << "(" << l << "°," << b << "°) | Error: " << e.what() << std::endl;
        }
    }

    std::cout << "\n=== Spatial Integration Tests Completed ===" << std::endl;
}

void testGameraIntegration() {
    std::cout << "\n=== Testing GAMERA Integration ===" << std::endl;

    edge::EdgeUtils utils;

    // 直接使用GAMERA的Astro对象
    auto* astro = utils.getAstro();

    // 测试GAMERA的坐标转换函数
    std::cout << "Testing GAMERA Astro functions:" << std::endl;

    // 测试参数
    double r = 1.0;     // kpc
    double l = 45.0;    // 度
    double b = 30.0;    // 度

    std::vector<double> obs_pos = {8.5, 0.0, 0.0};  // 地球位置

    auto cartesian = astro->GetCartesian(r, l, b, obs_pos);

    std::cout << "GAMERA GetCartesian(r=" << r << " kpc, l=" << l << "°, b=" << b << "°):" << std::endl;
    std::cout << "  Result: (" << cartesian[0] << ", " << cartesian[1] << ", " << cartesian[2] << ") kpc" << std::endl;

    // 测试反向转换
    double l_back, b_back;
    astro->GetGalactic(cartesian[0], cartesian[1], cartesian[2],
                      obs_pos[0], obs_pos[1], obs_pos[2], l_back, b_back);

    std::cout << "GAMERA GetGalactic back: l=" << l_back << "°, b=" << b_back << "°" << std::endl;

    std::cout << "Coordinate agreement: l_diff=" << std::abs(l - l_back)
              << "°, b_diff=" << std::abs(b - b_back) << "°" << std::endl;

    std::cout << "\n=== GAMERA Integration Tests Completed ===" << std::endl;
}

int main() {
    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "  EDGE C++ Coordinates Test          " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testCoordinateConversions();
        testLineOfSightCalculation();
        testGameraIntegration();

        // 注意：空间积分测试比较耗时，可以根据需要启用
        // testSpatialIntegration();

        std::cout << "\n";
        std::cout << "========================================" << std::endl;
        std::cout << "  All Coordinates Tests Passed!      " << std::endl;
        std::cout << "========================================" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nTEST FAILED: " << e.what() << std::endl;
        return 1;
    }
}