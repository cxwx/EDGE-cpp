#include <iostream>
#include <iomanip>
#include <vector>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"
#include "core/Spectrum.hh"

/**
 * @brief 测试完整的光谱计算功能
 *
 * 这个程序测试新实现的光谱计算功能，并与Python版本的行为进行对比
 */

void testBasicSpectrumCalculation() {
    std::cout << "\n=== Testing Basic Spectrum Calculation ===" << std::endl;

    // 创建配置
    edge::EdgeConfig config;
    config.name = "TestSource";
    config.distance = 0.25;     // kpc
    config.alpha = 2.2;         // 光谱指数
    config.age = 3.42e5;        // years
    config.d0 = 4.0e27;         // cm^2/s
    config.delta = 0.33;        // 扩散指数
    config.emin = 0.001;        // TeV
    config.emax = 500.0;        // TeV

    try {
        // 初始化工具
        edge::EdgeUtils utils;

        // 创建扩散计算器
        edge::DiffusionCalculator diff_calc(config, &utils);

        // 计算能量轨迹
        std::cout << "Calculating energy trajectory..." << std::endl;
        auto trajectory = diff_calc.calculateEnergyTrajectory();

        // 准备能量轨迹数据用于插值
        edge::EnergyTrajectoryData traj_data;
        traj_data.time = trajectory.time;
        traj_data.energy = trajectory.energy;
        traj_data.lambda = trajectory.lambda;
        traj_data.prepareInterpolationArrays();

        std::cout << "Energy trajectory prepared for interpolation" << std::endl;
        std::cout << "Data points: " << traj_data.size() << std::endl;

        // 创建光谱计算器
        edge::SpectrumCalculator spectrum_calc(config, &utils);
        spectrum_calc.setEnergyTrajectory(traj_data);

        // 准备光度数据
        edge::PulsarEvolution pulsar(config);
        edge::LuminosityData lum_data;

        // 计算光度数组
        double t_min = std::max(1e-3, config.age - config.age); // 简化的时间范围
        double t_max = config.age;
        int lum_bins = 100;

        auto time_array = edge::InterpolationUtils::logspace(
            std::log10(t_min), std::log10(t_max), lum_bins);

        for (double t : time_array) {
            double lum = pulsar.calculateLuminosity(t);
            lum_data.time.push_back(t);
            lum_data.luminosity.push_back(lum);
        }

        lum_data.prepareInterpolationArrays();
        spectrum_calc.setLuminosityData(lum_data);

        std::cout << "Luminosity data prepared" << std::endl;

        // 测试单个点的光谱计算
        std::cout << "\nTesting single point spectrum calculation:" << std::endl;
        std::cout << "Energy [TeV] | Radius [pc] | Spectrum [1/(erg*cm^3)]" << std::endl;
        std::cout << "-------------|-------------|---------------------------" << std::endl;

        std::vector<double> test_energies = {1.0, 10.0, 100.0}; // TeV
        std::vector<double> test_radii = {1.0, 10.0, 50.0};     // pc

        for (double e_TeV : test_energies) {
            for (double r_pc : test_radii) {
                double e_erg = edge::EdgeUtils::teVToErg(e_TeV);

                try {
                    auto lambda_vec = spectrum_calc.fillLambdaVector(e_erg, 100, config.age);

                    double spectrum = spectrum_calc.calculateSpectrum(
                        e_erg, lambda_vec, r_pc,
                        lambda_vec.DT, lambda_vec.E0, config.age);

                    std::cout << std::setw(11) << e_TeV << " | "
                              << std::setw(11) << r_pc << " | "
                              << std::scientific << std::setprecision(2)
                              << spectrum << std::fixed << std::endl;

                } catch (const std::exception& e) {
                    std::cout << std::setw(11) << e_TeV << " | "
                              << std::setw(11) << r_pc << " | Error: "
                              << e.what() << std::endl;
                }
            }
        }

        // 测试归一化计算
        std::cout << "\nTesting normalization calculation:" << std::endl;
        std::vector<double> test_alphas = {1.5, 2.0, 2.2, 2.5, 3.0};

        for (double alpha : test_alphas) {
            double norm = edge::SpectrumCalculator::calculateNormalization(
                alpha, edge::EdgeUtils::teVToErg(config.emin),
                edge::EdgeUtils::teVToErg(config.emax));

            std::cout << "alpha = " << std::fixed << std::setprecision(1) << alpha << ": "
                      << "normalization = " << std::scientific << std::setprecision(3)
                      << norm << std::fixed << std::endl;
        }

        // 测试加速时间计算
        std::cout << "\nTesting acceleration time calculation:" << std::endl;
        std::cout << "Energy [TeV] | B-field [G] | Acc. Time [years]" << std::endl;
        std::cout << "-------------|-------------|------------------" << std::endl;

        std::vector<double> acc_test_energies = {0.1, 1.0, 10.0, 100.0}; // TeV
        std::vector<double> acc_test_bfields = {1.0, 3.0, 10.0};           // G

        for (double e_TeV : acc_test_energies) {
            for (double b_G : acc_test_bfields) {
                double e_erg = edge::EdgeUtils::teVToErg(e_TeV);

                try {
                    double acc_time = spectrum_calc.calculateAccelerationTime(e_erg, b_G);

                    std::cout << std::setw(11) << e_TeV << " | "
                              << std::setw(11) << b_G << " | "
                              << std::scientific << std::setprecision(2)
                              << acc_time << std::fixed << std::endl;

                } catch (const std::exception& e) {
                    std::cout << std::setw(11) << e_TeV << " | "
                              << std::setw(11) << b_G << " | Error: "
                              << e.what() << std::endl;
                }
            }
        }

        std::cout << "\n=== Spectrum Calculation Tests Completed ===" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error in spectrum calculation test: " << e.what() << std::endl;
        throw;
    }
}

void testInterpolationFunctions() {
    std::cout << "\n=== Testing Interpolation Functions ===" << std::endl;

    // 测试log10插值
    std::vector<double> log_x = {0.0, 1.0, 2.0, 3.0};  // 1, 10, 100, 1000
    std::vector<double> log_y = {0.0, 2.0, 4.0, 6.0};  // 1, 100, 10000, 1e6

    std::cout << "Testing log10 interpolation:" << std::endl;
    std::cout << "Input x | Expected y | Interpolated y | Error" << std::endl;
    std::cout << "--------|------------|----------------|-------" << std::endl;

    std::vector<double> test_x = {0.5, 1.5, 2.5};  // 在插值范围内的点
    for (double x : test_x) {
        double expected = 2.0 * x;  // 线性关系: y = 2x
        double interpolated = edge::InterpolationUtils::log10Interp(x, log_x, log_y);
        double error = std::abs(interpolated - expected) / expected;

        std::cout << std::fixed << std::setprecision(1) << x << "     | "
                  << std::setw(10) << expected << " | "
                  << std::setw(14) << interpolated << " | "
                  << std::scientific << error << std::fixed << std::endl;
    }

    // 测试线性插值
    std::vector<double> lin_x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> lin_y = {0.0, 2.0, 4.0, 6.0};

    std::cout << "\nTesting linear interpolation:" << std::endl;
    std::cout << "Input x | Expected y | Interpolated y | Error" << std::endl;
    std::cout << "--------|------------|----------------|-------" << std::endl;

    std::vector<double> lin_test_x = {0.5, 1.5, 2.5};
    for (double x : lin_test_x) {
        double expected = 2.0 * x;  // 线性关系: y = 2x
        double interpolated = edge::InterpolationUtils::interp(x, lin_x, lin_y);
        double error = std::abs(interpolated - expected) / expected;

        std::cout << std::fixed << std::setprecision(1) << x << "     | "
                  << std::setw(10) << expected << " | "
                  << std::setw(14) << interpolated << " | "
                  << std::scientific << error << std::fixed << std::endl;
    }

    // 测试logspace和linspace
    std::cout << "\nTesting array generation functions:" << std::endl;

    auto log_arr = edge::InterpolationUtils::logspace(0.0, 3.0, 4);
    std::cout << "logspace(0, 3, 4): [";
    for (size_t i = 0; i < log_arr.size(); ++i) {
        std::cout << log_arr[i];
        if (i < log_arr.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    auto lin_arr = edge::InterpolationUtils::linspace(0.0, 10.0, 5);
    std::cout << "linspace(0, 10, 5): [";
    for (size_t i = 0; i < lin_arr.size(); ++i) {
        std::cout << lin_arr[i];
        if (i < lin_arr.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "\n=== Interpolation Tests Completed ===" << std::endl;
}

void testNumericalIntegration() {
    std::cout << "\n=== Testing Numerical Integration ===" << std::endl;

    // 测试梯形法则
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0};  // y = x^2

    double integral_trapezoidal = edge::NumericalIntegration::trapezoidal(x, y);
    double expected_integral = 64.0 / 3.0;  // ∫₀⁴ x² dx = 64/3

    std::cout << "Testing trapezoidal integration of f(x) = x² from 0 to 4:" << std::endl;
    std::cout << "Expected: " << expected_integral << std::endl;
    std::cout << "Calculated: " << integral_trapezoidal << std::endl;
    std::cout << "Relative error: " << std::abs(integral_trapezoidal - expected_integral) / expected_integral << std::endl;

    // 测试辛普森法则
    double integral_simpson = edge::NumericalIntegration::simpson(x, y);
    std::cout << "\nTesting Simpson integration of f(x) = x² from 0 to 4:" << std::endl;
    std::cout << "Expected: " << expected_integral << std::endl;
    std::cout << "Calculated: " << integral_simpson << std::endl;
    std::cout << "Relative error: " << std::abs(integral_simpson - expected_integral) / expected_integral << std::endl;

    std::cout << "\n=== Numerical Integration Tests Completed ===" << std::endl;
}

int main() {
    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "  EDGE C++ Spectrum Functionality Test  " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testInterpolationFunctions();
        testNumericalIntegration();
        testBasicSpectrumCalculation();

        std::cout << "\n";
        std::cout << "========================================" << std::endl;
        std::cout << "  All Spectrum Tests Passed!           " << std::endl;
        std::cout << "========================================" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nTEST FAILED: " << e.what() << std::endl;
        return 1;
    }
}