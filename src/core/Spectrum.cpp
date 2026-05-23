#include "core/Spectrum.hh"
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <gsl/gsl_integration.h>

namespace edge {

// SpectrumCalculator 实现

SpectrumCalculator::SpectrumCalculator(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

double SpectrumCalculator::calculateAccelerationTime(double energy, double b_field) const {
    // 计算粒子加速到给定能量所需的最小时间
    // 对应Python中的acc_time函数

    // 动量计算: p = (E - m_e) / c
    double momentum = (energy - m_e) / c_speed;

    // 回旋半径: r = p * c / (e * B)
    double gyroradius = momentum * c_speed / (el_charge * b_field);

    // 加速时间计算（这里使用简化模型）
    // 实际实现可能需要更复杂的加速模型
    double acc_time_sec = gyroradius / c_speed;  // 简化的特征时间尺度

    // 转换为年
    return acc_time_sec / yr_to_sec;
}

SpectrumCalculator::LambdaVectorResult
SpectrumCalculator::fillLambdaVector(double energy_now, int bins, double age) {
    // 对应Python中的FillLambdaVector函数

    LambdaVectorResult result;

    // 检查能量轨迹数据是否已准备好
    if (!trajectory_data_.isValid()) {
        throw std::runtime_error("Energy trajectory data not set or invalid");
    }

    // 插值得到当前能量对应的时间 DT
    double log_energy_now = std::log10(energy_now);
    double log_DT = InterpolationUtils::log10Interp(
        log_energy_now,
        trajectory_data_.traj_inv_log_energy,
        trajectory_data_.traj_inv_log_time
    );

    double DT = std::pow(10.0, log_DT);

    // 确定初始能量 E0
    double E0;
    if (DT < age) {
        // 如果插值时间小于年龄，可以使用EMAX作为初始能量
        E0 = EdgeUtils::teVToErg(config_.emax);
    } else {
        // 否则需要根据时间差计算初始能量
        double log_time_diff = std::log10(DT - age);
        double log_E0 = InterpolationUtils::log10Interp(
            log_time_diff,
            trajectory_data_.traj_cont_log_time,
            trajectory_data_.traj_cont_log_energy
        );
        E0 = std::pow(10.0, log_E0);
    }

    // 创建能量数组（从energy_now到E0）
    auto energies = InterpolationUtils::logspace(
        std::log10(energy_now),
        std::log10(E0),
        bins
    );

    result.energies = energies;
    result.DT = DT;
    result.E0 = E0;

    // 计算每个能量点的lambda值
    std::vector<double> lambdas;
    lambdas.reserve(bins);

    // 首先计算energy_now对应的lambda值
    double log_lambda_now = InterpolationUtils::log10Interp(
        log_energy_now,
        trajectory_data_.lambda_cont_log_energy,
        trajectory_data_.lambda_cont_log_lambda
    );
    double lambda_now = std::pow(10.0, log_lambda_now);

    // 对每个能量点计算lambda差值
    for (size_t i = 1; i < energies.size(); ++i) {
        double e = energies[i];
        double log_e = std::log10(e);

        // 计算该能量对应的lambda
        double log_lambda_e = InterpolationUtils::log10Interp(
            log_e,
            trajectory_data_.lambda_cont_log_energy,
            trajectory_data_.lambda_cont_log_lambda
        );
        double lambda_e = std::pow(10.0, log_lambda_e);

        // lambda差值: lambda(energy_now) - lambda(e)
        double lambda_diff = lambda_now - lambda_e;

        lambdas.push_back(lambda_diff);
    }

    // 填充结果（移除第一个能量点，因为没有对应的lambda）
    result.energies.erase(result.energies.begin());
    result.lambdas = lambdas;

    // 准备对数数组
    result.log_energies.clear();
    result.log_lambdas.clear();

    result.log_energies.reserve(result.energies.size());
    result.log_lambdas.reserve(result.lambdas.size());

    for (size_t i = 0; i < result.energies.size(); ++i) {
        result.log_energies.push_back(std::log10(result.energies[i]));
        // lambda可能是负数或零，需要处理
        if (result.lambdas[i] > 0) {
            result.log_lambdas.push_back(std::log10(result.lambdas[i]));
        } else {
            result.log_lambdas.push_back(-300.0); // log10(很小的数)
        }
    }

    return result;
}

double SpectrumCalculator::calculateSpectrum(double energy,
                                            const LambdaVectorResult& lambda_vector,
                                            double radius, double DT, double E0, double age) {
    // 对应Python中的Spectrum函数

    // 计算最小加速时间
    double t_min = calculateAccelerationTime(E0, config_.bfield);

    if (DT <= t_min) {
        return 0.0;
    }

    // 转换半径到厘米
    double radius_cm = EdgeUtils::pcToCm(radius);

    // 计算光谱归一化
    double emin = EdgeUtils::teVToErg(config_.emin);
    double emax = EdgeUtils::teVToErg(config_.emax);
    double norm = calculateNormalization(config_.alpha, emin, emax);

    // 准备时间数组
    double t_min_log = std::log10(std::max(1e-3, age - DT));
    double t_max_log = std::log10(age - t_min);

    auto time_array = InterpolationUtils::logspace(t_min_log, t_max_log, 2000);

    // 准备时间数组T2 = T - (age - t_min - DT)
    std::vector<double> T2_array;
    T2_array.reserve(time_array.size());

    double time_offset = age - t_min - DT;
    for (double t : time_array) {
        T2_array.push_back(t - time_offset);
    }

    // 计算每个时间点的初始能量e0
    std::vector<double> e0_array;
    e0_array.reserve(T2_array.size());

    for (double T2 : T2_array) {
        if (T2 > 0) {
            double log_T2 = std::log10(T2);
            double log_e0 = InterpolationUtils::log10Interp(
                log_T2,
                trajectory_data_.traj_cont_log_time,
                trajectory_data_.traj_cont_log_energy
            );
            e0_array.push_back(std::pow(10.0, log_e0));
        } else {
            e0_array.push_back(emin); // 如果T2 <= 0，使用最小能量
        }
    }

    // 计算每个时间点的lambda值
    std::vector<double> lambda_array;
    lambda_array.reserve(e0_array.size());

    for (double e0 : e0_array) {
        if (e0 > emin) {
            double log_e0 = std::log10(e0);
            double log_lambda = InterpolationUtils::log10Interp(
                log_e0,
                lambda_vector.log_energies,
                lambda_vector.log_lambdas
            );
            lambda_array.push_back(std::pow(10.0, log_lambda));
        } else {
            lambda_array.push_back(1e-300); // 很小的lambda
        }
    }

    // 计算每个时间点的光度Q
    std::vector<double> Q_array;
    Q_array.reserve(time_array.size());

    if (!luminosity_data_.isValid()) {
        // 如果没有光度数据，使用简化的幂律模型
        double Q0 = 1.0 / norm; // 简化的归一化光度
        (void)Q0; // Suppress unused variable warning
        for (double t : time_array) {
            double tau = config_.tau0;
            double n = config_.brind;
            double luminosity = config_.edot * std::pow(1.0 + t / tau, -(n + 1.0) / (n - 1.0));
            Q_array.push_back(config_.mu * luminosity / norm);
        }
    } else {
        // 使用插值光度数据
        for (double t : time_array) {
            double log_t = std::log10(t);
            double log_lum = InterpolationUtils::log10Interp(
                log_t,
                luminosity_data_.log_time,
                luminosity_data_.log_luminosity
            );
            double luminosity = std::pow(10.0, log_lum);
            Q_array.push_back(luminosity / norm);
        }
    }

    // 计算被积函数值
    std::vector<double> integrand;
    integrand.reserve(time_array.size());

    for (size_t i = 0; i < time_array.size(); ++i) {
        double e0 = e0_array[i];
        double lambda = lambda_array[i];
        double Q = Q_array[i];

        // 检查数值有效性
        if (lambda <= 0 || e0 <= 0 || std::isnan(lambda) || std::isnan(e0)) {
            integrand.push_back(0.0);
            continue;
        }

        // 计算被积函数
        double spatial_term = std::exp(-radius_cm * radius_cm / (4.0 * lambda));
        double energy_term = std::pow(e0, -config_.alpha) * e0 * e0 / (energy * energy);
        double geometric_term = 1.0 / std::pow(4.0 * M_PI * lambda, 1.5);

        double value = Q * energy_term * spatial_term * geometric_term;

        // 检查数值有效性
        if (std::isnan(value) || std::isinf(value)) {
            value = 0.0;
        }

        integrand.push_back(value);
    }

    // 数值积分
    double integral_result = integrateOverTime(time_array, integrand);

    return integral_result;
}

std::vector<std::vector<double>> SpectrumCalculator::calculateSpectrumGrid(
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    double age) {

    std::vector<std::vector<double>> result;
    result.reserve(energies.size());

    for (double energy : energies) {
        std::vector<double> row;
        row.reserve(radii.size());

        // 对每个能量计算一次lambda向量
        auto lambda_vec = fillLambdaVector(energy, 100, age);

        for (double radius : radii) {
            double spectrum = calculateSpectrum(
                energy, lambda_vec, radius,
                lambda_vec.DT, lambda_vec.E0, age
            );
            row.push_back(spectrum);
        }

        result.push_back(row);

        // 输出进度
        static int counter = 0;
        if (++counter % 10 == 0) {
            std::cout << "Progress: " << counter << "/" << energies.size() << " energy bins" << std::endl;
        }
    }

    return result;
}

double SpectrumCalculator::calculateNormalization(double alpha, double emin, double emax) {
    // 计算注入光谱的归一化因子
    if (std::abs(alpha - 2.0) < 1e-6) {
        return std::log(emax / emin);
    } else {
        return (1.0 / (alpha - 2.0)) * (std::pow(emin, -alpha + 2.0) - std::pow(emax, -alpha + 2.0));
    }
}

void SpectrumCalculator::setEnergyTrajectory(const EnergyTrajectoryData& trajectory) {
    trajectory_data_ = trajectory;
    trajectory_data_.prepareInterpolationArrays();
}

void SpectrumCalculator::setLuminosityData(const LuminosityData& luminosity_data) {
    luminosity_data_ = luminosity_data;
    luminosity_data_.prepareInterpolationArrays();
}

double SpectrumCalculator::integrateOverTime(const std::vector<double>& time_points,
                                             const std::vector<double>& values) const {
    // 使用梯形法则进行积分
    return NumericalIntegration::trapezoidal(time_points, values);
}

std::vector<double> SpectrumCalculator::calculateLuminosityArray(double t_min, double t_max, int bins) const {
    // 计算时间相关的光度数组
    auto time_array = InterpolationUtils::logspace(std::log10(t_min), std::log10(t_max), bins);

    std::vector<double> luminosity_array;
    luminosity_array.reserve(bins);

    double tau = config_.tau0;
    double n = config_.brind;

    for (double t : time_array) {
        double luminosity = config_.edot * std::pow(1.0 + t / tau, -(n + 1.0) / (n - 1.0));
        luminosity_array.push_back(config_.mu * luminosity);
    }

    return luminosity_array;
}

// NumericalIntegration 实现

double NumericalIntegration::trapezoidal(const std::vector<double>& x_array,
                                        const std::vector<double>& y_array) {
    if (x_array.size() != y_array.size() || x_array.size() < 2) {
        throw std::invalid_argument("Arrays must have same size and at least 2 elements");
    }

    double integral = 0.0;

    for (size_t i = 1; i < x_array.size(); ++i) {
        double dx = x_array[i] - x_array[i - 1];
        double avg_y = 0.5 * (y_array[i] + y_array[i - 1]);
        integral += dx * avg_y;
    }

    return integral;
}

double NumericalIntegration::simpson(const std::vector<double>& x_array,
                                     const std::vector<double>& y_array) {
    if (x_array.size() != y_array.size() || x_array.size() < 3) {
        throw std::invalid_argument("Arrays must have same size and at least 3 elements");
    }

    double integral = 0.0;

    // 辛普森法则要求等间距点，这里使用简化版本
    for (size_t i = 1; i < x_array.size() - 1; i += 2) {
        double h1 = x_array[i] - x_array[i - 1];
        double h2 = x_array[i + 1] - x_array[i];

        if (std::abs(h1 - h2) > 1e-6 * h1) {
            // 如果间距不相等，回退到梯形法则
            return trapezoidal(x_array, y_array);
        }

        integral += (h1 / 3.0) * (y_array[i - 1] + 4.0 * y_array[i] + y_array[i + 1]);
    }

    // 处理剩余的点
    if (x_array.size() % 2 == 0) {
        size_t n = x_array.size() - 1;
        double h = x_array[n] - x_array[n - 1];
        integral += (h / 2.0) * (y_array[n - 1] + y_array[n]);
    }

    return integral;
}

double NumericalIntegration::gslIntegration(double (*func)(double, void*), void* params,
                                            double a, double b, double tolerance) {
    // 使用GSL进行自适应积分
    gsl_function F;
    F.function = func;
    F.params = params;

    double result, error;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    int status = gsl_integration_qags(&F, a, b, 0, tolerance, 1000, w, &result, &error);

    gsl_integration_workspace_free(w);

    if (status != GSL_SUCCESS) {
        throw std::runtime_error("GSL integration failed");
    }

    return result;
}

double NumericalIntegration::integratePairs(const std::vector<std::pair<double, double>>& pairs,
                                           double x_min, double x_max) {
    if (pairs.empty()) {
        return 0.0;
    }

    // 找到在积分范围内的数据点
    std::vector<double> x_filtered, y_filtered;

    for (const auto& pair : pairs) {
        if (pair.first >= x_min && pair.first <= x_max) {
            x_filtered.push_back(pair.first);
            y_filtered.push_back(pair.second);
        }
    }

    if (x_filtered.empty()) {
        return 0.0;
    }

    return trapezoidal(x_filtered, y_filtered);
}

} // namespace edge