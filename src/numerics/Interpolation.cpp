#include "numerics/Interpolation.hh"
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace edge {

// InterpolationUtils 实现

double InterpolationUtils::log10Interp(double log_x, const std::vector<double>& log_x_array,
                                       const std::vector<double>& log_y_array) {
    if (log_x_array.size() != log_y_array.size()) {
        throw std::invalid_argument("log_x_array and log_y_array must have the same size");
    }

    if (log_x_array.empty()) {
        throw std::invalid_argument("Input arrays cannot be empty");
    }

    // 检查是否超出范围
    if (log_x <= log_x_array.front()) {
        return log_y_array.front();
    }
    if (log_x >= log_x_array.back()) {
        return log_y_array.back();
    }

    // 找到插值区间
    size_t idx = findNearestIndex(log_x, log_x_array);

    // 确保idx在有效范围内且log_x_array[idx] <= log_x
    while (idx < log_x_array.size() - 1 && log_x_array[idx + 1] <= log_x) {
        idx++;
    }
    while (idx > 0 && log_x_array[idx] > log_x) {
        idx--;
    }

    // 线性插值
    double x1 = log_x_array[idx];
    double x2 = log_x_array[idx + 1];
    double y1 = log_y_array[idx];
    double y2 = log_y_array[idx + 1];

    if (x2 == x1) {
        return y1;
    }

    double t = (log_x - x1) / (x2 - x1);
    return y1 + t * (y2 - y1);
}

double InterpolationUtils::interp(double x, const std::vector<double>& x_array,
                                  const std::vector<double>& y_array) {
    if (x_array.size() != y_array.size()) {
        throw std::invalid_argument("x_array and y_array must have the same size");
    }

    if (x_array.empty()) {
        throw std::invalid_argument("Input arrays cannot be empty");
    }

    // 检查是否超出范围
    if (x <= x_array.front()) {
        return y_array.front();
    }
    if (x >= x_array.back()) {
        return y_array.back();
    }

    // 找到插值区间
    size_t idx = findNearestIndex(x, x_array);

    // 确保idx在有效范围内且x_array[idx] <= x
    while (idx < x_array.size() - 1 && x_array[idx + 1] <= x) {
        idx++;
    }
    while (idx > 0 && x_array[idx] > x) {
        idx--;
    }

    // 线性插值
    double x1 = x_array[idx];
    double x2 = x_array[idx + 1];
    double y1 = y_array[idx];
    double y2 = y_array[idx + 1];

    if (x2 == x1) {
        return y1;
    }

    double t = (x - x1) / (x2 - x1);
    return y1 + t * (y2 - y1);
}

std::vector<double> InterpolationUtils::logspace(double log10_min, double log10_max, int num_bins) {
    if (num_bins <= 0) {
        throw std::invalid_argument("num_bins must be positive");
    }

    std::vector<double> result;
    result.reserve(num_bins);

    double delta = (log10_max - log10_min) / (num_bins - 1);

    for (int i = 0; i < num_bins; ++i) {
        double log_val = log10_min + i * delta;
        result.push_back(std::pow(10.0, log_val));
    }

    return result;
}

std::vector<double> InterpolationUtils::linspace(double min, double max, int num_bins) {
    if (num_bins <= 0) {
        throw std::invalid_argument("num_bins must be positive");
    }

    std::vector<double> result;
    result.reserve(num_bins);

    double delta = (max - min) / (num_bins - 1);

    for (int i = 0; i < num_bins; ++i) {
        result.push_back(min + i * delta);
    }

    return result;
}

size_t InterpolationUtils::findNearestIndex(double value, const std::vector<double>& array) {
    if (array.empty()) {
        throw std::invalid_argument("Array cannot be empty");
    }

    // 使用std::lower_bound找到第一个>=value的元素
    auto it = std::lower_bound(array.begin(), array.end(), value);

    if (it == array.begin()) {
        return 0;
    }

    if (it == array.end()) {
        return array.size() - 1;
    }

    // 检查哪个元素更接近
    size_t idx1 = std::distance(array.begin(), it) - 1;
    size_t idx2 = std::distance(array.begin(), it);

    if (std::abs(array[idx1] - value) <= std::abs(array[idx2] - value)) {
        return idx1;
    } else {
        return idx2;
    }
}

std::vector<size_t> InterpolationUtils::argsort(const std::vector<double>& array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&array](size_t i1, size_t i2) {
                  return array[i1] < array[i2];
              });

    return indices;
}

// SplineInterpolator 实现

InterpolationUtils::SplineInterpolator::SplineInterpolator(
    const std::vector<double>& x, const std::vector<double>& y,
    const gsl_interp_type* interp_type)
    : x_(x), y_(y) {

    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y must have the same size");
    }

    if (x.size() < 3) {
        throw std::invalid_argument("Need at least 3 points for spline interpolation");
    }

    // 检查x是否单调
    for (size_t i = 1; i < x.size(); ++i) {
        if (x[i] <= x[i-1]) {
            throw std::invalid_argument("x must be strictly monotonic");
        }
    }

    // 创建GSL插值对象
    spline_ = gsl_spline_alloc(interp_type, x.size());
    accel_ = gsl_interp_accel_alloc();

    // 初始化样条
    gsl_spline_init(spline_, x.data(), y.data(), x.size());
}

InterpolationUtils::SplineInterpolator::~SplineInterpolator() {
    if (spline_) {
        gsl_spline_free(spline_);
    }
    if (accel_) {
        gsl_interp_accel_free(accel_);
    }
}

double InterpolationUtils::SplineInterpolator::eval(double x) const {
    return gsl_spline_eval(spline_, x, accel_);
}

double InterpolationUtils::SplineInterpolator::eval_deriv(double x) const {
    return gsl_spline_eval_deriv(spline_, x, accel_);
}

double InterpolationUtils::SplineInterpolator::eval_deriv2(double x) const {
    return gsl_spline_eval_deriv2(spline_, x, accel_);
}

double InterpolationUtils::SplineInterpolator::eval_integ(double a, double b) const {
    return gsl_spline_eval_integ(spline_, a, b, accel_);
}

// EnergyTrajectoryData 实现

void EnergyTrajectoryData::prepareInterpolationArrays() {
    if (time.empty() || energy.empty() || lambda.empty()) {
        return;
    }

    // 准备对数数组
    log_time.clear();
    log_energy.clear();
    log_lambda.clear();

    log_time.reserve(time.size());
    log_energy.reserve(energy.size());
    log_lambda.reserve(lambda.size());

    for (size_t i = 0; i < time.size(); ++i) {
        log_time.push_back(std::log10(time[i]));
        log_energy.push_back(std::log10(energy[i]));
        log_lambda.push_back(std::log10(lambda[i]));
    }

    // 准备时间->能量的对数数组（ETRAJCONT）
    traj_cont_log_time = log_time;
    traj_cont_log_energy = log_energy;

    // 准备能量->时间的对数数组（ETRAJCONTINVERSE）
    // 需要按能量排序
    auto sorted_indices = InterpolationUtils::argsort(energy);

    traj_inv_log_energy.clear();
    traj_inv_log_time.clear();

    traj_inv_log_energy.reserve(energy.size());
    traj_inv_log_time.reserve(energy.size());

    for (size_t idx : sorted_indices) {
        traj_inv_log_energy.push_back(log_energy[idx]);
        traj_inv_log_time.push_back(log_time[idx]);
    }

    // 准备能量->lambda的对数数组（LAMBCONT）
    // 同样需要按能量排序
    lambda_cont_log_energy.clear();
    lambda_cont_log_lambda.clear();

    lambda_cont_log_energy.reserve(energy.size());
    lambda_cont_log_lambda.reserve(lambda.size());

    for (size_t idx : sorted_indices) {
        lambda_cont_log_energy.push_back(log_energy[idx]);
        lambda_cont_log_lambda.push_back(log_lambda[idx]);
    }
}

void EnergyTrajectoryData::clear() {
    time.clear();
    energy.clear();
    lambda.clear();
    log_time.clear();
    log_energy.clear();
    log_lambda.clear();
    traj_cont_log_time.clear();
    traj_cont_log_energy.clear();
    traj_inv_log_energy.clear();
    traj_inv_log_time.clear();
    lambda_cont_log_energy.clear();
    lambda_cont_log_lambda.clear();
}

// LuminosityData 实现

void LuminosityData::prepareInterpolationArrays() {
    if (time.empty() || luminosity.empty()) {
        return;
    }

    log_time.clear();
    log_luminosity.clear();

    log_time.reserve(time.size());
    log_luminosity.reserve(luminosity.size());

    for (size_t i = 0; i < time.size(); ++i) {
        log_time.push_back(std::log10(time[i]));
        log_luminosity.push_back(std::log10(luminosity[i]));
    }
}

void LuminosityData::clear() {
    time.clear();
    luminosity.clear();
    log_time.clear();
    log_luminosity.clear();
}

} // namespace edge