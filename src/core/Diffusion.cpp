#include "core/Diffusion.hh"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace edge {

// DiffusionCalculator 实现
DiffusionCalculator::DiffusionCalculator(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }

    std::cout << "DiffusionCalculator initialized with D0 = " << config_.d0
              << " cm^2/s, delta = " << config_.delta << std::endl;
}

double DiffusionCalculator::calculateDiffusionCoefficient(double energy) const {
    // 计算扩散系数 D(E) = D0 * (1 + E/E_star)^delta
    double E_star = calculateE_star();
    double energy_norm = 1.0 + energy / E_star;
    return config_.d0 * std::pow(energy_norm, config_.delta);
}

double DiffusionCalculator::calculateEnergyLossRate(double energy) const {
    if (config_.use_kn) {
        // 使用 GAMERA 的 Klein-Nishina 能量损失率
        return utils_->getParticles()->EnergyLossRate(energy);
    } else {
        // 使用简化的汤姆逊散射近似
        double gamma = utils_->calculateGamma(energy);
        return utils_->calculateThomsonLossRate(gamma, config_.tot_e_dens);
    }
}

DiffusionCalculator::EnergyTrajectory DiffusionCalculator::calculateEnergyTrajectory() {
    EnergyTrajectory traj;

    // 能量范围（erg）
    double E_max = EdgeUtils::teVToErg(config_.emax);
    double E_min = EdgeUtils::teVToErg(config_.emin);
    double E_star = calculateE_star();
    (void)E_star; // Suppress unused variable warning
    (void)E_min; // Suppress unused variable warning
    (void)E_max; // Suppress unused variable warning

    // 初始化
    double e = E_max;
    double t = 0.0;
    double diff_int = 0.0;  // 积分扩散系数

    std::cout << "Calculating energy trajectory from "
              << EdgeUtils::ergToTeV(E_max) << " TeV to "
              << EdgeUtils::ergToTeV(E_min) << " TeV" << std::endl;

    int steps = 0;
    const int max_steps = 100000;  // 安全限制

    while (e > E_min && steps < max_steps) {
        // 计算当前能量下的损失率
        double loss_rate = calculateEnergyLossRate(e);

        // 计算时间步长（能量的小分数）
        double dt = EdgeConstants::ENERGY_LOSS_FACTOR * e / loss_rate;

        // 更新能量和时间
        double energy_loss = dt * loss_rate;
        e -= energy_loss;
        t += dt / EdgeUtils::yrToSec(1.0);  // 转换为年

        // 计算扩散系数
        double D = calculateDiffusionCoefficient(e);

        // 积分扩散系数 ∫(D/E_dot) dE
        diff_int += dt * loss_rate * D / loss_rate;

        // 存储结果
        traj.time.push_back(t);
        traj.energy.push_back(e);
        traj.lambda.push_back(diff_int);
        traj.diff_coeff.push_back(D);

        steps++;

        // 每1000步输出进度
        if (steps % 1000 == 0) {
            std::cout << "Step " << steps << ": E = "
                      << EdgeUtils::ergToTeV(e) << " TeV, t = "
                      << t << " years" << std::endl;
        }
    }

    std::cout << "Energy trajectory calculation completed in " << steps << " steps" << std::endl;
    std::cout << "Final time: " << t << " years" << std::endl;

    trajectory_ = traj;
    return traj;
}

std::pair<std::vector<double>, std::vector<double>>
DiffusionCalculator::createEnergyRadiusGrid() {
    // 创建能量网格（对数间隔）
    std::vector<double> energies;
    double E_min = EdgeUtils::teVToErg(config_.emin);
    double E_max = EdgeUtils::teVToErg(config_.emax);

    double log_E_min = std::log10(E_min);
    double log_E_max = std::log10(E_max);
    double d_log_E = (log_E_max - log_E_min) / (config_.ebins - 1);

    for (int i = 0; i < config_.ebins; ++i) {
        double log_E = log_E_min + i * d_log_E;
        energies.push_back(std::pow(10.0, log_E));
    }

    // 创建半径网格（对数间隔，单位 pc）
    std::vector<double> radii;
    double r_min = 1e-3;  // pc
    double r_max = 1e3 * config_.distance * 20.0;  // pc
    double log_r_min = std::log10(r_min);
    double log_r_max = std::log10(r_max);
    double d_log_r = (log_r_max - log_r_min) / (config_.rbins - 1);

    for (int i = 0; i < config_.rbins; ++i) {
        double log_r = log_r_min + i * d_log_r;
        radii.push_back(std::pow(10.0, log_r));
    }

    std::cout << "Created energy grid: " << energies.size()
              << " bins from " << EdgeUtils::ergToTeV(energies.front())
              << " to " << EdgeUtils::ergToTeV(energies.back()) << " TeV" << std::endl;
    std::cout << "Created radius grid: " << radii.size()
              << " bins from " << radii.front() << " to "
              << radii.back() << " pc" << std::endl;

    return {energies, radii};
}

double DiffusionCalculator::calculateElectronDensity(double energy, double radius, double age) {
    // 计算电子数密度 n(E, r, t)
    // 使用扩散方程的格林函数解

    double radius_cm = EdgeUtils::pcToCm(radius);
    double age_sec = EdgeUtils::yrToSec(age);
    (void)age_sec; // Suppress unused variable warning

    // 查找对应的 λ 值（从轨迹中）
    auto it = std::lower_bound(trajectory_.energy.begin(), trajectory_.energy.end(), energy);
    if (it == trajectory_.energy.end()) {
        return 0.0;  // 能量超出范围
    }

    size_t idx = std::distance(trajectory_.energy.begin(), it);
    double lambda = trajectory_.lambda[idx];

    // 简化的格林函数解
    // n(E, r, t) ∝ Q(E0) * exp(-r²/4λ) / (4πλ)^(3/2)
    // 这里是简化版本，实际实现需要更复杂的积分

    double diffusion_length = 4.0 * lambda;
    double spatial_term = std::exp(-radius_cm * radius_cm / diffusion_length);
    double normalization = 1.0 / std::pow(4.0 * M_PI * lambda, 1.5);

    // 注入光谱 Q(E) ∝ E^(-α)
    double injection_spectrum = std::pow(energy, -config_.alpha);

    return injection_spectrum * spatial_term * normalization;
}

DiffusionCalculator::ElectronDensityGrid
DiffusionCalculator::calculateElectronDensityGrid(
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    double age) {

    ElectronDensityGrid grid(energies.size(), std::vector<double>(radii.size()));

    std::cout << "Calculating electron density grid: "
              << energies.size() << " x " << radii.size() << std::endl;

    // 确保轨迹已计算
    if (trajectory_.energy.empty()) {
        calculateEnergyTrajectory();
    }

    // 计算网格（这里可以并行化）
    for (size_t i = 0; i < energies.size(); ++i) {
        for (size_t j = 0; j < radii.size(); ++j) {
            grid[i][j] = calculateElectronDensity(energies[i], radii[j], age);
        }

        if ((i + 1) % 10 == 0) {
            std::cout << "Progress: " << (i + 1) << "/" << energies.size() << " energy bins" << std::endl;
        }
    }

    return grid;
}

double DiffusionCalculator::calculateE_star() const {
    return EdgeConstants::E_STAR_DEFAULT * TeV_to_erg;
}

double DiffusionCalculator::integrateLossFunction(double energy) {
    // 积分损失函数 ∫ D(E)/E_dot dE
    // 这个函数用于计算有效扩散长度

    double result = 0.0;
    double E = energy;
    double dE = energy * 0.001;  // 小步长

    while (E > EdgeUtils::teVToErg(config_.emin)) {
        double loss_rate = calculateEnergyLossRate(E);
        double diff_coeff = calculateDiffusionCoefficient(E);

        if (loss_rate > 0) {
            result += (diff_coeff / loss_rate) * dE;
        }

        E -= dE;
    }

    return result;
}

// ElectronSpectrum 实现
ElectronSpectrum::ElectronSpectrum(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

double ElectronSpectrum::calculateSpectrum(double energy, double radius, double age) {
    // 计算电子数密度谱
    // 这是一个简化实现，实际需要考虑完整的扩散-损失方程

    DiffusionCalculator diff_calc(config_, utils_);
    return diff_calc.calculateElectronDensity(energy, radius, age);
}

std::vector<std::vector<double>> ElectronSpectrum::calculateSpectrumGrid(
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    double age) {

    DiffusionCalculator diff_calc(config_, utils_);
    return diff_calc.calculateElectronDensityGrid(energies, radii, age);
}

double ElectronSpectrum::calculateInjectionSpectrum(double energy) const {
    // 注入光谱 Q(E) = Q_0 * E^(-α)
    // 归一化常数 Q_0 取决于总能量注入
    double E_min = EdgeUtils::teVToErg(config_.emin);
    double E_max = EdgeUtils::teVToErg(config_.emax);
    (void)E_min; // Suppress unused variable warning
    (void)E_max; // Suppress unused variable warning

    // 总能量 = μ * E_dot * τ
    // 这个简化实现只返回光谱形状
    return std::pow(energy, -config_.alpha);
}

// PulsarEvolution 实现
PulsarEvolution::PulsarEvolution(const EdgeConfig& config)
    : config_(config) {

    // 计算特征时间尺度
    tau_ = config_.tau0;
}

double PulsarEvolution::calculateLuminosity(double time) const {
    // 自旋向下光度 L(t) = L_0 / (1 + t/τ)^(n-1)
    // 其中 n 是衰减指数（brind）

    double L_0 = config_.edot;
    double n = config_.brind;
    double tau = config_.tau0;

    return L_0 / std::pow(1.0 + time / tau, n - 1.0);
}

double PulsarEvolution::calculateTotalEnergy(double age) const {
    // 总能量 = ∫ L(t) dt
    // 对于 L(t) = L_0 / (1 + t/τ)^(n-1)
    // 积分结果取决于 n

    double L_0 = config_.edot;
    double n = config_.brind;
    double tau = config_.tau0;

    if (n > 2.0) {
        // 收敛积分
        return (L_0 * tau / (n - 2.0)) *
               (1.0 - 1.0 / std::pow(1.0 + age / tau, n - 2.0));
    } else {
        // 发散情况，数值积分
        double total = 0.0;
        double dt = age / 1000.0;
        for (double t = 0; t < age; t += dt) {
            total += calculateLuminosity(t) * dt;
        }
        return total;
    }
}

double PulsarEvolution::calculateCharacteristicAge() const {
    // 特征年龄 τ = P / (2 * Ṗ)
    // 这里使用简化的计算

    if (config_.period <= 0 || config_.brind <= 1) {
        return config_.age;
    }

    // 使用初始周期计算
    return (config_.period * config_.period) /
           (2.0 * (config_.period - config_.p0)) *
           (1.0e6 / EdgeUtils::yrToSec(1.0));  // 转换为年
}

} // namespace edge