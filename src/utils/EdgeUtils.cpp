#include "utils/EdgeUtils.hh"
#include <iostream>
#include <cmath>
#include <cassert>

namespace edge {

// EdgeConfig 实现
bool EdgeConfig::validate() const {
    if (emax <= emin) {
        std::cerr << "Error: EMAX must be greater than EMIN" << std::endl;
        return false;
    }
    if (distance <= 0) {
        std::cerr << "Error: Distance must be positive" << std::endl;
        return false;
    }
    if (age <= 0) {
        std::cerr << "Error: Age must be positive" << std::endl;
        return false;
    }
    if (alpha <= 0) {
        std::cerr << "Error: Spectral index must be positive" << std::endl;
        return false;
    }
    if (ebins <= 0 || rbins <= 0) {
        std::cerr << "Error: Number of bins must be positive" << std::endl;
        return false;
    }
    return true;
}

void EdgeConfig::print() const {
    std::cout << "=== EDGE Configuration ===" << std::endl;
    std::cout << "Source name: " << name << std::endl;
    std::cout << "Profile file: " << profile_file << std::endl;
    std::cout << "Energy range: " << emin << " - " << emax << " TeV" << std::endl;
    std::cout << "Spectral index: " << alpha << std::endl;
    std::cout << "Distance: " << distance << " kpc" << std::endl;
    std::cout << "Age: " << age << " years" << std::endl;
    std::cout << "Diffusion coefficient D0: " << d0 << " cm^2/s" << std::endl;
    std::cout << "Diffusion index delta: " << delta << std::endl;
    std::cout << "Energy density: " << tot_e_dens << " eV/cm^3" << std::endl;
    std::cout << "Magnetic field: " << bfield << " G" << std::endl;
    std::cout << "Use Klein-Nishina: " << (use_kn ? "yes" : "no") << std::endl;
    std::cout << "Energy bins: " << ebins << std::endl;
    std::cout << "Radial bins: " << rbins << std::endl;
    std::cout << "=========================" << std::endl;
}

// EdgeUtils 实现
EdgeUtils::EdgeUtils()
    : initialized_(false) {

    // 初始化 GAMERA 对象
    initGAMERA();
}

void EdgeUtils::initGAMERA() {
    try {
        // GAMERA 对象在构造时自动初始化
        // 设置 GAMERA 对象之间的引用关系
        setupGAMERARelations();

        initialized_ = true;

        std::cout << "GAMERA objects initialized successfully" << std::endl;
        std::cout << "Using GAMERA C++ library for physics calculations" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error initializing GAMERA: " << e.what() << std::endl;
        initialized_ = false;
        throw;
    }
}

void EdgeUtils::setupGAMERARelations() {
    // GAMERA 对象之间的引用关系需要通过设置方法建立
    // 具体方法取决于 GAMERA 的 API
    // 这里暂时留空，如果需要可以查看 GAMERA 文档
}

double EdgeUtils::calculateGamma(double energy) {
    // 计算洛伦兹因子 γ = E/m_e
    return energy / m_e;  // m_e 是 GAMERA 定义的电子质量（erg）
}

double EdgeUtils::calculateThomsonLossRate(double gamma, double e_dens) {
    // 汤姆逊散射能量损失率
    // dE/dt = (4/3) * σ_T * c * γ^2 * energy_density
    double e_dens_erg_cm3 = e_dens * eV_to_erg;  // 转换 eV/cm^3 到 erg/cm^3
    return (4.0 / 3.0) * sigma_T * c_speed * gamma * gamma * e_dens_erg_cm3;
}

} // namespace edge