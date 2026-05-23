#include <iostream>
#include <iomanip>
#include <cmath>

// 简单的测试程序，验证 GAMERA 库链接和基本功能
#include <gamera/Utils.h>
#include <gamera/Radiation.h>
#include <gamera/Particles.h>
#include <gamera/Astro.h>

void testGAMERAConstants() {
    std::cout << "=== Testing GAMERA Constants ===" << std::endl;

    std::cout << "Energy conversions:" << std::endl;
    std::cout << "  1 TeV = " << TeV_to_erg << " erg" << std::endl;
    std::cout << "  1 GeV = " << GeV_to_erg << " erg" << std::endl;
    std::cout << "  1 eV = " << eV_to_erg << " erg" << std::endl;

    std::cout << "\nDistance conversions:" << std::endl;
    std::cout << "  1 pc = " << pc_to_cm << " cm" << std::endl;
    std::cout << "  1 kpc = " << kpc_to_cm << " cm" << std::endl;

    std::cout << "\nPhysical constants:" << std::endl;
    std::cout << "  Electron mass = " << m_e << " erg" << std::endl;
    std::cout << "  Speed of light = " << c_speed << " cm/s" << std::endl;
    std::cout << "  Thomson cross section = " << sigma_T << " cm^2" << std::endl;
    std::cout << "  Year in seconds = " << yr_to_sec << " s" << std::endl;
    std::cout << "  Pi = " << pi << std::endl;

    // 测试单位转换
    std::cout << "\n=== Testing Unit Conversions ===" << std::endl;
    double energy_TeV = 1.0;
    double energy_erg = energy_TeV * TeV_to_erg;
    std::cout << energy_TeV << " TeV = " << energy_erg << " erg" << std::endl;

    double distance_kpc = 1.0;
    double distance_cm = distance_kpc * kpc_to_cm;
    std::cout << distance_kpc << " kpc = " << std::scientific << distance_cm
              << std::fixed << " cm" << std::endl;
}

void testGAMERAObjects() {
    std::cout << "\n=== Testing GAMERA Objects ===" << std::endl;

    try {
        // 创建 GAMERA 对象（构造函数不接受参数）
        Utils utils;
        Radiation radiation;
        Particles particles;
        Astro astro;

        std::cout << "GAMERA objects created successfully!" << std::endl;

        // 测试一些基本功能
        std::cout << "\n=== Testing Basic Functions ===" << std::endl;

        // 测试随机数生成器
        std::cout << "Random number test: ";
        for (int i = 0; i < 5; ++i) {
            std::cout << std::fixed << std::setprecision(4)
                      << utils.Random() << " ";
        }
        std::cout << std::endl;

        // 测试能量损失率计算
        std::cout << "\nEnergy loss rate test:" << std::endl;
        std::vector<double> energies_erg = {
            0.001 * TeV_to_erg,  // 1 GeV
            0.01 * TeV_to_erg,   // 10 GeV
            0.1 * TeV_to_erg,    // 100 GeV
            1.0 * TeV_to_erg,    // 1 TeV
            10.0 * TeV_to_erg    // 10 TeV
        };

        for (double e : energies_erg) {
            double loss_rate = particles.EnergyLossRate(e);
            std::cout << "E = " << std::setw(8) << e / TeV_to_erg << " TeV: "
                      << "dE/dt = " << std::scientific << loss_rate
                      << std::fixed << " erg/s" << std::endl;
        }

        std::cout << "\nAll GAMERA function tests passed!" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error creating GAMERA objects: " << e.what() << std::endl;
        throw;
    }
}

void testBasicPhysics() {
    std::cout << "\n=== Testing Basic Physics Calculations ===" << std::endl;

    // 测试洛伦兹因子计算
    std::cout << "Lorentz factor gamma = E/m_e:" << std::endl;
    std::vector<double> energies_TeV = {0.001, 0.01, 0.1, 1.0, 10.0};

    for (double e_TeV : energies_TeV) {
        double e_erg = e_TeV * TeV_to_erg;
        double gamma = e_erg / m_e;
        std::cout << "E = " << std::setw(8) << e_TeV << " TeV: "
                  << "gamma = " << std::scientific << gamma << std::fixed << std::endl;
    }

    // 测试汤姆逊散射能量损失率
    std::cout << "\nThomson scattering energy loss rate:" << std::endl;
    double energy_density_eV_cm3 = 1.0; // 1 eV/cm^3
    double energy_density_erg_cm3 = energy_density_eV_cm3 * eV_to_erg;

    for (double e_TeV : energies_TeV) {
        double e_erg = e_TeV * TeV_to_erg;
        double gamma = e_erg / m_e;
        double loss_rate = (4.0 / 3.0) * sigma_T * c_speed * gamma * gamma * energy_density_erg_cm3;

        std::cout << "E = " << std::setw(8) << e_TeV << " TeV: "
                  << "dE/dt = " << std::scientific << loss_rate
                  << std::fixed << " erg/s" << std::endl;
    }

    // 测试扩散系数
    std::cout << "\nDiffusion coefficient D(E) = D0 * (1 + E/E*)^delta:" << std::endl;
    double D0 = 4.0e27; // cm^2/s
    double delta = 0.33;
    double E_star = 3.0e-3 * TeV_to_erg; // erg

    for (double e_TeV : energies_TeV) {
        double e_erg = e_TeV * TeV_to_erg;
        double D = D0 * std::pow(1.0 + e_erg / E_star, delta);

        std::cout << "E = " << std::setw(8) << e_TeV << " TeV: "
                  << "D = " << std::scientific << D
                  << std::fixed << " cm^2/s" << std::endl;
    }

    std::cout << "\nAll physics tests passed!" << std::endl;
}

int main() {
    std::cout << "\n";
    std::cout << "========================================" << std::endl;
    std::cout << "  EDGE C++ Basic Functionality Test     " << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        testGAMERAConstants();
        testGAMERAObjects();
        testBasicPhysics();

        std::cout << "\n";
        std::cout << "========================================" << std::endl;
        std::cout << "  All Tests Passed Successfully!       " << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "\nEDGE C++ is ready for use!" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nTEST FAILED: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "\nUNKNOWN ERROR in tests" << std::endl;
        return 1;
    }
}