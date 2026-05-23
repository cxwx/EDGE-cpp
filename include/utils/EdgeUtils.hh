#ifndef EDGE_UTILS_H
#define EDGE_UTILS_H

#include <string>
#include <vector>
#include <gamera/Utils.h>
#include <gamera/Radiation.h>
#include <gamera/Particles.h>
#include <gamera/Astro.h>

namespace edge {

// EDGE 全局常量（从 Python 版本迁移）
struct EdgeConstants {
    // 能量常量（默认值，可被命令行参数覆盖）
    static constexpr double EMAX_DEFAULT = 500.0;          // TeV
    static constexpr double EMIN_DEFAULT = 0.001;         // TeV
    static constexpr double E_STAR_DEFAULT = 3.0e-3;      // TeV（特征能量）

    // 扩散参数默认值
    static constexpr double D0_DEFAULT = 4.0e27;          // cm^2/s
    static constexpr double DELTA_DEFAULT = 0.33;         // 扩散指数
    static constexpr double DIST_DEFAULT = 0.25;          // kpc

    // 光谱参数默认值
    static constexpr double ALPHA_DEFAULT = 2.2;          // 光谱指数
    static constexpr double MU_DEFAULT = 0.5;             // 电子能量占比

    // 时间参数默认值
    static constexpr double AGE_DEFAULT = 3.42e5;         // years
    static constexpr double T0_DEFAULT = 1.2e4;           // years（初始自旋衰减时标）

    // 能量密度参数默认值
    static constexpr double TOT_E_DENS_DEFAULT = 1.06;    // eV/cm^3
    static constexpr double BCONT_DEFAULT = 3.0;          // G（磁场）

    // 脉冲星参数默认值
    static constexpr double EDOT_DEFAULT = 3.2e34;        // erg/s（自旋向下功率）
    static constexpr double BRIND_DEFAULT = 3.0;         // 衰减指数
    static constexpr double P_DEFAULT = 237.0;            // ms（脉冲星周期）
    static constexpr double P0_DEFAULT = 40.5;            // ms（初始周期）

    // 数值计算参数
    static constexpr double ENERGY_LOSS_FACTOR = 1.0e-3;  // 能量损失时间步长因子
};

// 命令行参数结构
struct EdgeConfig {
    // 源参数
    std::string name = "Source";
    std::string profile_file = "Data/GemingaProfile.dat";

    // 能量谱参数
    double alpha = EdgeConstants::ALPHA_DEFAULT;         // 光谱指数
    double emax = EdgeConstants::EMAX_DEFAULT;           // TeV
    double emin = EdgeConstants::EMIN_DEFAULT;           // TeV

    // 扩散参数
    double distance = EdgeConstants::DIST_DEFAULT;       // kpc
    double delta = EdgeConstants::DELTA_DEFAULT;         // 扩散指数
    double d0 = EdgeConstants::D0_DEFAULT;               // cm^2/s
    double size = 5.0;                                   // 源尺寸

    // 时间参数
    double age = EdgeConstants::AGE_DEFAULT;             // years
    double tau0 = EdgeConstants::T0_DEFAULT;             // years

    // 物理参数
    double mu = EdgeConstants::MU_DEFAULT;               // 电子能量占比
    double tot_e_dens = EdgeConstants::TOT_E_DENS_DEFAULT; // eV/cm^3
    double bfield = EdgeConstants::BCONT_DEFAULT;        // G
    bool use_kn = false;                                 // 使用 Klein-Nishina 近似

    // 脉冲星参数
    double edot = EdgeConstants::EDOT_DEFAULT;           // erg/s
    double brind = EdgeConstants::BRIND_DEFAULT;         // 衰减指数
    double period = EdgeConstants::P_DEFAULT;            // ms
    double p0 = EdgeConstants::P0_DEFAULT;               // ms

    // 数值参数
    int ebins = 100;                                     // 能量 bins
    int rbins = 400;                                     // 径向 bins

    // 运行时选项
    bool birth_period = false;
    bool all_pulsar = false;
    bool only_flux_earth = false;
    bool fig_eps = false;
    double timesupr = 0.0;        // 光度抑制时间 [years]

    // 验证配置有效性
    bool validate() const;
    void print() const;
};

// EDGE 工具类
class EdgeUtils {
public:
    // 构造函数
    EdgeUtils();

    // 初始化 GAMERA 对象
    void initGAMERA();

    // 获取 GAMERA 对象指针
    Utils* getUtils() { return &utils_; }
    Radiation* getRadiation() { return &radiation_; }
    Particles* getParticles() { return &particles_; }
    Astro* getAstro() { return &astro_; }

    // 单位转换（使用 GAMERA 定义的常量）
    static double teVToErg(double tev) { return tev * TeV_to_erg; }
    static double ergToTeV(double erg) { return erg / TeV_to_erg; }
    static double pcToCm(double pc) { return pc * pc_to_cm; }
    static double kpcToCm(double kpc) { return kpc * kpc_to_cm; }
    static double yrToSec(double yr) { return yr * yr_to_sec; }
    static double secToYr(double sec) { return sec / yr_to_sec; }

    // 物理计算辅助函数
    static double calculateGamma(double energy);          // 洛伦兹因子
    static double calculateThomsonLossRate(double gamma, double e_dens);

private:
    Utils utils_;
    Radiation radiation_;
    Particles particles_;
    Astro astro_;
    bool initialized_;

    // 设置 GAMERA 对象之间的引用关系
    void setupGAMERARelations();
};

} // namespace edge

#endif // EDGE_UTILS_H