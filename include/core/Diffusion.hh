#ifndef EDGE_DIFFUSION_H
#define EDGE_DIFFUSION_H

#include <vector>
#include <memory>
#include "utils/EdgeUtils.hh"

namespace edge {

/**
 * @brief 电子扩散计算类
 *
 * 这个类实现了电子从中央源向外扩散的物理计算，
 * 包括能量损失、扩散系数计算和粒子轨迹模拟。
 */
class DiffusionCalculator {
public:
    /**
     * @brief 构造函数
     * @param config EDGE 配置参数
     * @param utils EDGE 工具类指针
     */
    explicit DiffusionCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~DiffusionCalculator() = default;

    /**
     * @brief 计算扩散系数
     * @param energy 电子能量（erg）
     * @return 扩散系数（cm^2/s）
     */
    double calculateDiffusionCoefficient(double energy) const;

    /**
     * @brief 计算能量损失率
     * @param energy 电子能量（erg）
     * @return 能量损失率（erg/s）
     */
    double calculateEnergyLossRate(double energy) const;

    /**
     * @brief 计算电子能量轨迹
     *
     * 模拟电子从最大能量到最小能量的演化过程，考虑能量损失
     * 返回时间-能量关系和积分扩散系数
     */
    struct EnergyTrajectory {
        std::vector<double> time;        // 时间（years）
        std::vector<double> energy;      // 能量（erg）
        std::vector<double> lambda;      // 积分扩散系数
        std::vector<double> diff_coeff;  // 扩散系数
    };

    EnergyTrajectory calculateEnergyTrajectory();

    /**
     * @brief 创建二维电子能量-半径网格
     * @return 能量数组（erg）和半径数组（pc）
     */
    std::pair<std::vector<double>, std::vector<double>> createEnergyRadiusGrid();

    /**
     * @brief 计算电子密度分布
     *
     * 对于给定能量和半径，计算电子的数密度分布
     */
    double calculateElectronDensity(double energy, double radius, double age);

    /**
     * @brief 批量计算电子密度分布
     */
    using ElectronDensityGrid = std::vector<std::vector<double>>;
    ElectronDensityGrid calculateElectronDensityGrid(
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        double age);

    // 访问器
    const EnergyTrajectory& getTrajectory() const { return trajectory_; }
    const EdgeConfig& getConfig() const { return config_; }

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
    EnergyTrajectory trajectory_;

    // 内部辅助函数
    double calculateE_star() const;
    double integrateLossFunction(double energy);
};

/**
 * @brief 电子光谱计算类
 */
class ElectronSpectrum {
public:
    explicit ElectronSpectrum(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 计算电子数密度谱
     *
     * @param energy 电子能量（erg）
     * @param radius 距离源的距离（pc）
     * @param age 源年龄（years）
     * @return 电子数密度（1/(erg*cm^3)）
     */
    double calculateSpectrum(double energy, double radius, double age);

    /**
     * @brief 批量计算光谱
     */
    std::vector<std::vector<double>> calculateSpectrumGrid(
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        double age);

    /**
     * @brief 计算注入光谱
     */
    double calculateInjectionSpectrum(double energy) const;

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
};

/**
 * @brief 脉冲星光度和时间演化
 */
class PulsarEvolution {
public:
    explicit PulsarEvolution(const EdgeConfig& config);

    /**
     * @brief 计算自旋向下光度
     * @param time 时间（years）
     * @return 光度（erg/s）
     */
    double calculateLuminosity(double time) const;

    /**
     * @brief 计算总能量输出
     * @param age 年龄（years）
     * @return 总能量（erg）
     */
    double calculateTotalEnergy(double age) const;

    /**
     * @brief 从初始周期计算特征年龄
     */
    double calculateCharacteristicAge() const;

private:
    EdgeConfig config_;
    double tau_;  // 特征时间尺度（years）
};

} // namespace edge

#endif // EDGE_DIFFUSION_H