#ifndef EDGE_SPECTRUM_H
#define EDGE_SPECTRUM_H

#include <vector>
#include <memory>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"
#include "numerics/Interpolation.hh"

namespace edge {

/**
 * @brief 完整的光谱计算类
 *
 * 实现Python EDGE.py中的光谱计算功能，包括：
 * - FillLambdaVector: 填充lambda向量
 * - Spectrum: 主要的光谱计算函数
 * - acc_time: 加速时间计算
 */
class SpectrumCalculator {
public:
    /**
     * @brief 构造函数
     */
    SpectrumCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~SpectrumCalculator() = default;

    /**
     * @brief 计算加速时间
     *
     * 计算粒子加速到给定能量所需的最小时间
     *
     * @param energy 粒子能量（erg）
     * @param b_field 磁场强度（G）
     * @return 加速时间（years）
     */
    double calculateAccelerationTime(double energy, double b_field) const;

    /**
     * @brief 填充Lambda向量
     *
     * 对应Python中的FillLambdaVector函数
     * 计算从给定能量到初始能量的积分扩散系数
     *
     * @param energy_now 当前能量（erg）
     * @param bins 积分 bins 数量
     * @param age 源年龄（years）
     * @return LambdaVector结构
     */
    struct LambdaVectorResult {
        std::vector<double> energies;      // 能量数组（erg）
        std::vector<double> lambdas;       // lambda值数组
        std::vector<double> log_energies;  // 对数能量数组
        std::vector<double> log_lambdas;   // 对数lambda数组
        double DT;                         // 时间（years）
        double E0;                         // 初始能量（erg）
    };

    LambdaVectorResult fillLambdaVector(double energy_now, int bins, double age);

    /**
     * @brief 计算电子光谱
     *
     * 对应Python中的Spectrum函数
     * 计算在给定能量和半径处的电子微分数量
     *
     * @param energy 电子能量（erg）
     * @param lambda_vector lambda向量（来自fillLambdaVector）
     * @param radius 距离源的距离（pc）
     * @param DT 时间（years）
     * @param E0 初始能量（erg）
     * @param age 源年龄（years）
     * @return 电子数密度（1/(erg*cm^3)）
     */
    double calculateSpectrum(double energy,
                            const LambdaVectorResult& lambda_vector,
                            double radius, double DT, double E0, double age);

    /**
     * @brief 批量计算光谱
     *
     * @param energies 能量数组（erg）
     * @param radii 半径数组（pc）
     * @param age 源年龄（years）
     * @return 二维光谱网格 [energy_index][radius_index]
     */
    std::vector<std::vector<double>> calculateSpectrumGrid(
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        double age);

    /**
     * @brief 计算注入光谱归一化
     *
     * @param alpha 光谱指数
     * @param emin 最小能量（erg）
     * @param emax 最大能量（erg）
     * @return 归一化因子
     */
    static double calculateNormalization(double alpha, double emin, double emax);

    /**
     * @brief 设置能量轨迹数据
     *
     * @param trajectory 能量轨迹数据
     */
    void setEnergyTrajectory(const EnergyTrajectoryData& trajectory);

    /**
     * @brief 设置光度数据
     *
     * @param luminosity_data 光度数据
     */
    void setLuminosityData(const LuminosityData& luminosity_data);

    /**
     * @brief 获取能量轨迹数据
     */
    const EnergyTrajectoryData& getEnergyTrajectory() const { return trajectory_data_; }

    /**
     * @brief 获取光度数据
     */
    const LuminosityData& getLuminosityData() const { return luminosity_data_; }

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
    EnergyTrajectoryData trajectory_data_;
    LuminosityData luminosity_data_;

    // 内部辅助函数
    double integrateOverTime(const std::vector<double>& time_points,
                            const std::vector<double>& values) const;

    std::vector<double> calculateLuminosityArray(double t_min, double t_max, int bins) const;
};

/**
 * @brief 数值积分工具类
 *
 * 提供各种数值积分方法，兼容Python的数值积分行为
 */
class NumericalIntegration {
public:
    /**
     * @brief 梯形法则积分
     *
     * @param x_array x坐标数组（单调）
     * @param y_array y值数组
     * @return 积分结果
     */
    static double trapezoidal(const std::vector<double>& x_array,
                             const std::vector<double>& y_array);

    /**
     * @brief 辛普森法则积分
     *
     * @param x_array x坐标数组（单调，等间隔或接近等间隔）
     * @param y_array y值数组
     * @return 积分结果
     */
    static double simpson(const std::vector<double>& x_array,
                         const std::vector<double>& y_array);

    /**
     * @brief 使用GSL进行自适应积分
     *
     * @param func 被积函数
     * @param params 函数参数
     * @param a 积分下限
     * @param b 积分上限
     * @param tolerance 相对误差容忍度
     * @return 积分结果
     */
    static double gslIntegration(double (*func)(double, void*), void* params,
                                double a, double b, double tolerance = 1e-6);

    /**
     * @brief 向量对的积分
     *
     * 对应Python中的fu.Integrate(zip(T, val), t_min, t_max)
     *
     * @param pairs (x, y) 对的数组
     * @param x_min 积分下限
     * @param x_max 积分上限
     * @return 积分结果
     */
    static double integratePairs(const std::vector<std::pair<double, double>>& pairs,
                                double x_min, double x_max);
};

} // namespace edge

#endif // EDGE_SPECTRUM_H