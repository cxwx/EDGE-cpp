#ifndef EDGE_INTERPOLATION_H
#define EDGE_INTERPOLATION_H

#include <vector>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "utils/EdgeUtils.hh"

namespace edge {

/**
 * @brief 插值工具类
 *
 * 提供一维和二维插值功能，兼容Python/NumPy的插值行为
 */
class InterpolationUtils {
public:
    /**
     * @brief 对数空间的线性插值
     *
     * 相当于Python中的: 10**np.interp(np.log10(x), np.log10(x_array), np.log10(y_array))
     *
     * @param log_x 插值点的对数值
     * @param log_x_array 数据点x的对数值数组（必须单调）
     * @param log_y_array 数据点y的对数值数组
     * @return 插值结果
     */
    static double log10Interp(double log_x, const std::vector<double>& log_x_array,
                             const std::vector<double>& log_y_array);

    /**
     * @brief 常规线性插值
     *
     * 相当于Python中的: np.interp(x, x_array, y_array)
     *
     * @param x 插值点
     * @param x_array 数据点x数组（必须单调）
     * @param y_array 数据点y数组
     * @return 插值结果
     */
    static double interp(double x, const std::vector<double>& x_array,
                        const std::vector<double>& y_array);

    /**
     * @brief 创建对数间隔的数组
     *
     * 相当于Python中的: np.logspace(log10_min, log10_max, num_bins)
     *
     * @param log10_min 最小值的对数
     * @param log10_max 最大值的对数
     * @param num_bins 点的数量
     * @return 对数间隔的数组
     */
    static std::vector<double> logspace(double log10_min, double log10_max, int num_bins);

    /**
     * @brief 创建线性间隔的数组
     *
     * 相当于Python中的: np.linspace(min, max, num_bins)
     *
     * @param min 最小值
     * @param max 最大值
     * @param num_bins 点的数量
     * @return 线性间隔的数组
     */
    static std::vector<double> linspace(double min, double max, int num_bins);

    /**
     * @brief 使用GSL进行样条插值
     *
     * 创建样条插值对象，可以多次使用
     */
    class SplineInterpolator {
    public:
        SplineInterpolator(const std::vector<double>& x, const std::vector<double>& y,
                          const gsl_interp_type* interp_type = gsl_interp_linear);
        ~SplineInterpolator();

        double eval(double x) const;
        double eval_deriv(double x) const;
        double eval_deriv2(double x) const;
        double eval_integ(double a, double b) const;

    private:
        gsl_spline* spline_;
        gsl_interp_accel* accel_;
        std::vector<double> x_, y_;
    };

    /**
     * @brief 二维网格插值
     */
    template<typename T>
    static T bilinearInterp(double x, double y,
                           const std::vector<double>& x_array,
                           const std::vector<double>& y_array,
                           const std::vector<std::vector<T>>& data_array);

    /**
     * @brief 在数组中查找最近邻的索引
     */
    static size_t findNearestIndex(double value, const std::vector<double>& array);

    /**
     * @brief 数组排序并返回索引
     */
    static std::vector<size_t> argsort(const std::vector<double>& array);
};

/**
 * @brief 能量轨迹数据结构
 *
 * 存储能量-时间轨迹数据，用于插值
 */
struct EnergyTrajectoryData {
    std::vector<double> time;        // 时间数组（年）
    std::vector<double> energy;      // 能量数组（erg）
    std::vector<double> lambda;      // 积分扩散系数数组
    std::vector<double> log_time;    // 时间的对数（用于插值）
    std::vector<double> log_energy;  // 能量的对数（用于插值）
    std::vector<double> log_lambda;  // lambda的对数（用于插值）

    // 时间-能量和能量-时间对数数组（用于反向插值）
    std::vector<double> traj_cont_log_time;      // ETRAJCONT: 时间->能量
    std::vector<double> traj_cont_log_energy;    // ETRAJCONT: 时间->能量
    std::vector<double> traj_inv_log_energy;     // ETRAJCONTINVERSE: 能量->时间
    std::vector<double> traj_inv_log_time;       // ETRAJCONTINVERSE: 能量->时间
    std::vector<double> lambda_cont_log_energy;  // LAMBCONT: 能量->lambda
    std::vector<double> lambda_cont_log_lambda;  // LAMBCONT: 能量->lambda

    /**
     * @brief 从原始数据准备插值数组
     */
    void prepareInterpolationArrays();

    /**
     * @brief 清空所有数据
     */
    void clear();

    /**
     * @brief 获取数据点数量
     */
    size_t size() const { return time.size(); }

    /**
     * @brief 检查数据是否有效
     */
    bool isValid() const { return !time.empty() && time.size() == energy.size(); }
};

/**
 * @brief 光度演化数据结构
 *
 * 存储时间相关的光度数据
 */
struct LuminosityData {
    std::vector<double> time;           // 时间数组（年）
    std::vector<double> luminosity;     // 光度数组（erg/s）
    std::vector<double> log_time;       // 时间的对数
    std::vector<double> log_luminosity; // 光度的对数

    /**
     * @brief 准备插值数组
     */
    void prepareInterpolationArrays();

    /**
     * @brief 清空数据
     */
    void clear();

    /**
     * @brief 检查数据是否有效
     */
    bool isValid() const { return !time.empty() && time.size() == luminosity.size(); }
};

} // namespace edge

#endif // EDGE_INTERPOLATION_H