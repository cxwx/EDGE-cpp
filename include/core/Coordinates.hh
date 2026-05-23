#ifndef EDGE_COORDINATES_H
#define EDGE_COORDINATES_H

#include <vector>
#include <memory>
#include "utils/EdgeUtils.hh"

namespace edge {

/**
 * @brief 天文坐标转换类
 *
 * 提供银道坐标系统和笛卡尔坐标系统之间的转换
 * 基于GAMERA的Astro类功能
 */
class AstroCoordinates {
public:
    /**
     * @brief 银道坐标结构
     */
    struct GalacticCoords {
        double l;  // 银经（度）
        double b;  // 银纬（度）
        double distance;  // 距离（kpc）

        GalacticCoords() : l(0), b(0), distance(0) {}
        GalacticCoords(double longitude, double latitude, double dist)
            : l(longitude), b(latitude), distance(dist) {}
    };

    /**
     * @brief 笛卡尔坐标结构
     */
    struct CartesianCoords {
        double x, y, z;  // 笛卡尔坐标（kpc）

        CartesianCoords() : x(0), y(0), z(0) {}
        CartesianCoords(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

        // 计算到原点的距离
        double norm() const { return std::sqrt(x*x + y*y + z*z); }

        // 计算到另一点的距离
        double distanceTo(const CartesianCoords& other) const {
            double dx = x - other.x;
            double dy = y - other.y;
            double dz = z - other.z;
            return std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    };

    /**
     * @brief 构造函数
     */
    AstroCoordinates(EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~AstroCoordinates() = default;

    /**
     * @brief 设置银道参考点
     *
     * @param xyz_ref 参考点笛卡尔坐标（kpc）
     */
    void setGalacticReferencePoint(const CartesianCoords& xyz_ref);

    /**
     * @brief 银道坐标转笛卡尔坐标
     *
     * 相当于GAMERA的Astro::GetCartesian函数
     *
     * @param galactic 银道坐标
     * @param xyz_obs 观测者位置（kpc）
     * @return 笛卡尔坐标
     */
    CartesianCoords galacticToCartesian(const GalacticCoords& galactic,
                                       const CartesianCoords& xyz_obs);

    /**
     * @brief 笛卡尔坐标转银道坐标
     *
     * 相当于GAMERA的Astro::GetGalactic函数
     *
     * @param cartesian 笛卡尔坐标
     * @param xyz_ref 参考点坐标
     * @return 银道坐标
     */
    GalacticCoords cartesianToGalactic(const CartesianCoords& cartesian,
                                      const CartesianCoords& xyz_ref);

    /**
     * @brief 批量转换银道坐标到笛卡尔坐标
     */
    std::vector<CartesianCoords> galacticToCartesianBatch(
        const std::vector<GalacticCoords>& galactic_coords,
        const CartesianCoords& xyz_obs);

    /**
     * @brief 批量转换笛卡尔坐标到银道坐标
     */
    std::vector<GalacticCoords> cartesianToGalacticBatch(
        const std::vector<CartesianCoords>& cartesian_coords,
        const CartesianCoords& xyz_ref);

    /**
     * @brief 旋转坐标系
     *
     * 绕原点旋转坐标
     *
     * @param coords 要旋转的坐标
     * @param phi 绕x轴旋转角度（弧度）
     * @param theta 绕y轴旋转角度（弧度）
     * @param psi 绕z轴旋转角度（弧度）
     */
    static void rotateCoordinates(CartesianCoords& coords,
                                 double phi, double theta, double psi);

    /**
     * @brief 度数转弧度
     */
    static double degreesToRadians(double degrees) {
        return degrees * M_PI / 180.0;
    }

    /**
     * @brief 弧度转度数
     */
    static double radiansToDegrees(double radians) {
        return radians * 180.0 / M_PI;
    }

    /**
     * @brief 角度归一化到 [0, 360)
     */
    static double normalizeAngle(double angle_degrees) {
        double normalized = fmod(angle_degrees, 360.0);
        if (normalized < 0) normalized += 360.0;
        return normalized;
    }

    /**
     * @brief 获取GAMERA Astro对象（直接调用）
     */
    Astro* getGameraAstro() { return &astro_; }

private:
    EdgeUtils* utils_;
    Astro astro_;
    CartesianCoords reference_point_;
};

/**
 * @brief 视线积分计算类
 *
 * 实现沿视线的积分计算，用于计算地球位置的流量
 */
class LineOfSightCalculator {
public:
    /**
     * @brief 构造函数
     */
    LineOfSightCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~LineOfSightCalculator() = default;

    /**
     * @brief 生成视线点
     *
     * 沿给定银道坐标方向生成从地球出发的视线点
     *
     * @param l 银经（度）
     * @param b 银纬（度）
     * @param distance_values 距离数组（kpc）
     * @param earth_position 地球位置（kpc）
     * @return 视线点数组
     */
    std::vector<AstroCoordinates::CartesianCoords> generateLineOfSight(
        double l, double b,
        const std::vector<double>& distance_values,
        const AstroCoordinates::CartesianCoords& earth_position);

    /**
     * @brief 视线积分
     *
     * 对给定的二维密度数组进行视线积分
     * 相当于Python中的LineOfSightIntegration函数
     *
     * @param l 银经（度）
     * @param b 银纬（度）
     * @param density_2d 二维密度数组 [energy][radius]，单位 1/(erg*cm^3)
     * @param energies 能量数组（erg）
     * @param radii 半径数组（pc）
     * @param earth_position 地球位置（kpc）
     * @return 积分结果数组，每个能量一个值，单位 1/(erg*cm^2)
     */
    std::vector<double> lineOfSightIntegration(
        double l, double b,
        const std::vector<std::vector<double>>& density_2d,
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        const AstroCoordinates::CartesianCoords& earth_position);

    /**
     * @brief 体积积分
     *
     * 对给定的二维密度数组进行体积积分
     * 相当于Python中的LineOfSightVolumeIntegration函数
     *
     * @param l 银经（度）
     * @param b 银纬（度）
     * @param density_2d 二维密度数组
     * @param energies 能量数组（erg）
     * @param radii 半径数组（pc）
     * @param max_radius 最大积分半径（pc）
     * @return 体积积分结果
     */
    double lineOfSightVolumeIntegration(
        double l, double b,
        const std::vector<std::vector<double>>& density_2d,
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        double max_radius);

    /**
     * @brief 批量视线积分
     *
     * 对多个方向进行视线积分
     */
    std::vector<std::vector<double>> batchLineOfSightIntegration(
        const std::vector<double>& longitudes,
        const std::vector<double>& latitudes,
        const std::vector<std::vector<double>>& density_2d,
        const std::vector<double>& energies,
        const std::vector<double>& radii,
        const AstroCoordinates::CartesianCoords& earth_position);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
    AstroCoordinates coords_;
    int radial_bins_;

    // 内部辅助函数
    std::vector<double> generateRadialValues(double distance_kpc, int bins);
    size_t findRadiusIndex(double radius_pc, const std::vector<double>& radii);
};

/**
 * @brief 地球流量计算类
 *
 * 计算地球位置的电子/伽马射线流量
 */
class EarthFluxCalculator {
public:
    /**
     * @brief 构造函数
     */
    EarthFluxCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~EarthFluxCalculator() = default;

    /**
     * @brief 计算单个源的地球流量
     *
     * @param energies 能量数组（TeV）
     * @param source_position 源位置（银道坐标）
     * @param earth_position 地球位置（笛卡尔坐标，kpc）
     * @param age 源年龄（years）
     * @param density_2d 电子密度分布
     * @param radii 半径数组（pc）
     * @return 流量数组（1/(TeV*cm^2*s)）
     */
    std::vector<double> calculateSourceFlux(
        const std::vector<double>& energies,
        const AstroCoordinates::GalacticCoords& source_position,
        const AstroCoordinates::CartesianCoords& earth_position,
        double age,
        const std::vector<std::vector<double>>& density_2d,
        const std::vector<double>& radii);

    /**
     * @brief 计算所有脉冲星的总流量
     *
     * 相当于Python中的Flux_Earth_all_pulsars函数
     *
     * @param energies 能量数组（TeV）
     * @param pulsar_catalog 脉冲星目录
     * @param earth_position 地球位置
     * @return 总流量数组
     */
    std::vector<double> calculateTotalFlux(
        const std::vector<double>& energies,
        const std::vector<AstroCoordinates::GalacticCoords>& pulsar_catalog,
        const AstroCoordinates::CartesianCoords& earth_position);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
    LineOfSightCalculator los_calculator_;
    AstroCoordinates coords_;  // 添加坐标转换对象
};

} // namespace edge

#endif // EDGE_COORDINATES_H