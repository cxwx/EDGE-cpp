#ifndef EDGE_GAMMA_H
#define EDGE_GAMMA_H

#include <vector>
#include <memory>
#include "utils/EdgeUtils.hh"
#include "core/Diffusion.hh"

namespace edge {

/**
 * @brief 伽马射线计算类
 *
 * 计算伽马射线光谱、轮廓和地球位置的流量
 */
class GammaRayCalculator {
public:
    explicit GammaRayCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 计算伽马射线光谱
     *
     * @param electron_spectrum 电子谱密度
     * @param energies 能量数组（TeV）
     * @return 伽马射线光谱（1/(TeV*cm^2*s)）
     */
    std::vector<double> calculateGammaSpectrum(
        const std::vector<double>& energies,
        double radius,
        double age);

    /**
     * @brief 计算伽马射线表面亮度轮廓
     */
    struct GammaProfile {
        std::vector<double> radii;      // 半径（度）
        std::vector<double> surface_brightness; // 表面亮度
    };

    GammaProfile calculateGammaProfile(
        const std::vector<double>& radii_deg,
        double energy_TeV,
        double age);

    /**
     * @brief 计算总伽马射线通量
     */
    double calculateTotalGammaFlux(double energy_min_TeV, double energy_max_TeV);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
};

/**
 * @brief 视线积分计算
 */
class LineOfSightIntegration {
public:
    explicit LineOfSightIntegration(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 计算沿视线的积分
     *
     * @param l 银经（度）
     * @param b 银纬（度）
     * @param energy 能量（erg）
     * @param age 年龄（years）
     * @return 积分值（1/(erg*cm^2*s)）
     */
    double integrateLineOfSight(double l, double b, double energy, double age);

    /**
     * @brief 批量计算多个方向的积分
     */
    std::vector<std::vector<double>> integrateAllSightlines(
        const std::vector<double>& longitudes,
        const std::vector<double>& latitudes,
        double energy,
        double age);

    /**
     * @brief 计算体积积分
     */
    double integrateVolume(const std::vector<std::vector<double>>& density_map,
                          double max_radius_pc);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;

    // 内部坐标转换
    struct CartesianCoords {
        double x, y, z;
    };
    CartesianCoords sphericalToCartesian(double l, double b, double r) const;
};

/**
 * @brief 地球位置电子流量计算
 */
class EarthFluxCalculator {
public:
    explicit EarthFluxCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 计算地球位置的电子流量
     *
     * @param energies 能量数组（TeV）
     * @param source_params 源参数（距离、年龄等）
     * @return 电子流量（1/(TeV*cm^2*s)）
     */
    std::vector<double> calculateEarthFlux(
        const std::vector<double>& energies_TeV,
        const std::vector<double>& distances_kpc,
        const std::vector<double>& ages_years);

    /**
     * @brief 计算所有脉冲星对地球的总贡献
     */
    struct PulsarContribution {
        std::string name;
        double l;                      // 银经（度）
        double b;                      // 银纬（度）
        double distance;               // 距离
        double age;                    // 年龄
        std::vector<double> flux;      // 流量谱
    };

    std::vector<PulsarContribution> calculateAllPulsarContributions(
        const std::vector<double>& energies_TeV,
        const std::string& pulsar_catalog = "default");

    /**
     * @brief 计算总电子流量
     */
    std::vector<double> calculateTotalFlux(
        const std::vector<PulsarContribution>& contributions);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
    LineOfSightIntegration los_integrator_;

    // 均匀分布参数
    struct GalacticPulsarDistribution {
        double sn_rate;      // 超新星爆发率
        double scale_height; // 标高
        double scale_length; // 标长
    };

    GalacticPulsarDistribution getGalacticDistribution();
};

/**
 * @brief 伽马射线源光谱计算
 */
class GammaSourceSpectrum {
public:
    explicit GammaSourceSpectrum(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 计算源的伽马射线光谱
     *
     * 包括逆康普顿散射和韧致辐射
     */
    std::vector<double> calculateSourceSpectrum(
        const std::vector<double>& energies_TeV);

    /**
     * @brief 计算逆康普顿散射光谱
     */
    std::vector<double> calculateICSpectrum(
        const std::vector<double>& energies_TeV,
        const std::vector<double>& electron_spectrum);

    /**
     * @brief 计算韧致辐射光谱
     */
    std::vector<double> calculateBremsstrahlungSpectrum(
        const std::vector<double>& energies_TeV,
        const std::vector<double>& electron_spectrum);

private:
    EdgeConfig config_;
    EdgeUtils* utils_;
};

} // namespace edge

#endif // EDGE_GAMMA_H