#ifndef EDGE_MULTI_SOURCE_H
#define EDGE_MULTI_SOURCE_H

#include <vector>
#include <memory>
#include <string>
#include "utils/EdgeUtils.hh"
#include "core/Coordinates.hh"
#include "core/Spectrum.hh"

namespace edge {

/**
 * @brief 脉冲星数据结构
 *
 * 对应Python版本中的脉冲星目录数据
 */
struct PulsarData {
    std::string name;
    double age;                    // years
    double period;                 // ms
    double p0;                     // initial period (ms)
    double distance;               // kpc
    double l;                      // 银经（度）
    double b;                      // 银纬（度）
    double edot;                   // erg/s
    double alpha;                  // 光谱指数
    double characteristic_age;     // years

    // 笛卡尔坐标
    double x, y, z;                // kpc (银心坐标系)

    // 预计算的流量谱（可选）
    std::vector<double> flux_spectrum; // 1/(TeV*cm^2*s)
    bool has_flux_spectrum = false;

    // 计算到地球的距离
    double distanceToEarth(double earth_x = 8.5, double earth_y = 0.0, double earth_z = 0.0) const {
        double dx = x - earth_x;
        double dy = y - earth_y;
        double dz = z - earth_z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    // 获取银道坐标
    AstroCoordinates::GalacticCoords getGalacticCoords() const {
        return AstroCoordinates::GalacticCoords(l, b, distance);
    }
};

/**
 * @brief 脉冲星目录类
 *
 * 管理多个脉冲星的数据，对应Python的脉冲星分布功能
 */
class PulsarCatalog {
public:
    /**
     * @brief 构造函数
     */
    PulsarCatalog(EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~PulsarCatalog() = default;

    /**
     * @brief 从文件加载脉冲星目录
     *
     * @param filename 文件路径
     * @param format 文件格式 ('text', 'fits', 'csv')
     */
    void loadFromFile(const std::string& filename, const std::string& format = "text");

    /**
     * @brief 生成均匀分布的脉冲星
     *
     * 对应Python中的Homogeneus_distribution_pulsars函数
     *
     * @param age 最大年龄（years）
     * @param sn_rate 超新星爆发率（1/years）
     * @param spiralarms 是否考虑旋臂结构
     */
    void generateHomogeneousDistribution(double age, double sn_rate, bool spiralarms = true);

    /**
     * @brief 添加单个脉冲星
     */
    void addPulsar(const PulsarData& pulsar);

    /**
     * @brief 获取所有脉冲星
     */
    const std::vector<PulsarData>& getPulsars() const { return pulsars_; }

    /**
     * @brief 获取所有脉冲星（非const版本，用于修改）
     */
    std::vector<PulsarData>& getPulsars() { return pulsars_; }

    /**
     * @brief 脉冲星数量
     */
    size_t size() const { return pulsars_.size(); }

    /**
     * @brief 清空目录
     */
    void clear() { pulsars_.clear(); }

    /**
     * @brief 按距离过滤
     */
    std::vector<PulsarData> filterByDistance(double min_distance, double max_distance) const;

    /**
     * @brief 按年龄过滤
     */
    std::vector<PulsarData> filterByAge(double min_age, double max_age) const;

    /**
     * @brief 查找附近的脉冲星
     */
    std::vector<PulsarData> findNearbyPulsars(double max_distance, double earth_x = 8.5) const;

    /**
     * @brief 打印目录信息
     */
    void printInfo() const;

private:
    EdgeUtils* utils_;
    std::vector<PulsarData> pulsars_;

    // 内部辅助函数
    void calculateCartesianCoordinates();
    void generateLorimer06Distribution(double n_sources);
};

/**
 * @brief 多源流量计算器
 *
 * 对应Python中的Flux_Earth_all_pulsars函数
 */
class MultiSourceFluxCalculator {
public:
    /**
     * @brief 计算参数
     */
    struct CalculationParams {
        double Q0;                     // 常数注入率 [1/(erg*s)]
        double spectral_index;         // 稳态谱指数
        double min_distance_kpc;       // 最小距离（kpc）
        double max_age;                // 最大年龄（years）
        double sn_rate;                // 超新星爆发率 [1/years]

        // 地球位置
        double earth_x;                // kpc
        double earth_y;                // kpc
        double earth_z;                // kpc

        // 构造函数提供默认值
        CalculationParams()
            : Q0(5.0e32), spectral_index(-2.4), min_distance_kpc(1.0),
              max_age(1.0e7), sn_rate(2.0 / 100.0),
              earth_x(8.5), earth_y(0.0), earth_z(0.0) {}
    };

    /**
     * @brief 构造函数
     */
    MultiSourceFluxCalculator(const EdgeConfig& config, EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~MultiSourceFluxCalculator() = default;

    /**
     * @brief 计算所有脉冲星的总流量
     *
     * 对应Python中的Flux_Earth_all_pulsars函数
     *
     * @param energies 能量数组（TeV）
     * @param catalog 脉冲星目录
     * @param params 计算参数
     * @return 总流量数组 [1/(TeV*cm^2*s)]
     */
    std::vector<double> calculateTotalFlux(
        const std::vector<double>& energies,
        const PulsarCatalog& catalog,
        const CalculationParams& params = CalculationParams());

    /**
     * @brief 计算单个脉冲星的贡献
     *
     * @param energy 能量（erg）
     * @param pulsar 脉冲星数据
     * @param params 计算参数
     * @return 流量 [1/(erg*cm^2*s)]
     */
    double calculatePulsarContribution(
        double energy,
        const PulsarData& pulsar,
        const CalculationParams& params);

    /**
     * @brief 批量计算所有脉冲星的光谱
     */
    void calculateAllPulsarSpectra(
        const std::vector<double>& energies,
        PulsarCatalog& catalog,
        const CalculationParams& params);

    /**
     * @brief 分析脉冲星分布
     */
    struct DistributionAnalysis {
        int total_pulsars;
        int pulsars_within_1kpc;
        double average_distance;
        double total_flux;
        std::vector<PulsarData> nearby_pulsars;
    };

    DistributionAnalysis analyzeDistribution(
        const PulsarCatalog& catalog,
        double analysis_radius = 1.0); // kpc

private:
    EdgeConfig config_;
    EdgeUtils* utils_;

    // 稳态流量计算（Atoyan et al. 1995, Eq. 21）
    double calculateSteadyStateFlux(
        double energy_erg,
        double distance_cm,
        double pulsar_age,
        double D,
        double Q0,
        double spectral_index) const;

    // 计算冷却时间
    double calculateCoolingTime(double energy_erg, double loss_rate) const;
};

/**
 * @brief 数据文件管理类
 *
 * 处理数据输入输出，包括FITS文件
 */
class DataManager {
public:
    /**
     * @brief 构造函数
     */
    DataManager(EdgeUtils* utils);

    /**
     * @brief 析构函数
     */
    ~DataManager() = default;

    /**
     * @brief 读取轮廓数据文件
     *
     * 对应Python中的Data/GemingaProfile.dat文件
     */
    std::vector<std::pair<double, double>> readProfileData(const std::string& filename);

    /**
     * @brief 读取数据点文件
     */
    std::vector<std::vector<double>> readDataPoints(const std::string& filename);

    /**
     * @brief 写入结果到文件
     */
    void writeResults(const std::string& filename,
                     const std::vector<double>& energies,
                     const std::vector<double>& fluxes,
                     const std::string& format = "text");

    /**
     * @brief 写入网格数据
     */
    void writeGridData(const std::string& filename,
                      const std::vector<double>& x_coords,
                      const std::vector<double>& y_coords,
                      const std::vector<std::vector<double>>& data);

    /**
     * @brief 读取FITS格式的数据
     */
    std::vector<std::vector<double>> readFITS(const std::string& filename);

    /**
     * @brief 写入FITS格式的数据
     */
    void writeFITS(const std::string& filename,
                   const std::vector<std::string>& column_names,
                   const std::vector<std::vector<double>>& data);

private:
    EdgeUtils* utils_;

    // 内部辅助函数
    std::vector<std::string> splitString(const std::string& str, char delimiter) const;
    bool isComment(const std::string& line) const;
};

/**
 * @brief 完整的EDGE计算类
 *
 * 整合所有功能，提供完整的EDGE计算能力
 */
class EdgeCalculator {
public:
    /**
     * @brief 构造函数
     */
    EdgeCalculator(const EdgeConfig& config);

    /**
     * @brief 析构函数
     */
    ~EdgeCalculator() = default;

    /**
     * @brief 运行完整的EDGE计算
     *
     * @param energies 能量数组（TeV）
     * @return 计算结果
     */
    struct CalculationResult {
        // 能量网格
        std::vector<double> energies;          // TeV
        std::vector<double> radii;             // pc

        // 电子密度分布
        std::vector<std::vector<double>> electron_density; // [energy][radius]

        // 地球流量
        std::vector<double> earth_flux;        // 1/(TeV*cm^2*s)

        // 伽马射线光谱
        std::vector<double> gamma_spectrum;    // 1/(TeV*cm^2*s)

        // 脉冲星贡献
        std::vector<PulsarData> pulsars;
        std::vector<std::vector<double>> pulsar_fluxes; // [pulsar][energy]

        // 元数据
        double calculation_time;              // seconds
        std::vector<std::string> log_messages;
    };

    CalculationResult runCompleteCalculation(
        const std::vector<double>& energies);

    /**
     * @brief 设置脉冲星目录
     */
    void setPulsarCatalog(const PulsarCatalog& catalog);

    /**
     * @brief 获取配置
     */
    const EdgeConfig& getConfig() const { return config_; }

    /**
     * @brief 保存结果到文件
     */
    void saveResults(const CalculationResult& result, const std::string& output_dir);

private:
    EdgeConfig config_;
    EdgeUtils utils_;
    DataManager data_manager_;

    // 计算组件
    std::unique_ptr<DiffusionCalculator> diffusion_calc_;
    std::unique_ptr<SpectrumCalculator> spectrum_calc_;
    std::unique_ptr<MultiSourceFluxCalculator> multi_source_calc_;

    // 脉冲星目录（使用指针避免初始化问题）
    std::unique_ptr<PulsarCatalog> pulsar_catalog_;

    // 内部计算步骤
    void calculateElectronDensity(CalculationResult& result);
    void calculateEarthFlux(CalculationResult& result);
    void calculateGammaSpectrum(CalculationResult& result);
    void calculateMultiSourceContributions(CalculationResult& result);
};

} // namespace edge

#endif // EDGE_MULTI_SOURCE_H