#include "core/MultiSource.hh"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <iostream>

#ifdef HAS_FITS
#include <CCfits/CCfits>
#endif

namespace edge {

// PulsarCatalog 实现

PulsarCatalog::PulsarCatalog(EdgeUtils* utils)
    : utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

void PulsarCatalog::loadFromFile(const std::string& filename, const std::string& format) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    if (format == "text") {
        std::string line;
        while (std::getline(file, line)) {
            // 跳过注释和空行
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            PulsarData pulsar;

            // 解析数据：name age period p0 distance l b edot alpha
            if (iss >> pulsar.name >> pulsar.age >> pulsar.period >> pulsar.p0
                   >> pulsar.distance >> pulsar.l >> pulsar.b >> pulsar.edot >> pulsar.alpha) {

                // 计算特征年龄
                if (pulsar.period > 0 && pulsar.p0 > 0) {
                    pulsar.characteristic_age = (pulsar.period * 1e-3) /
                        (2 * M_PI * (pulsar.period - pulsar.p0) * 1e-3);
                    pulsar.characteristic_age *= 1.0 / (365.25 * 24 * 3600); // seconds to years
                } else {
                    pulsar.characteristic_age = pulsar.age;
                }

                pulsars_.push_back(pulsar);
            }
        }
    } else {
        throw std::invalid_argument("Unsupported file format: " + format);
    }

    calculateCartesianCoordinates();
    std::cout << "Loaded " << pulsars_.size() << " pulsars from " << filename << std::endl;
}

void PulsarCatalog::generateHomogeneousDistribution(double age, double sn_rate, bool spiralarms) {
    // 生成均匀分布的脉冲星
    double n_sources = age * sn_rate;

    std::cout << "Generating homogeneous pulsar distribution..." << std::endl;
    std::cout << "Number of sources: " << n_sources << std::endl;

    if (spiralarms) {
        generateLorimer06Distribution(n_sources);
    } else {
        // 简单的均匀分布
        pulsars_.reserve(static_cast<size_t>(n_sources));

        // 使用GAMERA的随机数生成器
        auto* gamera_utils = utils_->getUtils();

        for (size_t i = 0; i < static_cast<size_t>(n_sources); ++i) {
            PulsarData pulsar;

            // 随机位置（银盘）
            double r = 3.0 + 12.0 * gamera_utils->Random(); // 3-15 kpc
            double theta = 2 * M_PI * gamera_utils->Random();
            double z = (gamera_utils->Random() - 0.5) * 0.5; // ±0.25 kpc scale height

            pulsar.x = r * std::cos(theta);
            pulsar.y = r * std::sin(theta);
            pulsar.z = z;

            // 转换为银道坐标
            double d = std::sqrt(pulsar.x*pulsar.x + pulsar.y*pulsar.y + pulsar.z*pulsar.z);
            pulsar.distance = d;

            double l_rad = std::atan2(pulsar.y, pulsar.x);
            double b_rad = std::asin(pulsar.z / d);

            pulsar.l = l_rad * 180.0 / M_PI;
            pulsar.b = b_rad * 180.0 / M_PI;

            // 随机年龄（对数分布）
            double log_age = std::log10(1e3) + (std::log10(age) - std::log10(1e3)) * gamera_utils->Random();
            pulsar.age = std::pow(10.0, log_age);
            pulsar.characteristic_age = pulsar.age;

            // 典型参数
            pulsar.period = 0.1 + 0.5 * gamera_utils->Random(); // 0.1-0.6 s
            pulsar.p0 = pulsar.period * 0.9;
            pulsar.edot = 1e34 + 1e35 * gamera_utils->Random();
            pulsar.alpha = 2.0 + 0.5 * gamera_utils->Random();
            pulsar.name = "PSR_" + std::to_string(i);

            pulsars_.push_back(pulsar);
        }
    }

    std::cout << "Generated " << pulsars_.size() << " pulsars" << std::endl;
}

void PulsarCatalog::addPulsar(const PulsarData& pulsar) {
    pulsars_.push_back(pulsar);
    calculateCartesianCoordinates();
}

std::vector<PulsarData> PulsarCatalog::filterByDistance(double min_distance, double max_distance) const {
    std::vector<PulsarData> filtered;

    for (const auto& pulsar : pulsars_) {
        double distance = pulsar.distanceToEarth();
        if (distance >= min_distance && distance <= max_distance) {
            filtered.push_back(pulsar);
        }
    }

    return filtered;
}

std::vector<PulsarData> PulsarCatalog::filterByAge(double min_age, double max_age) const {
    std::vector<PulsarData> filtered;

    for (const auto& pulsar : pulsars_) {
        if (pulsar.age >= min_age && pulsar.age <= max_age) {
            filtered.push_back(pulsar);
        }
    }

    return filtered;
}

std::vector<PulsarData> PulsarCatalog::findNearbyPulsars(double max_distance, double earth_x) const {
    std::vector<PulsarData> nearby;

    for (const auto& pulsar : pulsars_) {
        double distance = pulsar.distanceToEarth(earth_x);
        if (distance <= max_distance) {
            nearby.push_back(pulsar);
        }
    }

    // 按距离排序
    std::sort(nearby.begin(), nearby.end(),
        [earth_x](const PulsarData& a, const PulsarData& b) {
            return a.distanceToEarth(earth_x) < b.distanceToEarth(earth_x);
        });

    return nearby;
}

void PulsarCatalog::printInfo() const {
    std::cout << "\n=== Pulsar Catalog Info ===" << std::endl;
    std::cout << "Total pulsars: " << pulsars_.size() << std::endl;

    if (!pulsars_.empty()) {
        // 统计信息
        double avg_distance = 0.0;
        double min_distance = 1e10;
        double max_distance = 0.0;

        for (const auto& pulsar : pulsars_) {
            double distance = pulsar.distanceToEarth();
            avg_distance += distance;
            min_distance = std::min(min_distance, distance);
            max_distance = std::max(max_distance, distance);
        }

        avg_distance /= pulsars_.size();

        std::cout << "Distance to Earth: min=" << min_distance << " kpc, "
                  << "avg=" << avg_distance << " kpc, "
                  << "max=" << max_distance << " kpc" << std::endl;

        // 年龄统计
        double avg_age = 0.0;
        for (const auto& pulsar : pulsars_) {
            avg_age += pulsar.age;
        }
        avg_age /= pulsars_.size();

        std::cout << "Average age: " << avg_age << " years" << std::endl;

        // 附近的脉冲星
        auto nearby = findNearbyPulsars(1.0);
        std::cout << "Pulsars within 1 kpc: " << nearby.size() << std::endl;

        if (!nearby.empty() && nearby.size() <= 10) {
            std::cout << "Nearby pulsars:" << std::endl;
            std::cout << "Distance [kpc] | x_Earth [kpc] | y_Earth [kpc] | Name" << std::endl;
            std::cout << "---------------|--------------|--------------|------" << std::endl;

            for (const auto& pulsar : nearby) {
                double distance = pulsar.distanceToEarth();
                std::cout << std::fixed << std::setprecision(3)
                          << std::setw(13) << distance << " | "
                          << std::setw(12) << (pulsar.x - 8.5) << " | "
                          << std::setw(12) << pulsar.y << " | "
                          << pulsar.name << std::endl;
            }
        }
    }

    std::cout << "=============================" << std::endl;
}

void PulsarCatalog::calculateCartesianCoordinates() {
    // 计算所有脉冲星的笛卡尔坐标
    for (auto& pulsar : pulsars_) {
        double l_rad = pulsar.l * M_PI / 180.0;
        double b_rad = pulsar.b * M_PI / 180.0;

        // 转换为笛卡尔坐标（相对于银心）
        double d = pulsar.distance;
        pulsar.x = d * std::cos(b_rad) * std::cos(l_rad);
        pulsar.y = d * std::cos(b_rad) * std::sin(l_rad);
        pulsar.z = d * std::sin(b_rad);
    }
}

void PulsarCatalog::generateLorimer06Distribution(double n_sources) {
    // 生成Lorimer et al. 2006的银河系脉冲星分布
    pulsars_.clear();
    pulsars_.reserve(static_cast<size_t>(n_sources));

    auto* gamera_utils = utils_->getUtils();

    // Lorimer 06分布参数
    double radial_scale = 4.5;  // kpc
    double vertical_scale = 0.33; // kpc

    for (size_t i = 0; i < static_cast<size_t>(n_sources); ++i) {
        PulsarData pulsar;

        // 径向分布（指数衰减）
        double u = gamera_utils->Random();
        double r = -radial_scale * std::log(1.0 - u * 0.95); // 限制在合理范围内

        // 方位角（可能考虑旋臂结构）
        double theta = 2 * M_PI * gamera_utils->Random();

        // 垂直分布（指数衰减）
        double v = gamera_utils->Random() - 0.5;
        double z = -vertical_scale * std::log(1.0 - 2.0 * std::abs(v));
        if (v < 0) z = -z;

        // 转换为笛卡尔坐标
        pulsar.x = r * std::cos(theta);
        pulsar.y = r * std::sin(theta);
        pulsar.z = z;

        // 银道坐标
        pulsar.distance = std::sqrt(r*r + z*z);
        double l_rad = std::atan2(pulsar.y, pulsar.x);
        double b_rad = std::asin(pulsar.z / pulsar.distance);

        pulsar.l = l_rad * 180.0 / M_PI;
        pulsar.b = b_rad * 180.0 / M_PI;

        // 年龄分布（对数均匀）
        double log_age = 3.0 + 4.0 * gamera_utils->Random(); // 1e3 to 1e7 years
        pulsar.age = std::pow(10.0, log_age);
        pulsar.characteristic_age = pulsar.age;

        // 周期分布（简化版：使用指数分布）
        double log_p_min = -1.0;  // 0.1 s
        double log_p_max = 0.0;   // 1.0 s
        double log_p = log_p_min + (log_p_max - log_p_min) * gamera_utils->Random();
        pulsar.period = std::pow(10.0, log_p); // seconds
        pulsar.p0 = pulsar.period * 0.95;

        // 自旋向下功率
        pulsar.edot = 1e34 * std::pow(10.0, 34.0 * (gamera_utils->Random() - 1.0));

        // 光谱指数
        pulsar.alpha = 2.0 + 0.4 * gamera_utils->Random();

        pulsar.name = "PSR_L06_" + std::to_string(i);

        pulsars_.push_back(pulsar);
    }
}

// MultiSourceFluxCalculator 实现

MultiSourceFluxCalculator::MultiSourceFluxCalculator(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

std::vector<double> MultiSourceFluxCalculator::calculateTotalFlux(
    const std::vector<double>& energies,
    const PulsarCatalog& catalog,
    const CalculationParams& params) {

    std::cout << "\n=== Calculating Total Flux from All Pulsars ===" << std::endl;
    std::cout << "Number of pulsars: " << catalog.size() << std::endl;
    std::cout << "Energy bins: " << energies.size() << std::endl;

    std::vector<double> total_flux(energies.size(), 0.0);

    // 对每个能量计算总流量 - 使用预计算的流量谱
    std::vector<double> flux_sum(energies.size(), 0.0);
    int total_contributing = 0;

    for (const auto& pulsar : catalog.getPulsars()) {
        if (!pulsar.has_flux_spectrum) continue;

        double distance_cm = pulsar.distanceToEarth(params.earth_x, params.earth_y, params.earth_z) * kpc_to_cm;

        // 只考虑距离大于min_distance的脉冲星
        if (distance_cm < params.min_distance_kpc * kpc_to_cm) continue;

        // 累加流量谱
        for (size_t i = 0; i < energies.size(); ++i) {
            flux_sum[i] += pulsar.flux_spectrum[i];
        }
        total_contributing++;
    }

    std::cout << "Contributing pulsars: " << total_contributing << std::endl;

    // 输出结果
    for (size_t i = 0; i < energies.size(); ++i) {
        total_flux[i] = flux_sum[i];
    }

    return total_flux;
}

double MultiSourceFluxCalculator::calculatePulsarContribution(
    double energy_erg,
    const PulsarData& pulsar,
    const CalculationParams& params) {

    // 计算扩散系数
    double D = config_.d0 * std::pow(1.0 + energy_erg / (3e-3 * TeV_to_erg), config_.delta);

    // 计算距离
    double distance_cm = pulsar.distanceToEarth(params.earth_x, params.earth_y, params.earth_z) * kpc_to_cm;

    // 计算冷却时间
    double gamma = energy_erg / m_e;
    double loss_rate = (4.0 / 3.0) * sigma_T * c_speed * gamma * gamma * 1.0; // assume 1 eV/cm^3
    double t_cool = energy_erg / loss_rate / yr_to_sec;

    // 考虑年龄限制
    double t_effective = std::min(pulsar.age, t_cool);

    // 使用稳态公式 (Atoyan et al. 1995, Eq. 21)
    double flux = params.Q0 * std::pow(energy_erg, params.spectral_index) /
                  (4.0 * M_PI * D * distance_cm) *
                  std::exp(-distance_cm * distance_cm / (4.0 * D * t_effective * yr_to_sec));

    return flux;
}

void MultiSourceFluxCalculator::calculateAllPulsarSpectra(
    const std::vector<double>& energies,
    PulsarCatalog& catalog,
    const CalculationParams& params) {

    std::cout << "\n=== Calculating All Pulsar Spectra ===" << std::endl;

    for (auto& pulsar : catalog.getPulsars()) {  // 注意：移除了const
        pulsar.flux_spectrum.clear();
        pulsar.flux_spectrum.reserve(energies.size());

        for (double energy_TeV : energies) {
            double energy_erg = energy_TeV * TeV_to_erg;
            double flux = calculatePulsarContribution(energy_erg, pulsar, params);
            pulsar.flux_spectrum.push_back(flux);
        }

        pulsar.has_flux_spectrum = true;
    }

    std::cout << "Calculated spectra for " << catalog.size() << " pulsars" << std::endl;
}

MultiSourceFluxCalculator::DistributionAnalysis
MultiSourceFluxCalculator::analyzeDistribution(
    const PulsarCatalog& catalog,
    double analysis_radius) {

    DistributionAnalysis analysis;

    analysis.total_pulsars = catalog.size();
    analysis.nearby_pulsars = catalog.findNearbyPulsars(analysis_radius);
    analysis.pulsars_within_1kpc = analysis.nearby_pulsars.size();

    if (!analysis.nearby_pulsars.empty()) {
        double total_distance = 0.0;
        for (const auto& pulsar : analysis.nearby_pulsars) {
            total_distance += pulsar.distanceToEarth();
        }
        analysis.average_distance = total_distance / analysis.nearby_pulsars.size();
    }

    return analysis;
}

double MultiSourceFluxCalculator::calculateSteadyStateFlux(
    double energy_erg,
    double distance_cm,
    double pulsar_age,
    double D,
    double Q0,
    double spectral_index) const {

    // 稳态解：Atoyan et al. 1995, Eq. 21
    double flux = Q0 * std::pow(energy_erg, spectral_index) / (4.0 * M_PI * D * distance_cm);

    // 扩散项
    double diffusion_term = std::exp(-distance_cm * distance_cm / (4.0 * D * pulsar_age * yr_to_sec));

    return flux * diffusion_term;
}

double MultiSourceFluxCalculator::calculateCoolingTime(double energy_erg, double loss_rate) const {
    return energy_erg / loss_rate / yr_to_sec; // years
}

// DataManager 实现

DataManager::DataManager(EdgeUtils* utils)
    : utils_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

std::vector<std::pair<double, double>> DataManager::readProfileData(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<double, double>> data;

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open profile file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (isComment(line)) continue;

        std::istringstream iss(line);
        double angle, flux;

        if (iss >> angle >> flux) {
            data.push_back({angle, flux});
        }
    }

    std::cout << "Read " << data.size() << " data points from " << filename << std::endl;
    return data;
}

std::vector<std::vector<double>> DataManager::readDataPoints(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open data file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (isComment(line)) continue;

        std::vector<std::string> parts = splitString(line, ' ');
        std::vector<double> row;

        for (const auto& part : parts) {
            if (!part.empty()) {
                try {
                    row.push_back(std::stod(part));
                } catch (const std::exception&) {
                    // 忽略无效数据
                }
            }
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    std::cout << "Read " << data.size() << " rows from " << filename << std::endl;
    return data;
}

void DataManager::writeResults(const std::string& filename,
                              const std::vector<double>& energies,
                              const std::vector<double>& fluxes,
                              const std::string& format) {

    if (format == "text") {
        std::ofstream file(filename);

        if (!file.is_open()) {
            throw std::runtime_error("Cannot create output file: " + filename);
        }

        file << "# Energy [TeV] Flux [1/(TeV*cm^2*s)]\n";

        for (size_t i = 0; i < energies.size(); ++i) {
            file << std::scientific << std::setprecision(6)
                 << energies[i] << " " << fluxes[i] << "\n";
        }

        std::cout << "Results written to " << filename << std::endl;
    } else {
        throw std::invalid_argument("Unsupported output format: " + format);
    }
}

void DataManager::writeGridData(const std::string& filename,
                               const std::vector<double>& x_coords,
                               const std::vector<double>& y_coords,
                               const std::vector<std::vector<double>>& data) {

    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot create grid file: " + filename);
    }

    // 写入网格数据（适合Python读取）
    file << "# Grid data: x[" << x_coords.size() << "], y[" << y_coords.size() << "]\n";

    // 写入x坐标
    file << "# x coordinates:\n";
    for (double x : x_coords) {
        file << x << " ";
    }
    file << "\n";

    // 写入y坐标
    file << "# y coordinates:\n";
    for (double y : y_coords) {
        file << y << " ";
    }
    file << "\n";

    // 写入数据
    file << "# data [y_index][x_index]:\n";
    for (const auto& row : data) {
        for (double value : row) {
            file << std::scientific << std::setprecision(6) << value << " ";
        }
        file << "\n";
    }

    std::cout << "Grid data written to " << filename << std::endl;
}

std::vector<std::vector<double>> DataManager::readFITS(const std::string& filename) {
    (void)filename; // Suppress unused parameter warning
#ifdef HAS_FITS
    try {
        CCfits::FITS file(filename, CCfits::Read, true);

        // 读取第一个HDU
        CCfits::PHDU& hdu = file.pHDU();

        // 获取数据
        std::valarray<double> contents;
        hdu.read(contents);

        // 重新整形为2D数组
        int axis1 = hdu.axis(0);
        int axis2 = (hdu.axes() > 1) ? hdu.axis(1) : 1;

        std::vector<std::vector<double>> data(axis2, std::vector<double>(axis1));

        for (int i = 0; i < axis2; ++i) {
            for (int j = 0; j < axis1; ++j) {
                data[i][j] = contents[i * axis1 + j];
            }
        }

        std::cout << "Read FITS file: " << axis2 << " x " << axis1 << " data" << std::endl;
        return data;

    } catch (CCfits::FitsException& e) {
        throw std::runtime_error("FITS error: " + std::string(e.message()));
    }
#else
    throw std::runtime_error("FITS support not enabled. Recompile with HAS_FITS=1.");
#endif
}

void DataManager::writeFITS(const std::string& filename,
                            const std::vector<std::string>& column_names,
                            const std::vector<std::vector<double>>& data) {
    (void)filename; // Suppress unused parameter warning
    (void)column_names; // Suppress unused parameter warning
    (void)data; // Suppress unused parameter warning
#ifdef HAS_FITS
    try {
        // 创建FITS文件
        long naxis = 2;
        long naxes[2] = {static_cast<long>(data[0].size()), static_cast<long>(data.size())};

        CCfits::FITS file(filename, USHORT_IMG, naxis, naxes);

        // 写入数据
        CCfits::Table* table = file.addTable("RESULTS", data.size(), column_names);

        for (size_t col = 0; col < data.size(); ++col) {
            std::vector<double> column_data = data[col];

            // 写入列数据
            table->column(col).write(column_data, 1, column_data.size());
        }

        file.close();

        std::cout << "FITS file written to " << filename << std::endl;

    } catch (CCfits::FitsException& e) {
        throw std::runtime_error("FITS error: " + std::string(e.message()));
    }
#else
    throw std::runtime_error("FITS support not enabled. Recompile with HAS_FITS=1.");
#endif
}

std::vector<std::string> DataManager::splitString(const std::string& str, char delimiter) const {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }

    return tokens;
}

bool DataManager::isComment(const std::string& line) const {
    // 跳过空行和注释行
    if (line.empty()) return true;

    size_t first_non_space = line.find_first_not_of(" \t");
    if (first_non_space == std::string::npos) return true;

    return line[first_non_space] == '#';
}

// EdgeCalculator 实现

EdgeCalculator::EdgeCalculator(const EdgeConfig& config)
    : config_(config), data_manager_(&utils_) {

    // 初始化计算组件
    diffusion_calc_ = std::make_unique<DiffusionCalculator>(config_, &utils_);
    spectrum_calc_ = std::make_unique<SpectrumCalculator>(config_, &utils_);
    multi_source_calc_ = std::make_unique<MultiSourceFluxCalculator>(config_, &utils_);
}

EdgeCalculator::CalculationResult EdgeCalculator::runCompleteCalculation(
    const std::vector<double>& energies) {

    std::cout << "\n";
    std::cout << "========================================================\n";
    std::cout << "  EDGE C++ Complete Calculation                        \n";
    std::cout << "========================================================\n";

    CalculationResult result;
    result.energies = energies;

    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // 1. 计算电子密度分布
        calculateElectronDensity(result);

        // 2. 计算地球流量
        calculateEarthFlux(result);

        // 3. 计算伽马射线光谱
        calculateGammaSpectrum(result);

        // 4. 计算多源贡献
        calculateMultiSourceContributions(result);

    } catch (const std::exception& e) {
        result.log_messages.push_back("ERROR: " + std::string(e.what()));
        throw;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    result.calculation_time = std::chrono::duration<double>(end_time - start_time).count();

    std::cout << "\n=== Calculation Summary ===" << std::endl;
    std::cout << "Calculation time: " << result.calculation_time << " seconds" << std::endl;
    std::cout << "Energy bins: " << result.energies.size() << std::endl;
    std::cout << "Radial bins: " << result.radii.size() << std::endl;
    std::cout << "Pulsars processed: " << result.pulsars.size() << std::endl;

    std::cout << "========================================================\n";

    return result;
}

void EdgeCalculator::setPulsarCatalog(const PulsarCatalog& catalog) {
    // 创建一个新的PulsarCatalog对象并复制数据
    pulsar_catalog_ = std::make_unique<PulsarCatalog>(&utils_);

    // 复制所有脉冲星数据
    for (const auto& pulsar : catalog.getPulsars()) {
        pulsar_catalog_->addPulsar(pulsar);
    }
}

void EdgeCalculator::saveResults(const CalculationResult& result, const std::string& output_dir) {
    // 创建输出目录
    std::string command = "mkdir -p " + output_dir;
    system(command.c_str());

    // 保存结果
    data_manager_.writeResults(output_dir + "/earth_flux.txt", result.energies, result.earth_flux);

    if (!result.electron_density.empty() && !result.radii.empty()) {
        data_manager_.writeGridData(output_dir + "/electron_density.txt",
                                   result.energies, result.radii, result.electron_density);
    }

    std::cout << "Results saved to " << output_dir << std::endl;
}

void EdgeCalculator::calculateElectronDensity(CalculationResult& result) {
    std::cout << "\n=== Step 1: Calculating Electron Density ===" << std::endl;

    // 创建网格
    auto [energies_erg, radii] = diffusion_calc_->createEnergyRadiusGrid();
    result.radii = radii;

    // 计算能量轨迹
    auto trajectory = diffusion_calc_->calculateEnergyTrajectory();

    // 准备插值数据
    EnergyTrajectoryData traj_data;
    traj_data.time = trajectory.time;
    traj_data.energy = trajectory.energy;
    traj_data.lambda = trajectory.lambda;
    traj_data.prepareInterpolationArrays();

    spectrum_calc_->setEnergyTrajectory(traj_data);

    // 计算电子密度
    result.electron_density = spectrum_calc_->calculateSpectrumGrid(energies_erg, radii, config_.age);

    std::cout << "Electron density calculation completed" << std::endl;
}

void EdgeCalculator::calculateEarthFlux(CalculationResult& result) {
    std::cout << "\n=== Step 2: Calculating Earth Flux ===" << std::endl;

    // 这里使用简化模型计算地球流量
    result.earth_flux.resize(result.energies.size(), 0.0);

    for (size_t i = 0; i < result.energies.size(); ++i) {
        double energy_TeV = result.energies[i];

        // 简化的流量估计
        result.earth_flux[i] = 1e-30 * std::pow(energy_TeV, -3.0); // 幂律谱
    }

    std::cout << "Earth flux calculation completed" << std::endl;
}

void EdgeCalculator::calculateGammaSpectrum(CalculationResult& result) {
    std::cout << "\n=== Step 3: Calculating Gamma Ray Spectrum ===" << std::endl;

    // 初始化伽马射线光谱
    result.gamma_spectrum.resize(result.energies.size(), 0.0);

    // 这里可以使用GAMERA的Radiation类来计算逆康普顿散射等过程
    // 暂时使用简化模型

    for (size_t i = 0; i < result.energies.size(); ++i) {
        double energy_TeV = result.energies[i];

        // 简化的伽马射线光谱
        result.gamma_spectrum[i] = 1e-32 * std::pow(energy_TeV, -2.7);
    }

    std::cout << "Gamma ray spectrum calculation completed" << std::endl;
}

void EdgeCalculator::calculateMultiSourceContributions(CalculationResult& result) {
    std::cout << "\n=== Step 4: Calculating Multi-Source Contributions ===" << std::endl;

    if (!pulsar_catalog_ || pulsar_catalog_->size() == 0) {
        std::cout << "No pulsar catalog provided, skipping multi-source calculation" << std::endl;
        return;
    }

    // 计算所有脉冲星的贡献
    MultiSourceFluxCalculator::CalculationParams params;
    pulsar_catalog_->printInfo();

    // 计算总流量（所有脉冲星的总和）
    auto total_flux = multi_source_calc_->calculateTotalFlux(result.energies, *pulsar_catalog_, params);

    // 将总流量添加到地球流量中（可选）
    if (result.earth_flux.empty()) {
        result.earth_flux = total_flux;
    } else {
        // 如果已有地球流量，则累加
        for (size_t i = 0; i < result.earth_flux.size() && i < total_flux.size(); ++i) {
            result.earth_flux[i] += total_flux[i];
        }
    }

    // 复制脉冲星数据
    for (const auto& pulsar : pulsar_catalog_->getPulsars()) {
        result.pulsars.push_back(pulsar);
    }

    std::cout << "Multi-source contribution calculation completed" << std::endl;
}

} // namespace edge