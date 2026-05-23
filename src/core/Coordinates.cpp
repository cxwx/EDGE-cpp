#include "core/Coordinates.hh"
#include "numerics/Interpolation.hh"
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace edge {

// AstroCoordinates 实现

AstroCoordinates::AstroCoordinates(EdgeUtils* utils)
    : utils_(utils), reference_point_(0.0, 0.0, 0.0) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

void AstroCoordinates::setGalacticReferencePoint(const CartesianCoords& xyz_ref) {
    reference_point_ = xyz_ref;

    // 设置GAMERA Astro对象的参考点
    std::vector<double> ref_vec = {xyz_ref.x, xyz_ref.y, xyz_ref.z};
    astro_.SetGalacticReferencePoint(ref_vec);
}

AstroCoordinates::CartesianCoords
AstroCoordinates::galacticToCartesian(const GalacticCoords& galactic,
                                     const CartesianCoords& xyz_obs) {
    // 使用GAMERA的Astro::GetCartesian函数
    std::vector<double> obs_coords = {xyz_obs.x, xyz_obs.y, xyz_obs.z};

    std::vector<double> result = astro_.GetCartesian(
        galactic.distance, galactic.l, galactic.b, obs_coords);

    return CartesianCoords(result[0], result[1], result[2]);
}

AstroCoordinates::GalacticCoords
AstroCoordinates::cartesianToGalactic(const CartesianCoords& cartesian,
                                      const CartesianCoords& xyz_ref) {
    // 使用GAMERA的Astro::GetGalactic函数
    double l, b;
    astro_.GetGalactic(cartesian.x, cartesian.y, cartesian.z,
                      xyz_ref.x, xyz_ref.y, xyz_ref.z, l, b);

    // 计算距离
    double distance = cartesian.norm();

    return GalacticCoords(l, b, distance);
}

std::vector<AstroCoordinates::CartesianCoords>
AstroCoordinates::galacticToCartesianBatch(
    const std::vector<GalacticCoords>& galactic_coords,
    const CartesianCoords& xyz_obs) {

    std::vector<CartesianCoords> result;
    result.reserve(galactic_coords.size());

    for (const auto& galactic : galactic_coords) {
        result.push_back(galacticToCartesian(galactic, xyz_obs));
    }

    return result;
}

std::vector<AstroCoordinates::GalacticCoords>
AstroCoordinates::cartesianToGalacticBatch(
    const std::vector<CartesianCoords>& cartesian_coords,
    const CartesianCoords& xyz_ref) {

    std::vector<GalacticCoords> result;
    result.reserve(cartesian_coords.size());

    for (const auto& cartesian : cartesian_coords) {
        result.push_back(cartesianToGalactic(cartesian, xyz_ref));
    }

    return result;
}

void AstroCoordinates::rotateCoordinates(CartesianCoords& coords,
                                        double phi, double theta, double psi) {
    // 实现简单的坐标旋转（使用欧拉角）
    double x = coords.x;
    double y = coords.y;
    double z = coords.z;

    // 绕z轴旋转psi
    double x1 = x * std::cos(psi) - y * std::sin(psi);
    double y1 = x * std::sin(psi) + y * std::cos(psi);
    double z1 = z;

    // 绕y轴旋转theta
    double x2 = x1 * std::cos(theta) + z1 * std::sin(theta);
    double y2 = y1;
    double z2 = -x1 * std::sin(theta) + z1 * std::cos(theta);

    // 绕x轴旋转phi
    coords.x = x2;
    coords.y = y2 * std::cos(phi) - z2 * std::sin(phi);
    coords.z = y2 * std::sin(phi) + z2 * std::cos(phi);
}

// LineOfSightCalculator 实现

LineOfSightCalculator::LineOfSightCalculator(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils), coords_(utils), radial_bins_(400) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

std::vector<AstroCoordinates::CartesianCoords>
LineOfSightCalculator::generateLineOfSight(
    double l, double b,
    const std::vector<double>& distance_values,
    const AstroCoordinates::CartesianCoords& earth_position) {

    std::vector<AstroCoordinates::CartesianCoords> line_of_sight;
    line_of_sight.reserve(distance_values.size());

    // 创建银道坐标
    AstroCoordinates::GalacticCoords galactic;
    galactic.l = l;
    galactic.b = b;

    // 对每个距离生成坐标点
    for (double distance : distance_values) {
        galactic.distance = distance;

        AstroCoordinates::CartesianCoords coords = coords_.galacticToCartesian(galactic, earth_position);
        line_of_sight.push_back(coords);
    }

    return line_of_sight;
}

std::vector<double> LineOfSightCalculator::lineOfSightIntegration(
    double l, double b,
    const std::vector<std::vector<double>>& density_2d,
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    const AstroCoordinates::CartesianCoords& earth_position) {

    // 生成径向距离数组（kpc）
    auto distance_values = generateRadialValues(config_.distance, radial_bins_);

    // 生成视线点
    auto line_of_sight = generateLineOfSight(l, b, distance_values, earth_position);

    // 对每个能量进行积分
    std::vector<double> results;
    results.reserve(energies.size());

    for (size_t i = 0; i < energies.size(); ++i) {
        std::vector<double> integrand;
        integrand.reserve(distance_values.size());

        // 对每个视线点计算被积函数
        for (const auto& coords : line_of_sight) {
            // 计算到源的距离
            double distance_pc = coords.distanceTo(earth_position) * 1000.0; // kpc to pc

            // 找到对应的半径索引
            size_t r_index = findRadiusIndex(distance_pc, radii);

            // 获取密度值
            double density = density_2d[i][r_index];
            integrand.push_back(density);
        }

        // 进行数值积分
        double integral = 0.0;
        for (size_t j = 1; j < distance_values.size(); ++j) {
            double dr_cm = (distance_values[j] - distance_values[j-1]) * kpc_to_cm;
            double avg_density = 0.5 * (integrand[j] + integrand[j-1]);
            integral += dr_cm * avg_density;
        }

        results.push_back(integral);
    }

    return results;
}

double LineOfSightCalculator::lineOfSightVolumeIntegration(
    double l, double b,
    const std::vector<std::vector<double>>& density_2d,
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    double max_radius) {

    (void)energies; // Suppress unused parameter warning
    // 生成径向距离数组
    auto distance_values = generateRadialValues(config_.distance, radial_bins_);

    // 生成视线点并计算密度
    auto line_of_sight = generateLineOfSight(l, b, distance_values,
                                            AstroCoordinates::CartesianCoords(config_.distance, 0.0, 0.0));

    double total_integral = 0.0;

    for (const auto& coords : line_of_sight) {
        double distance_pc = coords.norm() * 1000.0; // kpc to pc

        if (distance_pc > max_radius) continue;

        // 找到对应的半径索引
        size_t r_index = findRadiusIndex(distance_pc, radii);

        // 计算该点的平均密度（对所有能量）
        double avg_density = 0.0;
        for (const auto& energy_array : density_2d) {
            avg_density += energy_array[r_index];
        }
        avg_density /= density_2d.size();

        // 计算体积元
        double dr = distance_pc * 0.01; // 简化的步长
        double volume = 4.0 * M_PI * distance_pc * distance_pc * dr;

        total_integral += avg_density * volume;
    }

    return total_integral;
}

std::vector<std::vector<double>> LineOfSightCalculator::batchLineOfSightIntegration(
    const std::vector<double>& longitudes,
    const std::vector<double>& latitudes,
    const std::vector<std::vector<double>>& density_2d,
    const std::vector<double>& energies,
    const std::vector<double>& radii,
    const AstroCoordinates::CartesianCoords& earth_position) {

    std::vector<std::vector<double>> results;
    results.reserve(longitudes.size() * latitudes.size());

    // 对每个方向进行积分
    for (double l : longitudes) {
        for (double b : latitudes) {
            auto flux = lineOfSightIntegration(l, b, density_2d, energies, radii, earth_position);
            results.push_back(flux);
        }
    }

    return results;
}

std::vector<double> LineOfSightCalculator::generateRadialValues(double distance_kpc, int bins) {
    // 生成对数间隔的径向距离值（kpc）
    std::vector<double> rvals;

    // 生成从源到地球的精细网格
    auto log_rvals = InterpolationUtils::logspace(-6.0, std::log10(distance_kpc), bins);

    // 创建镜像数组（源的前后）
    std::vector<double> rvals_reversed;
    for (auto val : log_rvals) {
        rvals_reversed.push_back(distance_kpc - val);
    }
    std::reverse(rvals_reversed.begin(), rvals_reversed.end());

    // 合并数组
    rvals.insert(rvals.end(), rvals_reversed.begin(), rvals_reversed.end());
    rvals.insert(rvals.end(), log_rvals.begin(), log_rvals.end());

    return rvals;
}

size_t LineOfSightCalculator::findRadiusIndex(double radius_pc, const std::vector<double>& radii) {
    // 找到第一个大于radius_pc的半径索引
    auto it = std::upper_bound(radii.begin(), radii.end(), radius_pc);

    if (it == radii.begin()) {
        return 0;
    } else if (it == radii.end()) {
        return radii.size() - 1;
    } else {
        return std::distance(radii.begin(), it) - 1;
    }
}

// EarthFluxCalculator 实现

EarthFluxCalculator::EarthFluxCalculator(const EdgeConfig& config, EdgeUtils* utils)
    : config_(config), utils_(utils), los_calculator_(config, utils), coords_(utils) {

    if (!utils_) {
        throw std::invalid_argument("EdgeUtils pointer cannot be null");
    }
}

std::vector<double> EarthFluxCalculator::calculateSourceFlux(
    const std::vector<double>& energies,
    const AstroCoordinates::GalacticCoords& source_position,
    const AstroCoordinates::CartesianCoords& earth_position,
    double age,
    const std::vector<std::vector<double>>& density_2d,
    const std::vector<double>& radii) {

    (void)age; // Suppress unused parameter warning

    // 计算源到地球的方向
    AstroCoordinates::CartesianCoords source_cartesian = coords_.galacticToCartesian(source_position, earth_position);

    // 计算方向向量
    double dx = source_cartesian.x - earth_position.x;
    double dy = source_cartesian.y - earth_position.y;
    double dz = source_cartesian.z - earth_position.z;

    // 转换为银道坐标
    double l_deg = std::atan2(dy, dx) * 180.0 / M_PI;
    double dist_xy = std::sqrt(dx*dx + dy*dy);
    double b_deg = std::atan2(dz, dist_xy) * 180.0 / M_PI;

    // 进行视线积分
    auto flux = los_calculator_.lineOfSightIntegration(
        l_deg, b_deg, density_2d, energies, radii, earth_position);

    return flux;
}

std::vector<double> EarthFluxCalculator::calculateTotalFlux(
    const std::vector<double>& energies,
    const std::vector<AstroCoordinates::GalacticCoords>& pulsar_catalog,
    const AstroCoordinates::CartesianCoords& earth_position) {

    (void)earth_position; // Suppress unused parameter warning
    // 初始化总流量
    std::vector<double> total_flux(energies.size(), 0.0);

    // 这里需要预先计算电子密度分布，暂时留空
    // 实际实现需要调用SpectrumCalculator来获取密度分布

    // 对每个脉冲星计算贡献
    for (const auto& pulsar : pulsar_catalog) {
        (void)pulsar; // Suppress unused variable warning
        // 计算该脉冲星的贡献
        // 这里需要电子密度分布数据
        // 暂时跳过
    }

    return total_flux;
}

} // namespace edge