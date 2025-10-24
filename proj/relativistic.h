#ifndef	_RELATIVISTIC_H_
#define	_RELATIVISTIC_H_

#include "scene_env.h"
#include "def.h"

// 相対論的収差を適用したレイ方向の変換
Vector3d applyRelativisticAberration(const Vector3d& ray_dir, const Vector3d& velocity, double beta, double gamma) {

    if (beta < 1e-6) return ray_dir; // 速度が小さい場合は変換不要

    Vector3d v_norm = normalize(velocity);
    double cos_theta = dot(ray_dir,v_norm);

    // 相対論的収差の公式
    // cos(θ') = (cos(θ) - β) / (1 - β·cos(θ))
    double denominator = 1.0 - beta * cos_theta;
    if (fabs(denominator) < 1.0e-12) denominator = 1.0e-12;

    double cos_theta_prime = (cos_theta - beta) / denominator;

    // sin(θ') = sin(θ) / (γ(1 - β·cos(θ)))
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    double sin_theta_prime = sin_theta / (gamma * denominator);

    // レイ方向を速度方向成分と垂直成分に分解
    Vector3d parallel_component = v_norm * cos_theta;
    Vector3d perpendicular_component = ray_dir - parallel_component;
    double perp_length = perpendicular_component.length();

    // 垂直成分の方向を保持
    Vector3d perp_dir = (perp_length > 1e-10) ?
        perpendicular_component / perp_length :
        Vector3d(0, 0, 0);

    // 新しい方向を再構築
    Vector3d new_dir = v_norm * cos_theta_prime + perp_dir * sin_theta_prime;

    return normalize(new_dir);
}

// ドップラー因子を計算
double calculateDopplerFactor(const Vector3d& ray_dir, const Vector3d& velocity, double beta) {
    if (beta < 1e-6) return 1.0;

    Vector3d v_norm = normalize(velocity);
    double cos_theta = dot(ray_dir,v_norm);

    // ドップラー因子: sqrt((1-β)/(1+β)) / (1 - β·cos(θ))
    // または簡略化して: 1 / (γ(1 - β·cos(θ)))
    double denominator = 1.0 - beta * cos_theta;
    return 1.0 / denominator; // 相対論的ドップラー効果
}

// 可視光の波長範囲（ナノメートル）
const double WAVELENGTH_RED = 700.0;
const double WAVELENGTH_GREEN = 546.1;
const double WAVELENGTH_BLUE = 435.8;

// 波長からRGB強度を計算（簡易版）
// 波長が可視光範囲外の場合は減衰
double wavelengthToIntensity(double wavelength, double center, double width = 100.0) {
    double diff = wavelength - center;
    double intensity = std::exp(-(diff * diff) / (2.0 * width * width));

    // 可視光範囲外の減衰
    if (wavelength < 380.0 || wavelength > 750.0) {
        double fade = 1.0;
        if (wavelength < 380.0) {
            fade = std::exp(-(380.0 - wavelength) / 50.0);
        }
        else {
            fade = std::exp(-(wavelength - 750.0) / 50.0);
        }
        intensity *= fade;
    }

    return intensity;
}

// ドップラー効果をRGB色に適用（方法1: 波長シフト法）
Color applyDopplerShiftWavelength(const Color& original_color, double doppler_factor) {
    // 各RGB成分を対応する波長として扱い、シフト
    double shifted_r_wavelength = WAVELENGTH_RED / doppler_factor;
    double shifted_g_wavelength = WAVELENGTH_GREEN / doppler_factor;
    double shifted_b_wavelength = WAVELENGTH_BLUE / doppler_factor;

    // シフトした波長を新しいRGB値に変換
    Color shifted;

    // 赤成分: 元の赤がどの色に見えるか
    shifted.x += original_color.x * wavelengthToIntensity(shifted_r_wavelength, WAVELENGTH_RED);
    shifted.y += original_color.x * wavelengthToIntensity(shifted_r_wavelength, WAVELENGTH_GREEN);
    shifted.z += original_color.x * wavelengthToIntensity(shifted_r_wavelength, WAVELENGTH_BLUE);

    // 緑成分
    shifted.x += original_color.y * wavelengthToIntensity(shifted_g_wavelength, WAVELENGTH_RED);
    shifted.y += original_color.y * wavelengthToIntensity(shifted_g_wavelength, WAVELENGTH_GREEN);
    shifted.z += original_color.y * wavelengthToIntensity(shifted_g_wavelength, WAVELENGTH_BLUE);

    // 青成分
    shifted.x += original_color.z * wavelengthToIntensity(shifted_b_wavelength, WAVELENGTH_RED);
    shifted.y += original_color.z * wavelengthToIntensity(shifted_b_wavelength, WAVELENGTH_GREEN);
    shifted.z += original_color.z * wavelengthToIntensity(shifted_b_wavelength, WAVELENGTH_BLUE);

    return Color(Clamp(shifted.x,0,1), Clamp(shifted.y, 0, 1), Clamp(shifted.z, 0, 1));
}

// ドップラー効果をRGB色に適用（方法2: 簡易法）
// より高速だが物理的精度は低い
Color applyDopplerShiftSimple(const Color& original_color, double doppler_factor) {
    Color shifted = original_color;

    if (doppler_factor > 1.0) {
        // 青方偏移: 青を強調、赤を弱める
        double shift_amount = (doppler_factor - 1.0) * 2.0;
        shift_amount = std::min(shift_amount, 1.0);

        shifted.y += (1.0 - shifted.z) * shift_amount * 0.5;
        shifted.z += (1.0 - shifted.y) * shift_amount * 0.3;
        shifted.x *= (1.0 - shift_amount * 0.3);

    }
    else if (doppler_factor < 1.0) {
        // 赤方偏移: 赤を強調、青を弱める
        double shift_amount = (1.0 - doppler_factor) * 2.0;
        shift_amount = std::min(shift_amount, 1.0);

        shifted.x += (1.0 - shifted.x) * shift_amount * 0.5;
        shifted.y += (1.0 - shifted.y) * shift_amount * 0.2;
        shifted.z *= (1.0 - shift_amount * 0.3);
    }

    return Color(Clamp(shifted.x, 0, 1), Clamp(shifted.y, 0, 1), Clamp(shifted.z, 0, 1));
}

// ドップラー効果をRGB色に適用（方法3: 輝度変更も含む）
// 最も物理的に正確
Color applyDopplerShiftFull(const Color& original_color, double doppler_factor) {
    // ドップラー因子の4乗で輝度が変化（相対論的ビーミング効果）
    double intensity_factor = doppler_factor * doppler_factor *
        doppler_factor * doppler_factor;

    // 波長シフトを適用
    Color shifted = applyDopplerShiftWavelength(original_color, doppler_factor);

    // 輝度変更を適用
    shifted = shifted * intensity_factor;

    return Color(Clamp(shifted.x, 0, 1), Clamp(shifted.y, 0, 1), Clamp(shifted.z, 0, 1));
}
#endif
