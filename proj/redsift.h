// redshift_rgb.cpp
// Build: g++ -O2 -std=c++17 redshift_rgb.cpp -o redshift_rgb
// Usage:
//   ./redshift_rgb r Rin Gin Bin [--model 1|2] [--sigma 40]
// Example:
//   ./redshift_rgb 2.001 220 235 255 --model 2 --sigma 40
//
// Note: Units G=c=M=1. Valid for r>2. Performs linear-RGB math.

//#include <cmath>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <algorithm>
#include <iostream>
#include <tuple>

struct RGBd { double r,g,b; };

static inline double clamp01(double x){ return x<0?0:(x>1?1:x); }
static inline int toByte(double x){ return (int)std::lround(clamp01(x)*255.0); }

// sRGB (0..1) -> linear
static inline double srgb_to_linear_1(double c){
    // exact sRGB EOTF
    if (c <= 0.04045) return c/12.92;
    return std::pow((c+0.055)/1.055, 2.4);
}
static inline double linear_to_srgb_1(double c){
    // exact sRGB OETF
    if (c <= 0.0031308) return 12.92*c;
    return 1.055*std::pow(c, 1.0/2.4) - 0.055;
}

static inline RGBd srgb_to_linear_u8(int R, int G, int B){
    return {
        srgb_to_linear_1(R/255.0),
        srgb_to_linear_1(G/255.0),
        srgb_to_linear_1(B/255.0)
    };
}
static inline RGBd linear_to_srgb_u8(RGBd L){
    return {
        linear_to_srgb_1(L.r),
        linear_to_srgb_1(L.g),
        linear_to_srgb_1(L.b)
    };
}

// redshift factor g(r) = sqrt(1 - 2/r), r>2
static inline double redshift_g(double r){
    if (r <= 2.0){
        // guard: approach 2+ from outside
        const double eps = 1e-12;
        r = 2.0 + eps;
    }
    return std::sqrt(1.0 - 2.0/r);
}

// -------- Model 1: intensity only (color-preserving) --------
// L_out = g^3 * L_in  (in linear RGB)
static inline RGBd apply_model1(RGBd Lin, double g){
    double s = g*g*g;
    return { Lin.r*s, Lin.g*s, Lin.b*s };
}

inline Color Redsift(Color& rgb, double dist)
{
    double r = dist;
	int Rin = to_int(rgb.r/255);
	int Gin = to_int(rgb.g/255);
	int Bin = to_int(rgb.b/255);

    // sRGB -> linear
    RGBd Lin = srgb_to_linear_u8(Rin, Gin, Bin);

    double g = redshift_g(r);

    RGBd Lout_lin;
    Lout_lin = apply_model1(Lin, g);

    // simple tonemap (Reinhard) to avoid blowout, then linear->sRGB
    auto reinhard = [](double x){ return x / (1.0 + x); };
    Lout_lin.r = reinhard(Lout_lin.r);
    Lout_lin.g = reinhard(Lout_lin.g);
    Lout_lin.b = reinhard(Lout_lin.b);

    RGBd Lout_srgb = linear_to_srgb_u8(Lout_lin);

    int Rout = toByte(Lout_srgb.r);
    int Gout = toByte(Lout_srgb.g);
    int Bout = toByte(Lout_srgb.b);

	Color out = Color(Rout/255.0, Gout/255.0, Bout/255.0);
    return out;
}
