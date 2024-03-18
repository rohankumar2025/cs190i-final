#pragma once
#include "common.h"
#include "math_helpers.h"
#include <algorithm>
#include <cmath>
#include <iostream>
// #include <winuser.h>
using std::tuple;
using std::clamp;
using std::cout;
using std::endl;
using std::swap;

#define PI 3.141592653589
#define INV_PI 0.31830989
#define EXP 2.718281828459045

namespace muni {
struct Lambertian {
    Vec3f albedo;

    /** Evaluates the BRDF for the Lambertian material.
      \return The BRDF (fr) value.
  */
    Vec3f eval() const {
        return albedo * (float) INV_PI; 
     }

    /** Samples the BRDF for the Lambertian material.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
      float r = sqrt(std::max(0.0f, 1.0f - u[0] * u[0]));
      float phi = 2 * PI * u[1];

      Vec3f v = normalize(Vec3f(r * std::cos(phi), r * std::sin(phi), u[0]));
      Vec3f dir = from_local(v, normal);

      return {dir, 1.0/(2.0 * PI)};
    }
    /** Computes the PDF for the Lambertian material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
      if (dot(wo_world, normal) < 0.0f) return 0;
      return INV_PI / 4.0f;
    }

};

struct Microfacet {
    const Vec3f LOCAL_N{0.0f, 0.0f, 1.0f};
    float roughness;
    // refraction indices for RGB channels
    Vec3f n1;
    Vec3f n2;

    /** Computes the Fresnel term for the microfacet material.
      \param[in] wi The light incident direction in local space.
      \return The Fresnel term.
    */
    Vec3f F(Vec3f wi) const {
      Vec3f R0 = {pow((n1.x-n2.x) / (n1.x+n2.x), 2), pow((n1.y-n2.y) / (n1.y+n2.y), 2), pow((n1.z-n2.z) / (n1.z+n2.z), 2)};
      float cos_theta = dot(normalize(wi), LOCAL_N);
      Vec3f one_minus_R0 = Vec3f{1.0f} - R0;
      float one_minus_cos_theta_to_fifth = pow(1.0f - cos_theta, 5);
      Vec3f res = R0 + one_minus_cos_theta_to_fifth * one_minus_R0;
      return res;
    }
    /** Computes the Beckmann normal distribution function for the microfacet material.
      \param[in] h The half vector in local space.
      \return The normal distribution function.
    */
    float D(Vec3f h) const {
      float cos_theta_h = dot(h, LOCAL_N);
      float cos_theta_h_squared = cos_theta_h * cos_theta_h;
      float alpha_squared = pow(roughness, 2);
      float tan_theta_h_squared = (1.0f - cos_theta_h_squared) / cos_theta_h_squared;

      float exponent = -tan_theta_h_squared / alpha_squared;
      float denominator = PI * alpha_squared * cos_theta_h_squared * cos_theta_h_squared;

      return pow(EXP, exponent) / denominator;
    }


    /** Computes the shadowing-masking function for the microfacet material.
      \param[in] wo The outgoing direction in local space.
      \param[in] wi The light incident direction in local space.
      \return The shadowing-masking value.
    */
    float G(Vec3f wo, Vec3f wi) const {
        return G1(wo) * G1(wi);
    }

    float G1(Vec3f w) const {
      float cos_theta = dot(normalize(w), LOCAL_N);
      float tan_theta = sqrt(1 - cos_theta * cos_theta)/ cos_theta;
      float alpha = roughness;

      float a = 1 / (roughness * tan_theta);
      float lambda = (a < 1.6f) ? (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181 * a * a) : 0.0f;
      return 1 / (1 + lambda);
    }

    /** Evaluates the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] wi_world The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The BRDF (fr) value.
    */
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
      wo_world = normalize(wo_world);
      wi_world = normalize(wi_world);
      normal = normalize(normal);
      
      Vec3f wo = to_local(wo_world, normal);
      Vec3f wi = to_local(wi_world, normal);
      Vec3f h = normalize(0.5f * (wi + wo));

      Vec3f res = F(wi) * G(wo, wi) * D(h) / (4 * dot(normal, wo_world) * dot(normal, wi_world));

      // float theta = acos(dot(normal, h)) * 180 / PI;
      // spdlog::info("eval called with theta: {}", theta);
      // spdlog::info("F: {}, G: {}, D: {}, cos_theta_h: {}, res: {}", F(wi), G(wo, wi), D(h), dot(h, LOCAL_N), res);
      return res;
    }
    /** Computes the PDF for the microfacet material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
      normal = normalize(normal);
      Vec3f wh_world = normalize(0.5f * (wo_world + wi_world));

      float cos_theta_h = dot(wh_world, normal);
      float sin_theta_h = sqrt(1 - cos_theta_h * cos_theta_h);
      float tan_theta_h = sin_theta_h / cos_theta_h;
      float alpha_squared = pow(roughness, 2);

      float pdf_theta = 2.0f * sin_theta_h / (alpha_squared * pow(cos_theta_h, 3)) * pow(EXP, (-pow(tan_theta_h, 2) / alpha_squared));
      spdlog::info(pdf_theta / (2.0f * PI * sin_theta_h));
      return pdf_theta / (2.0f * PI * sin_theta_h);
    }

    /** Samples the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float>  sample(Vec3f wo_world, Vec3f normal, Vec2f u) const {
      float phi_h = 2.0f * PI * u[0];
      float pdf_phi = 1 / (2.0f * PI);
      float alpha_squared = pow(roughness, 2);

      float theta_h = atan(sqrt(-alpha_squared * log(1.0f - u[1])));
      float tan_theta_h = tan(theta_h), cos_theta_h = cos(theta_h), sin_theta_h = sin(theta_h);
  
      float denominator = alpha_squared * cos_theta_h * cos_theta_h * cos_theta_h;
      float exponent = -tan_theta_h * tan_theta_h / alpha_squared;

      float pdf_theta = 2.0f * sin_theta_h / denominator * pow(EXP, exponent); 
      //spdlog::info("theta_h: {}, denominator: {}, exp: {}, pdf: {}", theta_h, denominator, exponent, pdf_theta); 

      Vec3f wh_local = {sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h)};

      Vec3f wh_world = from_local(wh_local, normal);
      // spdlog::info("Sampled wh_world: {}, normal: {}", wh_world, normal);

      Vec3f wi_world = mirror_reflect(-wo_world, wh_world);
      
      // float pdf_theta = pdf(wo_world, wi_world, normal);
      float pdf_wi = pdf_theta * pdf_phi / (4.0f * sin_theta_h * dot(normalize(wo_world), normalize(wh_world)));
      // spdlog::info("pdf_wi: {}, pdf_theta: {}, pdf_phi: {}, cos: {}, sin_theta: {}", pdf_wi, pdf_theta, pdf_phi, dot(normalize(-wo_world), normalize(wh_world)), sin_theta_h);

      return {normalize(wi_world), pdf_wi};
    }

};

struct Dielectric {
  float eta, roughness;

  float Fresnel(float cos_thetaI) const {
    cos_thetaI = std::clamp(cos_thetaI, -1.0f, 1.0f);
    bool entering = cos_thetaI > 0.0f;
    float ei = eta, et = 1.0f;
    if (!entering)
      std::swap(ei, et);

    float sin_thetaT = ei / et * sqrt(std::max(0.0f, 1.0f - cos_thetaI * cos_thetaI));
    if (sin_thetaT >= 1.0f)
      return 1.0f;

    cos_thetaI = abs(cos_thetaI);
    float cos_thetaT = sqrt(std::max(0.0f, 1.0f - sin_thetaT * sin_thetaT));

    float Rs = (et * cos_thetaI - ei * cos_thetaT) / (et * cos_thetaI + ei * cos_thetaT);
    float Rp = (ei * cos_thetaI - et * cos_thetaT) / (ei * cos_thetaI + et * cos_thetaT);

    return 0.5f * (Rs * Rs + Rp * Rp);
  }

  float D(Vec3f h) const {   
    if (h.z < 0.0f) return 0.0f;

    float cos_thetaH_sq = h.z * h.z;
    float tan_thetaH_sq = (1.0f - cos_thetaH_sq) / cos_thetaH_sq;
    float alpha_sq = roughness * roughness;

    float sqrt_denom = alpha_sq + tan_thetaH_sq;

    return alpha_sq * INV_PI / (cos_thetaH_sq * cos_thetaH_sq * sqrt_denom * sqrt_denom);
  }

  float G(Vec3f wo, Vec3f wi, Vec3f h) const {
    return G1(wo, h) * G1(wi, h);
  }

  float G1(Vec3f w, Vec3f h) const {
    if (dot(w, h) * w.z <= 0.0f)
      return 0.0f;

    float tan_thetaW_sq = (1.0f - w.z * w.z) / (w.z * w.z);
    if (tan_thetaW_sq <= 0.0f)
      return 1.0f;

    float tan_thetaW = sqrt(tan_thetaW_sq);
    float root = roughness * tan_thetaW;

    return 2.0f / (1.0f + sqrt(1.0f + root * root));
  }

  float eval(Vec3f wo_world, Vec3f wi_world, Vec3f n) const {
    Vec3f wo = normalize(to_local(wo_world, n)), wi = normalize(to_local(wi_world, n));
    bool reflect = (wi.z * wo.z) > 0.0f;
    
    float etaT = (wo.z < 0.0f) ? 1.0f : eta;
    float etaI = (wo.z < 0.0f) ? eta : 1.0f;

    float reflect_coeff = (wi.z > 0.0f) ? 1.0f : -1.0f;
    Vec3f h = (reflect) ? normalize(wi + wo) * reflect_coeff : -normalize(wo * etaI + wi * etaT);

    float Dr = D(h);
    float Fr = Fresnel(dot(wo, h));
    float Gr = G(wo, wi, h);

    if (reflect) {
      float res = Fr * Dr * Gr / (4.0f * wi.z * wo.z);
      // spdlog::info("h: {}, Fr: {}, Dr: {}, Gr: {}, eval: {}", h, Fr, Dr, Gr, res);
    }

    float sqrt_denom = etaI * dot(wi, h) + etaT * dot(wo, h);

    float res = abs((1.0f - Fr) * Dr * Gr * etaT * etaT * dot(wi, h) * dot(wo, h) / (wi.z * wo.z * sqrt_denom * sqrt_denom)); 
    // spdlog::info("h: {}, Fr: {}, Dr: {}, Gr: {}, eval: {}", h, Fr, Dr, Gr, res);
    return res;
  }

  tuple<Vec3f, float> sample_normal(Vec2f u) const {
    float thetaM = atan(roughness * sqrt(u.x) / sqrt(1.0f - u.x));
    float phiM = 2 * PI * u.y;

    Vec3f m = normalize(Vec3f{sin(thetaM) * sin(phiM), sin(thetaM) * cos(phiM), cos(thetaM)});
    float pdf_m = D(m) * m.z;

    return {m, pdf_m};
  }


  tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f n, Vec3f u) const {
    Vec3f wo = normalize(to_local(wo_world, n));

    auto [m, pdf_m] = sample_normal(Vec2f{u.y, u.z});
    // spdlog::info("m: {}, pdf_m: {}", m, pdf_m);
  
    float cos_thetaO = dot(wo, m);
    float F = Fresnel(cos_thetaO);
    // spdlog::info("Fresnel_m: {}, Fresnel: {}", F, Fresnel(wo.z));

    if (u.x < F) {
      Vec3f wi = mirror_reflect(-wo, m);
      Vec3f wi_world = from_local(normalize(wi), n);
      return {wi_world, F};
    }

    bool entering = wo.z > 0.0f;
    float ei = eta, et = 1.0f;
    if (!entering) swap(ei, et);
    
    float sin_thetaO2 = 1.0f - cos_thetaO * cos_thetaO;
    float etaEff = ei / et;

    float sin_thetaT2 = etaEff * etaEff * sin_thetaO2;

    if (sin_thetaT2 >= 1.0f)
      return {Vec3f{}, 0.0f};

    float cos_thetaT = sqrt(1.0f - sin_thetaT2);

    if (entering)
      cos_thetaT *= -1;

    Vec3f wt = {etaEff * -wo.x, etaEff * -wo.y, cos_thetaT};

    // spdlog::info("wo: {}, wt: {}", wo, wt);

    Vec3f wt_world = from_local(normalize(wt), n);

    return {wt_world, (1.0f - F)};
  }



};


}  // namespace muni
