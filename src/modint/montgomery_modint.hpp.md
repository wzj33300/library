---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "Traceback (most recent call last):\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: common.hpp: line\
    \ -1: no such header\n"
  code: "#ifndef MONTGOMERY_MODINT_HPP\r\n#define MONTGOMERY_MODINT_HPP\r\n\r\n#include\
    \ \"common.hpp\"\r\n\r\n#ifdef LIB_DEBUG\r\n  #include <exception>\r\n#endif\r\
    \n#include <cstdint>\r\n#include <iostream>\r\n#include <type_traits>\r\n\r\n\
    LIB_BEGIN\r\n\r\ntemplate <std::uint32_t ModT>\r\nclass montgomery_modint30 {\r\
    \n  using i32 = std::int32_t;\r\n  using u32 = std::uint32_t;\r\n  using u64 =\
    \ std::uint64_t;\r\n\r\n  u32 v_;\r\n\r\n  static constexpr u32 get_r() {\r\n\
    \    u32 t = 2, iv = MOD * (t - MOD * MOD);\r\n    iv *= t - MOD * iv, iv *= t\
    \ - MOD * iv;\r\n    return iv * (MOD * iv - t);\r\n  }\r\n  static constexpr\
    \ u32 redc(u64 x) {\r\n    return (x + static_cast<u64>(static_cast<u32>(x) *\
    \ R) * MOD) >> 32;\r\n  }\r\n  static constexpr u32 norm(u32 x) { return x - (MOD\
    \ & -((MOD - 1 - x) >> 31)); }\r\n\r\n  enum : u32 { MOD = ModT, MOD2 = MOD *\
    \ 2, R = get_r(), R2 = -static_cast<u64>(MOD) % MOD };\r\n  enum : i32 { SMOD\
    \ = MOD };\r\n\r\n  static_assert(MOD & 1);\r\n  static_assert(-R * MOD == 1);\r\
    \n  static_assert((MOD >> 30) == 0);\r\n  static_assert(MOD != 1);\r\n\r\npublic:\r\
    \n  static constexpr u32 mod() { return MOD; }\r\n  static constexpr i32 smod()\
    \ { return SMOD; }\r\n  constexpr montgomery_modint30() : v_() {}\r\n  template\
    \ <typename Int, std::enable_if_t<std::is_integral_v<Int>, int> = 0>\r\n  constexpr\
    \ montgomery_modint30(Int v) : v_(redc(static_cast<u64>(v % SMOD + SMOD) * R2))\
    \ {}\r\n  constexpr u32 val() const { return norm(redc(v_)); }\r\n  constexpr\
    \ i32 sval() const { return norm(redc(v_)); }\r\n  template <typename IntT, std::enable_if_t<std::is_integral_v<IntT>,\
    \ int> = 0>\r\n  explicit constexpr operator IntT() const {\r\n    return static_cast<IntT>(val());\r\
    \n  }\r\n  constexpr montgomery_modint30 operator-() const {\r\n    montgomery_modint30\
    \ res;\r\n    res.v_ = (MOD2 & -(v_ != 0)) - v_;\r\n    return res;\r\n  }\r\n\
    \  constexpr montgomery_modint30 inv() const {\r\n    i32 x1 = 1, x3 = 0, a =\
    \ sval(), b = SMOD;\r\n    while (b != 0) {\r\n      i32 q = a / b, x1_old = x1,\
    \ a_old = a;\r\n      x1 = x3, x3 = x1_old - x3 * q, a = b, b = a_old - b * q;\r\
    \n    }\r\n#ifdef LIB_DEBUG\r\n    if (a != 1) throw std::runtime_error(\"modular\
    \ inverse error\");\r\n#endif\r\n    return montgomery_modint30(x1);\r\n  }\r\n\
    \  constexpr montgomery_modint30 &operator+=(const montgomery_modint30 &rhs) {\r\
    \n    v_ += rhs.v_ - MOD2, v_ += MOD2 & -(v_ >> 31);\r\n    return *this;\r\n\
    \  }\r\n  constexpr montgomery_modint30 &operator-=(const montgomery_modint30\
    \ &rhs) {\r\n    v_ -= rhs.v_, v_ += MOD2 & -(v_ >> 31);\r\n    return *this;\r\
    \n  }\r\n  constexpr montgomery_modint30 &operator*=(const montgomery_modint30\
    \ &rhs) {\r\n    v_ = redc(static_cast<u64>(v_) * rhs.v_);\r\n    return *this;\r\
    \n  }\r\n  constexpr montgomery_modint30 &operator/=(const montgomery_modint30\
    \ &rhs) {\r\n    return operator*=(rhs.inv());\r\n  }\r\n  constexpr montgomery_modint30\
    \ pow(u64 e) const {\r\n    for (montgomery_modint30 res(1u), x(*this);; x *=\
    \ x) {\r\n      if (e & 1) res *= x;\r\n      if ((e >>= 1) == 0) return res;\r\
    \n    }\r\n  }\r\n  constexpr void swap(montgomery_modint30 &rhs) {\r\n    u32\
    \ v = v_;\r\n    v_ = rhs.v_, rhs.v_ = v;\r\n  }\r\n  friend constexpr montgomery_modint30\
    \ operator+(const montgomery_modint30 &lhs,\r\n                              \
    \                   const montgomery_modint30 &rhs) {\r\n    return montgomery_modint30(lhs)\
    \ += rhs;\r\n  }\r\n  friend constexpr montgomery_modint30 operator-(const montgomery_modint30\
    \ &lhs,\r\n                                                 const montgomery_modint30\
    \ &rhs) {\r\n    return montgomery_modint30(lhs) -= rhs;\r\n  }\r\n  friend constexpr\
    \ montgomery_modint30 operator*(const montgomery_modint30 &lhs,\r\n          \
    \                                       const montgomery_modint30 &rhs) {\r\n\
    \    return montgomery_modint30(lhs) *= rhs;\r\n  }\r\n  friend constexpr montgomery_modint30\
    \ operator/(const montgomery_modint30 &lhs,\r\n                              \
    \                   const montgomery_modint30 &rhs) {\r\n    return montgomery_modint30(lhs)\
    \ /= rhs;\r\n  }\r\n  friend constexpr bool operator==(const montgomery_modint30\
    \ &lhs, const montgomery_modint30 &rhs) {\r\n    return norm(lhs.v_) == norm(rhs.v_);\r\
    \n  }\r\n  friend constexpr bool operator!=(const montgomery_modint30 &lhs, const\
    \ montgomery_modint30 &rhs) {\r\n    return norm(lhs.v_) != norm(rhs.v_);\r\n\
    \  }\r\n  friend std::istream &operator>>(std::istream &is, montgomery_modint30\
    \ &rhs) {\r\n    i32 x;\r\n    is >> x;\r\n    rhs = montgomery_modint30(x);\r\
    \n    return is;\r\n  }\r\n  friend std::ostream &operator<<(std::ostream &os,\
    \ const montgomery_modint30 &rhs) {\r\n    return os << rhs.val();\r\n  }\r\n\
    };\r\n\r\ntemplate <std::uint32_t MOD>\r\nusing mm30 = montgomery_modint30<MOD>;\r\
    \n\r\nLIB_END\r\n\r\n#endif"
  dependsOn:
  - src/common.hpp
  isVerificationFile: false
  path: src/modint/montgomery_modint.hpp
  requiredBy: []
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
documentation_of: src/modint/montgomery_modint.hpp
layout: document
title: Montgomery ModInt
---
