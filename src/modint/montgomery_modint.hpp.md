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
  code: "#ifndef MONTGOMERY_MODINT_HPP\n#define MONTGOMERY_MODINT_HPP\n\n#include\
    \ \"common.hpp\"\n\n#ifdef LIB_DEBUG\n  #include <exception>\n#endif\n#include\
    \ <cstdint>\n#include <iostream>\n#include <type_traits>\n\nLIB_BEGIN\n\ntemplate\
    \ <std::uint32_t ModT>\nclass montgomery_modint30 {\n  using i32 = std::int32_t;\n\
    \  using u32 = std::uint32_t;\n  using u64 = std::uint64_t;\n\n  u32 v_;\n\n \
    \ static constexpr u32 get_r() {\n    u32 t = 2, iv = MOD * (t - MOD * MOD);\n\
    \    iv *= t - MOD * iv, iv *= t - MOD * iv;\n    return iv * (MOD * iv - t);\n\
    \  }\n  static constexpr u32 redc(u64 x) {\n    return (x + static_cast<u64>(static_cast<u32>(x)\
    \ * R) * MOD) >> 32;\n  }\n  static constexpr u32 norm(u32 x) { return x - (MOD\
    \ & -((MOD - 1 - x) >> 31)); }\n\n  enum : u32 { MOD = ModT, MOD2 = MOD * 2, R\
    \ = get_r(), R2 = -static_cast<u64>(MOD) % MOD };\n  enum : i32 { SMOD = MOD };\n\
    \n  static_assert(MOD & 1);\n  static_assert(-R * MOD == 1);\n  static_assert((MOD\
    \ >> 30) == 0);\n  static_assert(MOD != 1);\n\npublic:\n  static constexpr u32\
    \ mod() { return MOD; }\n  static constexpr i32 smod() { return SMOD; }\n  constexpr\
    \ montgomery_modint30() : v_() {}\n  template <typename Int, std::enable_if_t<std::is_integral_v<Int>,\
    \ int> = 0>\n  constexpr montgomery_modint30(Int v) : v_(redc(static_cast<u64>(v\
    \ % SMOD + SMOD) * R2)) {}\n  constexpr u32 val() const { return norm(redc(v_));\
    \ }\n  constexpr i32 sval() const { return norm(redc(v_)); }\n  template <typename\
    \ IntT, std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  explicit constexpr\
    \ operator IntT() const {\n    return static_cast<IntT>(val());\n  }\n  constexpr\
    \ montgomery_modint30 operator-() const {\n    montgomery_modint30 res;\n    res.v_\
    \ = (MOD2 & -(v_ != 0)) - v_;\n    return res;\n  }\n  constexpr montgomery_modint30\
    \ inv() const {\n    i32 x1 = 1, x3 = 0, a = sval(), b = SMOD;\n    while (b !=\
    \ 0) {\n      i32 q = a / b, x1_old = x1, a_old = a;\n      x1 = x3, x3 = x1_old\
    \ - x3 * q, a = b, b = a_old - b * q;\n    }\n#ifdef LIB_DEBUG\n    if (a != 1)\
    \ throw std::runtime_error(\"modular inverse error\");\n#endif\n    return montgomery_modint30(x1);\n\
    \  }\n  constexpr montgomery_modint30 &operator+=(const montgomery_modint30 &rhs)\
    \ {\n    v_ += rhs.v_ - MOD2, v_ += MOD2 & -(v_ >> 31);\n    return *this;\n \
    \ }\n  constexpr montgomery_modint30 &operator-=(const montgomery_modint30 &rhs)\
    \ {\n    v_ -= rhs.v_, v_ += MOD2 & -(v_ >> 31);\n    return *this;\n  }\n  constexpr\
    \ montgomery_modint30 &operator*=(const montgomery_modint30 &rhs) {\n    v_ =\
    \ redc(static_cast<u64>(v_) * rhs.v_);\n    return *this;\n  }\n  constexpr montgomery_modint30\
    \ &operator/=(const montgomery_modint30 &rhs) {\n    return operator*=(rhs.inv());\n\
    \  }\n  constexpr montgomery_modint30 pow(u64 e) const {\n    for (montgomery_modint30\
    \ res(1u), x(*this);; x *= x) {\n      if (e & 1) res *= x;\n      if ((e >>=\
    \ 1) == 0) return res;\n    }\n  }\n  constexpr void swap(montgomery_modint30\
    \ &rhs) {\n    u32 v = v_;\n    v_ = rhs.v_, rhs.v_ = v;\n  }\n  friend constexpr\
    \ montgomery_modint30 operator+(const montgomery_modint30 &lhs,\n            \
    \                                     const montgomery_modint30 &rhs) {\n    return\
    \ montgomery_modint30(lhs) += rhs;\n  }\n  friend constexpr montgomery_modint30\
    \ operator-(const montgomery_modint30 &lhs,\n                                \
    \                 const montgomery_modint30 &rhs) {\n    return montgomery_modint30(lhs)\
    \ -= rhs;\n  }\n  friend constexpr montgomery_modint30 operator*(const montgomery_modint30\
    \ &lhs,\n                                                 const montgomery_modint30\
    \ &rhs) {\n    return montgomery_modint30(lhs) *= rhs;\n  }\n  friend constexpr\
    \ montgomery_modint30 operator/(const montgomery_modint30 &lhs,\n            \
    \                                     const montgomery_modint30 &rhs) {\n    return\
    \ montgomery_modint30(lhs) /= rhs;\n  }\n  friend constexpr bool operator==(const\
    \ montgomery_modint30 &lhs, const montgomery_modint30 &rhs) {\n    return norm(lhs.v_)\
    \ == norm(rhs.v_);\n  }\n  friend constexpr bool operator!=(const montgomery_modint30\
    \ &lhs, const montgomery_modint30 &rhs) {\n    return norm(lhs.v_) != norm(rhs.v_);\n\
    \  }\n  friend std::istream &operator>>(std::istream &is, montgomery_modint30\
    \ &rhs) {\n    i32 x;\n    is >> x;\n    rhs = montgomery_modint30(x);\n    return\
    \ is;\n  }\n  friend std::ostream &operator<<(std::ostream &os, const montgomery_modint30\
    \ &rhs) {\n    return os << rhs.val();\n  }\n};\n\ntemplate <std::uint32_t MOD>\n\
    using mm30 = montgomery_modint30<MOD>;\n\nLIB_END\n\n#endif"
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
