---
data:
  _extendedDependsOn:
  - icon: ':question:'
    path: common.hpp
    title: common.hpp
  - icon: ':question:'
    path: common.hpp
    title: common.hpp
  - icon: ':question:'
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: modint/long_montgomery_modint.hpp
    title: Long Montgomery ModInt
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/convolution.hpp\"\n\n\n\n#line 1 \"common.hpp\"\n\n\
    \n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n#define LIB_END }\n\
    #define LIB ::lib::\n\n\n#line 1 \"modint/long_montgomery_modint.hpp\"\n\n\n\n\
    #line 5 \"modint/long_montgomery_modint.hpp\"\n\n#ifdef LIB_DEBUG\n  #include\
    \ <stdexcept>\n#endif\n#include <cstdint>\n#include <iostream>\n#include <type_traits>\n\
    \nLIB_BEGIN\n\ntemplate <std::uint64_t ModT>\nclass montgomery_modint63 {\n  using\
    \ u32 = std::uint32_t;\n  using i64 = std::int64_t;\n  using u64 = std::uint64_t;\n\
    \n  u64 v_{};\n\n  static constexpr u64 get_r() {\n    u64 t = 2, iv = MOD * (t\
    \ - MOD * MOD);\n    iv *= t - MOD * iv, iv *= t - MOD * iv, iv *= t - MOD * iv;\n\
    \    return iv * (t - MOD * iv);\n  }\n  static constexpr u64 get_r2() {\n   \
    \ u64 iv = -u64(MOD) % MOD;\n    for (int i = 0; i != 64; ++i)\n      if ((iv\
    \ <<= 1) >= MOD) iv -= MOD;\n    return iv;\n  }\n  static constexpr u64 mul_high(u64\
    \ x, u64 y) {\n    u64 a = x >> 32, b = static_cast<u32>(x), c = y >> 32, d =\
    \ static_cast<u32>(y), ad = a * d,\n        bc = b * c;\n    return a * c + (ad\
    \ >> 32) + (bc >> 32) +\n           (((ad & 0xFFFFFFFF) + (bc & 0xFFFFFFFF) +\
    \ (b * d >> 32)) >> 32);\n  }\n  static constexpr u64 redc_mul(u64 x, u64 y) {\n\
    \    u64 res = mul_high(x, y) - mul_high(x * y * R, MOD);\n    return res + (MOD\
    \ & -(res >> 63));\n  }\n  static constexpr u64 norm(i64 x) { return x + (MOD\
    \ & -(x < 0)); }\n\n  static constexpr u64 MOD  = ModT;\n  static constexpr u64\
    \ R    = get_r();\n  static constexpr u64 R2   = get_r2();\n  static constexpr\
    \ i64 SMOD = static_cast<i64>(MOD);\n\n  static_assert(MOD & 1);\n  static_assert(R\
    \ * MOD == 1);\n  static_assert((MOD >> 63) == 0);\n  static_assert(MOD != 1);\n\
    \npublic:\n  static constexpr u64 mod() { return MOD; }\n  static constexpr i64\
    \ smod() { return SMOD; }\n  constexpr montgomery_modint63() {}\n  template <typename\
    \ IntT, std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  constexpr montgomery_modint63(IntT\
    \ v) : v_(redc_mul(norm(v % SMOD), R2)) {}\n  constexpr u64 val() const {\n  \
    \  u64 res = -mul_high(v_ * R, MOD);\n    return res + (MOD & -(res >> 63));\n\
    \  }\n  constexpr i64 sval() const { return val(); }\n  constexpr bool is_zero()\
    \ const { return v_ == 0; }\n  template <typename IntT, std::enable_if_t<std::is_integral_v<IntT>,\
    \ int> = 0>\n  explicit constexpr operator IntT() const {\n    return static_cast<IntT>(val());\n\
    \  }\n  constexpr montgomery_modint63 operator-() const {\n    montgomery_modint63\
    \ res;\n    res.v_ = (MOD & -(v_ != 0)) - v_;\n    return res;\n  }\n  constexpr\
    \ montgomery_modint63 inv() const {\n    i64 x1 = 1, x3 = 0, a = sval(), b = SMOD;\n\
    \    while (b != 0) {\n      i64 q = a / b, x1_old = x1, a_old = a;\n      x1\
    \ = x3, x3 = x1_old - x3 * q, a = b, b = a_old - b * q;\n    }\n#ifdef LIB_DEBUG\n\
    \    if (a != 1) throw std::runtime_error(\"modular inverse error\");\n#endif\n\
    \    return montgomery_modint63(x1);\n  }\n  constexpr montgomery_modint63 &operator+=(const\
    \ montgomery_modint63 &rhs) {\n    v_ += rhs.v_ - MOD, v_ += MOD & -(v_ >> 63);\n\
    \    return *this;\n  }\n  constexpr montgomery_modint63 &operator-=(const montgomery_modint63\
    \ &rhs) {\n    v_ -= rhs.v_, v_ += MOD & -(v_ >> 63);\n    return *this;\n  }\n\
    \  constexpr montgomery_modint63 &operator*=(const montgomery_modint63 &rhs) {\n\
    \    v_ = redc_mul(v_, rhs.v_);\n    return *this;\n  }\n  constexpr montgomery_modint63\
    \ &operator/=(const montgomery_modint63 &rhs) {\n    return operator*=(rhs.inv());\n\
    \  }\n  constexpr montgomery_modint63 pow(u64 e) const {\n    for (montgomery_modint63\
    \ res(1), x(*this);; x *= x) {\n      if (e & 1) res *= x;\n      if ((e >>= 1)\
    \ == 0) return res;\n    }\n  }\n  constexpr void swap(montgomery_modint63 &rhs)\
    \ {\n    auto v = v_;\n    v_ = rhs.v_, rhs.v_ = v;\n  }\n  friend constexpr montgomery_modint63\
    \ operator+(const montgomery_modint63 &lhs,\n                                \
    \                 const montgomery_modint63 &rhs) {\n    return montgomery_modint63(lhs)\
    \ += rhs;\n  }\n  friend constexpr montgomery_modint63 operator-(const montgomery_modint63\
    \ &lhs,\n                                                 const montgomery_modint63\
    \ &rhs) {\n    return montgomery_modint63(lhs) -= rhs;\n  }\n  friend constexpr\
    \ montgomery_modint63 operator*(const montgomery_modint63 &lhs,\n            \
    \                                     const montgomery_modint63 &rhs) {\n    return\
    \ montgomery_modint63(lhs) *= rhs;\n  }\n  friend constexpr montgomery_modint63\
    \ operator/(const montgomery_modint63 &lhs,\n                                \
    \                 const montgomery_modint63 &rhs) {\n    return montgomery_modint63(lhs)\
    \ /= rhs;\n  }\n  friend constexpr bool operator==(const montgomery_modint63 &lhs,\
    \ const montgomery_modint63 &rhs) {\n    return lhs.v_ == rhs.v_;\n  }\n  friend\
    \ constexpr bool operator!=(const montgomery_modint63 &lhs, const montgomery_modint63\
    \ &rhs) {\n    return lhs.v_ != rhs.v_;\n  }\n  friend std::istream &operator>>(std::istream\
    \ &is, montgomery_modint63 &rhs) {\n    i64 x;\n    is >> x;\n    rhs = montgomery_modint63(x);\n\
    \    return is;\n  }\n  friend std::ostream &operator<<(std::ostream &os, const\
    \ montgomery_modint63 &rhs) {\n    return os << rhs.val();\n  }\n};\n\ntemplate\
    \ <std::uint64_t ModT>\nusing mm63 = montgomery_modint63<ModT>;\n\nLIB_END\n\n\
    \n#line 1 \"math/radix2_ntt.hpp\"\n\n\n\n#line 5 \"math/radix2_ntt.hpp\"\n\n#include\
    \ <array>\n#include <cassert>\n#line 9 \"math/radix2_ntt.hpp\"\n#include <vector>\n\
    \nLIB_BEGIN\n\nnamespace detail {\n\ntemplate <typename IntT>\nconstexpr std::enable_if_t<std::is_integral_v<IntT>,\
    \ int> bsf(IntT v) {\n  if (static_cast<std::make_signed_t<IntT>>(v) <= 0) return\
    \ -1;\n  int res = 0;\n  for (; (v & 1) == 0; ++res) v >>= 1;\n  return res;\n\
    }\n\ntemplate <typename ModIntT>\nconstexpr ModIntT quadratic_nonresidue_prime()\
    \ {\n  auto mod = ModIntT::mod();\n  for (int i = 2;; ++i)\n    if (ModIntT(i).pow(mod\
    \ >> 1) == mod - 1) return ModIntT(i);\n}\n\ntemplate <typename ModIntT>\nconstexpr\
    \ ModIntT gen_of_sylow_2_subgroup() {\n  auto mod = ModIntT::mod();\n  return\
    \ quadratic_nonresidue_prime<ModIntT>().pow(mod >> bsf(mod - 1));\n}\n\ntemplate\
    \ <typename ModIntT>\nconstexpr std::array<ModIntT, bsf(ModIntT::mod() - 1) -\
    \ 1> root() {\n  std::array<ModIntT, bsf(ModIntT::mod() - 1) - 1> rt; // order(`rt[i]`)\
    \ = 2^(i + 2).\n  rt.back() = gen_of_sylow_2_subgroup<ModIntT>();\n  for (int\
    \ i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) rt[i] = rt[i + 1] * rt[i + 1];\n\
    \  return rt;\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> iroot() {\n  std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> irt;\n  irt.back() = gen_of_sylow_2_subgroup<ModIntT>().inv();\n \
    \ for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) irt[i] = irt[i + 1] *\
    \ irt[i + 1];\n  return irt;\n}\n\n} // namespace detail\n\n// Input:  integer\
    \ `n`.\n// Output: 2^(\u2308log_2(`n`)\u2309).\nint ntt_len(int n) {\n  --n;\n\
    \  n |= n >> 1;\n  n |= n >> 2;\n  n |= n >> 4;\n  n |= n >> 8;\n  return (n |\
    \ n >> 16) + 1;\n}\n\n// Input:           f(x) = `a[0]` + `a[1]`x + ... + `a[n\
    \ - 1]`x^(`n` - 1) where `n` is power of 2.\n// Output(inplace): reversed binary\
    \ permutation of [f(\u03B6^0), f(\u03B6), f(\u03B6^2), ..., f(\u03B6^(`n` - 1))].\n\
    template <typename IterT>\nvoid dft_n(IterT a, int n) {\n  assert((n & (n - 1))\
    \ == 0);\n  using T                  = typename std::iterator_traits<IterT>::value_type;\n\
    \  static constexpr auto rt = detail::root<T>();\n  static std::vector<T> root(1);\n\
    \  if (int s = static_cast<int>(root.size()); s << 1 < n) {\n    root.resize(n\
    \ >> 1);\n    for (int i = detail::bsf(s); (1 << i) < (n >> 1); ++i) {\n     \
    \ int j   = 1 << i;\n      root[j] = rt[i];\n      for (int k = j + 1; k < (j\
    \ << 1); ++k) root[k] = root[k - j] * root[j];\n    }\n  }\n  for (int j = 0,\
    \ l = n >> 1; j != l; ++j) {\n    T u(a[j]), v(a[j + l]);\n    a[j] = u + v, a[j\
    \ + l] = u - v;\n  }\n  for (int i = n >> 1; i >= 2; i >>= 1) {\n    for (int\
    \ j = 0, l = i >> 1; j != l; ++j) {\n      T u(a[j]), v(a[j + l]);\n      a[j]\
    \ = u + v, a[j + l] = u - v;\n    }\n    for (int j = i, l = i >> 1, m = 1; j\
    \ != n; j += i, ++m) {\n      for (int k = j; k != j + l; ++k) {\n        T u(a[k]),\
    \ v(a[k + l] * root[m]);\n        a[k] = u + v, a[k + l] = u - v;\n      }\n \
    \   }\n  }\n}\n\n// Input:           reversed binary permutation of [f(\u03B6\
    ^0), f(\u03B6), f(\u03B6^2), ..., f(\u03B6^(`n` - 1))].\n// Output(inplace): f(x)\
    \ = `a[0]` + `a[1]`x + ... + `a[n - 1]`x^(`n` - 1) where `n` is power of 2.\n\
    template <typename IterT>\nvoid idft_n(IterT a, int n) {\n  assert((n & (n - 1))\
    \ == 0);\n  using T                  = typename std::iterator_traits<IterT>::value_type;\n\
    \  static constexpr auto rt = detail::iroot<T>();\n  static std::vector<T> root(1);\n\
    \  if (int s = static_cast<int>(root.size()); s << 1 < n) {\n    root.resize(n\
    \ >> 1);\n    for (int i = detail::bsf(s); (1 << i) < (n >> 1); ++i) {\n     \
    \ int j   = 1 << i;\n      root[j] = rt[i];\n      for (int k = j + 1; k < (j\
    \ << 1); ++k) root[k] = root[k - j] * root[j];\n    }\n  }\n  for (int i = 2;\
    \ i < n; i <<= 1) {\n    for (int j = 0, l = i >> 1; j != l; ++j) {\n      T u(a[j]),\
    \ v(a[j + l]);\n      a[j] = u + v, a[j + l] = u - v;\n    }\n    for (int j =\
    \ i, l = i >> 1, m = 1; j != n; j += i, ++m) {\n      for (int k = j; k != j +\
    \ l; ++k) {\n        T u(a[k]), v(a[k + l]);\n        a[k] = u + v, a[k + l] =\
    \ (u - v) * root[m];\n      }\n    }\n  }\n  const T iv(T::mod() - T::mod() /\
    \ n);\n  for (int j = 0, l = n >> 1; j != l; ++j) {\n    T u(a[j] * iv), v(a[j\
    \ + l] * iv);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n}\n\ntemplate <typename\
    \ ContainerT>\nvoid dft(ContainerT &a) {\n  dft_n(a.begin(), a.size());\n}\n\n\
    template <typename ContainerT>\nvoid idft(ContainerT &a) {\n  idft_n(a.begin(),\
    \ a.size());\n}\n\nLIB_END\n\n\n#line 7 \"math/convolution.hpp\"\n\n#include <algorithm>\n\
    #line 12 \"math/convolution.hpp\"\n\nLIB_BEGIN\n\ntemplate <typename ModIntT>\n\
    std::vector<ModIntT> convolution(const std::vector<ModIntT> &lhs, const std::vector<ModIntT>\
    \ &rhs) {\n  int n = static_cast<int>(lhs.size()), m = static_cast<int>(rhs.size());\n\
    \  if (n == 0 || m == 0) return std::vector<ModIntT>{};\n  if (std::min(n, m)\
    \ <= 32) {\n    std::vector<ModIntT> res(n + m - 1);\n    for (int i = 0; i !=\
    \ n; ++i)\n      for (int j = 0; j != m; ++j) res[i + j] += lhs[i] * rhs[j];\n\
    \    return res;\n  }\n  int len = ntt_len(n + m - 1);\n  std::vector<ModIntT>\
    \ lhs_cpy(len), rhs_cpy(len);\n  std::copy_n(lhs.cbegin(), n, lhs_cpy.begin());\n\
    \  std::copy_n(rhs.cbegin(), m, rhs_cpy.begin());\n  dft_n(lhs_cpy.begin(), len),\
    \ dft_n(rhs_cpy.begin(), len);\n  for (int i = 0; i != len; ++i) lhs_cpy[i] *=\
    \ rhs_cpy[i];\n  idft_n(lhs_cpy.begin(), len);\n  lhs_cpy.resize(n + m - 1);\n\
    \  return lhs_cpy;\n}\n\ntemplate <typename IntT>\nstd::enable_if_t<std::is_integral_v<IntT>\
    \ && sizeof(IntT) <= sizeof(std::int32_t),\n                 std::vector<IntT>>\n\
    convolution_mod(const std::vector<IntT> &lhs, const std::vector<IntT> &rhs, const\
    \ IntT modular) {\n  using mint0 = mm63<0x3F9A000000000001>;\n  using mint1 =\
    \ mm63<0x3FC6000000000001>;\n  auto res0   = convolution(std::vector<mint0>(lhs.begin(),\
    \ lhs.end()),\n                          std::vector<mint0>(rhs.begin(), rhs.end()));\n\
    \  auto res1   = convolution(std::vector<mint1>(lhs.begin(), lhs.end()),\n   \
    \                       std::vector<mint1>(rhs.begin(), rhs.end()));\n  const\
    \ int n = res0.size();\n  std::vector<IntT> res(n);\n  //    a mod m_0 = a_0,\
    \ a mod m_1 = a_1\n  // -> a_0 + k_0m_0 = a_1 + k_1m_1\n  // -> a_0 - a_1 \u2261\
    \ k_1m_1 (mod m_0)\n  // -> k_1 \u2261 (a_0 - a_1) / m_1 (mod m_0)\n  static constexpr\
    \ mint0 im1_mod_m0(mint0(mint1::mod()).inv());\n  const IntT m1_mod_modular =\
    \ mint1::mod() % modular;\n  for (int i = 0; i != n; ++i) {\n    mint0 k1((res0[i]\
    \ - res1[i].val()) * im1_mod_m0);\n    res[i] = (k1.val() % modular * m1_mod_modular\
    \ + res1[i].val()) % modular;\n  }\n  return res;\n}\n\nLIB_END\n\n\n"
  code: "#ifndef CONVOLUTION_HPP\n#define CONVOLUTION_HPP\n\n#include \"../common.hpp\"\
    \n#include \"../modint/long_montgomery_modint.hpp\"\n#include \"radix2_ntt.hpp\"\
    \n\n#include <algorithm>\n#include <cstdint>\n#include <type_traits>\n#include\
    \ <vector>\n\nLIB_BEGIN\n\ntemplate <typename ModIntT>\nstd::vector<ModIntT> convolution(const\
    \ std::vector<ModIntT> &lhs, const std::vector<ModIntT> &rhs) {\n  int n = static_cast<int>(lhs.size()),\
    \ m = static_cast<int>(rhs.size());\n  if (n == 0 || m == 0) return std::vector<ModIntT>{};\n\
    \  if (std::min(n, m) <= 32) {\n    std::vector<ModIntT> res(n + m - 1);\n   \
    \ for (int i = 0; i != n; ++i)\n      for (int j = 0; j != m; ++j) res[i + j]\
    \ += lhs[i] * rhs[j];\n    return res;\n  }\n  int len = ntt_len(n + m - 1);\n\
    \  std::vector<ModIntT> lhs_cpy(len), rhs_cpy(len);\n  std::copy_n(lhs.cbegin(),\
    \ n, lhs_cpy.begin());\n  std::copy_n(rhs.cbegin(), m, rhs_cpy.begin());\n  dft_n(lhs_cpy.begin(),\
    \ len), dft_n(rhs_cpy.begin(), len);\n  for (int i = 0; i != len; ++i) lhs_cpy[i]\
    \ *= rhs_cpy[i];\n  idft_n(lhs_cpy.begin(), len);\n  lhs_cpy.resize(n + m - 1);\n\
    \  return lhs_cpy;\n}\n\ntemplate <typename IntT>\nstd::enable_if_t<std::is_integral_v<IntT>\
    \ && sizeof(IntT) <= sizeof(std::int32_t),\n                 std::vector<IntT>>\n\
    convolution_mod(const std::vector<IntT> &lhs, const std::vector<IntT> &rhs, const\
    \ IntT modular) {\n  using mint0 = mm63<0x3F9A000000000001>;\n  using mint1 =\
    \ mm63<0x3FC6000000000001>;\n  auto res0   = convolution(std::vector<mint0>(lhs.begin(),\
    \ lhs.end()),\n                          std::vector<mint0>(rhs.begin(), rhs.end()));\n\
    \  auto res1   = convolution(std::vector<mint1>(lhs.begin(), lhs.end()),\n   \
    \                       std::vector<mint1>(rhs.begin(), rhs.end()));\n  const\
    \ int n = res0.size();\n  std::vector<IntT> res(n);\n  //    a mod m_0 = a_0,\
    \ a mod m_1 = a_1\n  // -> a_0 + k_0m_0 = a_1 + k_1m_1\n  // -> a_0 - a_1 \u2261\
    \ k_1m_1 (mod m_0)\n  // -> k_1 \u2261 (a_0 - a_1) / m_1 (mod m_0)\n  static constexpr\
    \ mint0 im1_mod_m0(mint0(mint1::mod()).inv());\n  const IntT m1_mod_modular =\
    \ mint1::mod() % modular;\n  for (int i = 0; i != n; ++i) {\n    mint0 k1((res0[i]\
    \ - res1[i].val()) * im1_mod_m0);\n    res[i] = (k1.val() % modular * m1_mod_modular\
    \ + res1[i].val()) % modular;\n  }\n  return res;\n}\n\nLIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  - modint/long_montgomery_modint.hpp
  - common.hpp
  - math/radix2_ntt.hpp
  isVerificationFile: false
  path: math/convolution.hpp
  requiredBy: []
  timestamp: '2022-04-23 22:52:36+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
documentation_of: math/convolution.hpp
layout: document
title: Convolution
---
