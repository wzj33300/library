---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: common.hpp
    title: common.hpp
  - icon: ':heavy_check_mark:'
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.1.test.cpp
    title: remote_test/yosupo/math/convolution_mod.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/truncated_formal_power_series.hpp\"\n\n\n\n#line 1\
    \ \"common.hpp\"\n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n\
    #define LIB_END }\n#define LIB ::lib::\n\n\n#line 1 \"math/radix2_ntt.hpp\"\n\n\
    \n\n#line 5 \"math/radix2_ntt.hpp\"\n\n#include <array>\n#include <cassert>\n\
    #include <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\nnamespace detail {\n\
    \ntemplate <typename IntT>\nconstexpr std::enable_if_t<std::is_integral_v<IntT>,\
    \ int> bsf(IntT v) {\n  if (static_cast<std::make_signed_t<IntT>>(v) <= 0) return\
    \ -1;\n  int res = 0;\n  for (; (v & 1) == 0; ++res) v >>= 1;\n  return res;\n\
    }\n\ntemplate <typename ModIntT>\nconstexpr ModIntT quadratic_nonresidue_prime()\
    \ {\n  auto mod = ModIntT::mod();\n  for (int i = 2;; ++i)\n    if (ModIntT(i).pow(mod\
    \ >> 1) == mod - 1) return ModIntT(i);\n}\n\ntemplate <typename ModIntT>\nconstexpr\
    \ ModIntT gen_of_sylow2subgroup() {\n  auto mod = ModIntT::mod();\n  return quadratic_nonresidue_prime<ModIntT>().pow(mod\
    \ >> bsf(mod - 1));\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> root() {\n  std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> rt; // order(`rt[i]`) = 2^(i + 2).\n  rt.back() = gen_of_sylow2subgroup<ModIntT>();\n\
    \  for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) rt[i] = rt[i + 1] *\
    \ rt[i + 1];\n  return rt;\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> iroot() {\n  std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> irt;\n  irt.back() = gen_of_sylow2subgroup<ModIntT>().inv();\n  for\
    \ (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) irt[i] = irt[i + 1] * irt[i\
    \ + 1];\n  return irt;\n}\n\n} // namespace detail\n\n// Input:  integer `n`.\n\
    // Output: 2^(\u2308log_2(`n`)\u2309).\nint ntt_len(int n) {\n  --n;\n  n |= n\
    \ >> 1;\n  n |= n >> 2;\n  n |= n >> 4;\n  n |= n >> 8;\n  return (n | n >> 16)\
    \ + 1;\n}\n\n// Input:           f(x) = `a[0]` + `a[1]`x + ... + `a[n - 1]`x^(`n`\
    \ - 1) where `n` is power of 2.\n// Output(inplace): reversed binary permutation\
    \ of [f(\u03B6^0), f(\u03B6), f(\u03B6^2), ..., f(\u03B6^(`n` - 1))].\ntemplate\
    \ <typename IterT>\nvoid dft_n(IterT a, int n) {\n  assert((n & (n - 1)) == 0);\n\
    \  using T                  = typename std::iterator_traits<IterT>::value_type;\n\
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
    \ + l] * iv);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n}\n\nLIB_END\n\n\n#line\
    \ 6 \"math/truncated_formal_power_series.hpp\"\n\n#include <algorithm>\n#line\
    \ 9 \"math/truncated_formal_power_series.hpp\"\n#include <iostream>\n#include\
    \ <iterator>\n#line 13 \"math/truncated_formal_power_series.hpp\"\n\nLIB_BEGIN\n\
    \ntemplate <typename ModIntT>\nclass truncated_formal_power_series : public std::vector<ModIntT>\
    \ {\n  static_assert(std::is_same_v<typename std::vector<ModIntT>::value_type,\
    \ ModIntT>);\n\npublic:\n  using std::vector<ModIntT>::vector;\n\n  enum : int\
    \ { NEGATIVE_INFINITY = -1 };\n\n  // leading coefficient\n  ModIntT lc() const\
    \ {\n    int d = deg();\n    return d == NEGATIVE_INFINITY ? ModIntT() : this->operator[](d);\n\
    \  }\n  // degree\n  int deg() const {\n    // treat formal power series like\
    \ polynomials\n    int n = static_cast<int>(this->size());\n    while (n >= 0\
    \ && !this->operator[](n).is_zero()) --n;\n    return n == -1 ? NEGATIVE_INFINITY\
    \ : n;\n  }\n  // order\n  int ord() const;\n  bool is_zero() const { return deg()\
    \ == NEGATIVE_INFINITY; }\n  void shrink() { this->resize(deg() + 1); }\n  truncated_formal_power_series\
    \ operator-() {\n    truncated_formal_power_series res(*this);\n    for (auto\
    \ &&i : res) i = -i;\n    return res;\n  }\n\n  truncated_formal_power_series\
    \ &operator+=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) += rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator-=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) -= rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator*=(const truncated_formal_power_series &rhs);\n  truncated_formal_power_series\
    \ inv(int n) const;\n  truncated_formal_power_series log(int n) const;\n  truncated_formal_power_series\
    \ exp(int n) const;\n  truncated_formal_power_series div(const truncated_formal_power_series\
    \ &rhs, int n) const;\n\n  friend truncated_formal_power_series operator+(const\
    \ truncated_formal_power_series &lhs,\n                                      \
    \           const truncated_formal_power_series &rhs) {\n    return truncated_formal_power_series(lhs)\
    \ += rhs;\n  }\n  friend truncated_formal_power_series operator-(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) -= rhs;\n  }\n  friend\
    \ truncated_formal_power_series operator*(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) *= rhs;\n  }\n\n  friend\
    \ std::istream &operator>>(std::istream &lhs, truncated_formal_power_series &rhs)\
    \ {\n    for (auto &&i : rhs) lhs >> i;\n    return lhs;\n  }\n  friend std::ostream\
    \ &operator<<(std::ostream &lhs, const truncated_formal_power_series &rhs) {\n\
    \    int s = 0, e = static_cast<int>(rhs.size());\n    lhs << '[';\n    for (auto\
    \ &&i : rhs) {\n      lhs << i;\n      if (s >= 1) lhs << 'x';\n      if (s >\
    \ 1) lhs << '^' << s;\n      if (++s != e) lhs << \" + \";\n    }\n    return\
    \ lhs << ']';\n  }\n};\n\ntemplate <typename IterT>\ntruncated_formal_power_series(IterT,\
    \ IterT)\n    -> truncated_formal_power_series<typename std::iterator_traits<IterT>::value_type>;\n\
    \ntemplate <typename ModIntT>\nusing tfps = truncated_formal_power_series<ModIntT>;\n\
    \ntemplate <typename ModIntT>\ntfps<ModIntT> &tfps<ModIntT>::operator*=(const\
    \ tfps<ModIntT> &rhs) {\n  // 6E\n  int n = static_cast<int>(this->size()), m\
    \ = static_cast<int>(rhs.size());\n  if (n == 0 || m == 0) {\n    this->clear();\n\
    \    return *this;\n  }\n  if (std::min(n, m) <= 32) {\n    tfps<ModIntT> res(n\
    \ + m - 1);\n    for (int i = 0; i != n; ++i)\n      for (int j = 0; j != m; ++j)\
    \ res[i + j] += this->operator[](i) * rhs[j];\n    return this->operator=(res);\n\
    \  }\n  int len = ntt_len(n + m - 1);\n  tfps<ModIntT> rhs_cpy(len);\n  std::copy_n(rhs.cbegin(),\
    \ m, rhs_cpy.begin());\n  this->resize(len);\n  dft_n(this->begin(), len), dft_n(rhs_cpy.begin(),\
    \ len);\n  for (int i = 0; i != len; ++i) this->operator[](i) *= rhs_cpy[i];\n\
    \  idft_n(this->begin(), len);\n  this->resize(n + m - 1);\n  return *this;\n\
    }\n\ntemplate <typename ModIntT>\ntfps<ModIntT> tfps<ModIntT>::inv(int n) const\
    \ {\n  // 10E\n  assert(n > 0);\n  assert(!this->front().is_zero());\n  if (n\
    \ == 1) return tfps<ModIntT>{this->front().inv()};\n  int len = ntt_len(n);\n\
    \  tfps<ModIntT> res(len), temp0(len), temp1(len), cpy(len);\n  std::copy(this->cbegin(),\
    \ this->cend(), cpy.begin());\n  res.front() = this->front().inv();\n  for (int\
    \ i = 2; i <= len; i <<= 1) {\n    std::copy_n(cpy.cbegin(), i, temp0.begin());\n\
    \    dft_n(temp0.begin(), i); // 2E\n    std::copy_n(res.cbegin(), i, temp1.begin());\n\
    \    dft_n(temp1.begin(), i); // 2E\n    for (int j = 0; j != i; ++j) temp0[j]\
    \ *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n    std::fill_n(temp0.begin(),\
    \ i >> 1, ModIntT());\n    dft_n(temp0.begin(), i); // 2E\n    for (int j = 0;\
    \ j != i; ++j) temp0[j] *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n  \
    \  for (int j = i >> 1; j != i; ++j) res[j] = -temp0[j];\n  }\n  res.resize(n);\n\
    \  return res;\n}\n\nLIB_END\n\n\n"
  code: "#ifndef TRUNCATED_FORMAL_POWER_SERIES_HPP\n#define TRUNCATED_FORMAL_POWER_SERIES_HPP\n\
    \n#include \"../common.hpp\"\n#include \"radix2_ntt.hpp\"\n\n#include <algorithm>\n\
    #include <cassert>\n#include <iostream>\n#include <iterator>\n#include <type_traits>\n\
    #include <vector>\n\nLIB_BEGIN\n\ntemplate <typename ModIntT>\nclass truncated_formal_power_series\
    \ : public std::vector<ModIntT> {\n  static_assert(std::is_same_v<typename std::vector<ModIntT>::value_type,\
    \ ModIntT>);\n\npublic:\n  using std::vector<ModIntT>::vector;\n\n  enum : int\
    \ { NEGATIVE_INFINITY = -1 };\n\n  // leading coefficient\n  ModIntT lc() const\
    \ {\n    int d = deg();\n    return d == NEGATIVE_INFINITY ? ModIntT() : this->operator[](d);\n\
    \  }\n  // degree\n  int deg() const {\n    // treat formal power series like\
    \ polynomials\n    int n = static_cast<int>(this->size());\n    while (n >= 0\
    \ && !this->operator[](n).is_zero()) --n;\n    return n == -1 ? NEGATIVE_INFINITY\
    \ : n;\n  }\n  // order\n  int ord() const;\n  bool is_zero() const { return deg()\
    \ == NEGATIVE_INFINITY; }\n  void shrink() { this->resize(deg() + 1); }\n  truncated_formal_power_series\
    \ operator-() {\n    truncated_formal_power_series res(*this);\n    for (auto\
    \ &&i : res) i = -i;\n    return res;\n  }\n\n  truncated_formal_power_series\
    \ &operator+=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) += rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator-=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) -= rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator*=(const truncated_formal_power_series &rhs);\n  truncated_formal_power_series\
    \ inv(int n) const;\n  truncated_formal_power_series log(int n) const;\n  truncated_formal_power_series\
    \ exp(int n) const;\n  truncated_formal_power_series div(const truncated_formal_power_series\
    \ &rhs, int n) const;\n\n  friend truncated_formal_power_series operator+(const\
    \ truncated_formal_power_series &lhs,\n                                      \
    \           const truncated_formal_power_series &rhs) {\n    return truncated_formal_power_series(lhs)\
    \ += rhs;\n  }\n  friend truncated_formal_power_series operator-(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) -= rhs;\n  }\n  friend\
    \ truncated_formal_power_series operator*(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) *= rhs;\n  }\n\n  friend\
    \ std::istream &operator>>(std::istream &lhs, truncated_formal_power_series &rhs)\
    \ {\n    for (auto &&i : rhs) lhs >> i;\n    return lhs;\n  }\n  friend std::ostream\
    \ &operator<<(std::ostream &lhs, const truncated_formal_power_series &rhs) {\n\
    \    int s = 0, e = static_cast<int>(rhs.size());\n    lhs << '[';\n    for (auto\
    \ &&i : rhs) {\n      lhs << i;\n      if (s >= 1) lhs << 'x';\n      if (s >\
    \ 1) lhs << '^' << s;\n      if (++s != e) lhs << \" + \";\n    }\n    return\
    \ lhs << ']';\n  }\n};\n\ntemplate <typename IterT>\ntruncated_formal_power_series(IterT,\
    \ IterT)\n    -> truncated_formal_power_series<typename std::iterator_traits<IterT>::value_type>;\n\
    \ntemplate <typename ModIntT>\nusing tfps = truncated_formal_power_series<ModIntT>;\n\
    \ntemplate <typename ModIntT>\ntfps<ModIntT> &tfps<ModIntT>::operator*=(const\
    \ tfps<ModIntT> &rhs) {\n  // 6E\n  int n = static_cast<int>(this->size()), m\
    \ = static_cast<int>(rhs.size());\n  if (n == 0 || m == 0) {\n    this->clear();\n\
    \    return *this;\n  }\n  if (std::min(n, m) <= 32) {\n    tfps<ModIntT> res(n\
    \ + m - 1);\n    for (int i = 0; i != n; ++i)\n      for (int j = 0; j != m; ++j)\
    \ res[i + j] += this->operator[](i) * rhs[j];\n    return this->operator=(res);\n\
    \  }\n  int len = ntt_len(n + m - 1);\n  tfps<ModIntT> rhs_cpy(len);\n  std::copy_n(rhs.cbegin(),\
    \ m, rhs_cpy.begin());\n  this->resize(len);\n  dft_n(this->begin(), len), dft_n(rhs_cpy.begin(),\
    \ len);\n  for (int i = 0; i != len; ++i) this->operator[](i) *= rhs_cpy[i];\n\
    \  idft_n(this->begin(), len);\n  this->resize(n + m - 1);\n  return *this;\n\
    }\n\ntemplate <typename ModIntT>\ntfps<ModIntT> tfps<ModIntT>::inv(int n) const\
    \ {\n  // 10E\n  assert(n > 0);\n  assert(!this->front().is_zero());\n  if (n\
    \ == 1) return tfps<ModIntT>{this->front().inv()};\n  int len = ntt_len(n);\n\
    \  tfps<ModIntT> res(len), temp0(len), temp1(len), cpy(len);\n  std::copy(this->cbegin(),\
    \ this->cend(), cpy.begin());\n  res.front() = this->front().inv();\n  for (int\
    \ i = 2; i <= len; i <<= 1) {\n    std::copy_n(cpy.cbegin(), i, temp0.begin());\n\
    \    dft_n(temp0.begin(), i); // 2E\n    std::copy_n(res.cbegin(), i, temp1.begin());\n\
    \    dft_n(temp1.begin(), i); // 2E\n    for (int j = 0; j != i; ++j) temp0[j]\
    \ *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n    std::fill_n(temp0.begin(),\
    \ i >> 1, ModIntT());\n    dft_n(temp0.begin(), i); // 2E\n    for (int j = 0;\
    \ j != i; ++j) temp0[j] *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n  \
    \  for (int j = i >> 1; j != i; ++j) res[j] = -temp0[j];\n  }\n  res.resize(n);\n\
    \  return res;\n}\n\nLIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  - math/radix2_ntt.hpp
  isVerificationFile: false
  path: math/truncated_formal_power_series.hpp
  requiredBy: []
  timestamp: '2022-04-20 19:56:42+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.1.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
documentation_of: math/truncated_formal_power_series.hpp
layout: document
title: Truncated Formal Power Series
---
