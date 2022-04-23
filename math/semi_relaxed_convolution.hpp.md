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
    path: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/semi_relaxed_convolution.hpp\"\n\n\n\n#line 1 \"common.hpp\"\
    \n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n#define LIB_END\
    \ }\n#define LIB ::lib::\n\n\n#line 1 \"math/radix2_ntt.hpp\"\n\n\n\n#line 5 \"\
    math/radix2_ntt.hpp\"\n\n#include <array>\n#include <cassert>\n#include <type_traits>\n\
    #include <vector>\n\nLIB_BEGIN\n\nnamespace detail {\n\ntemplate <typename IntT>\n\
    constexpr std::enable_if_t<std::is_integral_v<IntT>, int> bsf(IntT v) {\n  if\
    \ (static_cast<std::make_signed_t<IntT>>(v) <= 0) return -1;\n  int res = 0;\n\
    \  for (; (v & 1) == 0; ++res) v >>= 1;\n  return res;\n}\n\ntemplate <typename\
    \ ModIntT>\nconstexpr ModIntT quadratic_nonresidue_prime() {\n  auto mod = ModIntT::mod();\n\
    \  for (int i = 2;; ++i)\n    if (ModIntT(i).pow(mod >> 1) == mod - 1) return\
    \ ModIntT(i);\n}\n\ntemplate <typename ModIntT>\nconstexpr ModIntT gen_of_sylow2subgroup()\
    \ {\n  auto mod = ModIntT::mod();\n  return quadratic_nonresidue_prime<ModIntT>().pow(mod\
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
    \ + l] * iv);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n}\n\ntemplate <typename\
    \ ContainerT>\nvoid dft(ContainerT &a) {\n  dft_n(a.begin(), a.size());\n}\n\n\
    template <typename ContainerT>\nvoid idft(ContainerT &a) {\n  idft_n(a.begin(),\
    \ a.size());\n}\n\nLIB_END\n\n\n#line 6 \"math/semi_relaxed_convolution.hpp\"\n\
    \n#include <algorithm>\n#line 9 \"math/semi_relaxed_convolution.hpp\"\n#include\
    \ <utility>\n#line 11 \"math/semi_relaxed_convolution.hpp\"\n\nLIB_BEGIN\n\ntemplate\
    \ <typename ModIntT, typename FnT>\nclass semi_relaxed_convolution {\n  std::vector<ModIntT>\
    \ fixed_A_{}, B_{}, c_{};\n  std::vector<std::vector<std::vector<ModIntT>>> dft_A_cache_{},\
    \ dft_B_cache_{};\n  int n_{};\n  FnT handle_;\n\n  enum : int { BASE_CASE_SIZE\
    \ = 32, LOG_BLOCK = 4, BLOCK = 1 << LOG_BLOCK, MASK = BLOCK - 1 };\n\n  static_assert((BASE_CASE_SIZE\
    \ & (BASE_CASE_SIZE - 1)) == 0);\n  static_assert(std::is_invocable_r_v<ModIntT,\
    \ FnT, int, const std::vector<ModIntT> &> ||\n                std::is_invocable_r_v<ModIntT,\
    \ FnT, int>);\n\npublic:\n  semi_relaxed_convolution(const std::vector<ModIntT>\
    \ &A, FnT &&handle)\n      : fixed_A_(A), c_(1024), handle_(std::forward<FnT>(handle))\
    \ {}\n\n  const std::vector<ModIntT> &get_multiplier() const { return B_; }\n\
    \  const std::vector<ModIntT> &get_multiplicand() const { return fixed_A_; }\n\
    \  semi_relaxed_convolution &await(int k) {\n    while (n_ < k) next();\n    return\
    \ *this;\n  }\n  ModIntT at(int k) {\n    while (n_ <= k) next();\n    return\
    \ c_[k];\n  }\n  ModIntT operator[](int k) { return at(k); }\n  ModIntT next();\n\
    };\n\ntemplate <typename ModIntT, typename FnT>\nModIntT semi_relaxed_convolution<ModIntT,\
    \ FnT>::next() {\n  {\n    // enlarge space\n    int len = ntt_len(n_ << 1 | 1);\n\
    \    if (static_cast<int>(c_.size()) < len) c_.resize(len);\n    if (static_cast<int>(fixed_A_.size())\
    \ < len) fixed_A_.resize(len);\n  }\n  if ((n_ & (BASE_CASE_SIZE - 1)) == 0) {\n\
    \    for (int t = n_ / BASE_CASE_SIZE, block_size = BASE_CASE_SIZE, lv = 0; t\
    \ != 0;\n         t >>= LOG_BLOCK, block_size <<= LOG_BLOCK, ++lv) {\n      if\
    \ (int i = t & MASK, block_size2 = block_size << 1, l = n_ - block_size; i !=\
    \ 0) {\n        if (block_size * i == n_) {\n          if (static_cast<int>(dft_A_cache_.size())\
    \ == lv) {\n            dft_A_cache_.emplace_back();\n            dft_B_cache_.emplace_back(BLOCK\
    \ - 1);\n          }\n          dft(dft_A_cache_[lv].emplace_back(fixed_A_.begin()\
    \ + (i - 1) * block_size,\n                                            fixed_A_.begin()\
    \ + (i + 1) * block_size));\n        }\n        auto &B_cache = dft_B_cache_[lv];\n\
    \        B_cache[i - 1].resize(block_size2);\n        std::fill_n(std::copy_n(B_.begin()\
    \ + l, block_size, B_cache[i - 1].begin()), block_size,\n                    ModIntT());\n\
    \        dft(B_cache[i - 1]);\n        std::vector<ModIntT> temp_sum(block_size2);\n\
    \        for (int j = 0; j != i; ++j)\n          for (int k = 0; k != block_size2;\
    \ ++k)\n            temp_sum[k] += dft_A_cache_[lv][i - 1 - j][k] * B_cache[j][k];\n\
    \        idft(temp_sum);\n        for (int j = block_size; j != block_size2; ++j)\
    \ c_[j + n_ - block_size] += temp_sum[j];\n        break;\n      }\n    }\n  }\n\
    \  for (int i = 0, l = n_ & ~(BASE_CASE_SIZE - 1); i < n_ - l; ++i)\n    c_[n_]\
    \ += fixed_A_[n_ - l - i] * B_[l + i];\n  if constexpr (std::is_invocable_r_v<ModIntT,\
    \ FnT, int, const std::vector<ModIntT> &>) {\n    c_[n_] += fixed_A_.front() *\
    \ B_.emplace_back(handle_(n_, c_));\n  } else {\n    c_[n_] += fixed_A_.front()\
    \ * B_.emplace_back(handle_(n_));\n  }\n  return c_[n_++];\n}\n\nLIB_END\n\n\n"
  code: "#ifndef SEMI_RELAXED_CONVOLUTION_HPP\n#define SEMI_RELAXED_CONVOLUTION_HPP\n\
    \n#include \"../common.hpp\"\n#include \"radix2_ntt.hpp\"\n\n#include <algorithm>\n\
    #include <type_traits>\n#include <utility>\n#include <vector>\n\nLIB_BEGIN\n\n\
    template <typename ModIntT, typename FnT>\nclass semi_relaxed_convolution {\n\
    \  std::vector<ModIntT> fixed_A_{}, B_{}, c_{};\n  std::vector<std::vector<std::vector<ModIntT>>>\
    \ dft_A_cache_{}, dft_B_cache_{};\n  int n_{};\n  FnT handle_;\n\n  enum : int\
    \ { BASE_CASE_SIZE = 32, LOG_BLOCK = 4, BLOCK = 1 << LOG_BLOCK, MASK = BLOCK -\
    \ 1 };\n\n  static_assert((BASE_CASE_SIZE & (BASE_CASE_SIZE - 1)) == 0);\n  static_assert(std::is_invocable_r_v<ModIntT,\
    \ FnT, int, const std::vector<ModIntT> &> ||\n                std::is_invocable_r_v<ModIntT,\
    \ FnT, int>);\n\npublic:\n  semi_relaxed_convolution(const std::vector<ModIntT>\
    \ &A, FnT &&handle)\n      : fixed_A_(A), c_(1024), handle_(std::forward<FnT>(handle))\
    \ {}\n\n  const std::vector<ModIntT> &get_multiplier() const { return B_; }\n\
    \  const std::vector<ModIntT> &get_multiplicand() const { return fixed_A_; }\n\
    \  semi_relaxed_convolution &await(int k) {\n    while (n_ < k) next();\n    return\
    \ *this;\n  }\n  ModIntT at(int k) {\n    while (n_ <= k) next();\n    return\
    \ c_[k];\n  }\n  ModIntT operator[](int k) { return at(k); }\n  ModIntT next();\n\
    };\n\ntemplate <typename ModIntT, typename FnT>\nModIntT semi_relaxed_convolution<ModIntT,\
    \ FnT>::next() {\n  {\n    // enlarge space\n    int len = ntt_len(n_ << 1 | 1);\n\
    \    if (static_cast<int>(c_.size()) < len) c_.resize(len);\n    if (static_cast<int>(fixed_A_.size())\
    \ < len) fixed_A_.resize(len);\n  }\n  if ((n_ & (BASE_CASE_SIZE - 1)) == 0) {\n\
    \    for (int t = n_ / BASE_CASE_SIZE, block_size = BASE_CASE_SIZE, lv = 0; t\
    \ != 0;\n         t >>= LOG_BLOCK, block_size <<= LOG_BLOCK, ++lv) {\n      if\
    \ (int i = t & MASK, block_size2 = block_size << 1, l = n_ - block_size; i !=\
    \ 0) {\n        if (block_size * i == n_) {\n          if (static_cast<int>(dft_A_cache_.size())\
    \ == lv) {\n            dft_A_cache_.emplace_back();\n            dft_B_cache_.emplace_back(BLOCK\
    \ - 1);\n          }\n          dft(dft_A_cache_[lv].emplace_back(fixed_A_.begin()\
    \ + (i - 1) * block_size,\n                                            fixed_A_.begin()\
    \ + (i + 1) * block_size));\n        }\n        auto &B_cache = dft_B_cache_[lv];\n\
    \        B_cache[i - 1].resize(block_size2);\n        std::fill_n(std::copy_n(B_.begin()\
    \ + l, block_size, B_cache[i - 1].begin()), block_size,\n                    ModIntT());\n\
    \        dft(B_cache[i - 1]);\n        std::vector<ModIntT> temp_sum(block_size2);\n\
    \        for (int j = 0; j != i; ++j)\n          for (int k = 0; k != block_size2;\
    \ ++k)\n            temp_sum[k] += dft_A_cache_[lv][i - 1 - j][k] * B_cache[j][k];\n\
    \        idft(temp_sum);\n        for (int j = block_size; j != block_size2; ++j)\
    \ c_[j + n_ - block_size] += temp_sum[j];\n        break;\n      }\n    }\n  }\n\
    \  for (int i = 0, l = n_ & ~(BASE_CASE_SIZE - 1); i < n_ - l; ++i)\n    c_[n_]\
    \ += fixed_A_[n_ - l - i] * B_[l + i];\n  if constexpr (std::is_invocable_r_v<ModIntT,\
    \ FnT, int, const std::vector<ModIntT> &>) {\n    c_[n_] += fixed_A_.front() *\
    \ B_.emplace_back(handle_(n_, c_));\n  } else {\n    c_[n_] += fixed_A_.front()\
    \ * B_.emplace_back(handle_(n_));\n  }\n  return c_[n_++];\n}\n\nLIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  - math/radix2_ntt.hpp
  isVerificationFile: false
  path: math/semi_relaxed_convolution.hpp
  requiredBy: []
  timestamp: '2022-04-23 15:00:57+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
documentation_of: math/semi_relaxed_convolution.hpp
layout: document
title: Semi-Relaxed Convolution
---
