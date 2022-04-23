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
    path: remote_test/yosupo/math/convolution_mod.3.test.cpp
    title: remote_test/yosupo/math/convolution_mod.3.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/relaxed_convolution.hpp\"\n\n\n\n#line 1 \"common.hpp\"\
    \n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n#define LIB_END\
    \ }\n#define LIB ::lib::\n\n\n#line 1 \"math/radix2_ntt.hpp\"\n\n\n\n#line 5 \"\
    math/radix2_ntt.hpp\"\n\n#include <array>\n#include <cassert>\n#include <type_traits>\n\
    #include <vector>\n\nLIB_BEGIN\n\nnamespace detail {\n\ntemplate <typename IntT>\n\
    constexpr std::enable_if_t<std::is_integral_v<IntT>, int> bsf(IntT v) {\n  if\
    \ (static_cast<std::make_signed_t<IntT>>(v) <= 0) return -1;\n  int res = 0;\n\
    \  for (; (v & 1) == 0; ++res) v >>= 1;\n  return res;\n}\n\ntemplate <typename\
    \ ModIntT>\nconstexpr ModIntT quadratic_nonresidue_prime() {\n  auto mod = ModIntT::mod();\n\
    \  for (int i = 2;; ++i)\n    if (ModIntT(i).pow(mod >> 1) == mod - 1) return\
    \ ModIntT(i);\n}\n\ntemplate <typename ModIntT>\nconstexpr ModIntT gen_of_sylow_2_subgroup()\
    \ {\n  auto mod = ModIntT::mod();\n  return quadratic_nonresidue_prime<ModIntT>().pow(mod\
    \ >> bsf(mod - 1));\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> root() {\n  std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> rt; // order(`rt[i]`) = 2^(i + 2).\n  rt.back() = gen_of_sylow_2_subgroup<ModIntT>();\n\
    \  for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) rt[i] = rt[i + 1] *\
    \ rt[i + 1];\n  return rt;\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
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
    \ a.size());\n}\n\nLIB_END\n\n\n#line 6 \"math/relaxed_convolution.hpp\"\n\n#include\
    \ <functional>\n#line 10 \"math/relaxed_convolution.hpp\"\n\nLIB_BEGIN\n\ntemplate\
    \ <typename ModIntT>\nclass relaxed_convolution {                       // O(n\
    \ log^2 n) impl\n  std::vector<ModIntT> a_{}, b_{}, c_{};          // `a_ * b_`\
    \ = `c_`\n  std::vector<std::vector<ModIntT>> ac_{}, bc_{}; // cached DFTs\n \
    \ std::function<ModIntT()> ha_{}, hb_{};          // handle for `a` and `b`\n\
    \  int n_{};                                       // counter\n\n  enum : int\
    \ { BASE_CASE_SIZE = 32 };\n\n  template <typename FnT>\n  static auto wrap(FnT\
    \ &&f, int &n, const std::vector<ModIntT> &c, std::vector<ModIntT> &e) {\n   \
    \ if constexpr (std::is_invocable_r_v<ModIntT, FnT, int, const std::vector<ModIntT>\
    \ &>) {\n      return std::bind(\n          [f](int n, const std::vector<ModIntT>\
    \ &c, std::vector<ModIntT> &e) mutable {\n            return ModIntT(e.emplace_back(f(n,\
    \ c)));\n          },\n          std::cref(n), std::cref(c), std::ref(e));\n \
    \   } else if constexpr (std::is_invocable_r_v<ModIntT, FnT, int>) {\n      return\
    \ std::bind(\n          [f](int n, std::vector<ModIntT> &e) mutable { return ModIntT(e.emplace_back(f(n)));\
    \ },\n          std::cref(n), std::ref(e));\n    } else if constexpr (std::is_invocable_r_v<ModIntT,\
    \ FnT>) {\n      return std::bind(\n          [f](std::vector<ModIntT> &e) mutable\
    \ { return ModIntT(e.emplace_back(f())); },\n          std::ref(e));\n    } else\
    \ {\n      throw;\n    }\n  }\n\npublic:\n  // `h0` multiplicand, `h1` multiplier\n\
    \  template <typename Fn0T, typename Fn1T>\n  relaxed_convolution(Fn0T &&h0, Fn1T\
    \ &&h1)\n      : c_(4), ha_(wrap(h0, n_, c_, a_)), hb_(wrap(h1, n_, c_, b_)) {}\n\
    \  const std::vector<ModIntT> &get_multiplicand() const { return a_; }\n  const\
    \ std::vector<ModIntT> &get_multiplier() const { return b_; }\n  ModIntT at(int\
    \ k) {\n    while (n_ <= k) next();\n    return c_[k];\n  }\n  ModIntT operator[](int\
    \ k) { return at(k); }\n  ModIntT next();\n};\n\ntemplate <typename ModIntT>\n\
    ModIntT relaxed_convolution<ModIntT>::next() {\n  {\n    // enlarge space\n  \
    \  int l = ntt_len(n_ << 1 | 1);\n    if (static_cast<int>(c_.size()) < l) c_.resize(l);\n\
    \  }\n  switch (n_) {\n  case 0: c_[0] = ha_() * hb_(); break;\n  case 1:\n  \
    \  c_[1] = ha_() * b_.front() + a_.front() * hb_();\n    c_[2] = a_[1] * b_[1];\n\
    \    break;\n  case 2:\n    c_[2] += ha_() * b_.front() + a_.front() * hb_();\n\
    \    c_[3] = a_[2] * b_[1] + a_[1] * b_[2];\n    break;\n  default:\n    if ((n_\
    \ & (n_ - 1)) == 0) {\n      int t0 = n_ >> 1, t1 = n_;\n      auto &c0 = ac_.emplace_back(a_.begin()\
    \ + t0, a_.begin() + t1);\n      auto &c1 = bc_.emplace_back(b_.begin() + t0,\
    \ b_.begin() + t1);\n      c0.resize(t1);\n      c1.resize(t1);\n      dft(c0),\
    \ dft(c1);\n      std::vector res(c0);\n      for (int i = 0; i < t1; ++i) res[i]\
    \ *= c1[i];\n      idft(res);\n      for (int i = 0; i < t1 - 1; ++i) c_[t1 +\
    \ i] += res[i];\n    }\n    c_[n_] += ha_() * b_.front() + a_.front() * hb_();\n\
    \    c_[n_ + 1] += a_[1] * b_.back() + a_.back() * b_[1];\n    for (int sft =\
    \ 1, offset = ntt_len(n_ + 1) >> 1, t = n_ + 1 - offset;\n         (t & 1) ==\
    \ 0 && 1 << sft < offset; ++sft, t >>= 1)\n      if (1 << sft <= BASE_CASE_SIZE)\
    \ {\n        for (int i = 0, m = n_ + 1 - (1 << sft); i != 1 << sft; ++i)\n  \
    \        for (int j = 0; j != 1 << sft; ++j)\n            c_[n_ + 1 + i + j] +=\
    \ a_[m + i] * b_[j + (1 << sft)] + a_[j + (1 << sft)] * b_[m + i];\n      } else\
    \ {\n        std::vector c0(a_.begin() + n_ + 1 - (1 << sft), a_.begin() + n_\
    \ + 1);\n        std::vector c1(b_.begin() + n_ + 1 - (1 << sft), b_.begin() +\
    \ n_ + 1);\n        c0.resize(2 << sft);\n        c1.resize(2 << sft);\n     \
    \   dft(c0), dft(c1);\n        for (int i = 0; i != 2 << sft; ++i)\n         \
    \ c0[i] = c0[i] * bc_[sft - 1][i] + c1[i] * ac_[sft - 1][i];\n        idft(c0);\n\
    \        for (int i = 0; i != (2 << sft) - 1; ++i) c_[n_ + 1 + i] += c0[i];\n\
    \      }\n  }\n  return c_[n_++];\n}\n\nLIB_END\n\n\n"
  code: "#ifndef RELAXED_CONVOLUTION_HPP\n#define RELAXED_CONVOLUTION_HPP\n\n#include\
    \ \"../common.hpp\"\n#include \"radix2_ntt.hpp\"\n\n#include <functional>\n#include\
    \ <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\ntemplate <typename ModIntT>\n\
    class relaxed_convolution {                       // O(n log^2 n) impl\n  std::vector<ModIntT>\
    \ a_{}, b_{}, c_{};          // `a_ * b_` = `c_`\n  std::vector<std::vector<ModIntT>>\
    \ ac_{}, bc_{}; // cached DFTs\n  std::function<ModIntT()> ha_{}, hb_{};     \
    \     // handle for `a` and `b`\n  int n_{};                                 \
    \      // counter\n\n  enum : int { BASE_CASE_SIZE = 32 };\n\n  template <typename\
    \ FnT>\n  static auto wrap(FnT &&f, int &n, const std::vector<ModIntT> &c, std::vector<ModIntT>\
    \ &e) {\n    if constexpr (std::is_invocable_r_v<ModIntT, FnT, int, const std::vector<ModIntT>\
    \ &>) {\n      return std::bind(\n          [f](int n, const std::vector<ModIntT>\
    \ &c, std::vector<ModIntT> &e) mutable {\n            return ModIntT(e.emplace_back(f(n,\
    \ c)));\n          },\n          std::cref(n), std::cref(c), std::ref(e));\n \
    \   } else if constexpr (std::is_invocable_r_v<ModIntT, FnT, int>) {\n      return\
    \ std::bind(\n          [f](int n, std::vector<ModIntT> &e) mutable { return ModIntT(e.emplace_back(f(n)));\
    \ },\n          std::cref(n), std::ref(e));\n    } else if constexpr (std::is_invocable_r_v<ModIntT,\
    \ FnT>) {\n      return std::bind(\n          [f](std::vector<ModIntT> &e) mutable\
    \ { return ModIntT(e.emplace_back(f())); },\n          std::ref(e));\n    } else\
    \ {\n      throw;\n    }\n  }\n\npublic:\n  // `h0` multiplicand, `h1` multiplier\n\
    \  template <typename Fn0T, typename Fn1T>\n  relaxed_convolution(Fn0T &&h0, Fn1T\
    \ &&h1)\n      : c_(4), ha_(wrap(h0, n_, c_, a_)), hb_(wrap(h1, n_, c_, b_)) {}\n\
    \  const std::vector<ModIntT> &get_multiplicand() const { return a_; }\n  const\
    \ std::vector<ModIntT> &get_multiplier() const { return b_; }\n  ModIntT at(int\
    \ k) {\n    while (n_ <= k) next();\n    return c_[k];\n  }\n  ModIntT operator[](int\
    \ k) { return at(k); }\n  ModIntT next();\n};\n\ntemplate <typename ModIntT>\n\
    ModIntT relaxed_convolution<ModIntT>::next() {\n  {\n    // enlarge space\n  \
    \  int l = ntt_len(n_ << 1 | 1);\n    if (static_cast<int>(c_.size()) < l) c_.resize(l);\n\
    \  }\n  switch (n_) {\n  case 0: c_[0] = ha_() * hb_(); break;\n  case 1:\n  \
    \  c_[1] = ha_() * b_.front() + a_.front() * hb_();\n    c_[2] = a_[1] * b_[1];\n\
    \    break;\n  case 2:\n    c_[2] += ha_() * b_.front() + a_.front() * hb_();\n\
    \    c_[3] = a_[2] * b_[1] + a_[1] * b_[2];\n    break;\n  default:\n    if ((n_\
    \ & (n_ - 1)) == 0) {\n      int t0 = n_ >> 1, t1 = n_;\n      auto &c0 = ac_.emplace_back(a_.begin()\
    \ + t0, a_.begin() + t1);\n      auto &c1 = bc_.emplace_back(b_.begin() + t0,\
    \ b_.begin() + t1);\n      c0.resize(t1);\n      c1.resize(t1);\n      dft(c0),\
    \ dft(c1);\n      std::vector res(c0);\n      for (int i = 0; i < t1; ++i) res[i]\
    \ *= c1[i];\n      idft(res);\n      for (int i = 0; i < t1 - 1; ++i) c_[t1 +\
    \ i] += res[i];\n    }\n    c_[n_] += ha_() * b_.front() + a_.front() * hb_();\n\
    \    c_[n_ + 1] += a_[1] * b_.back() + a_.back() * b_[1];\n    for (int sft =\
    \ 1, offset = ntt_len(n_ + 1) >> 1, t = n_ + 1 - offset;\n         (t & 1) ==\
    \ 0 && 1 << sft < offset; ++sft, t >>= 1)\n      if (1 << sft <= BASE_CASE_SIZE)\
    \ {\n        for (int i = 0, m = n_ + 1 - (1 << sft); i != 1 << sft; ++i)\n  \
    \        for (int j = 0; j != 1 << sft; ++j)\n            c_[n_ + 1 + i + j] +=\
    \ a_[m + i] * b_[j + (1 << sft)] + a_[j + (1 << sft)] * b_[m + i];\n      } else\
    \ {\n        std::vector c0(a_.begin() + n_ + 1 - (1 << sft), a_.begin() + n_\
    \ + 1);\n        std::vector c1(b_.begin() + n_ + 1 - (1 << sft), b_.begin() +\
    \ n_ + 1);\n        c0.resize(2 << sft);\n        c1.resize(2 << sft);\n     \
    \   dft(c0), dft(c1);\n        for (int i = 0; i != 2 << sft; ++i)\n         \
    \ c0[i] = c0[i] * bc_[sft - 1][i] + c1[i] * ac_[sft - 1][i];\n        idft(c0);\n\
    \        for (int i = 0; i != (2 << sft) - 1; ++i) c_[n_ + 1 + i] += c0[i];\n\
    \      }\n  }\n  return c_[n_++];\n}\n\nLIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  - math/radix2_ntt.hpp
  isVerificationFile: false
  path: math/relaxed_convolution.hpp
  requiredBy: []
  timestamp: '2022-04-23 22:52:36+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.3.test.cpp
documentation_of: math/relaxed_convolution.hpp
layout: document
title: Relaxed Convolution
---
