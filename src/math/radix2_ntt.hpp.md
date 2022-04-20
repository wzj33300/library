---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: src/math/convolution.hpp
    title: Convolution
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
  code: "#ifndef RADIX2_NTT_HPP\r\n#define RADIX2_NTT_HPP\r\n\r\n#include \"common.hpp\"\
    \r\n\r\n#include <array>\r\n#include <cassert>\r\n#include <type_traits>\r\n#include\
    \ <vector>\r\n\r\nLIB_BEGIN\r\n\r\nnamespace detail {\r\n\r\ntemplate <typename\
    \ IntT>\r\nconstexpr std::enable_if_t<std::is_integral_v<IntT>, int> bsf(IntT\
    \ v) {\r\n  if (static_cast<std::make_signed_t<IntT>>(v) <= 0) return -1;\r\n\
    \  int res = 0;\r\n  for (; (v & 1) == 0; ++res) v >>= 1;\r\n  return res;\r\n\
    }\r\n\r\ntemplate <typename ModIntT>\r\nconstexpr ModIntT quadratic_nonresidue_prime()\
    \ {\r\n  auto mod = ModIntT::mod();\r\n  for (int i = 2;; ++i)\r\n    if (ModIntT(i).pow(mod\
    \ >> 1) == mod - 1) return ModIntT(i);\r\n}\r\n\r\ntemplate <typename ModIntT>\r\
    \nconstexpr ModIntT gen_of_sylow2subgroup() {\r\n  auto mod = ModIntT::mod();\r\
    \n  return quadratic_nonresidue_prime<ModIntT>().pow(mod >> bsf(mod - 1));\r\n\
    }\r\n\r\ntemplate <typename ModIntT>\r\nconstexpr std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> root() {\r\n  std::array<ModIntT, bsf(ModIntT::mod() - 1) - 1> rt;\
    \ // order(`rt[i]`) = 2^(i + 2).\r\n  rt.back() = gen_of_sylow2subgroup<ModIntT>();\r\
    \n  for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) rt[i] = rt[i + 1] *\
    \ rt[i + 1];\r\n  return rt;\r\n}\r\n\r\ntemplate <typename ModIntT>\r\nconstexpr\
    \ std::array<ModIntT, bsf(ModIntT::mod() - 1) - 1> iroot() {\r\n  std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> irt;\r\n  irt.back() = gen_of_sylow2subgroup<ModIntT>().inv();\r\
    \n  for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) irt[i] = irt[i + 1]\
    \ * irt[i + 1];\r\n  return irt;\r\n}\r\n\r\n} // namespace detail\r\n\r\n// Input:\
    \  integer `n`.\r\n// Output: 2^(\u2308log_2(`n`)\u2309).\r\nint ntt_len(int n)\
    \ {\r\n  --n;\r\n  n |= n >> 1;\r\n  n |= n >> 2;\r\n  n |= n >> 4;\r\n  n |=\
    \ n >> 8;\r\n  return (n | n >> 16) + 1;\r\n}\r\n\r\n// Input:           f(x)\
    \ = `a[0]` + `a[1]`x + ... + `a[n - 1]`x^(`n` - 1) where `n` is power of 2.\r\n\
    // Output(inplace): reversed binary permutation of [f(\u03B6^0), f(\u03B6), f(\u03B6\
    ^2), ..., f(\u03B6^(`n` - 1))].\r\ntemplate <typename IterT>\r\nvoid dft_n(IterT\
    \ a, int n) {\r\n  assert((n & (n - 1)) == 0);\r\n  using T                  =\
    \ typename std::iterator_traits<IterT>::value_type;\r\n  static constexpr auto\
    \ rt = detail::root<T>();\r\n  static std::vector<T> root(1);\r\n  if (int s =\
    \ static_cast<int>(root.size()); s << 1 < n) {\r\n    root.resize(n >> 1);\r\n\
    \    for (int i = detail::bsf(s); (1 << i) < (n >> 1); ++i) {\r\n      int j \
    \  = 1 << i;\r\n      root[j] = rt[i];\r\n      for (int k = j + 1; k < (j <<\
    \ 1); ++k) root[k] = root[k - j] * root[j];\r\n    }\r\n  }\r\n  for (int j =\
    \ 0, l = n >> 1; j != l; ++j) {\r\n    T u(a[j]), v(a[j + l]);\r\n    a[j] = u\
    \ + v, a[j + l] = u - v;\r\n  }\r\n  for (int i = n >> 1; i >= 2; i >>= 1) {\r\
    \n    for (int j = 0, l = i >> 1; j != l; ++j) {\r\n      T u(a[j]), v(a[j + l]);\r\
    \n      a[j] = u + v, a[j + l] = u - v;\r\n    }\r\n    for (int j = i, l = i\
    \ >> 1, m = 1; j != n; j += i, ++m) {\r\n      for (int k = j; k != j + l; ++k)\
    \ {\r\n        T u(a[k]), v(a[k + l] * root[m]);\r\n        a[k] = u + v, a[k\
    \ + l] = u - v;\r\n      }\r\n    }\r\n  }\r\n}\r\n\r\n// Input:           reversed\
    \ binary permutation of [f(\u03B6^0), f(\u03B6), f(\u03B6^2), ..., f(\u03B6^(`n`\
    \ - 1))].\r\n// Output(inplace): f(x) = `a[0]` + `a[1]`x + ... + `a[n - 1]`x^(`n`\
    \ - 1) where `n` is power of 2.\r\ntemplate <typename IterT>\r\nvoid idft_n(IterT\
    \ a, int n) {\r\n  assert((n & (n - 1)) == 0);\r\n  using T                  =\
    \ typename std::iterator_traits<IterT>::value_type;\r\n  static constexpr auto\
    \ rt = detail::iroot<T>();\r\n  static std::vector<T> root(1);\r\n  if (int s\
    \ = static_cast<int>(root.size()); s << 1 < n) {\r\n    root.resize(n >> 1);\r\
    \n    for (int i = detail::bsf(s); (1 << i) < (n >> 1); ++i) {\r\n      int j\
    \   = 1 << i;\r\n      root[j] = rt[i];\r\n      for (int k = j + 1; k < (j <<\
    \ 1); ++k) root[k] = root[k - j] * root[j];\r\n    }\r\n  }\r\n  for (int i =\
    \ 2; i < n; i <<= 1) {\r\n    for (int j = 0, l = i >> 1; j != l; ++j) {\r\n \
    \     T u(a[j]), v(a[j + l]);\r\n      a[j] = u + v, a[j + l] = u - v;\r\n   \
    \ }\r\n    for (int j = i, l = i >> 1, m = 1; j != n; j += i, ++m) {\r\n     \
    \ for (int k = j; k != j + l; ++k) {\r\n        T u(a[k]), v(a[k + l]);\r\n  \
    \      a[k] = u + v, a[k + l] = (u - v) * root[m];\r\n      }\r\n    }\r\n  }\r\
    \n  const T iv(T::mod() - T::mod() / n);\r\n  for (int j = 0, l = n >> 1; j !=\
    \ l; ++j) {\r\n    T u(a[j] * iv), v(a[j + l] * iv);\r\n    a[j] = u + v, a[j\
    \ + l] = u - v;\r\n  }\r\n}\r\n\r\nLIB_END\r\n\r\n#endif"
  dependsOn:
  - src/common.hpp
  isVerificationFile: false
  path: src/math/radix2_ntt.hpp
  requiredBy:
  - src/math/convolution.hpp
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
documentation_of: src/math/radix2_ntt.hpp
layout: document
title: Radix-2 NTT
---
