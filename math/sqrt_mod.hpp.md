---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: common.hpp
    title: common.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/sqrt_mod.0.test.cpp
    title: remote_test/yosupo/math/sqrt_mod.0.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/sqrt_mod.hpp\"\n\n\n\n#line 1 \"common.hpp\"\n\n\n\n\
    #define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n#define LIB_END }\n#define\
    \ LIB ::lib::\n\n\n#line 5 \"math/sqrt_mod.hpp\"\n\n#include <random>\n#include\
    \ <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\ntemplate <typename ModIntT>\n\
    std::vector<ModIntT> sqrt_mod_prime(ModIntT a) {\n  // Bostan-Mori's algorithm\n\
    \  if (a.is_zero()) return {a};\n  const auto p = ModIntT::mod();\n  if (a.pow(p\
    \ >> 1) == -1) return {};\n  if ((p & 3) == 3) {\n    ModIntT b(a.pow((p + 1)\
    \ >> 2));\n    return {b, -b};\n  }\n  std::mt19937 gen(std::random_device{}());\n\
    \  std::uniform_int_distribution<std::remove_cv_t<decltype(p)>> dis(2, p - 1);\n\
    \  ModIntT t;\n  do { t = dis(gen); } while ((t * t - 4 * a).pow(p >> 1) != -1);\n\
    \  ModIntT k0(1), k1, k2(-t), k3(a);\n  for (auto e = (p + 1) >> 1;;) {\n    //\
    \ clang-format off\n    if (e & 1) k0 = k1 - k0 * k2, k1 *= k3;\n    else k1 =\
    \ k0 * k3 - k1 * k2;\n    // clang-format on\n    if ((e >>= 1) == 0) return {k0,\
    \ -k0};\n    k2 = k3 + k3 - k2 * k2, k3 *= k3;\n  }\n}\n\nLIB_END\n\n\n"
  code: "#ifndef SQRT_MOD_HPP\n#define SQRT_MOD_HPP\n\n#include \"../common.hpp\"\n\
    \n#include <random>\n#include <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\n\
    template <typename ModIntT>\nstd::vector<ModIntT> sqrt_mod_prime(ModIntT a) {\n\
    \  // Bostan-Mori's algorithm\n  if (a.is_zero()) return {a};\n  const auto p\
    \ = ModIntT::mod();\n  if (a.pow(p >> 1) == -1) return {};\n  if ((p & 3) == 3)\
    \ {\n    ModIntT b(a.pow((p + 1) >> 2));\n    return {b, -b};\n  }\n  std::mt19937\
    \ gen(std::random_device{}());\n  std::uniform_int_distribution<std::remove_cv_t<decltype(p)>>\
    \ dis(2, p - 1);\n  ModIntT t;\n  do { t = dis(gen); } while ((t * t - 4 * a).pow(p\
    \ >> 1) != -1);\n  ModIntT k0(1), k1, k2(-t), k3(a);\n  for (auto e = (p + 1)\
    \ >> 1;;) {\n    // clang-format off\n    if (e & 1) k0 = k1 - k0 * k2, k1 *=\
    \ k3;\n    else k1 = k0 * k3 - k1 * k2;\n    // clang-format on\n    if ((e >>=\
    \ 1) == 0) return {k0, -k0};\n    k2 = k3 + k3 - k2 * k2, k3 *= k3;\n  }\n}\n\n\
    LIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  isVerificationFile: false
  path: math/sqrt_mod.hpp
  requiredBy: []
  timestamp: '2022-04-26 19:23:58+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/sqrt_mod.0.test.cpp
documentation_of: math/sqrt_mod.hpp
layout: document
title: Square Roots in Finite Fields
---
