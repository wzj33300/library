---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: datastructure/disjoint_set.hpp
    title: Disjoint Set
  - icon: ':heavy_check_mark:'
    path: datastructure/weighted_disjoint_set.hpp
    title: Weighted Disjoint Set
  - icon: ':heavy_check_mark:'
    path: math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: math/czt.hpp
    title: Chirp Z-transform (Bluestein's algorithm)
  - icon: ':heavy_check_mark:'
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: math/relaxed_convolution.hpp
    title: Relaxed Convolution
  - icon: ':heavy_check_mark:'
    path: math/semi_relaxed_convolution.hpp
    title: Semi-Relaxed Convolution
  - icon: ':heavy_check_mark:'
    path: math/truncated_formal_power_series.hpp
    title: Truncated Formal Power Series
  - icon: ':heavy_check_mark:'
    path: modint/long_montgomery_modint.hpp
    title: Long Montgomery ModInt
  - icon: ':heavy_check_mark:'
    path: modint/montgomery_modint.hpp
    title: Montgomery ModInt
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/aizu/datastructure/weighted_union_find.0.test.cpp
    title: remote_test/aizu/datastructure/weighted_union_find.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/datastructure/union_find.0.test.cpp
    title: remote_test/yosupo/datastructure/union_find.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.1.test.cpp
    title: remote_test/yosupo/math/convolution_mod.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.1.test.cpp
    title: remote_test/yosupo/math/convolution_mod.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.2.test.cpp
    title: remote_test/yosupo/math/convolution_mod.2.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.2.test.cpp
    title: remote_test/yosupo/math/convolution_mod.2.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.3.test.cpp
    title: remote_test/yosupo/math/convolution_mod.3.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.3.test.cpp
    title: remote_test/yosupo/math/convolution_mod.3.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: '#line 1 "common.hpp"




    #define LIB_DEBUG


    #define LIB_BEGIN namespace lib {

    #define LIB_END }

    #define LIB ::lib::



    '
  code: '#ifndef COMMON_HPP

    #define COMMON_HPP


    #define LIB_DEBUG


    #define LIB_BEGIN namespace lib {

    #define LIB_END }

    #define LIB ::lib::


    #endif'
  dependsOn: []
  isVerificationFile: false
  path: common.hpp
  requiredBy:
  - math/semi_relaxed_convolution.hpp
  - math/radix2_ntt.hpp
  - math/relaxed_convolution.hpp
  - math/convolution.hpp
  - math/convolution.hpp
  - math/truncated_formal_power_series.hpp
  - math/czt.hpp
  - modint/montgomery_modint.hpp
  - modint/long_montgomery_modint.hpp
  - datastructure/weighted_disjoint_set.hpp
  - datastructure/disjoint_set.hpp
  timestamp: '2022-04-20 11:49:11+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/aizu/datastructure/weighted_union_find.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.1.test.cpp
  - remote_test/yosupo/math/convolution_mod.1.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.0.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.1.test.cpp
  - remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.3.test.cpp
  - remote_test/yosupo/math/convolution_mod.3.test.cpp
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.2.test.cpp
  - remote_test/yosupo/math/convolution_mod.2.test.cpp
  - remote_test/yosupo/datastructure/union_find.0.test.cpp
documentation_of: common.hpp
layout: document
redirect_from:
- /library/common.hpp
- /library/common.hpp.html
title: common.hpp
---
