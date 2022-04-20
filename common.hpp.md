---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: datastructure/disjoint_set.hpp
    title: Disjoint Set
  - icon: ':heavy_check_mark:'
    path: math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: modint/montgomery_modint.hpp
    title: Montgomery ModInt
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/datastructure/union_find.0.test.cpp
    title: remote_test/yosupo/datastructure/union_find.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
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
  - math/radix2_ntt.hpp
  - math/convolution.hpp
  - modint/montgomery_modint.hpp
  - datastructure/disjoint_set.hpp
  timestamp: '2022-04-20 11:49:11+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
  - remote_test/yosupo/datastructure/union_find.0.test.cpp
documentation_of: common.hpp
layout: document
redirect_from:
- /library/common.hpp
- /library/common.hpp.html
title: common.hpp
---
