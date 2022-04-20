---
data:
  _extendedDependsOn: []
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: src/datastructure/disjoint_set.hpp
    title: Disjoint Set
  - icon: ':heavy_check_mark:'
    path: src/math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: src/math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: src/modint/montgomery_modint.hpp
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
  bundledCode: "#line 1 \"src/common.hpp\"\n\n\n\r\n#define LIB_DEBUG\r\n\r\n#define\
    \ LIB_BEGIN namespace lib {\r\n#define LIB_END }\r\n#define LIB ::lib::\r\n\r\n\
    \n"
  code: "#ifndef COMMON_HPP\r\n#define COMMON_HPP\r\n\r\n#define LIB_DEBUG\r\n\r\n\
    #define LIB_BEGIN namespace lib {\r\n#define LIB_END }\r\n#define LIB ::lib::\r\
    \n\r\n#endif"
  dependsOn: []
  isVerificationFile: false
  path: src/common.hpp
  requiredBy:
  - src/math/radix2_ntt.hpp
  - src/math/convolution.hpp
  - src/modint/montgomery_modint.hpp
  - src/datastructure/disjoint_set.hpp
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
  - remote_test/yosupo/datastructure/union_find.0.test.cpp
documentation_of: src/common.hpp
layout: document
redirect_from:
- /library/src/common.hpp
- /library/src/common.hpp.html
title: src/common.hpp
---
