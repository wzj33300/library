---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  - icon: ':heavy_check_mark:'
    path: src/math/radix2_ntt.hpp
    title: Radix-2 NTT
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
  code: "#ifndef CONVOLUTION_HPP\r\n#define CONVOLUTION_HPP\r\n\r\n#include \"common.hpp\"\
    \r\n#include \"radix2_ntt.hpp\"\r\n\r\n#include <algorithm>\r\n#include <vector>\r\
    \n\r\nLIB_BEGIN\r\n\r\ntemplate <typename ModIntT>\r\nstd::vector<ModIntT> convolution(const\
    \ std::vector<ModIntT> &lhs, std::vector<ModIntT> &rhs) {\r\n  int n = static_cast<int>(lhs.size()),\
    \ m = static_cast<int>(rhs.size());\r\n  if (n == 0 || m == 0) return std::vector<ModIntT>{};\r\
    \n  if (std::min(n, m) <= 32) {\r\n    std::vector<ModIntT> res(n + m - 1);\r\n\
    \    for (int i = 0; i != n; ++i)\r\n      for (int j = 0; j != m; ++j) res[i\
    \ + j] += lhs[i] * rhs[j];\r\n    return res;\r\n  }\r\n  int len = ntt_len(n\
    \ + m - 1);\r\n  std::vector<ModIntT> lhs_cpy(len), rhs_cpy(len);\r\n  std::copy_n(lhs.cbegin(),\
    \ n, lhs_cpy.begin());\r\n  std::copy_n(rhs.cbegin(), m, rhs_cpy.begin());\r\n\
    \  dft_n(lhs_cpy.begin(), len), dft_n(rhs_cpy.begin(), len);\r\n  for (int i =\
    \ 0; i != len; ++i) lhs_cpy[i] *= rhs_cpy[i];\r\n  idft_n(lhs_cpy.begin(), len);\r\
    \n  lhs_cpy.resize(n + m - 1);\r\n  return lhs_cpy;\r\n}\r\n\r\nLIB_END\r\n\r\n\
    #endif"
  dependsOn:
  - src/common.hpp
  - src/math/radix2_ntt.hpp
  isVerificationFile: false
  path: src/math/convolution.hpp
  requiredBy: []
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
documentation_of: src/math/convolution.hpp
layout: document
title: Convolution
---
