---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  - icon: ':heavy_check_mark:'
    path: src/math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: src/math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: src/modint/montgomery_modint.hpp
    title: Montgomery ModInt
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/convolution_mod
    links:
    - https://judge.yosupo.jp/problem/convolution_mod
  bundledCode: "Traceback (most recent call last):\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: math/convolution.hpp:\
    \ line -1: no such header\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/convolution_mod\"\r\n\r\
    \n#include \"math/convolution.hpp\"\r\n#include \"modint/montgomery_modint.hpp\"\
    \r\n\r\n#include <iostream>\r\n#include <iterator>\r\n\r\nint main() {\r\n#ifdef\
    \ LOCAL\r\n  std::freopen(\"in\", \"r\", stdin), std::freopen(\"out\", \"w\",\
    \ stdout);\r\n#endif\r\n  std::ios::sync_with_stdio(false);\r\n  std::cin.tie(nullptr);\r\
    \n  int n, m;\r\n  std::cin >> n >> m;\r\n  using mint = lib::mm30<998244353>;\r\
    \n  std::vector<mint> a, b;\r\n  std::copy_n(std::istream_iterator<mint>(std::cin),\
    \ n, std::back_inserter(a));\r\n  std::copy_n(std::istream_iterator<mint>(std::cin),\
    \ m, std::back_inserter(b));\r\n  auto res = lib::convolution(a, b);\r\n  std::copy_n(res.begin(),\
    \ n + m - 1, std::ostream_iterator<mint>(std::cout, \" \"));\r\n  return 0;\r\n\
    }"
  dependsOn:
  - src/math/convolution.hpp
  - src/common.hpp
  - src/math/radix2_ntt.hpp
  - src/modint/montgomery_modint.hpp
  isVerificationFile: true
  path: remote_test/yosupo/math/convolution_mod.0.test.cpp
  requiredBy: []
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: remote_test/yosupo/math/convolution_mod.0.test.cpp
layout: document
redirect_from:
- /verify/remote_test/yosupo/math/convolution_mod.0.test.cpp
- /verify/remote_test/yosupo/math/convolution_mod.0.test.cpp.html
title: remote_test/yosupo/math/convolution_mod.0.test.cpp
---
