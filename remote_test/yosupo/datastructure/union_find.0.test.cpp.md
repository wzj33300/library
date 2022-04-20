---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  - icon: ':heavy_check_mark:'
    path: src/datastructure/disjoint_set.hpp
    title: Disjoint Set
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/unionfind
    links:
    - https://judge.yosupo.jp/problem/unionfind
  bundledCode: "Traceback (most recent call last):\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/documentation/build.py\"\
    , line 71, in _render_source_code_stat\n    bundled_code = language.bundle(stat.path,\
    \ basedir=basedir, options={'include_paths': [basedir]}).decode()\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus.py\"\
    , line 187, in bundle\n    bundler.update(path)\n  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 401, in update\n    self.update(self._resolve(pathlib.Path(included), included_from=path))\n\
    \  File \"/opt/hostedtoolcache/Python/3.10.4/x64/lib/python3.10/site-packages/onlinejudge_verify/languages/cplusplus_bundle.py\"\
    , line 260, in _resolve\n    raise BundleErrorAt(path, -1, \"no such header\"\
    )\nonlinejudge_verify.languages.cplusplus_bundle.BundleErrorAt: datastructure/disjoint_set.hpp:\
    \ line -1: no such header\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/unionfind\"\n\n#include\
    \ \"datastructure/disjoint_set.hpp\"\n\n#include <iostream>\n\nint main() {\n\
    #ifdef LOCAL\n  std::freopen(\"in\", \"r\", stdin), std::freopen(\"out\", \"w\"\
    , stdout);\n#endif\n  std::ios::sync_with_stdio(false);\n  std::cin.tie(nullptr);\n\
    \  int n, q;\n  std::cin >> n >> q;\n  lib::disjoint_set ds(n);\n  while (q--)\
    \ {\n    int cmd, u, v;\n    std::cin >> cmd >> u >> v;\n    if (cmd == 0) {\n\
    \      ds.unite(u, v);\n    } else {\n      std::cout << static_cast<int>(ds.is_same(u,\
    \ v)) << '\\n';\n    }\n  }\n  return 0;\n}"
  dependsOn:
  - src/datastructure/disjoint_set.hpp
  - src/common.hpp
  isVerificationFile: true
  path: remote_test/yosupo/datastructure/union_find.0.test.cpp
  requiredBy: []
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: remote_test/yosupo/datastructure/union_find.0.test.cpp
layout: document
redirect_from:
- /verify/remote_test/yosupo/datastructure/union_find.0.test.cpp
- /verify/remote_test/yosupo/datastructure/union_find.0.test.cpp.html
title: remote_test/yosupo/datastructure/union_find.0.test.cpp
---
