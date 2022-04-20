---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: src/common.hpp
    title: src/common.hpp
  _extendedRequiredBy: []
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/datastructure/union_find.0.test.cpp
    title: remote_test/yosupo/datastructure/union_find.0.test.cpp
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
  code: "#ifndef DISJOINT_SET_HPP\r\n#define DISJOINT_SET_HPP\r\n\r\n#include \"common.hpp\"\
    \r\n\r\n#include <numeric>\r\n#include <vector>\r\n\r\nLIB_BEGIN\r\n\r\nclass\
    \ disjoint_set {\r\n  mutable std::vector<int> p_; // parent\r\n  std::vector<int>\
    \ s_;         // size\r\n\r\npublic:\r\n  disjoint_set() = default;\r\n  disjoint_set(int\
    \ n) : p_(n), s_(n, 1) { std::iota(p_.begin(), p_.end(), 0); }\r\n  void make_set(int\
    \ n) {\r\n    p_.resize(n);\r\n    s_.assign(n, 1);\r\n    std::iota(p_.begin(),\
    \ p_.end(), 0);\r\n  }\r\n  int find(int u) const {\r\n    // path havling\r\n\
    \    while (p_[u] != p_[p_[u]]) u = p_[u] = p_[p_[u]];\r\n    return p_[u];\r\n\
    \  }\r\n  bool is_same(int u, int v) const { return find(u) == find(v); }\r\n\
    \  int unite(int u, int v) {\r\n    u = find(u), v = find(v);\r\n    if (u ==\
    \ v) return u;\r\n    if (s_[u] < s_[v]) {\r\n      s_[v] += s_[u];\r\n      return\
    \ p_[u] = v;\r\n    } else {\r\n      s_[u] += s_[v];\r\n      return p_[v] =\
    \ u;\r\n    }\r\n  }\r\n  int get_component_size(int u) const { return s_[find(u)];\
    \ }\r\n};\r\n\r\nLIB_END\r\n\r\n#endif"
  dependsOn:
  - src/common.hpp
  isVerificationFile: false
  path: src/datastructure/disjoint_set.hpp
  requiredBy: []
  timestamp: '2022-04-20 11:11:22+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/datastructure/union_find.0.test.cpp
documentation_of: src/datastructure/disjoint_set.hpp
layout: document
title: Disjoint Set
---
