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
  code: "#ifndef DISJOINT_SET_HPP\n#define DISJOINT_SET_HPP\n\n#include \"common.hpp\"\
    \n\n#include <numeric>\n#include <vector>\n\nLIB_BEGIN\n\nclass disjoint_set {\n\
    \  mutable std::vector<int> p_; // parent\n  std::vector<int> s_;         // size\n\
    \npublic:\n  disjoint_set() = default;\n  disjoint_set(int n) : p_(n), s_(n, 1)\
    \ { std::iota(p_.begin(), p_.end(), 0); }\n  void make_set(int n) {\n    p_.resize(n);\n\
    \    s_.assign(n, 1);\n    std::iota(p_.begin(), p_.end(), 0);\n  }\n  int find(int\
    \ u) const {\n    // path havling\n    while (p_[u] != p_[p_[u]]) u = p_[u] =\
    \ p_[p_[u]];\n    return p_[u];\n  }\n  bool is_same(int u, int v) const { return\
    \ find(u) == find(v); }\n  int unite(int u, int v) {\n    u = find(u), v = find(v);\n\
    \    if (u == v) return u;\n    if (s_[u] < s_[v]) {\n      s_[v] += s_[u];\n\
    \      return p_[u] = v;\n    } else {\n      s_[u] += s_[v];\n      return p_[v]\
    \ = u;\n    }\n  }\n  int get_component_size(int u) const { return s_[find(u)];\
    \ }\n};\n\nLIB_END\n\n#endif"
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
