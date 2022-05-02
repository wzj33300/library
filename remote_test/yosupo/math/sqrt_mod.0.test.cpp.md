---
data:
  _extendedDependsOn:
  - icon: ':question:'
    path: common.hpp
    title: common.hpp
  - icon: ':question:'
    path: common.hpp
    title: common.hpp
  - icon: ':question:'
    path: math/sqrt_mod.hpp
    title: Square Roots in Finite Fields
  - icon: ':heavy_check_mark:'
    path: modint/runtime_montgomery_modint.hpp
    title: Runtime Montgomery ModInt
  _extendedRequiredBy: []
  _extendedVerifiedWith: []
  _isVerificationFailed: false
  _pathExtension: cpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    '*NOT_SPECIAL_COMMENTS*': ''
    PROBLEM: https://judge.yosupo.jp/problem/sqrt_mod
    links:
    - https://judge.yosupo.jp/problem/sqrt_mod
  bundledCode: "#line 1 \"remote_test/yosupo/math/sqrt_mod.0.test.cpp\"\n#define PROBLEM\
    \ \"https://judge.yosupo.jp/problem/sqrt_mod\"\n\n#line 1 \"math/sqrt_mod.hpp\"\
    \n\n\n\n#line 1 \"common.hpp\"\n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace\
    \ lib {\n#define LIB_END }\n#define LIB ::lib::\n\n\n#line 5 \"math/sqrt_mod.hpp\"\
    \n\n#include <random>\n#include <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\
    \ntemplate <typename ModIntT>\nstd::vector<ModIntT> sqrt_mod_prime(ModIntT a)\
    \ {\n  // Bostan-Mori's algorithm\n  if (a.is_zero()) return {a};\n  const auto\
    \ p = ModIntT::mod();\n  if (a.pow(p >> 1) == -1) return {};\n  if ((p & 3) ==\
    \ 3) {\n    ModIntT b(a.pow((p + 1) >> 2));\n    return {b, -b};\n  }\n  std::mt19937\
    \ gen(std::random_device{}());\n  std::uniform_int_distribution<std::remove_cv_t<decltype(p)>>\
    \ dis(2, p - 1);\n  ModIntT t;\n  do { t = dis(gen); } while ((t * t - 4 * a).pow(p\
    \ >> 1) != -1);\n  ModIntT k0(1), k1, k2(-t), k3(a);\n  for (auto e = (p + 1)\
    \ >> 1;;) {\n    // clang-format off\n    if (e & 1) k0 = k1 - k0 * k2, k1 *=\
    \ k3;\n    else k1 = k0 * k3 - k1 * k2;\n    // clang-format on\n    if ((e >>=\
    \ 1) == 0) return {k0, -k0};\n    k2 = k3 + k3 - k2 * k2, k3 *= k3;\n  }\n}\n\n\
    LIB_END\n\n\n#line 1 \"modint/runtime_montgomery_modint.hpp\"\n\n\n\n#line 5 \"\
    modint/runtime_montgomery_modint.hpp\"\n\n#ifdef LIB_DEBUG\n  #include <stdexcept>\n\
    #endif\n#include <cstdint>\n#include <iostream>\n#line 12 \"modint/runtime_montgomery_modint.hpp\"\
    \n\nLIB_BEGIN\n\ntemplate <int /* IdT */>\nclass runtime_montgomery_modint30 {\n\
    \  using i32 = std::int32_t;\n  using u32 = std::uint32_t;\n  using u64 = std::uint64_t;\n\
    \n  u32 v_{};\n\n  static u32 redc(u64 x) { return (x + static_cast<u64>(static_cast<u32>(x)\
    \ * R) * MOD) >> 32; }\n  static u32 norm(u32 x) { return x - (MOD & -((MOD -\
    \ 1 - x) >> 31)); }\n\n  static u32 MOD, MOD2, R, R2;\n  static i32 SMOD;\n\n\
    public:\n  static bool set_mod(u32 m) {\n    if ((m & 1) == 0 || m == 1 || m >>\
    \ 30 != 0) return false;\n    MOD = m, MOD2 = MOD << 1;\n    {\n      // compute\
    \ R\n      u32 t = 2, iv = MOD * (t - MOD * MOD);\n      iv *= t - MOD * iv, iv\
    \ *= t - MOD * iv;\n      R = iv * (MOD * iv - t);\n    }\n    // compute R2\n\
    \    R2   = -static_cast<u64>(MOD) % MOD;\n    SMOD = static_cast<i32>(MOD);\n\
    \    return true;\n  }\n  static u32 mod() { return MOD; }\n  static i32 smod()\
    \ { return SMOD; }\n  runtime_montgomery_modint30() {}\n  template <typename IntT,\
    \ std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  runtime_montgomery_modint30(IntT\
    \ v) : v_(redc(static_cast<u64>(v % SMOD + SMOD) * R2)) {}\n  u32 val() const\
    \ { return norm(redc(v_)); }\n  i32 sval() const { return norm(redc(v_)); }\n\
    \  bool is_zero() const { return v_ == 0 || v_ == MOD; }\n  template <typename\
    \ IntT, std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  explicit operator\
    \ IntT() const {\n    return static_cast<IntT>(val());\n  }\n  runtime_montgomery_modint30\
    \ operator-() const {\n    runtime_montgomery_modint30 res;\n    res.v_ = (MOD2\
    \ & -(v_ != 0)) - v_;\n    return res;\n  }\n  runtime_montgomery_modint30 inv()\
    \ const {\n    i32 x1 = 1, x3 = 0, a = sval(), b = SMOD;\n    while (b != 0) {\n\
    \      i32 q = a / b, x1_old = x1, a_old = a;\n      x1 = x3, x3 = x1_old - x3\
    \ * q, a = b, b = a_old - b * q;\n    }\n#ifdef LIB_DEBUG\n    if (a != 1) throw\
    \ std::runtime_error(\"modular inverse error\");\n#endif\n    return runtime_montgomery_modint30(x1);\n\
    \  }\n  runtime_montgomery_modint30 &operator+=(const runtime_montgomery_modint30\
    \ &rhs) {\n    v_ += rhs.v_ - MOD2, v_ += MOD2 & -(v_ >> 31);\n    return *this;\n\
    \  }\n  runtime_montgomery_modint30 &operator-=(const runtime_montgomery_modint30\
    \ &rhs) {\n    v_ -= rhs.v_, v_ += MOD2 & -(v_ >> 31);\n    return *this;\n  }\n\
    \  runtime_montgomery_modint30 &operator*=(const runtime_montgomery_modint30 &rhs)\
    \ {\n    v_ = redc(static_cast<u64>(v_) * rhs.v_);\n    return *this;\n  }\n \
    \ runtime_montgomery_modint30 &operator/=(const runtime_montgomery_modint30 &rhs)\
    \ {\n    return operator*=(rhs.inv());\n  }\n  runtime_montgomery_modint30 pow(u64\
    \ e) const {\n    for (runtime_montgomery_modint30 res(1), x(*this);; x *= x)\
    \ {\n      if (e & 1) res *= x;\n      if ((e >>= 1) == 0) return res;\n    }\n\
    \  }\n  void swap(runtime_montgomery_modint30 &rhs) {\n    auto v = v_;\n    v_\
    \ = rhs.v_, rhs.v_ = v;\n  }\n  friend runtime_montgomery_modint30 operator+(const\
    \ runtime_montgomery_modint30 &lhs,\n                                        \
    \       const runtime_montgomery_modint30 &rhs) {\n    return runtime_montgomery_modint30(lhs)\
    \ += rhs;\n  }\n  friend runtime_montgomery_modint30 operator-(const runtime_montgomery_modint30\
    \ &lhs,\n                                               const runtime_montgomery_modint30\
    \ &rhs) {\n    return runtime_montgomery_modint30(lhs) -= rhs;\n  }\n  friend\
    \ runtime_montgomery_modint30 operator*(const runtime_montgomery_modint30 &lhs,\n\
    \                                               const runtime_montgomery_modint30\
    \ &rhs) {\n    return runtime_montgomery_modint30(lhs) *= rhs;\n  }\n  friend\
    \ runtime_montgomery_modint30 operator/(const runtime_montgomery_modint30 &lhs,\n\
    \                                               const runtime_montgomery_modint30\
    \ &rhs) {\n    return runtime_montgomery_modint30(lhs) /= rhs;\n  }\n  friend\
    \ bool operator==(const runtime_montgomery_modint30 &lhs,\n                  \
    \       const runtime_montgomery_modint30 &rhs) {\n    return norm(lhs.v_) ==\
    \ norm(rhs.v_);\n  }\n  friend bool operator!=(const runtime_montgomery_modint30\
    \ &lhs,\n                         const runtime_montgomery_modint30 &rhs) {\n\
    \    return norm(lhs.v_) != norm(rhs.v_);\n  }\n  friend std::istream &operator>>(std::istream\
    \ &is, runtime_montgomery_modint30 &rhs) {\n    i32 x;\n    is >> x;\n    rhs\
    \ = runtime_montgomery_modint30(x);\n    return is;\n  }\n  friend std::ostream\
    \ &operator<<(std::ostream &os, const runtime_montgomery_modint30 &rhs) {\n  \
    \  return os << rhs.val();\n  }\n};\n\ntemplate <int IdT>\ntypename runtime_montgomery_modint30<IdT>::u32\
    \ runtime_montgomery_modint30<IdT>::MOD;\ntemplate <int IdT>\ntypename runtime_montgomery_modint30<IdT>::u32\
    \ runtime_montgomery_modint30<IdT>::MOD2;\ntemplate <int IdT>\ntypename runtime_montgomery_modint30<IdT>::u32\
    \ runtime_montgomery_modint30<IdT>::R;\ntemplate <int IdT>\ntypename runtime_montgomery_modint30<IdT>::u32\
    \ runtime_montgomery_modint30<IdT>::R2;\ntemplate <int IdT>\ntypename runtime_montgomery_modint30<IdT>::i32\
    \ runtime_montgomery_modint30<IdT>::SMOD;\n\ntemplate <int IdT>\nusing rmm30 =\
    \ runtime_montgomery_modint30<IdT>;\n\nLIB_END\n\n\n#line 5 \"remote_test/yosupo/math/sqrt_mod.0.test.cpp\"\
    \n\n#include <cassert>\n#line 8 \"remote_test/yosupo/math/sqrt_mod.0.test.cpp\"\
    \n\nint main() {\n#ifdef LOCAL\n  std::freopen(\"in\", \"r\", stdin), std::freopen(\"\
    out\", \"w\", stdout);\n#endif\n  std::ios::sync_with_stdio(false);\n  std::cin.tie(nullptr);\n\
    \  int t;\n  std::cin >> t;\n  while (t--) {\n    int k, p;\n    std::cin >> k\
    \ >> p;\n    if (p == 2) {\n      std::cout << k << '\\n';\n    } else {\n   \
    \   bool ok = lib::rmm30<-1>::set_mod(p);\n      assert(ok);\n      auto res =\
    \ lib::sqrt_mod_prime<lib::rmm30<-1>>(k);\n      if (res.empty()) {\n        std::cout\
    \ << \"-1\\n\";\n      } else {\n        std::cout << res.front() << '\\n';\n\
    \      }\n    }\n  }\n  return 0;\n}\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/sqrt_mod\"\n\n#include\
    \ \"math/sqrt_mod.hpp\"\n#include \"modint/runtime_montgomery_modint.hpp\"\n\n\
    #include <cassert>\n#include <iostream>\n\nint main() {\n#ifdef LOCAL\n  std::freopen(\"\
    in\", \"r\", stdin), std::freopen(\"out\", \"w\", stdout);\n#endif\n  std::ios::sync_with_stdio(false);\n\
    \  std::cin.tie(nullptr);\n  int t;\n  std::cin >> t;\n  while (t--) {\n    int\
    \ k, p;\n    std::cin >> k >> p;\n    if (p == 2) {\n      std::cout << k << '\\\
    n';\n    } else {\n      bool ok = lib::rmm30<-1>::set_mod(p);\n      assert(ok);\n\
    \      auto res = lib::sqrt_mod_prime<lib::rmm30<-1>>(k);\n      if (res.empty())\
    \ {\n        std::cout << \"-1\\n\";\n      } else {\n        std::cout << res.front()\
    \ << '\\n';\n      }\n    }\n  }\n  return 0;\n}"
  dependsOn:
  - math/sqrt_mod.hpp
  - common.hpp
  - modint/runtime_montgomery_modint.hpp
  - common.hpp
  isVerificationFile: true
  path: remote_test/yosupo/math/sqrt_mod.0.test.cpp
  requiredBy: []
  timestamp: '2022-04-26 21:07:24+08:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: remote_test/yosupo/math/sqrt_mod.0.test.cpp
layout: document
redirect_from:
- /verify/remote_test/yosupo/math/sqrt_mod.0.test.cpp
- /verify/remote_test/yosupo/math/sqrt_mod.0.test.cpp.html
title: remote_test/yosupo/math/sqrt_mod.0.test.cpp
---
