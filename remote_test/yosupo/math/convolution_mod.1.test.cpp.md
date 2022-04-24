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
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  - icon: ':heavy_check_mark:'
    path: math/truncated_formal_power_series.hpp
    title: Truncated Formal Power Series
  - icon: ':question:'
    path: modint/montgomery_modint.hpp
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
  bundledCode: "#line 1 \"remote_test/yosupo/math/convolution_mod.1.test.cpp\"\n#define\
    \ PROBLEM \"https://judge.yosupo.jp/problem/convolution_mod\"\n\n#line 1 \"math/truncated_formal_power_series.hpp\"\
    \n\n\n\n#line 1 \"common.hpp\"\n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace\
    \ lib {\n#define LIB_END }\n#define LIB ::lib::\n\n\n#line 1 \"math/radix2_ntt.hpp\"\
    \n\n\n\n#line 5 \"math/radix2_ntt.hpp\"\n\n#include <array>\n#include <cassert>\n\
    #include <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\nnamespace detail {\n\
    \ntemplate <typename IntT>\nconstexpr std::enable_if_t<std::is_integral_v<IntT>,\
    \ int> bsf(IntT v) {\n  if (static_cast<std::make_signed_t<IntT>>(v) <= 0) return\
    \ -1;\n  int res = 0;\n  for (; (v & 1) == 0; ++res) v >>= 1;\n  return res;\n\
    }\n\ntemplate <typename ModIntT>\nconstexpr ModIntT quadratic_nonresidue_prime()\
    \ {\n  auto mod = ModIntT::mod();\n  for (int i = 2;; ++i)\n    if (ModIntT(i).pow(mod\
    \ >> 1) == mod - 1) return ModIntT(i);\n}\n\ntemplate <typename ModIntT>\nconstexpr\
    \ ModIntT gen_of_sylow_2_subgroup() {\n  auto mod = ModIntT::mod();\n  return\
    \ quadratic_nonresidue_prime<ModIntT>().pow(mod >> bsf(mod - 1));\n}\n\ntemplate\
    \ <typename ModIntT>\nconstexpr std::array<ModIntT, bsf(ModIntT::mod() - 1) -\
    \ 1> root() {\n  std::array<ModIntT, bsf(ModIntT::mod() - 1) - 1> rt; // order(`rt[i]`)\
    \ = 2^(i + 2).\n  rt.back() = gen_of_sylow_2_subgroup<ModIntT>();\n  for (int\
    \ i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) rt[i] = rt[i + 1] * rt[i + 1];\n\
    \  return rt;\n}\n\ntemplate <typename ModIntT>\nconstexpr std::array<ModIntT,\
    \ bsf(ModIntT::mod() - 1) - 1> iroot() {\n  std::array<ModIntT, bsf(ModIntT::mod()\
    \ - 1) - 1> irt;\n  irt.back() = gen_of_sylow_2_subgroup<ModIntT>().inv();\n \
    \ for (int i = bsf(ModIntT::mod() - 1) - 3; i >= 0; --i) irt[i] = irt[i + 1] *\
    \ irt[i + 1];\n  return irt;\n}\n\n} // namespace detail\n\n// Input:  integer\
    \ `n`.\n// Output: 2^(\u2308log_2(`n`)\u2309).\nint ntt_len(int n) {\n  --n;\n\
    \  n |= n >> 1, n |= n >> 2, n |= n >> 4, n |= n >> 8;\n  return (n | n >> 16)\
    \ + 1;\n}\n\n// Input:           f(x) = `a[0]` + `a[1]`x + ... + `a[n - 1]`x^(`n`\
    \ - 1) where `n` is power of 2.\n// Output(inplace): reversed binary permutation\
    \ of [f(\u03B6^0), f(\u03B6), f(\u03B6^2), ..., f(\u03B6^(`n` - 1))].\ntemplate\
    \ <typename IterT>\nvoid dft_n(IterT a, int n) {\n  assert((n & (n - 1)) == 0);\n\
    \  using T                  = typename std::iterator_traits<IterT>::value_type;\n\
    \  static constexpr auto rt = detail::root<T>();\n  static std::vector<T> root(1);\n\
    \  if (int s = static_cast<int>(root.size()); s << 1 < n) {\n    root.resize(n\
    \ >> 1);\n    for (int i = detail::bsf(s), j; 1 << i < n >> 1; ++i) {\n      root[j\
    \ = 1 << i] = rt[i];\n      for (int k = j + 1; k < j << 1; ++k) root[k] = root[k\
    \ - j] * root[j];\n    }\n  }\n  for (int j = 0, l = n >> 1; j != l; ++j) {\n\
    \    T u(a[j]), v(a[j + l]);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n  for\
    \ (int i = n >> 1; i >= 2; i >>= 1) {\n    for (int j = 0, l = i >> 1; j != l;\
    \ ++j) {\n      T u(a[j]), v(a[j + l]);\n      a[j] = u + v, a[j + l] = u - v;\n\
    \    }\n    for (int j = i, l = i >> 1, m = 1; j != n; j += i, ++m)\n      for\
    \ (int k = j; k != j + l; ++k) {\n        T u(a[k]), v(a[k + l] * root[m]);\n\
    \        a[k] = u + v, a[k + l] = u - v;\n      }\n  }\n}\n\n// Input:       \
    \    reversed binary permutation of [f(\u03B6^0), f(\u03B6), f(\u03B6^2), ...,\
    \ f(\u03B6^(`n` - 1))].\n// Output(inplace): f(x) = `a[0]` + `a[1]`x + ... + `a[n\
    \ - 1]`x^(`n` - 1) where `n` is power of 2.\ntemplate <typename IterT>\nvoid idft_n(IterT\
    \ a, int n) {\n  assert((n & (n - 1)) == 0);\n  using T                  = typename\
    \ std::iterator_traits<IterT>::value_type;\n  static constexpr auto rt = detail::iroot<T>();\n\
    \  static std::vector<T> root(1);\n  if (int s = static_cast<int>(root.size());\
    \ s << 1 < n) {\n    root.resize(n >> 1);\n    for (int i = detail::bsf(s), j;\
    \ 1 << i < n >> 1; ++i) {\n      root[j = 1 << i] = rt[i];\n      for (int k =\
    \ j + 1; k < j << 1; ++k) root[k] = root[k - j] * root[j];\n    }\n  }\n  for\
    \ (int i = 2; i < n; i <<= 1) {\n    for (int j = 0, l = i >> 1; j != l; ++j)\
    \ {\n      T u(a[j]), v(a[j + l]);\n      a[j] = u + v, a[j + l] = u - v;\n  \
    \  }\n    for (int j = i, l = i >> 1, m = 1; j != n; j += i, ++m)\n      for (int\
    \ k = j; k != j + l; ++k) {\n        T u(a[k]), v(a[k + l]);\n        a[k] = u\
    \ + v, a[k + l] = (u - v) * root[m];\n      }\n  }\n  const T iv(T::mod() - T::mod()\
    \ / n);\n  for (int j = 0, l = n >> 1; j != l; ++j) {\n    T u(a[j] * iv), v(a[j\
    \ + l] * iv);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n}\n\n// clang-format\
    \ off\ntemplate <typename ContainerT> void dft(ContainerT &&a) { dft_n(a.begin(),\
    \ a.size()); }\ntemplate <typename ContainerT> void idft(ContainerT &&a) { idft_n(a.begin(),\
    \ a.size()); }\ntemplate <typename IterT> void dft(IterT beg, IterT end) { dft_n(beg,\
    \ end - beg); }\ntemplate <typename IterT> void idft(IterT beg, IterT end) { idft_n(beg,\
    \ end - beg); }\n// clang-format on\n\ntemplate <typename T>\nvoid dft_doubling(const\
    \ std::vector<T> &a, std::vector<T> &dft_a) {\n  static constexpr auto rt = detail::root<T>();\n\
    \  int as = static_cast<int>(a.size()), n = static_cast<int>(dft_a.size());\n\
    \  // `dft_a` = `dft_n`(`a` mod (x^n - 1))\n  // doubling `dft_a` is just computing\
    \ dft_n((`a` mod (x^n + 1))(\u03B6^(2n))).\n  dft_a.resize(n << 1);\n  auto it\
    \ = dft_a.begin() + n;\n  for (int i = 0, is_even = 0, j; i != as; ++i) {\n  \
    \  if ((j = i & (n - 1)) == 0) is_even ^= 1;\n    it[j] += is_even ? a[i] : -a[i];\n\
    \  }\n  T r(n == 1 ? T(-1) : rt[detail::bsf(n) - 1]), v(1);\n  for (int i = 0;\
    \ i != n; ++i) it[i] *= v, v *= r;\n  dft_n(it, n);\n}\n\nLIB_END\n\n\n#line 6\
    \ \"math/truncated_formal_power_series.hpp\"\n\n#include <algorithm>\n#line 9\
    \ \"math/truncated_formal_power_series.hpp\"\n#include <iostream>\n#include <iterator>\n\
    #line 13 \"math/truncated_formal_power_series.hpp\"\n\nLIB_BEGIN\n\ntemplate <typename\
    \ ModIntT>\nclass truncated_formal_power_series : public std::vector<ModIntT>\
    \ {\n  static_assert(std::is_same_v<typename std::vector<ModIntT>::value_type,\
    \ ModIntT>);\n\npublic:\n  using std::vector<ModIntT>::vector;\n\n  enum : int\
    \ { NEGATIVE_INFINITY = -1 };\n\n  // leading coefficient\n  ModIntT lc() const\
    \ {\n    int d = deg();\n    return d == NEGATIVE_INFINITY ? ModIntT() : this->operator[](d);\n\
    \  }\n  // degree\n  int deg() const {\n    // treat formal power series like\
    \ polynomials\n    int n = static_cast<int>(this->size());\n    while (n >= 0\
    \ && !this->operator[](n).is_zero()) --n;\n    return n == -1 ? NEGATIVE_INFINITY\
    \ : n;\n  }\n  // order\n  int ord() const;\n  bool is_zero() const { return deg()\
    \ == NEGATIVE_INFINITY; }\n  void shrink() { this->resize(deg() + 1); }\n  truncated_formal_power_series\
    \ operator-() {\n    truncated_formal_power_series res(*this);\n    for (auto\
    \ &&i : res) i = -i;\n    return res;\n  }\n\n  truncated_formal_power_series\
    \ &operator+=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) += rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator-=(const truncated_formal_power_series &rhs) {\n    if (this->size()\
    \ < rhs.size()) this->resize(rhs.size());\n    for (int i = 0, e = static_cast<int>(rhs.size());\
    \ i != e; ++i) this->operator[](i) -= rhs[i];\n    return *this;\n  }\n  truncated_formal_power_series\
    \ &operator*=(const truncated_formal_power_series &rhs);\n  truncated_formal_power_series\
    \ inv(int n) const;\n  truncated_formal_power_series log(int n) const;\n  truncated_formal_power_series\
    \ exp(int n) const;\n  truncated_formal_power_series div(const truncated_formal_power_series\
    \ &rhs, int n) const;\n\n  friend truncated_formal_power_series operator+(const\
    \ truncated_formal_power_series &lhs,\n                                      \
    \           const truncated_formal_power_series &rhs) {\n    return truncated_formal_power_series(lhs)\
    \ += rhs;\n  }\n  friend truncated_formal_power_series operator-(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) -= rhs;\n  }\n  friend\
    \ truncated_formal_power_series operator*(const truncated_formal_power_series\
    \ &lhs,\n                                                 const truncated_formal_power_series\
    \ &rhs) {\n    return truncated_formal_power_series(lhs) *= rhs;\n  }\n\n  friend\
    \ std::istream &operator>>(std::istream &lhs, truncated_formal_power_series &rhs)\
    \ {\n    for (auto &&i : rhs) lhs >> i;\n    return lhs;\n  }\n  friend std::ostream\
    \ &operator<<(std::ostream &lhs, const truncated_formal_power_series &rhs) {\n\
    \    int s = 0, e = static_cast<int>(rhs.size());\n    lhs << '[';\n    for (auto\
    \ &&i : rhs) {\n      lhs << i;\n      if (s >= 1) lhs << 'x';\n      if (s >\
    \ 1) lhs << \"^(\" << s << ')';\n      if (++s != e) lhs << \" + \";\n    }\n\
    \    return lhs << ']';\n  }\n};\n\ntemplate <typename IterT>\ntruncated_formal_power_series(IterT,\
    \ IterT)\n    -> truncated_formal_power_series<typename std::iterator_traits<IterT>::value_type>;\n\
    \ntemplate <typename ModIntT>\nusing tfps = truncated_formal_power_series<ModIntT>;\n\
    \ntemplate <typename ModIntT>\ntfps<ModIntT> &tfps<ModIntT>::operator*=(const\
    \ tfps<ModIntT> &rhs) {\n  // 6E\n  int n = static_cast<int>(this->size()), m\
    \ = static_cast<int>(rhs.size());\n  if (n == 0 || m == 0) {\n    this->clear();\n\
    \    return *this;\n  }\n  if (std::min(n, m) <= 32) {\n    tfps<ModIntT> res(n\
    \ + m - 1);\n    for (int i = 0; i != n; ++i)\n      for (int j = 0; j != m; ++j)\
    \ res[i + j] += this->operator[](i) * rhs[j];\n    return this->operator=(res);\n\
    \  }\n  int len = ntt_len(n + m - 1);\n  tfps<ModIntT> rhs_cpy(len);\n  std::copy_n(rhs.cbegin(),\
    \ m, rhs_cpy.begin());\n  this->resize(len);\n  dft_n(this->begin(), len), dft_n(rhs_cpy.begin(),\
    \ len);\n  for (int i = 0; i != len; ++i) this->operator[](i) *= rhs_cpy[i];\n\
    \  idft_n(this->begin(), len);\n  this->resize(n + m - 1);\n  return *this;\n\
    }\n\ntemplate <typename ModIntT>\ntfps<ModIntT> tfps<ModIntT>::inv(int n) const\
    \ {\n  // 10E\n  assert(n >= 0);\n  if (n == 0) return tfps<ModIntT>{};\n  assert(!this->front().is_zero());\n\
    \  if (n == 1) return tfps<ModIntT>{this->front().inv()};\n  int len = ntt_len(n);\n\
    \  tfps<ModIntT> res(len), temp0(len), temp1(len), cpy(len);\n  std::copy(this->cbegin(),\
    \ this->cend(), cpy.begin());\n  res.front() = this->front().inv();\n  for (int\
    \ i = 2; i <= len; i <<= 1) {\n    std::copy_n(cpy.cbegin(), i, temp0.begin());\n\
    \    dft_n(temp0.begin(), i); // 2E\n    std::copy_n(res.cbegin(), i, temp1.begin());\n\
    \    dft_n(temp1.begin(), i); // 2E\n    for (int j = 0; j != i; ++j) temp0[j]\
    \ *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n    std::fill_n(temp0.begin(),\
    \ i >> 1, ModIntT());\n    dft_n(temp0.begin(), i); // 2E\n    for (int j = 0;\
    \ j != i; ++j) temp0[j] *= temp1[j];\n    idft_n(temp0.begin(), i); // 2E\n  \
    \  for (int j = i >> 1; j != i; ++j) res[j] = -temp0[j];\n  }\n  res.resize(n);\n\
    \  return res;\n}\n\nLIB_END\n\n\n#line 1 \"modint/montgomery_modint.hpp\"\n\n\
    \n\n#line 5 \"modint/montgomery_modint.hpp\"\n\n#ifdef LIB_DEBUG\n  #include <stdexcept>\n\
    #endif\n#include <cstdint>\n#line 12 \"modint/montgomery_modint.hpp\"\n\nLIB_BEGIN\n\
    \ntemplate <std::uint32_t ModT>\nclass montgomery_modint30 {\n  using i32 = std::int32_t;\n\
    \  using u32 = std::uint32_t;\n  using u64 = std::uint64_t;\n\n  u32 v_{};\n\n\
    \  static constexpr u32 get_r() {\n    u32 t = 2, iv = MOD * (t - MOD * MOD);\n\
    \    iv *= t - MOD * iv, iv *= t - MOD * iv;\n    return iv * (MOD * iv - t);\n\
    \  }\n  static constexpr u32 redc(u64 x) {\n    return (x + static_cast<u64>(static_cast<u32>(x)\
    \ * R) * MOD) >> 32;\n  }\n  static constexpr u32 norm(u32 x) { return x - (MOD\
    \ & -((MOD - 1 - x) >> 31)); }\n\n  static constexpr u32 MOD  = ModT;\n  static\
    \ constexpr u32 MOD2 = MOD << 1;\n  static constexpr u32 R    = get_r();\n  static\
    \ constexpr u32 R2   = -static_cast<u64>(MOD) % MOD;\n  static constexpr i32 SMOD\
    \ = static_cast<i32>(MOD);\n\n  static_assert(MOD & 1);\n  static_assert(-R *\
    \ MOD == 1);\n  static_assert((MOD >> 30) == 0);\n  static_assert(MOD != 1);\n\
    \npublic:\n  static constexpr u32 mod() { return MOD; }\n  static constexpr i32\
    \ smod() { return SMOD; }\n  constexpr montgomery_modint30() {}\n  template <typename\
    \ IntT, std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  constexpr montgomery_modint30(IntT\
    \ v) : v_(redc(static_cast<u64>(v % SMOD + SMOD) * R2)) {}\n  constexpr u32 val()\
    \ const { return norm(redc(v_)); }\n  constexpr i32 sval() const { return norm(redc(v_));\
    \ }\n  constexpr bool is_zero() const { return v_ == 0 || v_ == MOD; }\n  template\
    \ <typename IntT, std::enable_if_t<std::is_integral_v<IntT>, int> = 0>\n  explicit\
    \ constexpr operator IntT() const {\n    return static_cast<IntT>(val());\n  }\n\
    \  constexpr montgomery_modint30 operator-() const {\n    montgomery_modint30\
    \ res;\n    res.v_ = (MOD2 & -(v_ != 0)) - v_;\n    return res;\n  }\n  constexpr\
    \ montgomery_modint30 inv() const {\n    i32 x1 = 1, x3 = 0, a = sval(), b = SMOD;\n\
    \    while (b != 0) {\n      i32 q = a / b, x1_old = x1, a_old = a;\n      x1\
    \ = x3, x3 = x1_old - x3 * q, a = b, b = a_old - b * q;\n    }\n#ifdef LIB_DEBUG\n\
    \    if (a != 1) throw std::runtime_error(\"modular inverse error\");\n#endif\n\
    \    return montgomery_modint30(x1);\n  }\n  constexpr montgomery_modint30 &operator+=(const\
    \ montgomery_modint30 &rhs) {\n    v_ += rhs.v_ - MOD2, v_ += MOD2 & -(v_ >> 31);\n\
    \    return *this;\n  }\n  constexpr montgomery_modint30 &operator-=(const montgomery_modint30\
    \ &rhs) {\n    v_ -= rhs.v_, v_ += MOD2 & -(v_ >> 31);\n    return *this;\n  }\n\
    \  constexpr montgomery_modint30 &operator*=(const montgomery_modint30 &rhs) {\n\
    \    v_ = redc(static_cast<u64>(v_) * rhs.v_);\n    return *this;\n  }\n  constexpr\
    \ montgomery_modint30 &operator/=(const montgomery_modint30 &rhs) {\n    return\
    \ operator*=(rhs.inv());\n  }\n  constexpr montgomery_modint30 pow(u64 e) const\
    \ {\n    for (montgomery_modint30 res(1), x(*this);; x *= x) {\n      if (e &\
    \ 1) res *= x;\n      if ((e >>= 1) == 0) return res;\n    }\n  }\n  constexpr\
    \ void swap(montgomery_modint30 &rhs) {\n    auto v = v_;\n    v_ = rhs.v_, rhs.v_\
    \ = v;\n  }\n  friend constexpr montgomery_modint30 operator+(const montgomery_modint30\
    \ &lhs,\n                                                 const montgomery_modint30\
    \ &rhs) {\n    return montgomery_modint30(lhs) += rhs;\n  }\n  friend constexpr\
    \ montgomery_modint30 operator-(const montgomery_modint30 &lhs,\n            \
    \                                     const montgomery_modint30 &rhs) {\n    return\
    \ montgomery_modint30(lhs) -= rhs;\n  }\n  friend constexpr montgomery_modint30\
    \ operator*(const montgomery_modint30 &lhs,\n                                \
    \                 const montgomery_modint30 &rhs) {\n    return montgomery_modint30(lhs)\
    \ *= rhs;\n  }\n  friend constexpr montgomery_modint30 operator/(const montgomery_modint30\
    \ &lhs,\n                                                 const montgomery_modint30\
    \ &rhs) {\n    return montgomery_modint30(lhs) /= rhs;\n  }\n  friend constexpr\
    \ bool operator==(const montgomery_modint30 &lhs, const montgomery_modint30 &rhs)\
    \ {\n    return norm(lhs.v_) == norm(rhs.v_);\n  }\n  friend constexpr bool operator!=(const\
    \ montgomery_modint30 &lhs, const montgomery_modint30 &rhs) {\n    return norm(lhs.v_)\
    \ != norm(rhs.v_);\n  }\n  friend std::istream &operator>>(std::istream &is, montgomery_modint30\
    \ &rhs) {\n    i32 x;\n    is >> x;\n    rhs = montgomery_modint30(x);\n    return\
    \ is;\n  }\n  friend std::ostream &operator<<(std::ostream &os, const montgomery_modint30\
    \ &rhs) {\n    return os << rhs.val();\n  }\n};\n\ntemplate <std::uint32_t ModT>\n\
    using mm30 = montgomery_modint30<ModT>;\n\nLIB_END\n\n\n#line 5 \"remote_test/yosupo/math/convolution_mod.1.test.cpp\"\
    \n\n#line 8 \"remote_test/yosupo/math/convolution_mod.1.test.cpp\"\n\nint main()\
    \ {\n#ifdef LOCAL\n  std::freopen(\"in\", \"r\", stdin), std::freopen(\"out\"\
    , \"w\", stdout);\n#endif\n  std::ios::sync_with_stdio(false);\n  std::cin.tie(nullptr);\n\
    \  int n, m;\n  std::cin >> n >> m;\n  using mint = lib::mm30<998244353>;\n  lib::tfps<mint>\
    \ a, b;\n  std::copy_n(std::istream_iterator<mint>(std::cin), n, std::back_inserter(a));\n\
    \  std::copy_n(std::istream_iterator<mint>(std::cin), m, std::back_inserter(b));\n\
    \  auto ab = a * b;\n  std::copy(ab.cbegin(), ab.cend(), std::ostream_iterator<mint>(std::cout,\
    \ \" \"));\n  return 0;\n}\n"
  code: "#define PROBLEM \"https://judge.yosupo.jp/problem/convolution_mod\"\n\n#include\
    \ \"math/truncated_formal_power_series.hpp\"\n#include \"modint/montgomery_modint.hpp\"\
    \n\n#include <iostream>\n#include <iterator>\n\nint main() {\n#ifdef LOCAL\n \
    \ std::freopen(\"in\", \"r\", stdin), std::freopen(\"out\", \"w\", stdout);\n\
    #endif\n  std::ios::sync_with_stdio(false);\n  std::cin.tie(nullptr);\n  int n,\
    \ m;\n  std::cin >> n >> m;\n  using mint = lib::mm30<998244353>;\n  lib::tfps<mint>\
    \ a, b;\n  std::copy_n(std::istream_iterator<mint>(std::cin), n, std::back_inserter(a));\n\
    \  std::copy_n(std::istream_iterator<mint>(std::cin), m, std::back_inserter(b));\n\
    \  auto ab = a * b;\n  std::copy(ab.cbegin(), ab.cend(), std::ostream_iterator<mint>(std::cout,\
    \ \" \"));\n  return 0;\n}"
  dependsOn:
  - math/truncated_formal_power_series.hpp
  - common.hpp
  - math/radix2_ntt.hpp
  - modint/montgomery_modint.hpp
  - common.hpp
  isVerificationFile: true
  path: remote_test/yosupo/math/convolution_mod.1.test.cpp
  requiredBy: []
  timestamp: '2022-04-24 22:23:52+08:00'
  verificationStatus: TEST_ACCEPTED
  verifiedWith: []
documentation_of: remote_test/yosupo/math/convolution_mod.1.test.cpp
layout: document
redirect_from:
- /verify/remote_test/yosupo/math/convolution_mod.1.test.cpp
- /verify/remote_test/yosupo/math/convolution_mod.1.test.cpp.html
title: remote_test/yosupo/math/convolution_mod.1.test.cpp
---
