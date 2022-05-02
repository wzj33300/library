---
data:
  _extendedDependsOn:
  - icon: ':heavy_check_mark:'
    path: common.hpp
    title: common.hpp
  - icon: ':heavy_check_mark:'
    path: math/radix2_ntt.hpp
    title: Radix-2 NTT
  _extendedRequiredBy:
  - icon: ':heavy_check_mark:'
    path: math/convolution.hpp
    title: Convolution
  - icon: ':heavy_check_mark:'
    path: math/polynomial.hpp
    title: Polynomial
  - icon: ':heavy_check_mark:'
    path: math/truncated_formal_power_series.hpp
    title: Truncated Formal Power Series
  _extendedVerifiedWith:
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod.1.test.cpp
    title: remote_test/yosupo/math/convolution_mod.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
    title: remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/division_of_polynomials.0.test.cpp
    title: remote_test/yosupo/math/division_of_polynomials.0.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/exp_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/exp_of_formal_power_series.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/inv_of_formal_power_series.2.test.cpp
    title: remote_test/yosupo/math/inv_of_formal_power_series.2.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/log_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/log_of_formal_power_series.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/pow_of_formal_power_series.1.test.cpp
    title: remote_test/yosupo/math/pow_of_formal_power_series.1.test.cpp
  - icon: ':heavy_check_mark:'
    path: remote_test/yosupo/math/sqrt_of_formal_power_series.0.test.cpp
    title: remote_test/yosupo/math/sqrt_of_formal_power_series.0.test.cpp
  _isVerificationFailed: false
  _pathExtension: hpp
  _verificationStatusIcon: ':heavy_check_mark:'
  attributes:
    links: []
  bundledCode: "#line 1 \"math/truncated_fourier_transform.hpp\"\n\n\n\n#line 1 \"\
    common.hpp\"\n\n\n\n#define LIB_DEBUG\n\n#define LIB_BEGIN namespace lib {\n#define\
    \ LIB_END }\n#define LIB ::lib::\n\n\n#line 1 \"math/radix2_ntt.hpp\"\n\n\n\n\
    #line 5 \"math/radix2_ntt.hpp\"\n\n#include <algorithm>\n#include <array>\n#include\
    \ <cassert>\n#include <type_traits>\n#include <vector>\n\nLIB_BEGIN\n\nnamespace\
    \ detail {\n\ntemplate <typename IntT>\nconstexpr std::enable_if_t<std::is_integral_v<IntT>,\
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
    \ end - beg); }\n// clang-format on\n\ntemplate <typename ModIntT>\nvoid dft_doubling(const\
    \ std::vector<ModIntT> &a, std::vector<ModIntT> &dft_a) {\n  static constexpr\
    \ auto rt = detail::root<ModIntT>();\n  int as = static_cast<int>(a.size()), n\
    \ = static_cast<int>(dft_a.size());\n  // `dft_a` = dft_n(`a` mod (x^n - 1))\n\
    \  // doubling `dft_a` is just computing dft_n((`a` mod (x^n + 1))(\u03B6^(2n))).\n\
    \  dft_a.resize(n << 1);\n  auto it = dft_a.begin() + n;\n  for (int i = 0, is_even\
    \ = 0, j; i != as; ++i) {\n    if ((j = i & (n - 1)) == 0) is_even ^= 1;\n   \
    \ it[j] += is_even ? a[i] : -a[i];\n  }\n  ModIntT r(n == 1 ? ModIntT(-1) : rt[detail::bsf(n)\
    \ - 1]), v(1);\n  for (int i = 0; i != n; ++i) it[i] *= v, v *= r;\n  dft_n(it,\
    \ n);\n}\n\ntemplate <typename ModIntT>\nvoid dft_doubling(std::vector<ModIntT>\
    \ &dft_a) {\n  static constexpr auto rt = detail::root<ModIntT>();\n  int n  \
    \                  = static_cast<int>(dft_a.size());\n  dft_a.resize(n << 1);\n\
    \  auto it = dft_a.begin() + n;\n  std::copy_n(dft_a.cbegin(), n, it);\n  idft_n(it,\
    \ n);\n  ModIntT r(n == 1 ? ModIntT(-1) : rt[detail::bsf(n) - 1]), v(1);\n  for\
    \ (int i = 0; i != n; ++i) it[i] *= v, v *= r;\n  dft_n(it, n);\n}\n\nLIB_END\n\
    \n\n#line 6 \"math/truncated_fourier_transform.hpp\"\n\n#line 8 \"math/truncated_fourier_transform.hpp\"\
    \n#include <utility>\n#line 10 \"math/truncated_fourier_transform.hpp\"\n\nLIB_BEGIN\n\
    \ntemplate <typename ContainerT>\nvoid tft(ContainerT &&a) {\n  using Container\
    \          = std::remove_cv_t<std::remove_reference_t<ContainerT>>;\n  using T\
    \                  = typename Container::value_type;\n  static constexpr auto\
    \ rt = detail::root<T>();\n  static std::vector<T> root(1);\n  const int n = static_cast<int>(a.size());\n\
    \  if ((n & (n - 1)) == 0) return dft(std::forward<ContainerT>(a));\n  const int\
    \ len = ntt_len(n);\n  if (int s = static_cast<int>(root.size()); s << 1 < len)\
    \ {\n    root.resize(len >> 1);\n    for (int i = detail::bsf(s), j; 1 << i <\
    \ len >> 1; ++i) {\n      root[j = 1 << i] = rt[i];\n      for (int k = j + 1;\
    \ k < j << 1; ++k) root[k] = root[k - j] * root[j];\n    }\n  }\n  a.resize(len);\n\
    \  for (int j = 0, l = len >> 1; j != l; ++j) {\n    T u(a[j]), v(a[j + l]);\n\
    \    a[j] = u + v, a[j + l] = u - v;\n  }\n  for (int i = len >> 1; i >= 2; i\
    \ >>= 1) {\n    for (int j = 0, l = i >> 1; j != l; ++j) {\n      T u(a[j]), v(a[j\
    \ + l]);\n      a[j] = u + v, a[j + l] = u - v;\n    }\n    for (int j = i, l\
    \ = i >> 1, m = 1; j < n && j != len; j += i, ++m)\n      for (int k = j; k !=\
    \ j + l; ++k) {\n        T u(a[k]), v(a[k + l] * root[m]);\n        a[k] = u +\
    \ v, a[k + l] = u - v;\n      }\n  }\n  a.resize(n);\n}\n\ntemplate <typename\
    \ ContainerT>\nvoid itft(ContainerT &&a) {\n  using Container           = std::remove_cv_t<std::remove_reference_t<ContainerT>>;\n\
    \  using T                   = typename Container::value_type;\n  static constexpr\
    \ auto rt  = detail::root<T>();\n  static constexpr auto irt = detail::iroot<T>();\n\
    \  static std::vector<T> root{T(1)}, iroot{T(1)};\n  const int n = static_cast<int>(a.size());\n\
    \  if ((n & (n - 1)) == 0) return idft(std::forward<ContainerT>(a));\n  const\
    \ int len = ntt_len(n);\n  if (int s = static_cast<int>(root.size()); s << 1 <\
    \ len) {\n    root.resize(len >> 1);\n    iroot.resize(len >> 1);\n    for (int\
    \ i = detail::bsf(s), j; 1 << i < len >> 1; ++i) {\n      root[j = 1 << i] = rt[i],\
    \ iroot[j] = irt[i];\n      for (int k = j + 1; k < j << 1; ++k)\n        root[k]\
    \ = root[k - j] * root[j], iroot[k] = iroot[k - j] * iroot[j];\n    }\n  }\n \
    \ a.resize(len);\n  struct itft_rec {\n    itft_rec(Container &a) : a_(a), i2_(T(2).inv())\
    \ {}\n    // [`head`, `tail`), [`tail`, `last`)\n    void run(int head, int tail,\
    \ int last) {\n      if (head >= tail) return;\n      if (int mid = (last - head)\
    \ / 2 + head, len = mid - head; mid <= tail) {\n        // pull up [`head`, `mid`)\n\
    \        T i2p(1);\n        for (int i = 1; i != len; i <<= 1, i2p *= i2_)\n \
    \         for (int j = head, m = head / (i << 1); j != mid; j += i << 1, ++m)\n\
    \            for (int k = j; k != j + i; ++k) {\n              T u(a_[k]), v(a_[k\
    \ + i]);\n              a_[k] = u + v, a_[k + i] = (u - v) * iroot[m];\n     \
    \       }\n        for (int i = head; i != mid; ++i) a_[i] *= i2p;\n        //\
    \ push down [`tail`, `last`)\n        T r(root[head / (len << 1)] * 2);\n    \
    \    for (int i = tail; i != last; ++i) a_[i] = a_[i - len] - a_[i] * r;\n   \
    \     run(mid, tail, last);\n        // pull up [`head`, `last`)\n        r =\
    \ iroot[head / (len << 1)] * i2_;\n        for (int i = head; i != mid; ++i) {\n\
    \          T u(a_[i]), v(a_[i + len]);\n          a_[i] = (u + v) * i2_, a_[i\
    \ + len] = (u - v) * r;\n        }\n      } else {\n        T r(root[head / (len\
    \ << 1)]);\n        // push down [`tail`, `mid`)\n        for (int i = tail; i\
    \ != mid; ++i) a_[i] += a_[i + len] * r;\n        run(head, tail, mid);\n    \
    \    // pull up [`head`, `mid`)\n        for (int i = head; i != mid; ++i) a_[i]\
    \ -= a_[i + len] * r;\n      }\n    }\n    Container &a_;\n    const T i2_;\n\
    \  } rec(a);\n  rec.run(0, n, len);\n  a.resize(n);\n}\n\nLIB_END\n\n\n"
  code: "#ifndef TRUNCATED_FOURIER_TRANSFORM_HPP\n#define TRUNCATED_FOURIER_TRANSFORM_HPP\n\
    \n#include \"../common.hpp\"\n#include \"radix2_ntt.hpp\"\n\n#include <type_traits>\n\
    #include <utility>\n#include <vector>\n\nLIB_BEGIN\n\ntemplate <typename ContainerT>\n\
    void tft(ContainerT &&a) {\n  using Container          = std::remove_cv_t<std::remove_reference_t<ContainerT>>;\n\
    \  using T                  = typename Container::value_type;\n  static constexpr\
    \ auto rt = detail::root<T>();\n  static std::vector<T> root(1);\n  const int\
    \ n = static_cast<int>(a.size());\n  if ((n & (n - 1)) == 0) return dft(std::forward<ContainerT>(a));\n\
    \  const int len = ntt_len(n);\n  if (int s = static_cast<int>(root.size()); s\
    \ << 1 < len) {\n    root.resize(len >> 1);\n    for (int i = detail::bsf(s),\
    \ j; 1 << i < len >> 1; ++i) {\n      root[j = 1 << i] = rt[i];\n      for (int\
    \ k = j + 1; k < j << 1; ++k) root[k] = root[k - j] * root[j];\n    }\n  }\n \
    \ a.resize(len);\n  for (int j = 0, l = len >> 1; j != l; ++j) {\n    T u(a[j]),\
    \ v(a[j + l]);\n    a[j] = u + v, a[j + l] = u - v;\n  }\n  for (int i = len >>\
    \ 1; i >= 2; i >>= 1) {\n    for (int j = 0, l = i >> 1; j != l; ++j) {\n    \
    \  T u(a[j]), v(a[j + l]);\n      a[j] = u + v, a[j + l] = u - v;\n    }\n   \
    \ for (int j = i, l = i >> 1, m = 1; j < n && j != len; j += i, ++m)\n      for\
    \ (int k = j; k != j + l; ++k) {\n        T u(a[k]), v(a[k + l] * root[m]);\n\
    \        a[k] = u + v, a[k + l] = u - v;\n      }\n  }\n  a.resize(n);\n}\n\n\
    template <typename ContainerT>\nvoid itft(ContainerT &&a) {\n  using Container\
    \           = std::remove_cv_t<std::remove_reference_t<ContainerT>>;\n  using\
    \ T                   = typename Container::value_type;\n  static constexpr auto\
    \ rt  = detail::root<T>();\n  static constexpr auto irt = detail::iroot<T>();\n\
    \  static std::vector<T> root{T(1)}, iroot{T(1)};\n  const int n = static_cast<int>(a.size());\n\
    \  if ((n & (n - 1)) == 0) return idft(std::forward<ContainerT>(a));\n  const\
    \ int len = ntt_len(n);\n  if (int s = static_cast<int>(root.size()); s << 1 <\
    \ len) {\n    root.resize(len >> 1);\n    iroot.resize(len >> 1);\n    for (int\
    \ i = detail::bsf(s), j; 1 << i < len >> 1; ++i) {\n      root[j = 1 << i] = rt[i],\
    \ iroot[j] = irt[i];\n      for (int k = j + 1; k < j << 1; ++k)\n        root[k]\
    \ = root[k - j] * root[j], iroot[k] = iroot[k - j] * iroot[j];\n    }\n  }\n \
    \ a.resize(len);\n  struct itft_rec {\n    itft_rec(Container &a) : a_(a), i2_(T(2).inv())\
    \ {}\n    // [`head`, `tail`), [`tail`, `last`)\n    void run(int head, int tail,\
    \ int last) {\n      if (head >= tail) return;\n      if (int mid = (last - head)\
    \ / 2 + head, len = mid - head; mid <= tail) {\n        // pull up [`head`, `mid`)\n\
    \        T i2p(1);\n        for (int i = 1; i != len; i <<= 1, i2p *= i2_)\n \
    \         for (int j = head, m = head / (i << 1); j != mid; j += i << 1, ++m)\n\
    \            for (int k = j; k != j + i; ++k) {\n              T u(a_[k]), v(a_[k\
    \ + i]);\n              a_[k] = u + v, a_[k + i] = (u - v) * iroot[m];\n     \
    \       }\n        for (int i = head; i != mid; ++i) a_[i] *= i2p;\n        //\
    \ push down [`tail`, `last`)\n        T r(root[head / (len << 1)] * 2);\n    \
    \    for (int i = tail; i != last; ++i) a_[i] = a_[i - len] - a_[i] * r;\n   \
    \     run(mid, tail, last);\n        // pull up [`head`, `last`)\n        r =\
    \ iroot[head / (len << 1)] * i2_;\n        for (int i = head; i != mid; ++i) {\n\
    \          T u(a_[i]), v(a_[i + len]);\n          a_[i] = (u + v) * i2_, a_[i\
    \ + len] = (u - v) * r;\n        }\n      } else {\n        T r(root[head / (len\
    \ << 1)]);\n        // push down [`tail`, `mid`)\n        for (int i = tail; i\
    \ != mid; ++i) a_[i] += a_[i + len] * r;\n        run(head, tail, mid);\n    \
    \    // pull up [`head`, `mid`)\n        for (int i = head; i != mid; ++i) a_[i]\
    \ -= a_[i + len] * r;\n      }\n    }\n    Container &a_;\n    const T i2_;\n\
    \  } rec(a);\n  rec.run(0, n, len);\n  a.resize(n);\n}\n\nLIB_END\n\n#endif"
  dependsOn:
  - common.hpp
  - math/radix2_ntt.hpp
  isVerificationFile: false
  path: math/truncated_fourier_transform.hpp
  requiredBy:
  - math/polynomial.hpp
  - math/convolution.hpp
  - math/truncated_formal_power_series.hpp
  timestamp: '2022-05-01 20:49:46+08:00'
  verificationStatus: LIBRARY_ALL_AC
  verifiedWith:
  - remote_test/yosupo/math/convolution_mod.1.test.cpp
  - remote_test/yosupo/math/pow_of_formal_power_series.1.test.cpp
  - remote_test/yosupo/math/convolution_mod_1000000007.0.test.cpp
  - remote_test/yosupo/math/convolution_mod.0.test.cpp
  - remote_test/yosupo/math/exp_of_formal_power_series.1.test.cpp
  - remote_test/yosupo/math/inv_of_formal_power_series.2.test.cpp
  - remote_test/yosupo/math/division_of_polynomials.0.test.cpp
  - remote_test/yosupo/math/log_of_formal_power_series.1.test.cpp
  - remote_test/yosupo/math/sqrt_of_formal_power_series.0.test.cpp
documentation_of: math/truncated_fourier_transform.hpp
layout: document
title: Truncated Fourier Transform
---

A very interesting algorithm. The recursion steps are elegant.

## Reference

1. J. van der Hoeven. The Truncated Fourier Transform and Applications.
2. Li Zhang. Implementation Techniques for the Truncated Fourier Transform.