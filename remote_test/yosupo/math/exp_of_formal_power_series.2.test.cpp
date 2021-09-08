#define PROBLEM "https://judge.yosupo.jp/problem/exp_of_formal_power_series"

#include <iostream>
#include <vector>

#include "math/formal_power_series/prime_binomial.hpp"
#include "math/formal_power_series/semi_relaxed_convolution.hpp"
#include "modint/Montgomery_modint.hpp"

int main() {
#ifdef LOCAL
  std::freopen("in", "r", stdin), std::freopen("out", "w", stdout);
#endif
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);
  using mint = lib::MontModInt<998244353>;
  int n;
  std::cin >> n;
  std::vector<mint> A(n), B;
  for (auto &i : A) std::cin >> i;
  for (int i = 1; i < n; ++i) A[i - 1] = A[i] * i;
  lib::PrimeBinomial<mint> bi(n);
  lib::SemiRelaxedConvolution<mint> rc(A, B, [&bi](int idx, const std::vector<mint> &contri) {
    if (idx == 0) return mint(1);
    return contri[idx - 1] * bi.inv_unsafe(idx);
  });
  while (n--) rc.next();
  for (auto i : B) std::cout << i << ' ';
  return 0;
}