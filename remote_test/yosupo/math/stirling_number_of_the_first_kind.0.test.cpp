#define PROBLEM "https://judge.yosupo.jp/problem/stirling_number_of_the_first_kind"

#include <iostream>

#include "math/famous_sequence/Stirling_numbers.hpp"
#include "math/formal_power_series/polynomial.hpp"
#include "modint/Montgomery_modint.hpp"

int main() {
#ifdef LOCAL
  std::freopen("in", "r", stdin), std::freopen("out", "w", stdout);
#endif
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);
  using poly = lib::Poly<lib::MontModInt<998244353>>;
  int n;
  std::cin >> n;
  poly res;
  lib::Stirling1st_row(n, res);
  for (int i = 0; i <= n; ++i) std::cout << ((n - i & 1) ? (-res[i]) : res[i]) << ' ';
  return 0;
}