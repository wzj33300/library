#ifndef FPS_COMPOSITION_HPP
#define FPS_COMPOSITION_HPP

#include "../common.hpp"
#include "taylor_shift.hpp"
#include "truncated_formal_power_series.hpp"

#include <algorithm>
#include <utility>
#include <vector>

LIB_BEGIN

// returns f(g) mod x^n
// reference: noshi91's blog: https://noshi91.hatenablog.com/entry/2024/03/16/224034
// for detailed definition of this code check pseudocode written by me:
// https://codeforces.com/blog/entry/127674
template <typename ModIntT>
tfps<ModIntT> composition(const tfps<ModIntT> &f, const tfps<ModIntT> &g, int n) {
  if (g.empty() || n <= 0) return {};
  struct composition_rec {
    composition_rec(tfps<ModIntT> &&P) : P_(std::move(P)) {}
    tfps<ModIntT> run(const tfps<ModIntT> Q, int d, int n) {
      // [0,n] => [y^0]Q, [n+1,2n+1] => [y^1]Q, ...
      assert(static_cast<int>(Q.size()) == (d + 1) * (n + 1));
      if (n == 0) {
        tfps<ModIntT> res(d);
        std::copy_n(P_.cbegin(), std::min(P_.size(), res.size()), res.rbegin());
        return res;
      }
      // let y=x^(2n+2) => [0,2n+2) = [y^0]Q, ...
      // y^0[0,2n+2), y^1[2n+2,4n+4), ..., y^(2d)[2d(2n+2),(2d+1)(2n+2)-1)
      const int len = ntt_len((d << 1 | 1) * ((n + 1) << 1) - 1);
      tfps<ModIntT> dftQ(len), VV(len >> 1), V((d << 1 | 1) * ((n >> 1) + 1));
      for (int i = 0; i <= d; ++i)
        std::copy_n(Q.cbegin() + i * (n + 1), n + 1, dftQ.begin() + i * ((n + 1) << 1));
      dft(dftQ);
      // apply dft trick from Bostan&Mori's paper
      for (int i = 0; i != len; i += 2) VV[i >> 1] = dftQ[i] * dftQ[i + 1];
      idft(VV);
      for (int i = 0; i <= d << 1; ++i)
        for (int j = 0; j <= n >> 1; ++j) V[i * ((n >> 1) + 1) + j] = VV[i * (n + 1) + j];
      const auto TT = run(std::move(V), d << 1, n >> 1);
      VV.assign(len >> 1, ModIntT());
      for (int i = 0; i != d << 1; ++i)
        for (int j = 0; j <= n >> 1; ++j) VV[i * (n + 1) + j] = TT[i * ((n >> 1) + 1) + j];
      dft(VV);
      auto &&T = dftQ;
      // apply dft trick from Bostan&Mori's paper
      for (int i = 0; i != len; i += 2) {
        auto l = VV[i >> 1] * dftQ[i + 1], r = VV[i >> 1] * dftQ[i];
        T[i] = l, T[i + 1] = r;
      }
      idft(T);
      auto &&U = VV;
      U.assign(d * (n + 1), ModIntT());
      // [y^(-d+1..0)]T => [y^(d..2d-1)]T
      for (int i = 0; i != d; ++i)
        for (int j = 0; j <= n; ++j) U[i * (n + 1) + j] = T[(i + d) * ((n + 1) << 1) + j];
      return U;
    }

  private:
    const tfps<ModIntT> P_;
  } a(taylor_shift(f, g.front()));
  tfps<ModIntT> Q(n << 1); // [0,n)=1, [n,2n)=-g
  Q.front()   = ModIntT(1);
  const int s = static_cast<int>(g.size());
  for (int i = n + 1; i - n < s && i != n << 1; ++i) Q[i] = -g[i - n];
  return a.run(std::move(Q), 1, n - 1);
}

LIB_END

#endif