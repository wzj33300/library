#pragma once

#include "fft.hpp"
#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

// returns coefficients generated by closure
// closure: gen(index, current_product)
template <typename Tp, typename Closure>
inline std::enable_if_t<std::is_invocable_r_v<Tp, Closure, int, const std::vector<Tp> &>,
                        std::vector<Tp>>
semi_relaxed_convolution(const std::vector<Tp> &A, Closure gen, int n) {
    enum { BaseCaseSize = 32 };
    static_assert((BaseCaseSize & (BaseCaseSize - 1)) == 0);

    static const int Block[]     = {16, 16, 16, 16, 16};
    static const int BlockSize[] = {
        BaseCaseSize,
        BaseCaseSize * Block[0],
        BaseCaseSize * Block[0] * Block[1],
        BaseCaseSize * Block[0] * Block[1] * Block[2],
        BaseCaseSize * Block[0] * Block[1] * Block[2] * Block[3],
        BaseCaseSize * Block[0] * Block[1] * Block[2] * Block[3] * Block[4],
    };

    // returns (which_block, level)
    auto blockinfo = [](int ind) {
        int i = ind / BaseCaseSize, lv = 0;
        while ((i & (Block[lv] - 1)) == 0) i /= Block[lv++];
        return std::make_pair(i & (Block[lv] - 1), lv);
    };

    std::vector<Tp> B(n), AB(n);
    std::vector<std::vector<std::vector<Tp>>> dftA, dftB;

    for (int i = 0; i < n; ++i) {
        const int s = i & (BaseCaseSize - 1);

        // blocked contribution
        if (i >= BaseCaseSize && s == 0) {
            const auto [j, lv]  = blockinfo(i);
            const int blocksize = BlockSize[lv];

            if (blocksize * j == i) {
                if ((int)dftA.size() == lv) {
                    dftA.emplace_back();
                    dftB.emplace_back(Block[lv] - 1);
                }
                if ((j - 1) * blocksize < (int)A.size()) {
                    dftA[lv]
                        .emplace_back(A.begin() + (j - 1) * blocksize,
                                      A.begin() + std::min((j + 1) * blocksize, (int)A.size()))
                        .resize(blocksize * 2);
                    fft(dftA[lv][j - 1]);
                } else {
                    dftA[lv].emplace_back(blocksize * 2);
                }
            }

            dftB[lv][j - 1].resize(blocksize * 2);
            std::copy_n(B.begin() + (i - blocksize), blocksize, dftB[lv][j - 1].begin());
            std::fill_n(dftB[lv][j - 1].begin() + blocksize, blocksize, Tp());
            fft(dftB[lv][j - 1]);

            // middle product
            std::vector<Tp> mp(blocksize * 2);
            for (int k = 0; k < j; ++k)
                for (int l = 0; l < blocksize * 2; ++l)
                    mp[l] += dftA[lv][j - 1 - k][l] * dftB[lv][k][l];
            inv_fft(mp);

            for (int k = 0; k < blocksize && i + k < n; ++k) AB[i + k] += mp[k + blocksize];
        }

        // basecase contribution
        for (int j = std::max(i - s, i - (int)A.size() + 1); j < i; ++j) AB[i] += A[i - j] * B[j];
        B[i] = gen(i, AB);
        if (!A.empty()) AB[i] += A[0] * B[i];
    }

    return B;
}
