---
title: Arbitrary Modulo Square Matrix
documentation_of: ./arbitrary_modulo_square_matrix.hpp
---

## 任意模数方阵的 Gauss 消元

在矩阵 Gauss 消元的过程中，我们想利用第 $i$ 行去消去第 $j$ 行使得第 $j$ 行第 $i$ 个元素为零，此时因为没有办法使用乘法逆元（模数非素数），需要进行类似于辗转相除的过程，为了使得这个过程尽可能短，最好选择一列中“数值”非零且最小的当做主元。

因为实际上只是对一行中的主元和另一行同一列的两个元素进行辗转相除，那么将其表为线性组合后辗转相除计算最后的线性组合后再应用于这两行即可，见 <https://github.com/hly1204/library/issues/13> 。

同样的我们也可以使用这种方法使得任意模数的方阵变为上 Hessenberg 矩阵，而通过上 Hessenberg 矩阵求出特征多项式是没有除法的，沿用之前的方法即可。