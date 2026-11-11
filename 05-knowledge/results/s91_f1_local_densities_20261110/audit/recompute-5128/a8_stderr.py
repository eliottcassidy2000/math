# Auditor recompute-5128: empirical noise level of the window excess E at large primes (the theorem quotes a
# standard error 0.022 at p = 65537), from the distribution of the class sums over ALL classes of the window
# [2^24, 2^27); also the word-pair weighting (s + 1) of the counts.
import sys
import math
import numpy as np
from a0_defs import psi_table32
from a7_misses import rep_counts


def main():
    B = 27
    a = 24
    psi = psi_table32(((1 << B) + 1) // 3 + 4)
    r0, r1, sat = rep_counts(psi, B)
    i0, i1 = 1 << (a - 1), 1 << (B - 1)
    N = 2 * np.arange(i0, i1, dtype=np.int64) + 1
    s0 = r0[i0:i1].astype(np.float64)
    s1 = r1[i0:i1].astype(np.float64)
    Z = np.load("census_28.npz")
    all_p = [int(p) for p in Z["all_p"]]
    wb = list(range(a + 1, B + 1))
    print(f"== window [2^{a}, 2^{B}): class sums of r_sigma over all classes; distinguished class is M = N + [sigma = 0] = 0 (mod p) ==")
    print("   p     sheet  E(dist. class)   sd of E over the other classes   z      E_nd      z(E_nd)   min / max E over other classes   Poisson sd sqrt(p/T)")
    for p in [251, 257, 263, 677, 691, 1039, 2729, 2731, 2741, 8191, 65521, 65537, 65539]:
        j = all_p.index(p)
        for sheet in (0, 1):
            M = N + (1 if sheet == 0 else 0)
            w = s0 if sheet == 0 else s1
            C = np.bincount((M % p).astype(np.int64), weights=w, minlength=p)
            T = C.sum()
            others = (T - C[0]) / (p - 1)
            E = C / others - 1
            sd = np.std(E[1:])
            dsm = int(Z["Dsm"][sheet, wb, j].sum())
            dlg = int(Z["Dlg"][sheet, wb, j].sum())
            dbo = int(Z["Dboth"][sheet, wb, j].sum())
            ddg = int(Z["Ddiag"][sheet, wb, j].sum())
            anyd = dsm + dlg - dbo + ddg
            assert abs(C[0] - int(Z["dist"][sheet, wb, j].sum())) < 0.5
            End = (C[0] - anyd) / others - 1  # others unchanged: divisible triples are all in the distinguished class
            print(f" {p:6d}   {sheet}    {E[0]:+.4f}          {sd:.4f}                      {E[0] / sd:+6.2f}   {End:+.4f}   {End / sd:+6.2f}     {E[1:].min():+.4f} / {E[1:].max():+.4f}            {math.sqrt(p / T):.4f}")
    print()
    print("== word-pair weighting: sum over s of (s + 1) * #{triples with padding s, N < 2^b} ==")
    tots = Z["tots"]
    cs = np.cumsum(tots, axis=1)
    Bz = tots.shape[1] - 1
    sw = np.arange(tots.shape[0]) + 1
    Wc = (cs * sw[:, None]).sum(axis=0)
    for b in range(20, Bz + 1):
        X = 2.0 ** b
        slope = (Wc[b] / X - Wc[b - 1] / (X / 2)) / math.log(2)
        print(f"   b = {b}: weighted count {int(Wc[b])}; /(X ln X) = {Wc[b] / (X * math.log(X)):.6f}; two-point slope = {slope:.6f} (pi^2/12 = {math.pi ** 2 / 12:.6f}); share of sheet 0 = {cs[0, b] / Wc[b]:.5f} (limit 0.19641)")


if __name__ == "__main__":
    main()
