#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adversarial_two_ursell_kps.py

Scrutinize the worker's SUBTLE claims:
 (A) c_2 = alpha_2 - alpha_1^2/2 = -|E| - alpha_1/2 = b_2 - alpha_1/2.
 (B) The "two Ursell families" distinction: is the GRAPHICAL b_k (Mayer connected
     subgraph) really DIFFERENT from the analytic Taylor cumulant c_k of log I?
     The worker says b_2 = -|E| but c_2 = -|E| - alpha_1/2.
 (C) Is the worker's claim that "log I(Omega,z) = sum b_k z^k" honest? The worker's
     OWN script comment (line 20) literally writes "log I(Omega,z)= sum b_k z^k"
     for the GRAPHICAL b_k. But analytically log I has Taylor coeffs c_k, NOT b_k.
     So which is the true Taylor expansion of log I?

I test on several conflict graphs by directly Taylor-expanding log(I(Omega,z)).
"""
from fractions import Fraction


def taylor_log_poly(alpha, K):
    """Taylor coefficients d_k of log( sum_k alpha_k z^k ), alpha_0=1."""
    a = [Fraction(c) for c in alpha] + [Fraction(0)] * (K + 2)
    assert a[0] == 1
    d = [Fraction(0)] * (K + 1)
    # log(1+u) where I = 1 + (a1 z + a2 z^2 + ...). Use recurrence:
    # I = exp(L), L = sum d_k z^k => I' = L' I => m a_m = sum_{k=1..m} k d_k a_{m-k}
    for m in range(1, K + 1):
        s = Fraction(m) * a[m]
        for k in range(1, m):
            s -= Fraction(k) * d[k] * a[m - k]
        d[m] = s / Fraction(m)
    return d[1:K + 1]


def main():
    print("Testing the 'two Ursell families' claim on explicit conflict graphs.\n")
    # Build a few alpha vectors directly:
    # Case 1: K3 (3 pairwise-adjacent): alpha = [1,3,0]; E=3, alpha1=3, alpha2=0.
    # Case 2: path P3 (3 vtx, 2 edges): alpha = [1,3,1]; E=2.
    # Case 3: edgeless 3 vtx: alpha=[1,3,3,1]; E=0.
    # Case 4: Paley T_7 conflict: alpha=[1,80,7]; E = C(80,2)-7 = 3160-7=3153.
    cases = {
        "K3 (E=3)": ([1, 3, 0], 3),
        "P3 path (E=2)": ([1, 3, 1], 2),
        "edgeless3 (E=0)": ([1, 3, 3, 1], 0),
        "PaleyT7 (E=3153)": ([1, 80, 7], 3160 - 7),
        "single edge (E=1)": ([1, 2, 0], 1),
        "C4 (E=4)": ([1, 4, 2], 4),  # 4-cycle: alpha2 = 2 (two opposite pairs)
    }
    print(f"{'graph':22s} {'alpha1':>6s} {'alpha2':>6s} {'E':>6s} "
          f"{'b2=-E':>8s} {'c2(Taylor)':>12s} {'b2-a1/2':>10s} {'match?':>7s}")
    for name, (alpha, E) in cases.items():
        a1 = alpha[1] if len(alpha) > 1 else 0
        a2 = alpha[2] if len(alpha) > 2 else 0
        d = taylor_log_poly(alpha, 3)
        c2 = d[1]  # coefficient of z^2 in log I
        b2 = Fraction(-E)
        predicted_c2 = b2 - Fraction(a1, 2)
        also = Fraction(a2) - Fraction(a1 * a1, 2)
        print(f"{name:22s} {a1:6d} {a2:6d} {E:6d} {str(b2):>8s} "
              f"{str(c2):>12s} {str(predicted_c2):>10s} "
              f"{str(c2 == predicted_c2 == also):>7s}")
    print()
    print("So: the TAYLOR coefficient of z^2 in log I is c2 = alpha2 - alpha1^2/2,")
    print("which EQUALS b2 - alpha1/2 (since b2=-E=alpha2-C(alpha1,2)=alpha2-alpha1^2/2+alpha1/2).")
    print("=> The literal Taylor log I has c_k, NOT the graphical b_k. The worker's")
    print("   claim that these are DIFFERENT series is mathematically correct.")
    print("   The script's docstring line 20 'log I = sum b_k z^k' is loose/wrong as")
    print("   a Taylor statement, but the worker's RETURNED text corrects it explicitly.")


if __name__ == "__main__":
    main()
