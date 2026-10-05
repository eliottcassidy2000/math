#!/usr/bin/env python3
"""
Exact pair-coherence census of the q-towers (opus, 2026-10-05).

For the uniform-start q-tower, f_N(0) = sum_a w_a e(R Phi_q(a)), w_a = 2^-|a|, Phi_q(a) = sum_i q^(i-1) 2^-S_i (mod 1),
and E_R |f_N(0)|^2 = sum_a w_a^2 = 3^-N exactly for Haar R (Fourier note 4k: Phi_q is injective on paths).
The frequencies of 2-adic valuation v see only the binary digits of Phi beyond position v:
    E_{R'} |f_N at frequency 2^v R'|^2 = D_v := sum_{a,b} w_a w_b 1[ 2^v (Phi(a) - Phi(b)) in Z ],
the weighted mass of path PAIRS whose value difference has denominator dividing 2^v ("coherent pairs": their relative
phase takes at most 2^v values as R varies).  D_0 = D_1 = 3^-N (injectivity; denominators >= 4), D_v -> 1 as v -> inf.
This script computes D_v exactly (Fractions, valuations <= cap) for several multipliers q, i.e. the profile of the
mean square over the dyadic shells of the frequency, which is where the towers can differ while all sharing 3^-N.
Usage: python collatz_pascal_tower_pair_coherence_20261005.py [N=6] [cap=10]
"""
import sys, math, time
from fractions import Fraction
from itertools import product
from collections import defaultdict

def census(q, N, cap):
    # enumerate paths; group by Phi mod 1 truncated to v digits for v = 0..V
    V = 2 * N + cap + 2
    # value as an integer numerator over 2^V... Phi has denominator <= 2^(N*cap); use exact Fractions mod 1
    masses = [defaultdict(float) for _ in range(V + 1)]
    deep = [defaultdict(float) for _ in range(V + 1)]     # classes of frac(2^v Phi): digits beyond v agree
    tot_w = 0.0
    for cs in product(range(1, cap + 1), repeat=N):
        S = 0; phi = Fraction(0)
        for i, c in enumerate(cs):
            S += c
            phi += Fraction(q ** i, 2 ** S)
        phi = phi - (phi.numerator // phi.denominator)      # mod 1, in [0,1)
        w = 2.0 ** (-S)
        tot_w += w
        # class of phi modulo 2^-v: floor(phi * 2^v)
        num, den = phi.numerator, phi.denominator          # den = 2^e
        e = den.bit_length() - 1
        for v in range(V + 1):
            if v >= e:
                key = num << (v - e)
            else:
                key = num >> (e - v)
            masses[v][key] += w
            # deep class: frac(2^v Phi) = ((num << v) mod den) / den, as a reduced pair
            if v >= e:
                dkey = (0, 1)
            else:
                nn = (num << v) % den
                g = math.gcd(nn, den)
                dkey = (nn // g, den // g)
            deep[v][dkey] += w
    C = [sum(m * m for m in masses[v].values()) for v in range(V + 1)]
    D = [sum(m * m for m in deep[v].values()) for v in range(V + 1)]
    return C, D, tot_w

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    out = []
    def P(s): print(s, flush=True); out.append(s)
    P(f"pair-coherence census: N={N}, valuations <= {cap} (weight kept {1 - N * 2.0 ** -cap:.6f} of 1); 3^-N = {3.0 ** -N:.3e}")
    P("C_v = collision mass of the FIRST v binary digits of Phi mod 1 = mean square of f over the low frequencies R < 2^v (coarse-scale concentration of the law of Phi mod 1; uniform law: C_v = 2^-v);")
    P("D_v = collision mass of the digits BEYOND position v = pair mass with 2^v (Phi_a - Phi_b) integral = mean square over the frequencies 2^v R' (D_0 = D_1 = sum w_a^2, injectivity).")
    t0 = time.time()
    for q in (1, -1, 3, -3, 5, -5, 7, -7, 9, 11, 13, 15, 17):
        C, D, tw = census(q, N, cap)
        P(f"q={q:>3}: C_v/3^-N v=0..{len(C)-1}: " + " ".join(f"{d * 3.0 ** N:.3f}" for d in C) + f"   [{time.time()-t0:.0f}s]")
        P(f"       D_v/3^-N v=0..{len(D)-1}: " + " ".join(f"{d * 3.0 ** N:.3f}" for d in D))
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
