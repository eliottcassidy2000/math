#!/usr/bin/env python3
"""Lane twoapex - checks added by the FINISHING agent (the lane was interrupted after d1 reported one failed check).

usage: python z1_finish_checks.py [STATE_FILE]

(A) The two-block family F1 = (02, 23).  d1_growth_closed_forms.py tested
        H(0^p a'', b' 1^q) = x(b') x(a'') + (h(b') + 1)(C(a'') + 1)
    and found 4096 violations.  The identity misses one term.  Correct statement (proof below):
        H = x(b') x(a'') + (h(b') + 1)(C(a'') + 1) - [p = 0 or b'_p = 0] [q = 0 or a''_1 = 1].
    Proof.  Block 1 = s_1..s_p (letters (0, b'_i)), block 2 = s_(p+1)..s_n (letters (a''_j, 1)).  T1 = block 1 + v
    and T2 = block 2 + u are one-apex tournaments with counts x(b'), x(a'').  Every arc between T1 and T2 goes
    from T1 to T2 except u -> v.  A Hamiltonian path that does not use u -> v is a Hamiltonian path of T1 followed
    by one of T2 whose junction is not the (absent) arc v -> u: x(b') x(a'') - e1 e2, where e1 = number of paths of
    T1 ending at v = [p = 0 or s_p -> v] = [p = 0 or b'_p = 0] and e2 = number of paths of T2 starting at
    u = [q = 0 or u -> s_(p+1)] = [q = 0 or a''_1 = 1].  A path that uses u -> v is A1 B1 u v A2 B2 with (A1, A2)
    increasing runs partitioning block 1, (B1, B2) partitioning block 2 (one cannot return from block 2 or u to
    block 1 except through u -> v); the only conditions are: B1 empty or last(B1) -> u (a'' = 0 there), A2 empty or
    v -> first(A2) (b' = 1 there); all other junctions are arcs (block 1 -> block 2, block 1 -> u, v -> block 2).
    Choosing first(A2) = F with b'_F = 1 leaves 2^(p - F) splittings, so (A1, A2) has 1 + h(b') choices;
    choosing last(B1) = L with a''_L = 0 leaves 2^(L - 1), so (B1, B2) has 1 + C(a'') choices.  QED.
    Three bilinear terms: this is the seam rank 3 reported for (02, 23) by e2_two_block_templates.py and the form
    H = g(u) g(v) - e(u) e(v) + (1 + u)(1 + v) of THM-4875.
    Checked here on ALL instances with p, q <= 7.

(B) Exact first moment from the 17-dimensional representation (Python integers), n = 1..60:
        sum_{|w| = n} H(w) = 4^n (3^n + 3 * 2^n + 3) / 8.
    (Direct proof: for a uniform word the three runs R1, R2, R3 of an arrangement impose independent conditions
    of probability 1/2, 1/4, 1/2 when nonempty; R2 empty is allowed for one arrangement only; summing over the
    3-colourings by their set of nonempty runs gives 2[(3^n - 3 2^n + 3)/16 + (2^n - 2)/4 + 1/4] + (2^n - 2)/4 + 1.)
    Exact second moment through the Kronecker square (289 dimensions): E[H^2] / 9^n for n up to 60.

(C) The words 3^p 0^p with exact integers (d1 used int64 and overflowed for n >= 40): the constant
        (1 - H(3^p 0^p) / (2 * 3^n)) / (2/3)^(n/2).

(D) The initial state is fixed by the letter 0 (x(empty) A_0 = x(empty), exactly), so the Markov chain of prefix
    states modulo m has a self-loop at the start state: aperiodic on its (irreducible) orbit.

(E) Dispersion of the multiplicities (needs the census state of b1_census.py): for each dyadic window of odd N,
    prim(N) = number of primitive words (first letter not 0, last letter not 3) of length <= L with H = N:
    mean, variance within blocks of 256 consecutive odd N, variance / mean, lag-1 correlation, and the
    Chebyshev bound variance / mean^2 on the share of N with prim(N) = 0.
"""
import itertools
import json
import sys
from fractions import Fraction

import numpy as np

from ta_core import H36
from ta_fast import H24


def x_one_apex(a):
    n = len(a)
    if n == 0:
        return 1
    t = (1 if a[0] == 1 else 0) + (1 if a[-1] == 0 else 0)
    for k in range(n - 1):
        if a[k] == 0 and a[k + 1] == 1:
            t += 1
    for F in range(n):
        if a[F]:
            for L in range(F + 1, n):
                if a[L] == 0:
                    t += 1 << (L - F - 1)
    return t


def main():
    ok = True
    # ---- (A)
    bad_old = bad_new = cnt = 0
    PQ = 7
    for p in range(0, PQ + 1):
        for q in range(0, PQ + 1):
            for bp in itertools.product((0, 1), repeat=p):
                x1 = x_one_apex(bp)
                y1 = sum(bp[F] << (p - 1 - F) for F in range(p)) + 1
                e1 = 1 if (p == 0 or bp[-1] == 0) else 0
                for aq in itertools.product((0, 1), repeat=q):
                    word = [2 * t for t in bp] + [t + 2 for t in aq]
                    x2 = x_one_apex(aq)
                    y2 = sum((1 - aq[L]) << L for L in range(q)) + 1
                    e2 = 1 if (q == 0 or aq[0] == 1) else 0
                    h = H24(word)
                    cnt += 1
                    if h != x1 * x2 + y1 * y2:
                        bad_old += 1
                    if h != x1 * x2 + y1 * y2 - e1 * e2:
                        bad_new += 1
    print(f"(A) F1 = (02, 23), all instances with p, q <= {PQ}: {cnt}; two-term formula of d1 fails on {bad_old}; "
          f"three-term formula x x' + (h + 1)(C + 1) - e e' fails on {bad_new}")
    ok &= bad_new == 0

    # ---- (B)
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    A = d["A"]
    x0 = d["x0"]
    S = [[sum(A[l][i][j] for l in range(4)) for j in range(r)] for i in range(r)]
    v = list(x0)
    bad = 0
    NB = 60
    for n in range(1, NB + 1):
        v = [sum(v[i] * S[i][j] for i in range(r) if S[i][j]) for j in range(r)]
        if 8 * v[0] != 4 ** n * (3 ** n + 3 * 2 ** n + 3):
            bad += 1
    print(f"(B) sum over the 4^n words of H = 4^n (3^n + 3 2^n + 3)/8 for n = 1..{NB} (exact, 17-dim representation): "
          f"violations {bad}")
    ok &= bad == 0
    # Kronecker square: state z[i * r + j] = sum_w x_i(w) x_j(w)
    z = [x0[i] * x0[j] for i in range(r) for j in range(r)]
    # sparse lists of the four matrices
    sp = [[[(j, A[l][i][j]) for j in range(r) if A[l][i][j]] for i in range(r)] for l in range(4)]
    line = []
    prev_ratio = None
    for n in range(1, NB + 1):
        new = [0] * (r * r)
        for l in range(4):
            sl = sp[l]
            for i in range(r):
                for j in range(r):
                    c = z[i * r + j]
                    if c:
                        for (i2, c1) in sl[i]:
                            cc = c * c1
                            base = i2 * r
                            for (j2, c2) in sl[j]:
                                new[base + j2] += cc * c2
        z = new
        ratio = Fraction(z[0], 36 ** n)
        if n in (5, 9, 13, 20, 30, 40, 50, 60):
            line.append(f"n={n}: {float(ratio):.10f}")
        if n == 13:
            e13 = float(ratio)
        prev_ratio = ratio
    print("(B) E[H^2] / 9^n (exact rationals, shown as floats): " + "; ".join(line))
    print(f"    (d1 printed 0.032685 at n = 13; here {e13:.6f})")

    # ---- (C)
    for n in (20, 30, 40, 60, 80, 120):
        p = n // 2
        hv = H36([3] * p + [0] * (n - p))
        val = (1 - Fraction(hv, 2 * 3 ** n)) / Fraction(2, 3) ** p
        lb = 2 * 3 ** n - 4 * (3 ** p * 2 ** (n - p) + 2 ** p * 3 ** (n - p)) - 2 ** n
        print(f"(C) n = {n}: (1 - H(3^{p} 0^{p}) / (2 3^n)) / (2/3)^(n/2) = {float(val):.9f}; "
              f"2 3^n - 4 (3^p 2^q + 2^p 3^q) - 2^n <= H <= 2 3^n: {lb <= hv <= 2 * 3 ** n}")
        ok &= lb <= hv <= 2 * 3 ** n

    # ---- (D)
    fixed = [sum(x0[i] * A[0][i][j] for i in range(r)) for j in range(r)] == list(x0)
    print(f"(D) x(empty) A_0 = x(empty) exactly (self-loop of every chain modulo m at the start state): {fixed}")
    ok &= fixed

    # ---- (E)
    if len(sys.argv) > 1:
        st = np.load(sys.argv[1])
        last, prev, L = st["last"].astype(np.int64), st["prev"].astype(np.int64), int(st["upto"])
        prim = last - prev
        del last, prev
        print(f"(E) primitive words of length <= {L}; blocks of 256 consecutive odd N")
        print("    window          mean     var(in blocks)  var/mean  lag-1 corr   min    (min-mean)/sd   var/mean^2")
        B = 256
        for k in range(10, 26):
            lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
            seg = prim[lo:hi].astype(np.float64)
            nb = len(seg) // B
            blk = seg[: nb * B].reshape(nb, B)
            mu = blk.mean(axis=1, keepdims=True)
            res = blk - mu
            var = float((res ** 2).sum() / (nb * (B - 1)))
            mean = float(seg.mean())
            lag = float((res[:, 1:] * res[:, :-1]).sum() / (res ** 2).sum())
            mn = float(seg.min())
            print(f"    [2^{k:2d}, 2^{k + 1:2d}) {mean:9.2f} {var:14.2f} {var / mean:9.3f} {lag:10.4f} {int(mn):6d} "
                  f"{(mn - mean) / np.sqrt(var):10.2f} {var / mean ** 2:14.6f}")
    print("ALL CHECKS PASS" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
