#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_wsb_verify_gapaudit_kps-S9-wf.py   (kind-pasteur-2026-06-19-S9-wf)

V8 GAP AUDIT (efficient).  The claimed advance asserts that even AFTER the
support-6 floor (Lemma A), the ABSOLUTE bound is still lossy at the AP:
  sum_{supp>=6 relations} |K(n)|  ~ 1.2  >  delta_8 ~ 0.30,
so a SECOND-LAYER signed cancellation among the support-6 terms is still
required.  The full box range(-9,10)^7 is ~9e8 tuples (infeasible); we
restrict to the relation LATTICE of consec_8 and enumerate ONLY relations
that can have support>=6 cheaply, by iterating over the 2D solution lattice.

For consec_8, e=(1,2,3,4,5,6,7).  A relation n in Z^7 with sum_i n_i*(i+1)=0.
We enumerate all n with |n_i|<=L and sum n_i (i+1)=0, and split by the count
of nonzero non-7 coords.  L=6,7,8 are feasible (15^7 ~ 1.7e8 is borderline;
we use a meet-in-the-middle / partial-sum DFS to stay fast).
"""
from __future__ import annotations
import sys, itertools, math, cmath
from fractions import Fraction as F
from math import comb

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

TWO_PI_I = 2j * math.pi

def shat(n, j):
    if n == 0: return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)

SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))

def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7: total += x1 - x0
    return total

def enumerate_relations(coeffs, L):
    """All n in Z^d with |n_i|<=L and sum n_i*coeffs[i]=0, via DFS with pruning."""
    d = len(coeffs)
    # max remaining |contribution| from suffix
    suffix_max = [0] * (d + 1)
    for i in range(d - 1, -1, -1):
        suffix_max[i] = suffix_max[i + 1] + L * abs(coeffs[i])
    out = []
    cur = [0] * d
    def dfs(i, partial):
        if i == d:
            if partial == 0:
                out.append(tuple(cur))
            return
        # prune: remaining can reach [-suffix_max[i], suffix_max[i]]
        if abs(partial) > suffix_max[i]:
            return
        for v in range(-L, L + 1):
            cur[i] = v
            dfs(i + 1, partial + v * coeffs[i])
        cur[i] = 0
    dfs(0, 0)
    return out


def main():
    print("V8 GAP AUDIT (efficient) -- consec_8 relation lattice")
    print("=" * 70)
    E = list(range(8)); nz = [e for e in E if e != 0]  # 1..7
    d8 = measS7(E) - M7(8)
    print(f"  delta_8 = meas(S7(consec_8)) - M7(8) = {d8} = {float(d8):.6f}")
    for L in (5, 6, 7, 8):
        rels = enumerate_relations(nz, L)
        abs_all = 0.0; abs_supp6 = 0.0; signed = 0j
        nontrivial = 0
        for n in rels:
            if all(x == 0 for x in n): continue
            nontrivial += 1
            kv = Kk(list(n)); signed += kv; abs_all += abs(kv)
            nzc = sum(1 for x in n if x != 0 and x % 7 != 0)
            if nzc >= 6: abs_supp6 += abs(kv)
        print(f"\n  box |n_i|<={L}: {nontrivial} nonzero relations")
        print(f"    sum|K| all                = {abs_all:.5f}")
        print(f"    sum|K| support>=6 (floor) = {abs_supp6:.5f}")
        print(f"    signed sum (=corr_trunc)  = {signed.real:+.6f}  "
              f"(true corr={float(d8):.6f})")
        print(f"    floor removes {100*(1-abs_supp6/abs_all):.1f}% of abs mass")
        print(f"    abs(supp>=6) > delta_8?   = {abs_supp6 > float(d8)}  "
              f"(claim: YES => abs bound after floor STILL lossy)")
    print("\n  INTERPRETATION:")
    print("  - corr (signed) converges to delta_8=+0.3027 as L grows (the identity).")
    print("  - sum|K|(support>=6) is the absolute bound AFTER Lemma A's floor.")
    print("  - If sum|K|(supp>=6) STILL exceeds delta_8, the floor alone does NOT")
    print("    close the absolute bound -> a SECOND signed cancellation among the")
    print("    support-6 terms is required (the claimed remaining gap (b)).")


if __name__ == "__main__":
    main()
