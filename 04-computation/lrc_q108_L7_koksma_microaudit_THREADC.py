#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_koksma_microaudit_THREADC.py   (THREAD C audit, 2026-06-21)

ADVERSARIAL audit of the INTERMEDIATE steps of the elementary D<=14/p proof,
not just its conclusion. The proof claims:

  (A) Fix i. {qv} in I_i picks q sub-intervals of v of length 1/(7q). On the m-th
      one, {pv} sweeps an arc [a_m, a_m+L), L = p/(7q), with
        a_m = {pi/(7q) + pm/q}.
      Because gcd(p,q)=1, {pm/q : m=0..q-1} = {0,1/q,...,(q-1)/q}, so the a_m are
      EXACTLY equally spaced (gap 1/q, with a global shift).
  (B) mu(i,j) = (q/p)*(1/q) sum_m h_j(a_m), where h_j(a)=overlap([a,a+L), I_j).
  (C) h_j is a trapezoid in a with total variation Var(h_j)=2*min(L,1/7)=2/7 (since L>1/7).
  (D) Koksma on the q equally-spaced points: |(1/q)sum h_j(a_m) - int_0^1 h_j| <= Var(h_j)*D*
      with D* <= 1/q (star-discrepancy of q equally spaced points, with shift).
  (E) int_0^1 h_j da = L/7 = p/(49q).   [the average overlap of a length-L arc with a 1/7 bin]
  Combine: |mu(i,j)-1/49| = (q/p)|err| <= (q/p)(2/7)(1/q) = 2/(7p); sum 49 cells => 14/p.

We verify (A)-(E) EXACTLY (rational) for many (p,q) and report any deviation.
We also EXACTLY recompute the LHS |mu(i,j)-1/49| and compare against the per-cell
2/(7p) claim and the proof's reconstruction (q/p)(1/q)|sum h_j(a_m) - q*int h_j|.
"""
from fractions import Fraction as Fr
from math import gcd

P = 7

def frac(x):
    return x - (x.numerator // x.denominator)

def sec(yf):
    # sector index of a fractional value in [0,1)
    return int(P * yf)  # floor(7*yf)

# ---- exact mu via breakpoint enumeration (ground truth, same as the repo) ----
def mu_cells_exact(p, q):
    bp = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(0, P * f):
            bp.add(Fr(t, P * f))
    vs = sorted(bp)
    mu = {}
    for a, b in zip(vs, vs[1:]):
        mid = (a + b) / 2
        key = (sec(frac(q * mid)), sec(frac(p * mid)))
        mu[key] = mu.get(key, Fr(0)) + (b - a)
    return mu

# ---- the proof's overlap trapezoid h_j(a) = |[a,a+L) cap I_j| , I_j=[j/7,(j+1)/7) ----
def overlap_arc_bin(a, L, j):
    """exact measure of overlap of the circular arc [a, a+L) (length L<=1) with bin I_j=[j/7,(j+1)/7).
       a in [0,1). Handles wraparound."""
    # represent arc as up to two linear intervals on [0,1)
    a = frac(a)
    end = a + L
    segs = []
    if end <= 1:
        segs.append((a, end))
    else:
        segs.append((a, Fr(1)))
        segs.append((Fr(0), end - 1))
    lo = Fr(j, P); hi = Fr(j + 1, P)
    tot = Fr(0)
    for (s, e) in segs:
        tot += max(Fr(0), min(e, hi) - max(s, lo))
    return tot

# ---- proof reconstruction of mu(i,j): mu = (1/p)*sum_m h_j(a_m), with the (q/p)(1/q) bookkeeping
# Actually mu(i,j) = sum over the q v-subintervals of [length-of-{pv} in I_j within that subinterval].
# On the m-th subinterval v=(i+7m+s)/(7q), s in [0,1): pv = p(i+7m+s)/(7q) = a_m' + s*p/(7q),
# so {pv} sweeps arc starting at a_m={p(i+7m)/(7q)} of length L=p/(7q), parametrized by s in [0,1).
# The v-length contributed is (dv/ds)=1/(7q) per unit s, and the fraction of s with {pv} in I_j
# is (1/L)*overlap([a_m,a_m+L),I_j)... wait: as s goes 0->1, {pv} moves by L. So the s-measure with
# {pv} in I_j is overlap([a_m,a_m+L),I_j)/L * 1 ??? No: s->{pv} is linear with slope L over s in[0,1],
# i.e. {pv}=a_m+ s*L (mod 1). s-measure with {pv} in I_j = (1/L)*overlap. Times dv=ds/(7q).
# So contribution of subinterval m to mu(i,j) = (1/(7q)) * (1/L) * overlap = (1/(7q))*(7q/p)*overlap
#   = overlap/p.  Summing m: mu(i,j) = (1/p) sum_m h_j(a_m).
# The repo writeup says mu=(q/p)(1/q)sum h_j(a_m) = (1/p) sum h_j(a_m). SAME. Good — check it.
def mu_proof_reconstruction(p, q):
    L = Fr(p, P * q)
    mu = {}
    for i in range(P):
        for m in range(q):
            a_m = frac(Fr(p * (i + P * m), P * q))
            for j in range(P):
                ov = overlap_arc_bin(a_m, L, j)
                mu[(i, j)] = mu.get((i, j), Fr(0)) + ov / p
    return mu

# ---- step (A): a_m equally spaced with gap 1/q ----
def check_equal_spacing(p, q):
    for i in range(P):
        starts = sorted(frac(Fr(p * (i + P * m), P * q)) for m in range(q))
        # equally spaced gap 1/q?
        # sorted starts should be {c, c+1/q, ...} for some c; check consecutive diffs all equal mod
        diffs = set()
        for k in range(q):
            d = frac(starts[(k + 1) % q] - starts[k]) if k < q - 1 else frac(starts[0] + 1 - starts[q - 1])
            diffs.add(d)
        if diffs != {Fr(1, q)}:
            return False, (i, starts, diffs)
    return True, None

# ---- step (C): Var(h_j) over a in [0,1) is exactly 2/7 (when L>=1/7) ----
def check_var_hj(p, q):
    L = Fr(p, P * q)
    # sample h_j on a fine exact grid of breakpoints; trapezoid breakpoints occur where
    # a, a+L cross j/7,(j+1)/7 -> a in {j/7, (j+1)/7, j/7-L, (j+1)/7-L} mod 1
    results = {}
    for j in range(P):
        crit = set()
        for c in (Fr(j, P), Fr(j + 1, P)):
            crit.add(frac(c)); crit.add(frac(c - L))
        crit |= {Fr(0)}
        pts = sorted(crit)
        # evaluate h_j at midpoints of consecutive crit intervals AND at the crit points to get the TV
        # Total variation of a continuous piecewise-linear periodic function = sum |jumps in value| at vertices
        allpts = sorted(crit | {Fr(0), Fr(1)})
        vals = [(a, overlap_arc_bin(a, L, j)) for a in allpts]
        tv = Fr(0)
        for k in range(len(vals) - 1):
            tv += abs(vals[k + 1][1] - vals[k][1])
        results[j] = tv
    return results

# ---- step (E): int_0^1 h_j(a) da = L/7 ----
def check_integral_hj(p, q):
    L = Fr(p, P * q)
    # integral over a in [0,1) of overlap([a,a+L),I_j) = L * |I_j| = L/7 (Fubini)
    # compute exactly by integrating piecewise-linear
    out = {}
    for j in range(P):
        crit = {Fr(0), Fr(1)}
        for c in (Fr(j, P), Fr(j + 1, P)):
            crit.add(frac(c)); crit.add(frac(c - L))
        pts = sorted(crit)
        integral = Fr(0)
        for k in range(len(pts) - 1):
            a0, a1 = pts[k], pts[k + 1]
            mid = (a0 + a1) / 2
            # piecewise linear -> trapezoidal exact: average of endpoints * width
            h0 = overlap_arc_bin(a0, L, j); h1 = overlap_arc_bin(a1, L, j)
            integral += (h0 + h1) / 2 * (a1 - a0)
        out[j] = integral
    return out

def main():
    print("=" * 78)
    print("ADVERSARIAL micro-audit of the elementary D<=14/p Koksma proof steps")
    print("=" * 78)
    fail_eqsp = fail_var = fail_int = fail_recon = fail_percell = 0
    worst_percell = (0, 0, 0, 0, 0.0)  # p,q,i,j, |mu-1/49|*7p/2 ratio (should be <=1)
    checked = 0
    for q in range(1, 26):
        pmax = int(Fr(43, 20) * q)
        for p in range(q + 1, pmax + 1):
            if gcd(p, q) != 1:
                continue
            r = Fr(p, q)
            if not (Fr(1) < r <= Fr(43, 20)):
                continue
            checked += 1
            # (A) equal spacing
            ok, info = check_equal_spacing(p, q)
            if not ok:
                fail_eqsp += 1
                if fail_eqsp <= 3:
                    print(f"  EQSP FAIL p/q={p}/{q}: {info}")
            # (C) Var = 2/7
            vh = check_var_hj(p, q)
            for j, tv in vh.items():
                if tv != Fr(2, 7):
                    fail_var += 1
                    if fail_var <= 5:
                        print(f"  VAR!=2/7 p/q={p}/{q} j={j}: Var(h_j)={tv} ({float(tv):.5f})")
            # (E) integral = L/7
            ih = check_integral_hj(p, q)
            L = Fr(p, P * q)
            for j, integ in ih.items():
                if integ != L / 7:
                    fail_int += 1
                    if fail_int <= 5:
                        print(f"  INT!=L/7 p/q={p}/{q} j={j}: int={integ} vs L/7={L/7}")
            # (B) proof reconstruction of mu == exact mu
            mu_exact = mu_cells_exact(p, q)
            mu_recon = mu_proof_reconstruction(p, q)
            keys = set(mu_exact) | set(mu_recon)
            for k in keys:
                if mu_exact.get(k, Fr(0)) != mu_recon.get(k, Fr(0)):
                    fail_recon += 1
                    if fail_recon <= 5:
                        print(f"  RECON MISMATCH p/q={p}/{q} cell={k}: exact={mu_exact.get(k,Fr(0))} recon={mu_recon.get(k,Fr(0))}")
                    break
            # per-cell |mu-1/49| <= 2/(7p)
            inv = Fr(1, 49)
            for i in range(P):
                for j in range(P):
                    dev = abs(mu_exact.get((i, j), Fr(0)) - inv)
                    bound = Fr(2, 7 * p)
                    if dev > bound:
                        fail_percell += 1
                        if fail_percell <= 8:
                            print(f"  PERCELL VIOLATION p/q={p}/{q} cell=({i},{j}): |mu-1/49|={dev}({float(dev):.5f}) > 2/(7p)={bound}({float(bound):.5f})")
                    ratio = float(dev / bound) if bound != 0 else 0.0
                    if ratio > worst_percell[4]:
                        worst_percell = (p, q, i, j, ratio)
    print(f"\n  checked {checked} ratios p/q in (1,43/20], q<=25")
    print(f"  (A) equal-spacing failures        : {fail_eqsp}")
    print(f"  (C) Var(h_j)!=2/7 failures         : {fail_var}")
    print(f"  (E) integral!=L/7 failures         : {fail_int}")
    print(f"  (B) mu reconstruction mismatches  : {fail_recon}")
    print(f"  per-cell |mu-1/49|<=2/(7p) violations: {fail_percell}")
    print(f"  worst per-cell tightness ratio = {worst_percell[4]:.4f} at p/q={worst_percell[0]}/{worst_percell[1]} cell=({worst_percell[2]},{worst_percell[3]})")
    print("  (ratio<=1 means the per-cell bound 2/(7p) holds; ratio close to 1 means tight)")

if __name__ == "__main__":
    main()
