#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) SIGNED-QUOTIENT SHELL DECAY  (kind-pasteur-2026-06-19-S12)

ANGLE: make the "signed/coset quotient is the ruler" (HYP-2640) PRECISE and
QUANTITATIVE for the support-6 Fourier sum
    meas(S7(E)) = M7(k) + sum_{0!=n in Lambda(E)} K(n).
We enumerate the actual relation lattice Lambda(E) up to an infinity-norm bound L,
restrict to support-6 relations (THM-538 floor), and watch:

  (a) SIGNED partial sum  C_signed(L) = sum_{|n|_inf <= L} K(n)   -> should -> p0-M7
  (b) ABSOLUTE partial sum A(L) = sum |K(n)|                       (diverges, MISTAKE-078)
  (c) per-shell signed sum  s(L) = sum_{|n|_inf = L} K(n)          -- the decay we want
  (d) per-mod-7-coset signed sums within each shell                -- the cancellation engine

We then estimate the decay rate of |s(L)| and compare AP vs wide set.

K(n) factorization (THM-538 / HYP-2632):  with support exactly {1..6} coords on values v_j,
   K(n) = sum_{T subset {1..6}} (-1)^|T| prod_j chat_T(n_j),
   chat_T(m) = Fourier coeff of indicator of [0,1) minus union_{j in T}[j/7,(j+1)/7) at frequency m.
For m != 0 (and 7 not| m): chat_T(m) = - sum_{j in T} cbar(j,m) where
   cbar(j,m) = integral_{j/7}^{(j+1)/7} e^{-2 pi i m t} dt
             = e^{-2 pi i m (j+ .. )}... we just use the exact closed form below.
We use the EXACT character form (HYP-2632): for nonzero non-7 freq m,
   shat_T(m) := chat_T(m) = (1/(2 pi i m)) * (1 - e^{-2 pi i m/7}) * sum_{j in T} e^{-2 pi i m j /7}   [sign bookkeeping below]
We instead just compute chat_T(m) numerically from the definition to avoid sign errors,
then verify K reconstructs p0-M7 against the exact dist_p engine.
"""
import sys, itertools, cmath, math
from fractions import Fraction
from functools import reduce
from math import gcd

sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---------- exact measure engine (from lrc14_empty_sector_distribution_kps.py) ----------
def dist_p0(E):
    """exact p0 = meas(S7) as a Fraction."""
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p0 = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        N = sum(1 for j in range(1,7) if j not in hit)
        if N == 0:
            p0 += (hi-lo)
    return p0

# ---------- Fourier kernel ----------
def chat_T(T, m):
    """Fourier coeff (freq m) of f_T = indicator of [0,1) minus union_{j in T}[j/7,(j+1)/7).
       T subset of {1..6}.  Returns complex.  m integer."""
    if m == 0:
        return 1.0 - len(T)/7.0
    if m % 7 == 0:
        return 0.0  # THM-503 apex prime vanishing for the SECTOR pieces; full f_T also -> see note
    # f_T = 1 - sum_{j in T} 1_{[j/7,(j+1)/7)}
    # Fourier coeff of 1 over [0,1) at m!=0 is 0.
    # coeff of 1_{[a,b)} at freq m: integral_a^b e^{-2pi i m t} dt = (e^{-2pi i m a}-e^{-2pi i m b})/(2 pi i m)
    s = 0.0+0.0j
    for j in T:
        a = j/7.0; b = (j+1)/7.0
        c = (cmath.exp(-2j*math.pi*m*a) - cmath.exp(-2j*math.pi*m*b))/(2j*math.pi*m)
        s += c
    return -s

def K_of_n(vals):
    """K(n) for a support-6 relation with nonzero coordinate values 'vals' (the n_j on the 6 supported offsets).
       vals: list of 6 integers (none mult of 7, none zero).  Sum over T subset {1..6}."""
    total = 0.0+0.0j
    for r in range(7):
        for T in itertools.combinations(range(1,7), r):
            prod = 1.0+0.0j
            for v in vals:
                prod *= chat_T(T, v)
            total += ((-1)**r) * prod
    return total

# precompute chat table once we know which freqs appear -- but K_of_n recomputes; cache:
_chat_cache = {}
def chat_T_cached(T, m):
    key = (T, m)
    c = _chat_cache.get(key)
    if c is None:
        c = chat_T(T, m); _chat_cache[key] = c
    return c

def K_of_n_fast(vals):
    total = 0.0+0.0j
    for r in range(7):
        for T in itertools.combinations(range(1,7), r):
            prod = 1.0+0.0j
            for v in vals:
                prod *= chat_T_cached(T, v)
            total += ((-1)**r) * prod
    return total

# ---------- relation lattice enumeration ----------
def relations_support6(E, L):
    """All n in Z^k with sum_j n_j E[j] = 0, support EXACTLY 6 (6 nonzero coords),
       no nonzero coord a multiple of 7, |n|_inf <= L.
       E is the offset list (0 included).  Returns list of (vals6, idxs6, n_inf, residues mod 7).
       Optimization: enumerate 5 free coords, solve the 6th from sum n_j e_j = 0."""
    k = len(E)
    out = []
    for idxs in itertools.combinations(range(k), 6):
        es = [E[i] for i in idxs]
        # pick the coordinate with the LARGEST |e| as the dependent one (best chance of small n)
        dep = max(range(6), key=lambda i: abs(es[i]))
        e_dep = es[dep]
        if e_dep == 0:
            # offset 0 is in support: its coefficient is free; constraint doesn't involve it.
            # then constraint is on the other 5; the 0-offset coord ranges freely. Handle separately.
            free = [i for i in range(6) if i != dep]
            efree = [es[i] for i in free]
            ranges = [range(-L, L+1) for _ in range(5)]
            for vfree in itertools.product(*ranges):
                if any(c == 0 or c % 7 == 0 for c in vfree):
                    continue
                if sum(c*e for c, e in zip(vfree, efree)) != 0:
                    continue
                # dependent (the 0-offset) coord can be any nonzero non-7 in [-L,L]
                for vd in range(-L, L+1):
                    if vd == 0 or vd % 7 == 0:
                        continue
                    combo = [0]*6
                    for i, c in zip(free, vfree): combo[i] = c
                    combo[dep] = vd
                    combo = tuple(combo)
                    ninf = max(abs(c) for c in combo)
                    out.append((combo, idxs, ninf, tuple(c % 7 for c in combo)))
            continue
        free = [i for i in range(6) if i != dep]
        efree = [es[i] for i in free]
        ranges = [range(-L, L+1) for _ in range(5)]
        for vfree in itertools.product(*ranges):
            if any(c == 0 or c % 7 == 0 for c in vfree):
                continue
            s = sum(c*e for c, e in zip(vfree, efree))
            if s % e_dep != 0:
                continue
            vd = -s // e_dep
            if vd == 0 or vd % 7 == 0 or abs(vd) > L:
                continue
            combo = [0]*6
            for i, c in zip(free, vfree): combo[i] = c
            combo[dep] = vd
            combo = tuple(combo)
            ninf = max(abs(c) for c in combo)
            out.append((combo, idxs, ninf, tuple(c % 7 for c in combo)))
    return out

def analyze(E, L, label):
    print(f"\n=== {label}  E={E}  L={L} ===")
    p0 = dist_p0(E)
    # M7(k):
    k = len(E)
    M7 = sum(((-1)**t)*math.comb(6,t)*((1-t/7.0)**(k-1)) for t in range(7))
    target = float(p0) - M7
    print(f"  p0(exact)={float(p0):.6f}  M7={M7:.6f}  correction target p0-M7 = {target:.6f}")

    rels = relations_support6(E, L)
    # dedup: n and -n give conjugate K, both counted (lattice symmetric). Keep all.
    # group by shell (ninf) and by mod-7 coset (residue tuple up to global F7^* scaling? use raw residue)
    from collections import defaultdict
    shell_signed = defaultdict(complex)
    shell_abs = defaultdict(float)
    coset_signed = defaultdict(complex)   # key (shell, residue)
    total_signed = 0.0+0.0j
    total_abs = 0.0
    for combo, idxs, ninf, res in rels:
        Kv = K_of_n_fast(combo)
        shell_signed[ninf] += Kv
        shell_abs[ninf] += abs(Kv)
        coset_signed[(ninf,res)] += Kv
        total_signed += Kv
        total_abs += abs(Kv)
    print(f"  #support-6 relations |n|inf<=L : {len(rels)}")
    print(f"  SIGNED partial sum  Re = {total_signed.real:.6f}  (Im={total_signed.imag:.2e})")
    print(f"  ABSOLUTE partial sum    = {total_abs:.6f}   (abs/|signed| ratio = {total_abs/max(abs(total_signed.real),1e-12):.1f})")
    print(f"  recovered fraction of target: {total_signed.real/target if target!=0 else float('nan'):.4f}")
    print(f"  per-shell SIGNED s(L) and ABSOLUTE a(L):")
    cumS = 0.0
    for L0 in sorted(shell_signed):
        sgn = shell_signed[L0].real
        ab = shell_abs[L0]
        cumS += sgn
        print(f"    shell |n|inf={L0:2d}: signed={sgn:+.6f}  abs={ab:.6f}  abs/|sgn|={ab/max(abs(sgn),1e-12):6.1f}  cumSigned={cumS:+.6f}")
    return target, shell_signed, shell_abs, coset_signed

if __name__ == "__main__":
    print("LRC(14) SIGNED-QUOTIENT SHELL DECAY  (kps-S12)")
    print("Comparing AP (consec) vs wide odd-AP set at k=8, support-6 relation lattice.")
    # AP consec k=8
    AP8 = list(range(8))
    # wide set: the (5,0,5) packet from HYP-2640 default bank
    WIDE8 = [0,1,3,5,7,9,11,13]
    for L in [3, 5, 7]:
        analyze(AP8, L, f"AP8 (consec)")
    for L in [3, 5, 7]:
        analyze(WIDE8, L, f"WIDE8 odd-AP")
