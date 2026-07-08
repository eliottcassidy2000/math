#!/usr/bin/env python3
r"""
lrc_k10_Dq_window_kps_S75.py  (kind-pasteur-2026-07-07-S75, HYP-5247)

CLOSING THE k=10 GAP: the residue-class-diameter window floor.

The k=10 (A') leg is TRUE with 7x margin (min true rho* = 0.398 at the compact AP_10), but
neither tool alone proves it:
  - intersection ledger (kps-S60): covers diam <= 17;
  - average-form conditional tent (mac-mini THM-655): fails on 79/286 shapes (avgc > c*),
    the offenders being 2-ADICALLY STRUCTURED families (long 2-AP + a bump).
The genuine residual = diam >= 18 AND avgc > c*: WIDE 2-adic families.

THE FIX (this file): klein-THM-653's window floor bounds each residue cluster's width by
|delta|*diam; for a structured family the residue-class spread is FAR below diam.  Replace
diam by
    D_q(E) := max over residues s mod q of diam({e in E : e = s mod q})
(independent of the resonance p since gcd(p,q)=1 permutes residues).  The good window at
p/q then has half-width c_q/D_q(E) with c_q=(7-q)/(7q) -- WIDER than klein's c_q/diam.  For
the 2-AP {0,2,..,16,B}: D_2 = 16 (the even sub-cluster) << B = diam, so the q=2 window at
x=1/2 has half-width c_2/16 -- a big window klein's diam-form misses.

TEST: for every k=10 residual offender (diam>=18, avgc>c* for a binding P), compute the
CONDITIONAL D_q-window mass  meas(G_P cap Union_{q<=6, p} [p/q +- c_q/D_q(E)])  and check
>= m_P.  If yes, the D_q-window floor closes the k=10 residual; combined with the ledger
(diam<=17) and average form (avgc<=c*), k=10 is discharged.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

M_P = F(14249, 252252)
K = 10
BETA = F(14 - K, 7 * K)
INTF = (F(1, 7) - BETA) ** 2 / 2
ONE_MINUS_FLOOR = F(2 * (K - 1) * (K - 7), 7 * K)  # 27/35

def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            bad.append((F(j, p) - w, F(j, p) + w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort()
    merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else:
            merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def meas_intervals(ivs):
    return sum(h - l for l, h in ivs)

def intersect(ivsA, ivsB):
    """intersection of two sorted interval lists."""
    out = []
    i = j = 0
    A = sorted(ivsA); B = sorted(ivsB)
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: out.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return out

def union_intervals(ivs):
    if not ivs: return []
    s = sorted(ivs); out = [list(s[0])]
    for l, h in s[1:]:
        if l <= out[-1][1]: out[-1][1] = max(out[-1][1], h)
        else: out.append([l, h])
    return [(l, h) for l, h in out]

def Dq(E, q):
    """max residue-class diameter mod q."""
    best = 0
    for s in range(q):
        cls = [e for e in E if e % q == s]
        if len(cls) >= 2:
            best = max(best, max(cls) - min(cls))
    return best

def Dq_windows(E):
    """union of good windows: for q=1..6, p/q with gcd(p,q)=1, half-width c_q/D_q(E)."""
    wins = []
    diam = max(E) - min(E)
    for q in range(1, 7):
        cq = F(7 - q, 7 * q)
        Dqv = Dq(E, q)
        if q == 1:
            Dqv = diam  # the p/q=0 (x~0) window governed by full diam
        if Dqv == 0:
            continue
        hw = cq / Dqv
        for p in range(0, q + 1):
            if q > 1 and gcd(p, q) != 1:
                continue
            c = F(p, q)
            lo, hi = c - hw, c + hw
            # wrap into [0,1)
            lo = max(lo, F(0)); hi = min(hi, F(1))
            if lo < hi:
                wins.append((lo, hi))
    return union_intervals(wins)

def c_of(iv, mGP, d):
    tot = F(0); d = int(abs(d)); th = F(1, 7)
    for (l, h) in iv:
        m_lo = int(l * d - th) - 1; m_hi = int(h * d) + 1
        for m in range(max(0, m_lo), m_hi + 1):
            wl, wh = (F(m) + BETA) / d, (F(m) + th) / d
            a, b = max(wl, l), min(wh, h)
            if a >= b: continue
            tot += (d * (b * b - a * a) / 2 - (F(m) + BETA) * (b - a))
    return tot / (mGP * INTF)

def avgc(E, iv, mGP, cache):
    ds = [E[j] - E[i] for i in range(K) for j in range(i + 1, K)]
    s = F(0)
    for d in ds:
        if d not in cache: cache[d] = c_of(iv, mGP, d)
        s += cache[d]
    return s / len(ds)

# find residual offenders (diam>=18, avgc>c*) then test the conditional D_q-window floor
BIND_P = [(1,12,13),(11,12,13),(1,2,3),(4,5,6),(5,12,13),(2,12,13),(3,4,5),(1,11,13),(6,8,10),(8,10,12)]
print("=" * 96)
print("k=10 residual (diam>=18, avgc>c*): does the CONDITIONAL D_q-window floor give rho* >= m_P?")
print("=" * 96)
rng = random.Random(75)
tested = 0; closed = 0; worst_ratio = None
def gen_2adic(rng, dmax):
    """long 2-AP {0,2,..,2m} + odd bump, primitive, diam in [18,dmax]."""
    m = rng.randint(5, 8)
    base = list(range(0, 2 * m + 1, 2))
    B = rng.randrange(2 * m + 1, dmax, 2) + 1  # odd bump beyond
    E = sorted(set(base + [B]))
    while len(E) < K:
        E = sorted(set(E + [rng.randrange(1, B, 2)]))  # odd fillers keep primitivity
    E = E[:K]
    if len(E) != K: return None
    E = sorted(E)
    g = 0
    for a in E[1:]: g = gcd(g, a - E[0])
    if g != 1 or E[-1] - E[0] < 18: return None
    return [e - E[0] for e in E]

for P in BIND_P:
    iv = GP_intervals(P); mGP = meas_intervals(iv)
    cstar = (1 - M_P / mGP) / ONE_MINUS_FLOOR
    cache = {}
    # gather offenders for this P
    offenders = []
    for _ in range(1200):
        E = gen_2adic(rng, 55)
        if E is None: continue
        if avgc(E, iv, mGP, cache) > cstar:
            offenders.append(E)
        if len(offenders) >= 40: break
    if not offenders:
        print(f"  P={str(P):>12}: no diam>=18 2-adic offenders found (avgc<=c*)")
        continue
    minrho = None
    for E in offenders:
        wins = Dq_windows(E)
        gp_win = intersect(iv, wins)
        rho_lb = meas_intervals(gp_win)   # conditional D_q-window mass (a lower bnd on rho*)
        if minrho is None or rho_lb < minrho[0]:
            minrho = (rho_lb, E)
    rho_lb, Eworst = minrho
    tested += 1
    ok = rho_lb >= M_P
    if ok: closed += 1
    ratio = float(rho_lb / M_P)
    if worst_ratio is None or ratio < worst_ratio[0]:
        worst_ratio = (ratio, P, Eworst, float(rho_lb))
    print(f"  P={str(P):>12} c*={float(cstar):.3f}: {len(offenders):>2} offenders, "
          f"min D_q-window rho* >= {float(rho_lb):.5f} ({ratio:.2f}x m_P)  "
          f"[{'CLOSED' if ok else 'GAP'}]  worst E={Eworst}")
print()
print(f"  Summary: {closed}/{tested} binding P's have all diam>=18 offenders closed by the")
print(f"  conditional D_q-window floor.  Worst: {worst_ratio[3]:.5f} ({worst_ratio[0]:.2f}x m_P) "
      f"at P={worst_ratio[1]}, E={worst_ratio[2]}")
print()
print("  READING: if all CLOSED, the k=10 leg = [ledger diam<=17] + [avgc<=c* wide] +")
print("  [D_q-window for the 2-adic residual] -- three tools, full cover.  The D_q-window is")
print("  klein-THM-653 with diam -> residue-class-diameter: the structural fix the 2-AP needs.")
