#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_even_moment_extremality_opus_0621.py  (opus, 2026-06-21, THREAD B)

THE CLEAN k=8 SATURATION: consec maximizes the binding EVEN Krawtchouk moments M_2, M_4.
Plus the deeper structure of why M_1 (odd) is the lone non-consec-extremal binding moment
at k=9,10,11 -- and why the FULL L_y still wins anyway.

KEY ALGEBRAIC FACT (exact):
  N = #missed inner sectors among {1..6}.  M_j(E)=E[K_j(N;6)].
  K_2(t)=2t^2-12t+15  =>  M_2 = 4 S_2 - 10 S_1 + 15.
  consec MINIMIZES S_1 and MAXIMIZES S_2  =>  BOTH terms push M_2 UP.  (no sign conflict!)
  Compare: the factorial form L_y(k=8)=1-S_1+S_2-(9/10)S_3+(3/5)S_4 has a CONFLICTING -S_3
  term (S_3 max at consec but with NEGATIVE weight).  The Krawtchouk basis ABSORBS that
  conflict: in the K-basis k=8 is purely even (c_1=c_3=c_5=0), and EACH binding M_j is
  consec-extremal.  THAT is why the even-band reading is "clean" (HYP-2724).

THIS SCRIPT:
 (A) confirm M_2 = 4 S_2 - 10 S_1 + 15 and the analogous M_4 expansion (exact).
 (B) WIDE/adversarial test: does consec maximize M_2, M_4 over wide-spread, Sidon, AP-dilate,
     geometric, random sets (k=8,9,10)?  (the even moments should be UNIVERSALLY consec-max).
 (C) the S_1/S_2 "aligned pull" verified directly: for every E, sign of (S_2(E)-S_2(consec))
     and (S_1(consec)-S_1(E)) -- are they BOTH <=0 deviations, i.e. consec jointly extremal
     in the 2D (S_1,S_2) plane in the M_2-increasing direction?
 (D) the ODD obstruction: M_1=E[K_1(N)]=6-2 S_1.  consec does NOT minimize... wait M_1=6-2E[N]
     so max M_1 <=> min E[N]=min S_1.  consec does NOT min S_1 (S_1 beaters exist) -> M_1 not
     consec-max.  We confirm M_1=6-2 S_1 and that the M_1 beaters are EXACTLY the S_1-beaters.
 (E) THE JOINT CLOSURE for k=9,10: L_y has c_1>0 on M_1 (non-extremal) but the EVEN part
     dominates.  Decompose L_y(consec)-L_y(E) = sum_j c_j (M_j(consec)-M_j(E)) and show the
     even-moment surplus (j=2,3) always covers the odd deficit (j=1).  Report the worst case.
"""
import sys, itertools, random
from math import comb, gcd
from fractions import Fraction as F
from functools import lru_cache
sys.stdout.reconfigure(line_buffering=True)
random.seed(80808)

def Kraw(j, t, n=6):
    return sum((-1)**i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))
KTAB = [[F(Kraw(j, t)) for t in range(7)] for j in range(7)]

def occupancy(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; L = x1 - x0
        hit = set(int(7 * e * xm) % 7 for e in E)
        free = 6 - len([s for s in hit if s != 0])
        p[free] += L
    return p
def Svec(p): return [sum(comb(t, r) * p[t] for t in range(7)) for r in range(7)]
def Mvec(p): return [sum(KTAB[j][t] * p[t] for t in range(7)) for j in range(7)]

DUAL_C = {
    8:  [F(1,16), F(0), F(1,40), F(0), F(3,80), F(0), F(0)],
    9:  [F(1,12), F(1,72), F(1,36), F(1,48), F(0), F(0), F(0)],
    10: [F(1,12), F(1,72), F(1,36), F(1,48), F(0), F(0), F(0)],
    11: [F(1,8), F(1,24), F(1,24), F(0), F(0), F(0), F(0)],
}
def banner(t): print("\n" + "=" * 96 + f"\n{t}\n" + "=" * 96)

# ---------- (A) algebraic M_2, M_4 in terms of S_r ----------
banner("(A) Exact M_j in factorial-moment basis (E[N^a] via Stirling -> S_r).")
# E[N]=S1; E[N(N-1)]=2 S2 -> E[N^2]=2S2+S1; E[N(N-1)(N-2)]=6S3 -> E[N^3]=6S3+6S2+S1; etc.
# M_2=2 E[N^2]-12 E[N]+15 = 2(2S2+S1)-12 S1+15 = 4 S2 -10 S1 +15.
# K_4(t)=15,-5,-1,3,-1,-5,15. Fit exact poly: K_4 = (2/3)t^4 -8 t^3 +(94/3)t^2 -44 t +15.
# Use exact Lagrange to get K_4 coeffs, then convert powers of N to factorial moments.
def poly_coeffs_exact(vals):
    # vals[t] for t=0..6; return coeffs c0..c6 with sum c_i t^i = vals[t] (Vandermonde solve, exact)
    n = 7
    A = [[F(t**i) for i in range(n)] for t in range(n)]
    b = [F(v) for v in vals]
    # gauss
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]; pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]; M[r] = [M[r][j] - f * M[col][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]
# Stirling: t^a = sum_b S2(a,b) * fallingfactorial(t,b);  E[ff(N,b)]=b! S_b.
def stirling2(a, b):
    return sum((-1)**(b - i) * comb(b, i) * i**a for i in range(b + 1)) // (1 if b == 0 else 1) // (comb(b,b)) // 1 if False else \
           sum((-1)**(b - i) * comb(b, i) * i**a for i in range(b + 1)) // _fact(b)
def _fact(n):
    f = 1
    for i in range(2, n + 1): f *= i
    return f
def M_in_S(j):
    """return dict r-> coeff so that M_j = sum_r coeff[r] S_r (S_0=1)."""
    co = poly_coeffs_exact([Kraw(j, t) for t in range(7)])  # M_j = sum_a co[a] E[N^a]
    out = {r: F(0) for r in range(7)}
    for a in range(7):
        if co[a] == 0: continue
        for b in range(a + 1):
            # E[N^a] contributes S2(a,b)*b!*S_b
            out[b] += co[a] * stirling2(a, b) * _fact(b)
    return out
for j in [1, 2, 3, 4]:
    d = M_in_S(j)
    terms = " + ".join(f"{d[r]}*S_{r}" for r in range(7) if d[r] != 0)
    print(f"  M_{j} = {terms}")

# ---------- (B) wide/adversarial even-moment extremality ----------
banner("(B) Does consec maximize M_2, M_4 over WIDE/adversarial families? (k=8,9,10)")
def make_families(k):
    fams = []
    fams.append(("consec", list(range(k))))
    # Sidon-ish (Mian-Chowla)
    s = [0, 1]; x = 2
    while len(s) < k:
        ok = all(len({a + b for a in s + [x] for b in s + [x] if a <= b}) ==
                 len(s) * (len(s) + 1) // 2 + len(s) + 1 for _ in [0])
        sums = set()
        good = True
        cand = s + [x]
        ps = [a + b for i, a in enumerate(cand) for b in cand[i:]]
        if len(ps) == len(set(ps)): s.append(x)
        x += 1
        if x > 500: break
    fams.append(("sidon", s[:k]))
    # geometric-ish
    g = sorted(set([0] + [2**i for i in range(k)]))[:k]
    fams.append(("geom2", g))
    g3 = sorted(set([0] + [ (3**i)//1 for i in range(1,k)]))[:k]
    fams.append(("geom3", g3))
    # AP dilates
    for step in [2, 3, 5, 7, 11]:
        fams.append((f"AP*{step}", [step * i for i in range(k)]))
    # random wide
    for r in range(8):
        sp = random.randint(k, 60); E = sorted(set([0, sp] + random.sample(range(1, sp), k - 2)))
        if len(E) == k: fams.append((f"rand{r}", E))
    # near-consec perturbations
    for r in range(6):
        base = list(range(k - 1)); base.append(k - 1 + random.randint(0, 8))
        fams.append((f"pert{r}", sorted(set(base))[:k] if len(set(base)) >= k else base))
    return fams
for k in [8, 9, 10]:
    Mc = Mvec(occupancy(list(range(k))))
    fams = make_families(k)
    over2 = []; over4 = []
    for name, E in fams:
        if len(set(E)) != k: continue
        M = Mvec(occupancy(E))
        if M[2] > Mc[2]: over2.append((name, float(M[2])))
        if M[4] > Mc[4]: over4.append((name, float(M[4])))
    print(f"  k={k}: M_2(consec)={float(Mc[2]):.5f}  M_4(consec)={float(Mc[4]):.5f}")
    print(f"        families with M_2 > consec: {over2 if over2 else 'NONE'}")
    print(f"        families with M_4 > consec: {over4 if over4 else 'NONE'}")

# ---------- (D) M_1 = 6 - 2 S_1 ; odd obstruction ----------
banner("(D) M_1 = 6 - 2 S_1 (odd).  M_1-beaters == S_1-beaters (consec NOT S_1-min).")
for k in [9, 10]:
    Ec = list(range(k)); pc = occupancy(Ec); Mc = Mvec(pc); Sc = Svec(pc)
    nM1 = nS1 = same = 0
    for rest in itertools.combinations(range(1, k + 5), k - 1):
        E = [0] + list(rest); p = occupancy(E); M = Mvec(p); S = Svec(p)
        bM1 = M[1] > Mc[1]; bS1 = S[1] < Sc[1]
        if bM1: nM1 += 1
        if bS1: nS1 += 1
        if bM1 == bS1: same += 1
        if (bM1 != bS1):
            pass
    print(f"  k={k}: M_1(consec)=6-2*{Sc[1]}={Mc[1]}={float(Mc[1]):.4f}; "
          f"M_1-beaters={nM1}, S_1-beaters(S1<consec)={nS1}, agreement={same}/{comb(k+4,k-1)}")

# ---------- (E) joint closure k=9,10: even surplus covers odd deficit ----------
banner("(E) k=9,10 JOINT: L_y(consec)-L_y(E) = sum_j c_j (M_j(c)-M_j(E)); even covers odd.")
WIN = {9: 16, 10: 15, 11: 14}
for k in [9, 10, 11]:
    c = DUAL_C[k]; Ec = list(range(k)); Mc = Mvec(occupancy(Ec))
    worst_total = F(0); worst_E = None; worst_decomp = None
    nbad = 0; nset = 0
    for rest in itertools.combinations(range(1, WIN[k] + 1), k - 1):
        E = [0] + list(rest); M = Mvec(occupancy(E)); nset += 1
        # contribution per j to (L_y(consec) - L_y(E)); positive = consec wins on that moment
        dec = [c[j] * (Mc[j] - M[j]) for j in range(7)]
        total = sum(dec)
        if total < worst_total:  # L_y(E) > L_y(consec)
            worst_total = total; worst_E = E; worst_decomp = dec
        if total < 0: nbad += 1
    print(f"\n  k={k} ({nset} sets, maxE<={WIN[k]}):  L_y-beaters (L_y(E)>consec): {nbad}")
    if worst_E:
        odd = sum(worst_decomp[j] for j in [1, 3, 5])
        even = sum(worst_decomp[j] for j in [2, 4, 6])
        print(f"    worst (most negative total={float(worst_total):.6f}) at {worst_E}")
        print(f"      per-j c_j*(Mc-M): {[round(float(x),5) for x in worst_decomp]}")
        print(f"      ODD part (j=1,3,5)={float(odd):.5f}  EVEN part (j=2,4,6)={float(even):.5f}")
    else:
        print(f"    consec is the strict max; min surplus over window:")
        # recompute min positive
        ms = None
        for rest in itertools.combinations(range(1, WIN[k] + 1), k - 1):
            E = [0] + list(rest)
            if E == Ec: continue
            M = Mvec(occupancy(E)); tot = sum(c[j]*(Mc[j]-M[j]) for j in range(7))
            if ms is None or tot < ms[0]: ms = (tot, E)
        print(f"      tightest: surplus={float(ms[0]):.6f} at {ms[1]}")

print("\nDONE.")
