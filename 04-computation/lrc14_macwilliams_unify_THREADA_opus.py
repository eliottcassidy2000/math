#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (opus, 2026-06-21): UNIFY corr=weight-enumerator (HYP-2724/kps) with
measS7=P(N=0)=even-Krawtchouk read (HYP-2724/mac-mini) via MacWILLIAMS.

THREE EXACT OBJECTS, all over Q (exact Fractions):

 (1) DEPTH LAW   pi_E(h) = meas{x in [0,1): exactly h of the 7 residues
     {floor(7 frac(e_i x)) : e_i in E} are HIT}.   h = 7 - N, N=#empty sectors.
     measS7(E) = pi_E(7) = P(N=0) = P(all 7 sectors hit).

 (2) COVERAGE / Bonferroni weights  g_J(h) = sum_{j=0}^{J} (-1)^j C(7-h, j).
     S_J(E) = <g_J, pi_E>  (J-th IE truncation).  measS7 = S_6 = <g_6, pi_E>,
     g_6(h)=[h=7].  These g_J ARE the partial alternating sums of the binary
     KRAWTCHOUK polynomials in the variable N = 7-h:
          full row g_6(h) = sum_{j=0}^{6}(-1)^j C(N, j)  with N = 7-h
                          = (-1)^N * [N=0]  ... no: it = [N==0] exactly
     (verified below).  The even partials g_2,g_4 are convex-in-N; odd g_1,g_3,g_5
     are not.  measS7 is the TERMINAL even partial.

 (3) RELATION CODE  Lambda(E) = { n in Z^k : sum_i n_i e_i = 0 }, a rank-(k-1)
     lattice/code [k, k-1, d], d = min support among NONZERO-offset relations.
     corr(E) = measS7(E) - iid_k = sum_{0 != n in Lambda(E)} K(n)   (mac HYP-2719).
     consec = ANTI-MDS (d=2: 2*e_1 = 1*e_2 i.e. 2*1=1*2 -> relation support-2).
     Sidon/arc = MDS (d>=3).
     We build the WEIGHT ENUMERATOR W_E(x,y) = sum_{n in Lambda, ||n||<=B}
        x^{k - supp(n)} y^{supp(n)} restricted to a bounded box (the lattice is
        infinite; we truncate at a height box and read the LOW-WEIGHT shells,
        which are the ones that matter for corr).

THE UNIFICATION TO BUILD/TEST:
  Krawtchouk is the MacWilliams transform kernel.  So:
   - corr = sum over the relation code of a Krawtchouk-weighted shell count,
   - measS7 = <g_6, pi_E> = a Krawtchouk read of the depth law,
  should be MacWilliams-DUAL.  consec's LOW-WEIGHT-heavy relation enumerator
  <=> convex-dominant occupancy N <=> top even-Krawtchouk band maximal.

This script: exact pi_E, exact g_J/Krawtchouk identity check, exact corr via
depth law, relation-code low-support shell counts, and the cross-tabulation that
makes the MacWilliams link concrete.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce

# ---------- depth law (exact) ----------
def depth_law(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7*abs(e)+1):
            bps.add(F(a, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    law = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*xm
            v = v - (v.numerator // v.denominator)   # frac in [0,1)
            hit.add((v.numerator*7)//v.denominator)
        h = len(hit)
        law[h] += hi-lo
    return law

def measS7(E):
    return depth_law(E)[7]

# ---------- Krawtchouk / coverage weights ----------
def gJ(h, J):
    return sum((-1)**j * comb(7-h, j) for j in range(J+1))

# Binary Krawtchouk K_j^{(m)}(w) = sum_t (-1)^t C(w,t) C(m-w, j-t)
def krawtchouk(j, w, m):
    return sum((-1)**t * comb(w, t) * comb(m-w, j-t) for t in range(j+1))

# ---------- iid baseline ----------
def stirling2(n, k):
    return sum((-1)**(k-i) * comb(k, i) * i**n for i in range(k+1)) // 1 if False else \
           sum((-1)**(k-i) * comb(k, i) * i**n for i in range(k+1)) // 1
def stirling2_exact(n, k):
    # S(n,k) integer
    s = 0
    for i in range(k+1):
        s += (-1)**(k-i) * comb(k, i) * i**n
    return s // 1  # exact via factorial division below
def iid_k(k):
    # iid_k = 7! * S(k,7) / 7^k   (prob that k iid uniform sectors cover all 7)
    Snk = 0
    for i in range(8):
        Snk += (-1)**(7-i) * comb(7, i) * i**k
    Snk = Snk  # = 7! * S(k,7)
    return F(Snk, 7**k)

def corr(E):
    k = len(E)
    return measS7(E) - iid_k(k)

# ---------- relation code Lambda(E): low-support shells in a height box ----------
def relation_shells(E, B=3):
    """Count NONZERO-offset relations sum n_i e_i = 0 with |n_i|<=B, by support size.
       Returns dict supp -> count (excluding the all-zero vector and trivial e_0=0 support-1)."""
    E = list(E)
    k = len(E)
    # the e_0=0 coordinate: a relation supported ONLY on the 0-offset (n_i e_i with e_i=0)
    # contributes trivially.  We follow the prompt: exclude relations whose support is
    # ONLY the zero-offset coordinate(s).  d = min support among relations touching a NONZERO e_i.
    zero_idx = [i for i,e in enumerate(E) if e == 0]
    nz_idx   = [i for i,e in enumerate(E) if e != 0]
    from collections import Counter
    cnt = Counter()
    ranges = [range(-B, B+1)]*k
    for n in itertools.product(*ranges):
        if all(v == 0 for v in n): continue
        if sum(n[i]*E[i] for i in range(k)) != 0: continue
        supp_nz = sum(1 for i in nz_idx if n[i] != 0)
        if supp_nz == 0:   # supported only on zero-offset coords -> trivial
            continue
        cnt[supp_nz] += 1
    return dict(cnt)

def min_distance(E, B=3):
    sh = relation_shells(E, B)
    return min(sh) if sh else None

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0], 0) == 1

def consec(k):
    return list(range(k))

if __name__ == "__main__":
    print("="*78)
    print("PART 0 — KRAWTCHOUK IDENTITY: g_J(h) = partial Krawtchouk in N=7-h")
    print("="*78)
    print(" Claim: g_6(h) = sum_{j=0}^{6}(-1)^j C(7-h,j) = [N==0], N=7-h.")
    print(" And the FULL alternating row equals the binary Krawtchouk evaluation")
    print(" that the indicator [N=0] = (1/2^6) sum_w (number-of-weight-w) ... (MacWilliams).")
    ok = True
    for h in range(8):
        N = 7-h
        g6 = gJ(h,6)
        ind = 1 if N==0 else 0
        # full alternating sum sum_{j=0}^{N}(-1)^j C(N,j) = [N==0]
        alt = sum((-1)**j*comb(N,j) for j in range(N+1))
        mark = "" if (g6==ind==alt) else "  MISMATCH"
        if g6!=ind: ok=False
        print(f"   h={h} N={N}: g_6={g6}  [N=0]={ind}  sum(-1)^j C(N,j)={alt}{mark}")
    print(f"  -> g_6 = [N=0] identity holds: {ok}  [VERIFIED]\n")

    print("="*78)
    print("PART 1 — baseline: depth law, measS7, corr, relation min-distance d")
    print("="*78)
    for k in [8,9,10]:
        C = consec(k)
        pl = depth_law(C)
        m = pl[7]; ik = iid_k(k); cr = m - ik
        d = min_distance(C, B=2)
        sh = relation_shells(C, B=2)
        print(f" k={k} consec={C}")
        print(f"   measS7={m} = {float(m):.6f}   iid_k={float(ik):.6f}   corr={float(cr):+.6f}")
        print(f"   relation min-support d (box B=2) = {d}   shells(supp->count)={dict(sorted(sh.items()))}")
    print()

# ============================================================================
# PART 2 — the Weyl/Fourier expansion of measS7 and the carrier kernel K(n).
# ============================================================================
# measS7(E) = integral_0^1 prod_{j=1}^{6} (sector j is HIT at x) ... is an
# inclusion-exclusion over the 7 sectors of products of indicators.  The CLEAN
# Fourier route (mac HYP-2719): each sector indicator 1_{sector s hit by some e_i x}
# expands in characters e(m * e_i x); integrating over x in [0,1) kills all
# frequencies except those with a RELATION sum_i n_i e_i = 0.  So
#     measS7(E) = iid_k + sum_{0 != n in Lambda(E)} K(n),
# where K(n) is a FIXED kernel (independent of E) depending only on the integer
# vector n (its support pattern and the residues n_i mod 7 and the sign/magnitude
# structure of the sector-IE Fourier coefficients).
#
# We do NOT need the closed form of K(n) to test the MacWilliams link; we need
# the empirical fact that corr is a SHELL SUM over Lambda(E) graded by support
# (=weight in the code sense).  We verify the support-graded structure exactly:
# build corr for many shapes, and the relation-shell weight enumerator for the
# same shapes, then test the weight-enumerator <-> corr / depth-law relationship.

def relation_shells_signed(E, B=2):
    """Full shell census incl. residue mod 7 of each nonzero n_i.
       Returns: dict supp -> count, and dict supp -> sum over relations of a
       proxy 'low-frequency weight' = prod over nonzero coords of 1/|n_i| (the
       Fourier coeff of a sector indicator at freq m decays ~1/m)."""
    E=list(E); k=len(E)
    nz_idx=[i for i,e in enumerate(E) if e!=0]
    from collections import Counter
    cnt=Counter(); fw=Counter()
    for n in itertools.product(*([range(-B,B+1)]*k)):
        if all(v==0 for v in n): continue
        if sum(n[i]*E[i] for i in range(k))!=0: continue
        supp=sum(1 for i in nz_idx if n[i]!=0)
        if supp==0: continue
        cnt[supp]+=1
        w=F(1)
        for i in nz_idx:
            if n[i]!=0: w*= F(1, abs(n[i]))
        fw[supp]+=w
    return dict(cnt), {s:fw[s] for s in fw}

# ============================================================================
# PART 3 — THE MACWILLIAMS CROSS-TAB.
# Build a bank; for each shape record (corr, depth-law moments S_2,S_4,measS7,
# and relation weight-enumerator A_2,A_3,A_4 low-support counts + freq-weighted).
# Test: which relation-code statistic best PREDICTS corr / the even bands?
# ============================================================================
def even_band(E, J):
    pl = depth_law(E)
    return sum(pl[h]*gJ(h,J) for h in range(8))

def pearson(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    vx=sum((x-mx)**2 for x in xs); vy=sum((y-my)**2 for y in ys)
    if vx==0 or vy==0: return float('nan')
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/((vx*vy)**0.5)

if __name__ == "__main__":
    print("="*78)
    print("PART 3 — MacWilliams cross-tab: relation weight-enumerator vs depth-law bands")
    print("="*78)
    k=8; W=11; B=2
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E)]
    print(f" k={k} span<= {W} -> {len(bank)} primitive shapes, relation box B={B}")
    recs=[]
    for E in bank:
        cr=float(corr(E))
        pl=depth_law(E)
        s2=float(sum(pl[h]*gJ(h,2) for h in range(8)))
        s4=float(sum(pl[h]*gJ(h,4) for h in range(8)))
        m7=float(pl[7])
        cnt,fw=relation_shells_signed(E,B)
        A2=cnt.get(2,0); A3=cnt.get(3,0); A4=cnt.get(4,0)
        fw2=float(fw.get(2,0)); fw3=float(fw.get(3,0)); fw4=float(fw.get(4,0))
        dmin=min(cnt) if cnt else 99
        recs.append(dict(E=E,corr=cr,s2=s2,s4=s4,m7=m7,A2=A2,A3=A3,A4=A4,
                         fw2=fw2,fw3=fw3,fw4=fw4,dmin=dmin))
    Cm7=max(recs,key=lambda r:r['m7'])
    print(f"  measS7 maximizer = {Cm7['E']} (consec={list(range(k))}?  {list(Cm7['E'])==list(range(k))})  measS7={Cm7['m7']:.5f}")
    corrs=[r['corr'] for r in recs]
    print("\n  Pearson corr( corr(E) , relation statistic ):")
    for key in ['A2','A3','A4','fw2','fw3','fw4']:
        print(f"     {key:>4s}: {pearson(corrs,[r[key] for r in recs]):+.4f}")
    print("\n  Pearson corr( measS7 , relation statistic ):")
    for key in ['A2','A3','A4','fw2','fw3','fw4']:
        print(f"     {key:>4s}: {pearson([r['m7'] for r in recs],[r[key] for r in recs]):+.4f}")
    # Does consec have minimal d (anti-MDS)?  And maximal low-support counts?
    print("\n  ANTI-MDS check: distribution of relation min-distance d over the bank,")
    print("  and whether the measS7-max shape sits at minimal d with max low-support mass.")
    from collections import Counter
    dd=Counter(r['dmin'] for r in recs)
    print(f"     d-histogram (box B={B}): {dict(sorted(dd.items()))}")
    by_m7=sorted(recs,key=lambda r:-r['m7'])
    print("     top-5 measS7 shapes:   E, measS7, d, A2,A3,A4")
    for r in by_m7[:5]:
        print(f"       {list(r['E'])} m7={r['m7']:.4f} d={r['dmin']} A=({r['A2']},{r['A3']},{r['A4']})")
    print("     bottom-5 measS7 shapes:")
    for r in by_m7[-5:]:
        print(f"       {list(r['E'])} m7={r['m7']:.4f} d={r['dmin']} A=({r['A2']},{r['A3']},{r['A4']})")
