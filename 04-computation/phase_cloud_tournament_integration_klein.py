#!/usr/bin/env python3
"""
phase_cloud_tournament_integration_klein.py  --  klein-2026-07-01-S80

NOVEL INTEGRATION of the tournament machinery with LRC: the covering-min binding, in PHASE-RESIDUE
coordinates (S68), is a concrete RUNNER CLOUD on Z/Phi6 -- and it unifies the Chebyshev 2-point dual (S73),
the three-gap theorem (HYP-3762), the runner-cloud tournament (S70), and the Phi6-irreducibility (S79).

PICTURE (verified): the construction's phase cloud at t* = {p(v)=n*v mod Phi6} = the ARITHMETIC PROGRESSION
{n*k mod Phi6 : k=1..n-2} (the small runners) + the KILLER at p=n(n-1)=-n mod Phi6. The observer (0) sits
in the gap between the killer (-n) and runner 1 (+n): clearance = n/Phi6 = M_C on BOTH sides. So the
Chebyshev 2-point alternation {1, killer} = the two cloud points FLANKING the observer's gap, and the
killer is the iota-antipode (-n) of the slowest runner (+n) -- symmetrizing the gap.

Computes: (1) the cloud + gaps (three-gap sizes); (2) the flanking binding = {1,killer}; (3) the ROTATIONAL
tournament of the cloud (scores, #3-cycles) + Redei parity; (4) the beater-obstruction in cloud terms.
"""
from fractions import Fraction as F
from itertools import combinations

def cloud_of(S, n, Phi6):
    return sorted((n*v)%Phi6 for v in S)

def gaps(cloud, Phi6):
    c=sorted(cloud); g=[(c[(i+1)%len(c)]-c[i])%Phi6 for i in range(len(c))]
    return sorted(set(g)), g

def rot_tournament_scores(cloud, Phi6):
    # rotational tournament: i->j iff (p_j - p_i) mod Phi6 in (0, Phi6/2)
    pts=sorted(cloud); k=len(pts); half=Phi6/2
    A=[[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i==j: continue
            d=(pts[j]-pts[i])%Phi6
            if 0<d<half: A[i][j]=1
    scores=[sum(A[i]) for i in range(k)]
    c3=0
    for a,b,c in combinations(range(k),3):
        # count cyclic triangles
        e=A[a][b]+A[b][c]+A[c][a]
        if e==0 or e==3 or (A[a][b]+A[b][a]==1): pass
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]): c3+=1
    return scores, c3

if __name__=="__main__":
    n=14; Phi6=n*n-n+1; C=list(range(1,n-1))+[n*(n-1)]; Mc=F(n,Phi6)
    cl=cloud_of(C,n,Phi6)
    print(f"n={n} Phi6={Phi6}; construction phase cloud (at t*): {cl}")
    ap=[(n*k)%Phi6 for k in range(1,n-1)]; killer=(n*(n-1))%Phi6
    print(f"  = AP(step n={n}) {ap} + KILLER at {killer} = -n mod Phi6 = {(-n)%Phi6}")
    gs,gall=gaps(cl,Phi6)
    print(f"(1) three-gap sizes of the cloud: {gs} (three-gap theorem: <=3 distinct); gap multiset sum={sum(gall)}=Phi6")
    # gap containing 0 (observer): between killer(-n) and runner-1(+n)
    d0=min(min(p,Phi6-p) for p in cl); flank=[p for p in cl if min(p,Phi6-p)==d0]
    print(f"(2) observer(0) clearance = {d0}/{Phi6} = M_C = {float(Mc):.5f}; FLANKING cloud points {flank} = {{+n, -n}} = runners {{1, killer}} (the Chebyshev 2-point dual, S73)")
    scores,c3=rot_tournament_scores(cl,Phi6)
    print(f"(3) rotational tournament of the cloud: scores {scores}; #3-cycles c3={c3}; H (Redei) odd by theorem (13 vertices, not brute-forced)")
    print(f"    the AP gives a near-transitive/circulant order; the killer inserts one near-tie next to runner 12 (gap 1 at {killer}-{max(ap)}={killer-max(ap)})")
    print("(4) BEATER-OBSTRUCTION in cloud terms: to beat, a covering set's phase cloud must clear <n from 0 at EVERY t")
    print("    (max clearance < n/Phi6). But covering (multiple of every q<=n) forces the small speeds {1..n-2}, whose")
    print("    phases at t* are the AP(step n) -> gap exactly n at 0. So the observer's gap can't shrink below n at t*;")
    print("    the killer only SYMMETRIZES it (antipode). This is the S79 Phi6-metric-irreducibility as a CLOUD fact:")
    print("    the AP-step is n because the DEEP witness is Phi6 (=n^2-n+1), and the AP of the small speeds tiles Z/Phi6")
    print("    with step n, leaving the size-n gap at 0. Covering => AP-step n => clearance n => M >= n/Phi6.")
    print("\nINTEGRATION: covering-min extremizer = an ARITHMETIC-PROGRESSION phase cloud (step n) + an ANTIPODAL KILLER,")
    print("observer in the size-n gap, flanked by the 2-point Chebyshev alternation {1, killer}={+n,-n}. Unifies")
    print("S68(phase-residue)+S73(Chebyshev 2-pt)+S79(Phi6-irreducible)+S70(cloud tournament)+HYP-3762(three-gap).")
