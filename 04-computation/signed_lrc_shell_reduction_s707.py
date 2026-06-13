#!/usr/bin/env python3
"""
signed_lrc_shell_reduction_s707.py
monad-explorer-2026-06-06-S707

ANGLE: structural reduction -- how the SIGNED-LRC additive/shell structure
(mod C=2n-1) transfers to the regular (unsigned) LRC gap (which lives mod n,
THM-403).  Builds on THM-401 (shell modulus 2n-1), THM-403 (n-th-root witness),
THM-414 / HYP-2272 (additive face = multiplicative energy; E_+, r_+(0)).

We test FOUR things rigorously:

 (A) SYMMETRIC-CLOSURE IDENTITY.  sym(AP) = {0,+-1,...,+-(n-1)} has exactly
     2n-1 elements and equals all of Z/(2n-1).  The THM-401 summand shells
     P_a={a,-a} are exactly the orbits of the sign-gauge group {+-1} on Z/C.
     => the modulus 2n-1 IS the size of the symmetric closure of the AP.

 (B) GLOBAL ADDITIVE-ENERGY MAX.  Is the AP {1,...,n-1} the GLOBAL maximizer of
     E_+ = sum_s r_+(s)^2 (r_+(s)=#{i<j: a_i+a_j=s mod C}) over ALL (n-1)-subsets
     of (Z/C)\{0}?  If yes, HYP-2272's bound "tight => E_+ <= E_+(AP)" is a
     GENERAL fact (true for every config), so its real content is the EQUALITY
     case.  We also characterize the maximizers: are they exactly the dilates
     u*AP, u in (Z/C)^x  (the regular orbit)?

 (C) TIGHT-CONFIG SEARCH (exact rational M).  Exhaustively find tight configs
     (M(S)=1/n) in an integer box, compute their (E_+, r_+(0)) mod C, confirm
     the AP energy-dominates the tight set, locate sporadics, and TEST whether
     (E_+, r_+(0)) is a complete invariant of the tight set (HYP-2272 C1).

 (D) THE REDUCTION: relate shell-structure (mod 2n-1, the SIGNED face) to
     witness-structure (mod n, the UNSIGNED gap).  Does a shell-collision force
     anything about the n-clock?  Is "tight + shell-perfect" => dilate of AP?
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd
import numpy as np

def shells(C):
    """Summand shells P_a={a,C-a} = sign-gauge orbits on (Z/C)\{0}."""
    seen=set(); out=[]
    for a in range(1,C):
        if a in seen: continue
        b=(C-a)%C
        out.append(tuple(sorted({a,b})))
        seen.add(a); seen.add(b)
    return out

def r_plus(resset, C):
    """r_+(s)=#{i<j: a_i+a_j = s mod C} as a dict."""
    rp={}
    a=list(resset)
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            s=(a[i]+a[j])%C
            rp[s]=rp.get(s,0)+1
    return rp

def E_plus(resset,C):
    rp=r_plus(resset,C)
    return sum(v*v for v in rp.values()), rp.get(0,0)

def units(C):
    return [u for u in range(1,C) if gcd(u,C)==1]

def dilate_orbit(resset,C):
    """orbit of resset under (Z/C)^x (multiplicative dilation)."""
    orb=set()
    for u in units(C):
        orb.add(frozenset((u*x)%C for x in resset))
    return orb

# ---------- exact M(S) = max_t min_i ||v_i t|| for integer speeds -----------
def frac_dist(x):
    """||x|| = distance from rational x to nearest integer, as Fraction."""
    f = x - (x.numerator//x.denominator)   # fractional part in [0,1)
    return min(f, 1-f)

def M_approx(speeds, ts):
    """float approx of M over a precomputed candidate-t array ts (numpy)."""
    v=np.array(sorted(set(speeds)),dtype=np.float64)
    P=np.outer(ts,v)          # |ts| x |v|
    fr=P-np.floor(P)
    d=np.minimum(fr,1-fr)
    mn=d.min(axis=1)
    return mn.max()

def cand_ts(B):
    """candidate breakpoints t in (0,1/2] for any speeds in [1,B]:
       m/(2v) (kinks) and m/(s) for sums/diffs s up to 2B."""
    s=set()
    for den in range(1,2*B+1):
        m=1
        while m*2<=den:        # t=m/den <= 1/2
            s.add((m,den)); m+=1
    return np.array(sorted({Fraction(m,d) for (m,d) in s}),dtype=np.float64)

def M_exact(speeds):
    """Exact M(S). Optimum of piecewise-linear min_i||v_i t|| on (0,1/2] is at a
    breakpoint: kinks t=m/(2 v_i), crossings t=m/(v_i +- v_j).  We evaluate the
    min at all such candidate t in (0,1/2] (||.|| is even & 1-periodic, so [0,1/2]
    suffices) and take the max.  Exact rational arithmetic -> rigorous."""
    S=sorted(set(speeds))
    cands=set()
    half=Fraction(1,2)
    def add(den):
        if den<=0: return
        m=1
        while True:
            t=Fraction(m,den)
            if t>half: break
            cands.add(t); m+=1
    for v in S: add(2*v)
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            add(S[i]+S[j]); add(S[j]-S[i])
    best=Fraction(0)
    for t in cands:
        mn=min(frac_dist(v*t) for v in S)
        if mn>best: best=mn
    return best

# --------------------------------------------------------------------------
def partA(nmax=9):
    print("="*70)
    print("(A) SYMMETRIC-CLOSURE IDENTITY: sym(AP)=Z/(2n-1); shells=gauge orbits")
    print("="*70)
    for n in range(3,nmax+1):
        C=2*n-1
        AP=set(range(1,n))           # {1,...,n-1}
        sym={0}|{x%C for x in AP}|{(-x)%C for x in AP}
        full=set(range(C))
        sh=shells(C)
        print(f" n={n:2d} C={C:2d}: |sym(AP)|={len(sym):2d} (=C? {len(sym)==C}); "
              f"sym==Z/C: {sym==full}; #shells={len(sh)} (=n-1? {len(sh)==n-1})")
    print()

def partB(nmax=9):
    print("="*70)
    print("(B) GLOBAL ADDITIVE-ENERGY MAX over ALL (n-1)-subsets of (Z/C)\\{0}")
    print("="*70)
    for n in range(4,nmax+1):
        C=2*n-1; k=n-1
        AP=frozenset(range(1,n))
        Eap,_=E_plus(AP,C)
        best=-1; maxers=[]
        total=0
        for combo in combinations(range(1,C),k):   # nonzero residues
            total+=1
            E,_=E_plus(combo,C)
            if E>best:
                best=E; maxers=[frozenset(combo)]
            elif E==best:
                maxers.append(frozenset(combo))
        ap_orbit=dilate_orbit(AP,C)
        ap_is_max = (Eap==best)
        maxers_set=set(maxers)
        maxers_eq_orbit = (maxers_set==ap_orbit)
        # are all maximizers dilates of AP?
        all_dilates = maxers_set.issubset(ap_orbit)
        print(f" n={n:2d} C={C:2d} k={k}: #subsets={total}, E_+(AP)={Eap}, "
              f"globalmax={best}, AP is max? {ap_is_max}")
        print(f"          #maximizers={len(maxers)}, |AP dilate-orbit|={len(ap_orbit)}, "
              f"maximizers==AP-orbit? {maxers_eq_orbit}, maximizers ⊆ AP-orbit? {all_dilates}")
        if not all_dilates:
            extra=[sorted(m) for m in maxers_set-ap_orbit][:3]
            print(f"          NON-AP maximizers exist, e.g.: {extra}")
    print()

def partC():
    print("="*70)
    print("(C) TIGHT-CONFIG SEARCH (exact M): additive invariants of the worry-set")
    print("="*70)
    # exhaustive over integer configs with speeds in [1,B], gcd=1, size k=n-1.
    for n in [5,6,7]:
        C=2*n-1; k=n-1; thr=Fraction(1,n); thrf=1.0/n
        B = 3*n           # search box; covers sporadics (n=7 (1,5,6,11,16,17), max 17)
        ts=cand_ts(B)     # shared float candidate set (rigorous superset of breakpoints)
        tight=[]
        for combo in combinations(range(1,B+1),k):
            if gcd_list(combo)!=1: continue
            Ma=M_approx(combo,ts)
            if abs(Ma-thrf)>1e-7:     # float prefilter: far from 1/n -> not tight
                continue
            if M_exact(combo)!=thr:    # exact confirmation
                continue
            E,r0=E_plus([x%C for x in combo],C)
            tight.append((combo,E,r0))
        # dedup by residue-set up to dilation for reporting
        Eap,_=E_plus(set(range(1,n)),C)
        print(f" n={n:2d} C={C:2d} thr=1/{n}, box[1,{B}]: #tight integer configs={len(tight)}")
        # group by (E,r0)
        from collections import Counter
        sig=Counter((E,r0) for (_,E,r0) in tight)
        print(f"     E_+(AP)={Eap}; distinct (E_+, r_+(0)) signatures among tight:")
        for (E,r0),cnt in sorted(sig.items(),reverse=True):
            flag=" <-- AP energy" if E==Eap else ""
            print(f"        (E_+={E:4d}, r_+(0)={r0}) : {cnt} configs{flag}")
        maxE=max(E for (_,E,_) in tight)
        print(f"     max E_+ among TIGHT = {maxE}; equals E_+(AP)? {maxE==Eap}; "
              f"AP is tight-energy-max? {maxE==Eap}")
        # examples of distinct configs at each signature
        seen=set();
        for (combo,E,r0) in sorted(tight,key=lambda z:(-z[1],z[2])):
            key=(E,r0)
            if key in seen: continue
            seen.add(key)
            print(f"        example (E_+={E},r0={r0}): speeds={combo} -> res mod {C}={sorted(x%C for x in combo)}")
    print()

def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

def partD():
    print("="*70)
    print("(D) REDUCTION: shell clock (mod 2n-1, SIGNED) vs witness clock (mod n)")
    print("="*70)
    # The n=14 floor configs (from repo): AP, 2*AP, V*.
    n=14; C=2*n-1
    AP=list(range(1,14))                      # {1,...,13}
    AP2=[2*x for x in AP]                      # 2,4,...,26
    V =[1,2,3,4,5,6,7,8,9,10,11,13,24]
    for name,S in [("AP",AP),("2*AP",AP2),("V*",V)]:
        res=[x%C for x in S]
        E,r0=E_plus(res,C)
        sh=shells(C)
        # shell occupancy
        occ={}
        for x in res:
            for P in sh:
                if x%C in P:
                    occ[P]=occ.get(P,0)+1
        empties=[P for P in sh if occ.get(P,0)==0]
        doubles=[P for P in sh if occ.get(P,0)>=2]
        # witness clock mod n: residues, collisions, zeros
        wn=[x%n for x in S]
        from collections import Counter
        wc=Counter(wn)
        zero_n = wc.get(0,0)
        coll_n = sum(c-1 for c in wc.values() if c>1)
        print(f" {name:5s}: res mod {C}={sorted(res)}")
        print(f"        E_+={E}, r_+(0)=shell-collisions={r0}; empty shells={[P for P in empties]}; doubled shells={[P for P in doubles]}")
        print(f"        mod n={n}: residues={sorted(wn)}; #(==0)={zero_n}; #witness-collisions={coll_n}")
    print()
    print(" Interpretation:")
    print("  * shell-perfect (r_+(0)=0) <=> one runner per summand shell <=> sym-closure")
    print("    is a transversal; AP & 2*AP are shell-perfect, V* has the 3+24=27 collision.")
    print("  * BUT all three are tight (gap 1/14) -- so shell structure does NOT change M;")
    print("    it refines the tight set.  The witness clock mod n only forbids residue 0;")
    print("    it is BLIND to the shell collision.  => signed face is strictly finer.")

if __name__=="__main__":
    partA(9)
    partB(9)
    partC()
    partD()
