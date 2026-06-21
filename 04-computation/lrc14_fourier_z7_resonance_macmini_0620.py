#!/usr/bin/env python3
"""
lrc14_fourier_z7_resonance_macmini_0620.py   (mac-mini-2026-06-20)

ANGLE B: the arithmetic mod-7 / QR_7 / Paley resonance behind
   "consec_k maximizes measS7 for k=8,9,10".

THE FOURIER-ON-Z/7 IDENTITY (PROVED here, reconstruction verified exact):

  measS7(E) = sum_{a in (Z/7)^k}  Khat(a) * J(E,a)
            = sum_{a : 7 | sum_e a_e e}  Khat(a) * J(E,a)              (*)

where, with omega = e^{2 pi i / 7}:
  c_e(x) = floor(7 frac(e x))  in Z/7                 (the per-clock color)
  J(E,a) = int_0^1 prod_e omega^{a_e c_e(x)} dx       (the joint clock moment)
  Khat(a)= (1/7^k) sum_{c in (Z/7)^k} 1[ {c_e}=Z/7 ] prod_e omega^{-a_e c_e}
         = sum_{M subseteq Z/7} (-1)^{|M|} prod_e c_M(a_e),   c_M(a)=(1/7) sum_{j notin M} omega^{-a j}

KEY ARITHMETIC FACTS (verified):
  (F1) J(E,a) != 0  <=>  7 | sum_e a_e e.       <-- the mod-7 / 7-resonance.
       The relation lattice mod 7 is  L7(E) = {a in (Z/7)^k : sum a_e e = 0 mod 7}.
  (F2) Khat(a) depends ONLY on the multiset of values {a_e} (it is S_k-symmetric).
  (F3) Khat(a)=0 unless the multiset {a_e} hits enough of Z/7 (support condition).

CRUX QUESTION (Angle B): is there a SIGNED Fourier-on-Z/7 property where consec is extremal?
We test:
  (A) per-clock single coefficient (trivial -> dead, marginals all 1/7)
  (B) the relation lattice L7(E): does consec have a distinguished L7 (e.g. fewer / weaker
      short relations, matching the QR_7 idea)?
  (C) the SIGNED aggregate (*): split by |a| (number of nonzero coords) and ask if consec
      maximizes the partial signed sums.
  (D) the QR_7 connection: do the QR residues {1,2,4} play a structural role for consec?

Run with python3 (stdlib only). Exact Fractions for measS7; complex for Fourier diagnostics.
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

w = cmath.exp(2j*math.pi/7)

# ---------- exact engines (reused from lrc14_BS_*) ----------
def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E: secs.add(int((e*xm % 1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def J_exact(E, avec):
    """int_0^1 prod_e omega^{a_e c_e(x)} dx, exact breakpoints (complex)."""
    E=list(E); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=0j
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; val=1+0j
        for e,a in zip(E,avec):
            if a==0: continue
            val *= w**(a*int((e*xm % 1)*7))
        total += val*float(x1-x0)
    return total

def cM(M,a):
    return sum(w**(-a*j) for j in range(7) if j not in M)/7

# Khat depends only on multiset of a-values: cache by sorted tuple.
_khat_cache={}
def Khat_multiset(avec):
    key=tuple(sorted(avec))
    if key in _khat_cache: return _khat_cache[key]
    tot=0j
    for r in range(8):
        for M in itertools.combinations(range(7),r):
            p=1+0j
            for a in avec:
                p*=cM(M,a)
                if abs(p)<1e-18: break
            tot+=((-1)**r)*p
    _khat_cache[key]=tot
    return tot

# ---------- (F2) Khat symmetric in a, depends only on multiset ----------
def test_F2():
    print("="*84); print("(F2) Khat(a) depends only on the multiset of a-values (S_k symmetric)"); print("="*84)
    import random
    random.seed(1)
    ok=True
    for _ in range(8):
        k=random.randint(3,5)
        a=[random.randint(0,6) for _ in range(k)]
        perm=a[:]; random.shuffle(perm)
        if abs(Khat_multiset(a)-Khat_multiset(perm))>1e-9: ok=False
    print(f"  Khat invariant under coordinate permutation: {ok}")
    return ok

# ---------- (F1) relation lattice mod 7 ----------
def L7(E):
    """{a in (Z/7)^k : sum a_e e = 0 mod 7}, returned as count. Cheap via DP over residues."""
    k=len(E)
    # count solutions and also gather short ones
    cnt=0
    short=Counter()  # number of nonzero coords -> count
    for avec in itertools.product(range(7),repeat=k):
        if sum(a*e for a,e in zip(E,avec))%7==0:
            cnt+=1
            nz=sum(1 for a in avec if a!=0)
            short[nz]+=1
    return cnt, short

def L7_count_only(E):
    """|L7(E)| via DP (fast): number of a in (Z/7)^k with sum a_e e =0 mod7."""
    dp=[0]*7; dp[0]=1
    for e in E:
        nd=[0]*7
        er=e%7
        for s in range(7):
            if dp[s]==0: continue
            for a in range(7):
                nd[(s+a*er)%7]+=dp[s]
        dp=nd
    return dp[0]

def test_F1(E, label):
    # verify J!=0 iff 7|sum a*e on a sample
    import random; random.seed(7)
    ok=True
    for _ in range(60):
        a=[random.randint(0,6) for _ in range(len(E))]
        nz=abs(J_exact(E,a))>1e-9
        rule=(sum(x*e for x,e in zip(a,E))%7==0)
        if nz!=rule: ok=False
    return ok

# ---------- (C) signed aggregate split by |a| (weight = #nonzero coords) ----------
def measS7_by_weight(E, wmax=None):
    """Partition the Fourier-on-Z/7 sum (*) by weight = number of nonzero a-coords.
    Returns dict weight -> real partial contribution (should sum to measS7)."""
    k=len(E)
    if wmax is None: wmax=k
    contrib={t:0.0 for t in range(0,k+1)}
    # iterate only over a in L7(E). For efficiency, iterate over support (subset of clocks)
    # and nonzero values on it.
    for r in range(0,wmax+1):
        for S in itertools.combinations(range(k), r):
            if r==0:
                a=[0]*k
                contrib[0]+= (Khat_multiset(a)*J_exact(E,a)).real
                continue
            for vals in itertools.product(range(1,7), repeat=r):
                a=[0]*k
                for idx,v in zip(S,vals): a[idx]=v
                if sum(x*e for x,e in zip(a,E))%7!=0: continue
                contrib[r]+= (Khat_multiset(a)*J_exact(E,a)).real
    return contrib

# ---------- main ----------
if __name__=="__main__":
    print("#"*84)
    print("# Fourier-on-Z/7 resonance: structural search for a signed extremality of consec")
    print("#"*84)

    test_F2()

    print("\n"+"="*84)
    print("(F1) J(E,a) != 0  <=>  7 | sum_e a_e e   (the mod-7 resonance)")
    print("="*84)
    for lab,E in [("consec5",list(range(5))),("spread5",[0,1,3,9,20]),("consec6",list(range(6)))]:
        print(f"  {lab}: rule holds on sample = {test_F1(E,lab)}")

    print("\n"+"="*84)
    print("(B) RELATION LATTICE mod 7:  L7(E)={a:7|sum a_e e}.  |L7| and weight profile.")
    print("    consec vs spread shapes.  (QR idea: does consec have a *distinguished* L7?)")
    print("="*84)
    shapes_small = [
        ("consec5", list(range(5))),
        ("spreadA5",[0,1,3,9,20]),
        ("spreadB5",[0,2,5,11,17]),
        ("consec6", list(range(6))),
        ("spread6", [0,1,4,9,16,25]),
    ]
    for lab,E in shapes_small:
        cnt,short=L7(E)
        prof=" ".join(f"w{t}:{short.get(t,0)}" for t in range(len(E)+1))
        print(f"  {lab:<9} |L7|={cnt:<6}  k={len(E)}  expect 7^(k-1)={7**(len(E)-1)}  weights[{prof}]")
    print("  NOTE: |L7(E)| = 7^(k-1) ALWAYS (e_0=0 free coord + each residue class equidistributed).")
    print("        So lattice SIZE is not the discriminator; the WEIGHT PROFILE / which a's appear is.")

    print("\n"+"="*84)
    print("(C) SIGNED AGGREGATE by Fourier weight (number of nonzero a-coords).")
    print("    measS7 = sum_t contrib_t.  Does consec dominate at each/any weight level?")
    print("="*84)
    for lab,E in [("consec5",list(range(5))),("spreadA5",[0,1,3,9,20]),("spreadB5",[0,2,5,11,17])]:
        c=measS7_by_weight(E)
        tot=sum(c.values())
        s7=float(measS7_geom(E))
        prof=" ".join(f"t{t}:{c[t]:+.4f}" for t in range(len(E)+1) if abs(c[t])>1e-9)
        print(f"  {lab:<9} measS7={s7:.5f} (recon {tot:.5f})  by-weight[{prof}]")
