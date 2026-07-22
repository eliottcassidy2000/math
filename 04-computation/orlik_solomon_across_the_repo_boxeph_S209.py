#!/usr/bin/env python3
"""orlik_solomon_across_the_repo_boxeph_S209.py -- boxeph-2026-07-21-S209

Where else does Orlik-Solomon (hyperplane/toric arrangement) structure govern repo objects?
Abstract OS signature: (geometric/arithmetic lattice of flats) + (Mobius/inclusion-exclusion count)
+ (characteristic polynomial via finite-field point-count) + (complement cohomology = OS algebra)
+ (localization-at-a-flat product factorization). We verify FOUR arrangement types the repo meets:

  P1 BRAID  A_{n-1}  <->  tournaments / NC2 Vandermonde (S208): OS Poincare = prod(1+kt),
            Betti = Stirling-1st -> a GRADED cohomology invariant of tournament/ordering space.
  P2 SHI / deformed braid <-> LRC resonance {x_i - x_j = integer}: finite-field char polynomial.
  P3 TORIC arrangement (De Concini-Procesi) <-> LRC relation lattice {k.v=0} (THM-1820):
            good-set measure |G_delta| = toric-complement volume = arithmetic-Mobius sum
            = the repo's LRCMod mod-q ladder (finite-field method on the torus).
  P4 TIGHT = RELATION-RICHEST = max Mobius/Betti: the AP (1,2,3) maximizes arrangement richness.
"""
from math import comb, factorial, sin, pi, gcd
from fractions import Fraction as F
from itertools import product

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

# ---------------------------------------------------------------------------
sep("P1  BRAID OS Poincare = prod_{k=1}^{n-1}(1+kt) = Stirling-1st -> tournament cohomology")
# The complement M(A_{n-1}) = ordered-configuration space; its cohomology (OS algebra) has
# Poincare polynomial pi(t) = prod_{k=1}^{n-1}(1+kt). Betti b_i = c(n, n-i) (unsigned Stirling 1st
# kind = #permutations of n with n-i cycles). This is a GRADED invariant of tournament/ordering
# space, strictly finer than char_A (which is only the top data). Verify the identity.
def stirling1_unsigned(n):
    # row of |s(n,j)|, j=0..n ; c(n,j)=#perms with j cycles
    row=[0]*(n+1); row[0]=1 if n==0 else 0
    prev=[1]
    for m in range(1,n+1):
        cur=[0]*(m+1)
        for j in range(1,m+1):
            cur[j]=prev[j-1]+(m-1)*(prev[j] if j<len(prev) else 0)
        prev=cur
    return prev  # |s(n,j)| for j=0..n
def poincare_braid(n):
    # coefficients of prod_{k=1}^{n-1}(1+k t)
    p=[1]
    for k in range(1,n):
        p=[ (p[i] if i<len(p) else 0) + k*(p[i-1] if i-1>=0 else 0) for i in range(len(p)+1)]
    return p
for n in range(2,8):
    pc=poincare_braid(n)                          # b_0..b_{n-1}
    st=stirling1_unsigned(n)                       # |s(n,j)|, j=0..n ; b_i = |s(n, n-i)|
    betti_from_stirling=[st[n-i] for i in range(n)]
    print(f"  n={n}: Poincare {pc}  == Stirling-1st [c(n,n-i)] {betti_from_stirling}? {pc==betti_from_stirling}"
          f"  (top Betti={pc[-1]}=(n-1)!={factorial(n-1)})")
print("  => tournament/ordering space has cohomology pi(t)=prod(1+kt); Betti = Stirling 1st kind.")
print("     A NEW graded lens on the tournament continuum, finer than char_A.")

# ---------------------------------------------------------------------------
sep("P2  finite-field characteristic polynomial: braid (falling fact) & SHI (deformed braid)")
# Athanasiadis finite-field method: chi_A(q) = #points of F_q^n on NO hyperplane, for q >> 0.
# BRAID {x_i=x_j}: distinct coords -> chi=q(q-1)...(q-n+1); regions=(-1)^n chi(-1)=n!.
# SHI {x_i-x_j in {0,1}, i<j}: chi=q(q-n)^{n-1}; regions=(n+1)^{n-1} (parking functions).
def count_offall_braid(n,q):
    c=0
    for pt in product(range(q),repeat=n):
        if len(set(pt))==n: c+=1
    return c
def count_offall_shi(n,q):
    c=0
    for pt in product(range(q),repeat=n):
        ok=True
        for i in range(n):
            for j in range(i+1,n):
                d=(pt[i]-pt[j])%q
                if d==0 or d==1: ok=False;break     # x_i=x_j or x_i-x_j=1
            if not ok:break
        if ok:c+=1
    return c
for n in (3,4):
    for q in (11,13):                # q>>n primes
        cb=count_offall_braid(n,q); fb=1
        for i in range(n): fb*=(q-i)
        cs=count_offall_shi(n,q);  fs=q*(q-n)**(n-1)
        print(f"  n={n},q={q}: braid #off-all={cb}==q(q-1)..={fb}?{cb==fb} | shi #off-all={cs}==q(q-n)^(n-1)={fs}?{cs==fs}")
print(f"  braid regions=n! ; shi regions=(n+1)^(n-1): n=3 -> {factorial(3)}, {4**2}=16 (parking fns).")
print("  => LRC resonances x_i - x_j = integer ARE a deformed braid arrangement; the finite-field")
print("     point-count (Diophantine!) is the native tool -- matching LRC's number-theoretic flavor.")

# ---------------------------------------------------------------------------
sep("P3  LRC = TORIC arrangement: |G_delta| = complement volume = arithmetic-Mobius (relation lattice)")
# THM-1820 bridge: int_0^1 prod_j g(v_j t) dt = sum_{k: k.v=0} prod_j ghat(k_j).
# g = far-set indicator 1[||x||>=delta]: ghat(0)=1-2delta, ghat(k)=-sin(2 pi k delta)/(pi k).
# LHS = toric-arrangement complement volume (good-set measure); RHS = sum over the RELATION LATTICE
# {k.v=0} = the LAYERS of the toric arrangement (De Concini-Procesi). Verify they match.
def ghat(k,delta):
    if k==0: return 1-2*delta
    return -sin(2*pi*k*delta)/(pi*k)
def Gmeasure_direct(v,delta,M=200000):
    n=len(v); cnt=0
    for i in range(M):
        t=(i+0.5)/M
        if all( min((v[j]*t)%1,1-(v[j]*t)%1) >= delta for j in range(n)): cnt+=1
    return cnt/M
def Gmeasure_lattice(v,delta,K=60):
    n=len(v); s=0.0
    # enumerate k with |k_j|<=K and k.v=0
    if n==3:
        a,b,c=v
        for k1 in range(-K,K+1):
            for k2 in range(-K,K+1):
                num=-(a*k1+b*k2)
                if num % c: continue
                k3=num//c
                if abs(k3)>K: continue
                s+=ghat(k1,delta)*ghat(k2,delta)*ghat(k3,delta)
    return s
for v,delta in [((1,2,3),0.20),((1,2,3),0.15),((1,3,5),0.15)]:
    d=Gmeasure_direct(v,delta); l=Gmeasure_lattice(v,delta)
    print(f"  v={v} delta={delta}: |G| direct={d:.5f}  relation-lattice sum={l:.5f}  match? {abs(d-l)<3e-3}")
print("  => good-set measure = toric-complement volume = arithmetic-Mobius sum over layers {k.v=0}.")
print("     The repo's LRCMod-mod-q ladders (THM-1820) ARE the finite-field method on the torus.")

# ---------------------------------------------------------------------------
sep("P4  TIGHT = RELATION-RICHEST = max arrangement Mobius/Betti richness: AP (1,2,3) wins")
# THM-1820 (B2): the ONLY delta*=1/4 tight n=3 family is (1,2,3); it is the RELATION-RICHEST
# (most small relations = most toric-arrangement layers = highest arithmetic-Mobius / Betti mass).
def relation_richness(v,B=2):
    n=len(v); cnt=0
    for k in product(range(-B,B+1),repeat=n):
        if any(k) and sum(k[j]*v[j] for j in range(n))==0: cnt+=1
    return cnt//2   # count +/-k once
triples=[(a,b,c) for a in range(1,7) for b in range(a+1,8) for c in range(b+1,10) if gcd(gcd(a,b),c)==1]
scored=sorted(((relation_richness(t),t) for t in triples),reverse=True)
print("  relation richness N_R (small relations |k|<=2) over primitive n=3 triples, top 6:")
for nr,t in scored[:6]:
    print(f"     {t}: N_R={nr}" + ("   <- the AP / tight extremal" if t==(1,2,3) else ""))
maxnr=scored[0][0]; ap=relation_richness((1,2,3))
print(f"  (1,2,3) N_R={ap}; is it the max? {ap==maxnr}  (tight = relation-richest = max Betti/Mobius mass)")

sep("SUMMARY")
print("""  ONE structure -- an arrangement with a Mobius lattice of flats/layers -- recurs across the repo:
   BRAID  (tournaments/NC2, S208): OS cohomology pi=prod(1+kt), Stirling Betti = graded lens.
   SHI/deformed braid (LRC resonance x_i-x_j=integer): finite-field char polynomial (Diophantine).
   TORIC / De Concini-Procesi (LRC relation lattice {k.v=0}, THM-1820): good-set = complement volume
     = arithmetic-Mobius sum = LRCMod mod-q = finite-field method on the torus. Layer-localization is
     the toric analog of the S208 flat-factorization -> a candidate Wall-A tool.
   TIGHT = RELATION-RICHEST: the AP maximizes toric-arrangement Betti/Mobius mass (THM-1820 B2).
  LEVERAGE: (a) tournament cohomology beyond char_A; (b) the finite-field/arithmetic-Mobius method as
  the native engine for LRC characteristic quasi-polynomials; (c) toric layer-localization factoring
  |G_delta| near a resonance, mirroring the braid flat-factorization that gave HYP-8775.""")
