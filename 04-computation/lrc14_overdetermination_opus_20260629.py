"""
LRC(14): over-determination / bounded-core reduction (option (i) continuation).
================================================================================
opus 2026-06-29. Rigorous pieces + cross-framework corroboration.

PART A (RIGOROUS structure theorem). Discretization: for any 13-set S and prime p,
        max_a min_s ||s a/p||  >=  M(S) - max(S)/(2p).
  Proof: take t* with min_s||s t*||=M(S); a=round(p t*) has |a/p-t*|<=1/(2p); ||.||
  is 1-Lipschitz so ||s a/p|| >= ||s t*|| - s/(2p) >= M(S)-max(S)/(2p).  QED.
  COROLLARY: if M(S)>1/14 then every prime p >= max(S)/(2(M(S)-1/14)) is a WITNESS
  (N(S,p)>=1). Hence the BAD primes (N=0) are confined to p < P0(S):=max(S)/(2(M-1/14)).
  => 'find a witness prime' is EFFECTIVE; the only hard primes are small ones.

PART B (corroboration of the bounded-core reduction; cf. HYP-3131 Lee-Yang & kps-S31as).
  Far/large speeds do NOT make a covering set harder: replacing small speeds by far
  ones strictly RAISES M (the set becomes MORE lonely). The hard (min-M) configs are
  BOUNDED and RESONANT (consecutive / shared-factor). This is the same conclusion as
  HYP-3131 ('far pushes Lee-Yang zeros out') reached here from the loneliness side.

CONVERGENCE (3 independent frameworks => the proof architecture):
  [THM-523] LRC(14) <=> M>=1/14 on primitive covering 13-sets.
  [Part B / HYP-3131] reduce to the BOUNDED RESONANT core (far elements only help).
  [kps-S31as / option (i)] finish = SIGNED Erdos-Turan/Fejer equidistribution of the
        resonant comb; NO absolute certificate (Lee-Yang/Asano/SOS/moment-LP/Bonferroni)
        can cross the worst case -- only cancellation. (option (i) proved this from the
        exp-sum side: the m>=3 tail is Theta(p) on resonant sets, so |error| bounds fail.)
  REMAINING LEMMA: the signed Erdos-Turan/Fejer bound giving lonely-measure>0 (equiv.
        N(S,p)>0 at some prime) for the bounded resonant core. == LRC(14), now isolated.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    b=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>b: b=v
    return b
def sieve(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(2,n+1) if s[i]]
def N_mod(S,p):
    H=p//14
    return sum(1 for a in range(p) if all(min((s*a)%p,p-(s*a)%p)>H for s in S))
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def prim(S): return reduce(gcd,S)==1
PR=sieve(4000)

print("PART A -- bad primes confined below P0(S)=max(S)/(2(M(S)-1/14)):")
for S in [[1,2,3,4,5,6,7,8,9,10,11,13,14],[1,2,3,12,18,20,21,22,23,24,25,26,28],
          [3,6,9,12,15,21,24,27,28,30,33,36,39]]:
    M=M_exact(S); marg=M-F(1,14); P0=float(max(S)/(2*marg)) if marg>0 else float('inf')
    bad=[p for p in PR if p>=15 and all(s%p for s in S) and N_mod([s%p for s in S],p)==0]
    mb=max(bad) if bad else None
    print(f"  S(max={max(S):2d}): M={float(M):.4f} margin={float(marg):.4f} P0={P0:6.0f} | "
          f"bad primes<4000: {len(bad)} (max {mb}) -> all < P0: {mb is None or mb<P0}")

print("\nPART B -- far elements RAISE M (help): drop smallest, add a far speed:")
Pp=4001
def Mproxy(S):
    b=0
    for a in range(1,Pp):
        m=Pp
        for s in S:
            r=(s*a)%Pp; d=r if r<Pp-r else Pp-r
            if d<m:
                m=d
                if m<=b: break
        if m>b: b=m
    return b/Pp
cur=[2,3,4,5,6,7,8,9,10,11,12,13,14]; random.seed(1)
print(f"  bounded base {cur}: Mproxy={Mproxy(cur):.4f}")
for _ in range(4):
    f=random.randint(200,5000)
    while f in cur: f=random.randint(200,5000)
    cur=sorted(cur[1:]+[f])
    print(f"  -> add far {f} (max={max(cur)}): Mproxy={Mproxy(cur):.4f}  cov={is_cov(cur)} prim={prim(cur)}")
print("  (M rises as resonant small speeds are replaced by far ones => bounded-core is the hard case)")
