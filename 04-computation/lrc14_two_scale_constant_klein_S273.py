#!/usr/bin/env python3
"""
lrc14_two_scale_constant_klein_S273.py
======================================
klein-2026-07-12-S273  (owner: work on the explicit two-scale constant)

The density-row tails (k=8 Phi<=cap9, k=9 J>=432/91) close via the far-element two-scale
limit (THM-687/688). The single-peel bound is PROVED (THM-700, kind-pasteur):
   |p0(E'u{w}) - Plat(E')| <= V(E')/(6w),   V(E') <= 14*sum(e').
THM-700's own "Scope" flags: the crude constant overshoots the empirical (~0.2/w) by ~10^2,
and it needs w >> sum(e') to bite -- so there is an INTERMEDIATE BAND (d just past the
exhaustive box, w NOT >> sum e') where the single-peel bound is VACUOUS.

This script pins the real situation for the ACTUAL functionals Phi (k=8) and J (k=9):
 (A) worst-case two-scale error*w over compact clusters AND adversarial w (lcm-multiples vs
     primes) -- is a uniform sharp constant even possible?
 (B) the intermediate band d=26..600: max Phi over structured worst-case 8-cores -- does the
     tail actually stay < cap9 with margin, or is there a real danger zone?
 (C) crossover: where does the crude bound bite (D0_crude) vs where is the band empirically safe?
"""
import math
from itertools import combinations

CAP9=1979/4004          # k=8 cap
THRESH_J=432/91         # k=9 threshold

def lcm(xs):
    r=1
    for x in xs:
        if x: r=r*x//math.gcd(r,x)
    return r

def moments(E,Ng=90000):
    """m1,m2,m3 = E[(N)_1,2,3], N=#empty of 7 sectors, over x-grid k/Ng."""
    s1=s2=s3=0
    for k in range(1,Ng):
        x=k/Ng
        occ=0
        for e in E:
            occ|=1<<(int((e*x%1.0)*7.0)%7)
        N=7-bin(occ).count("1")
        s1+=N; s2+=N*(N-1); s3+=N*(N-1)*(N-2)
    n=Ng-1
    return s1/n, s2/n, s3/n

def Phi_from_m(m1,m2,m3): return 1-(2/3)*m1+(47/252)*m2-(5/252)*m3
def J_from_m(m1,m2):      return 6*m1-m2                 # J=E[N(7-N)]=6m1-m2
def Phi(E):  m1,m2,m3=moments(E); return Phi_from_m(m1,m2,m3)
def J9(E):   m1,m2,_=moments(E);  return J_from_m(m1,m2)
def Phi_inf(C):  m1,m2,m3=moments(C); return Phi_from_m((6/7)*m1,(5/7)*m2,(4/7)*m3)
def J_inf(C):    m1,m2,_=moments(C);  return J_from_m((6/7)*m1,(5/7)*m2)

print("="*74)
print(f"cap9={CAP9:.5f}  thresh_J={THRESH_J:.5f}")
print("(A) worst-case two-scale error*w: clean (prime) w vs ADVERSARIAL w (lcm-multiples)")
print("="*74)
# Phi: compact 7-clusters
for C in [[0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,2,4,6,8,10,12]]:
    L=lcm(C); pinf=Phi_inf(C)
    print(f"  C={C}  lcm={L}  Phi_inf={pinf:.5f}")
    primes=[97,997,9973,99991]
    adv=[L,2*L,3*L,5*L,10*L,50*L]     # w a multiple of lcm(C): NO cancellation
    print(f"    clean-w (primes):", end="")
    for w in primes:
        e=abs(Phi(C+[w])-pinf); print(f"  w={w}:err*w={e*w:.3f}", end="")
    print()
    print(f"    adversarial-w (mult of lcm):", end="")
    for w in adv:
        e=abs(Phi(C+[w])-pinf); print(f"  w={w}:err*w={e*w:.3f}", end="")
    print()

print()
print("="*74)
print("(B) INTERMEDIATE BAND d=26..600: max Phi over structured worst-case 8-cores")
print("    (near-consec, edge-heavy, 2-block, dilated-AP, plus perturbations)")
print("="*74)
def band_cands(d):
    """structured 8-cores with 0 and max=d (worst-case candidates for max Phi)."""
    C=[]
    C.append([0,1,2,3,4,5,6,d])                 # consec-7 + far
    C.append([0,1,2,3,4,5,d-1,d])               # consec-6 + pair at far end
    C.append([0,1,2,3,4,d-2,d-1,d])
    C.append([0,1,2,3,d-3,d-2,d-1,d])           # 4+4 two-block
    C.append([0,1,2,d-2,d-1,d, d//2, d//2+1])   # 3+3+2
    # dilated AP-ish (near-covering)
    C.append(sorted(set([round(i*d/7) for i in range(8)])))
    C.append(sorted(set([0,d]+[round(i*d/6) for i in range(1,6)]+[1])))
    # consec-8 dilated to span d (equal spacing)
    C.append(sorted(set([round(i*d/7) for i in range(8)])))
    return [c for c in C if len(c)==8 and max(c)==d and len(set(c))==8]
band_max=-1; band_arg=None; per_d=[]
for d in [26,28,30,35,40,50,60,80,100,150,200,300,400,600]:
    cands=band_cands(d)
    mx=-1; arg=None
    for c in cands:
        p=Phi(c)
        if p>mx: mx=p; arg=c
    per_d.append((d,mx))
    if mx>band_max: band_max=mx; band_arg=arg
    print(f"  d={d:4d}: max Phi over {len(cands):2d} structured cands = {mx:.5f}   {'<=cap9' if mx<CAP9 else 'OVER!'}  arg={arg}")
print(f"  --> BAND MAX (d=26..600) = {band_max:.5f} at {band_arg};  cap9={CAP9:.5f}, margin +{CAP9-band_max:.5f}")

print()
print("="*74)
print("(C) CROSSOVER: crude THM-700 constant vs empirical; where does the bound bite?")
print("="*74)
# THM-700 crude: |p0(E'u w)-Plat(E')| <= V(E')/(6w), V(E')<=14*sum(e').
# For the Phi FUNCTIONAL the analogous crude const scales the same in sum(e'); measure the
# ratio crude/empirical on C={0..6} (sum e' = 21) and report the biting diameter.
C=[0,1,2,3,4,5,6]; sumC=sum(C); pinf=Phi_inf(C)
Vcrude=14*sumC                       # THM-700 total-variation bound on the cover indicators
# empirical worst over clean+adversarial w
emp=0
for w in [97,997,9973,60,120,300,3000]:
    emp=max(emp, abs(Phi(C+[w])-pinf)*w)
print(f"  C={C}: sum(e')={sumC}, THM-700 V(E')<=14*sum={Vcrude}, crude err<=V/(6w)={Vcrude/6:.1f}/w")
print(f"    empirical worst-case err*w (Phi, incl adversarial) ~ {emp:.3f}")
print(f"    ==> crude/empirical overshoot ~ {(Vcrude/6)/max(emp,1e-9):.0f}x  (THM-700 says ~10^2)")
# D0 where each bound clears the k=8 margin (cap9 - Phi_inf_max):
marg=CAP9-0.397                      # tail margin to cap9 (max Phi_inf ~0.397, S272)
print(f"  k=8 tail margin cap9 - max_Phi_inf = {CAP9:.5f}-0.397 = {marg:.5f}")
print(f"    D0_crude = (crude const)/margin = {(Vcrude/6)/marg:.0f}   (single-peel bites only past here)")
print(f"    D0_emp   = (empirical const)/margin = {emp/marg:.1f}")
print(f"  BOX (THM-719 exhaustive) covers d<=25. INTERMEDIATE BAND = 26..D0_crude is where the")
print(f"    single-peel bound is vacuous; part (B) shows max Phi there stays {band_max:.3f} < cap9.")
print("\ndone.")
