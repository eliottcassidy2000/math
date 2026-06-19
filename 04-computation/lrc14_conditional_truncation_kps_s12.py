#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) CONDITIONAL TRUNCATION TEST (kps-S12 part 3)

Part 2 showed: K(n)=D7(n mod7)/prod(n_j) EXACTLY, but the inf-norm-truncated
sum over support-6 relations PLATEAUS at ~0.02 (only ~2-8% of the true
correction 0.30).  The sum sum_{n} D7(n mod7)/prod(n_j) is only CONDITIONALLY
convergent, so the inf-norm box is the wrong summation order.

THIS SCRIPT tests three things:
 (1) Is the correction really NOT recovered by support-6 alone?  -> compute the
     EXACT contribution of ALL supports (the identity is over the FULL lattice,
     and THM-538 says support<6 -> K=0, so support-6,7,...,k all contribute;
     but here k=8 so supports 6,7,8).  We add support-7 and support-8.
 (2) Does pairing n with -n (real part) + a SYMMETRIC product-weight order
     converge faster?
 (3) Does a "peel the largest coordinate" Abel/partial-summation order (the
     HYP-2610 stranger contraction) converge to the target?

Conclusion target: identify the CORRECT summation order under which the signed
coset sum converges, and measure its rate.
"""
import sys, itertools, math, cmath
from fractions import Fraction
from collections import defaultdict
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def dist_p0(E):
    E=sorted(set(E)); bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p0=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if sum(1 for j in range(1,7) if j not in hit)==0: p0+=(hi-lo)
    return p0

def chat(T,m):
    if m==0: return 1.0-len(T)/7.0
    if m%7==0: return 0.0
    s=0+0j
    for j in T:
        a=j/7.0;b=(j+1)/7.0
        s+=(cmath.exp(-2j*math.pi*m*a)-cmath.exp(-2j*math.pi*m*b))/(2j*math.pi*m)
    return -s

# Factored K (THM-538 verified): K(n) = (1/prod_{nonzero-non7} n_j) * D7(those residues),
# but ONLY when support>=6 (else 0). For support exactly s, the nonzero non-7 coords are the support.
def zeta(p): return cmath.exp(-2j*math.pi*p/7.0)
def A(r): return (1-cmath.exp(-2j*math.pi*r/7.0))/(2j*math.pi)
def h_T(T,r): return -A(r)*sum(zeta(r*j) for j in T)
_D7cache={}
def D7(res):
    res=tuple(res)
    c=_D7cache.get(res)
    if c is not None: return c
    tot=0+0j
    for cnt in range(7):
        for T in itertools.combinations(range(1,7),cnt):
            p=1+0j
            for r in res: p*=h_T(T,r)
            tot+=((-1)**cnt)*p
    _D7cache[res]=tot
    return tot

def K_full(nvec):
    eff=[v for v in nvec if v!=0 and v%7!=0]
    if len(eff)<6:   # THM-538 floor
        return 0+0j
    ip=1.0
    for v in eff: ip/=v
    return ip*D7(tuple(v%7 for v in eff))

def all_relations(E,L):
    """ALL n in Z^k (any support), sum n_j e_j=0, |n|inf<=L, n!=0, computed by
       enumerating k-1 free coords and solving the dependent (largest |e|)."""
    k=len(E); out=[]
    dep=max(range(k),key=lambda i:abs(E[i])); e_dep=E[dep]
    free=[i for i in range(k) if i!=dep]; efree=[E[i] for i in free]
    ranges=[range(-L,L+1) for _ in range(k-1)]
    for vf in itertools.product(*ranges):
        s=sum(c*e for c,e in zip(vf,efree))
        if e_dep!=0:
            if s%e_dep!=0: continue
            vd=-s//e_dep
            if abs(vd)>L: continue
            nvec=[0]*k
            for i,c in zip(free,vf): nvec[i]=c
            nvec[dep]=vd
            if all(x==0 for x in nvec): continue
            out.append(tuple(nvec))
        else:
            if s!=0: continue
            for vd in range(-L,L+1):
                nvec=[0]*k
                for i,c in zip(free,vf): nvec[i]=c
                nvec[dep]=vd
                if all(x==0 for x in nvec): continue
                out.append(tuple(nvec))
    return out

if __name__=="__main__":
    print("LRC(14) CONDITIONAL TRUNCATION TEST (kps-S12 part 3)\n")
    for E in ([0,1,2,3,4,5,6,7],[0,1,3,5,7,9,11,13]):
        k=len(E)
        p0=dist_p0(E)
        M7=sum(((-1)**t)*math.comb(6,t)*((1-t/7.0)**(k-1)) for t in range(7))
        target=float(p0)-M7
        print(f"=== E={E}  target p0-M7 = {target:.6f} ===")
        for L in [3,4,5]:
            rels=all_relations(E,L)
            # split by support (count nonzero non-7 coords)
            by_supp=defaultdict(complex)
            total=0+0j
            for nvec in rels:
                Kv=K_full(nvec)
                supp=sum(1 for v in nvec if v!=0 and v%7!=0)
                by_supp[supp]+=Kv
                total+=Kv
            supp_str=" ".join(f"s{s}={by_supp[s].real:+.5f}" for s in sorted(by_supp) if abs(by_supp[s])>1e-9)
            print(f"  L={L}: #rels={len(rels):8d} total corr={total.real:+.6f} frac={total.real/target:.4f}  [{supp_str}]")
        print()
