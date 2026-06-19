#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) COSET-QUOTIENT DECAY  (kind-pasteur-2026-06-19-S12, part 2)

S12 part 1 (lrc14_signed_shell_decay) showed: grouping K(n) by |n|_inf SHELLS does
NOT give signed decay -- the per-shell signed sum oscillates and the inf-norm-<=7
partial sum recovers only ~8% of the correction.  The dominant mass lives at LARGE
coordinates (reciprocal 1/prod n_j tail).  So inf-norm is the wrong ruler.

HYP-2632/2633 says the RIGHT layer is the mod-7 CHARACTER packet:
   K(n) depends on n only through (residues mod 7) up to a reciprocal envelope.
THM-538/HYP-2632 finite kernel:  K(n) = (envelope from |n|) * (finite char value depending on n mod 7).
Concretely (THM-538 LEMMA B): chat_T(n) = ŝ_T(n) with |ŝ_T(n)| = |sin(pi n/7)|/(pi |n|),
and the PHASE / sign is a function of (n mod 7).  So
   K(n) = [PURE mod-7 character factor]  *  prod_j 1/n_j-type reciprocal magnitude (with signs from n mod 7).

THIS SCRIPT tests the precise factorization:
   K(n) = C7(n mod 7) * R(n),     where R(n) = prod_j w(n_j),  w(m)= sin(pi m/7)/(pi m) (signed, real).
i.e. does K(n) split EXACTLY as a finite mod-7 character coefficient times a product of
1-D reciprocal weights?  If YES, then
   correction = sum_{n in Lambda(E)} C7(n mod 7) R(n)
              = sum_{cosets c mod 7} C7(c) * [ sum_{n in Lambda(E), n=c mod 7} R(n) ].
The inner sum is an ABSOLUTELY CONVERGENT reciprocal lattice sum (prod 1/n_j over a
codim>=? sublattice) -- and the OUTER finite sum carries the cancellation (the -108U+54U).
We verify the factorization, then compute the coset reciprocal sums and the resulting
signed quotient bound.
"""
import sys, itertools, math
from fractions import Fraction
from functools import reduce

sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---- exact p0 engine ----
def dist_p0(E):
    E = sorted(set(E))
    bps = set([Fraction(0), Fraction(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(Fraction(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p0 = Fraction(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*mid; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        N = sum(1 for j in range(1,7) if j not in hit)
        if N == 0: p0 += (hi-lo)
    return p0

# ---- 1-D real reciprocal weight w(m) = sin(pi m /7)/(pi m) (the |shat_T| magnitude, signed) ----
def w(m):
    return math.sin(math.pi*m/7.0)/(math.pi*m)

# ---- The PURE mod-7 finite character factor C7(residue tuple) ----
# Claim to test: chat_T(m) (m!=0, 7 not| m) = w(m) * phase7_T(m mod 7), where phase7 is a
# UNIMODULAR-ish factor depending only on (m mod 7) and T.  Then
#   K(n) = sum_T (-1)^|T| prod_j chat_T(n_j)
#        = sum_T (-1)^|T| prod_j [ w(n_j) * g_T(n_j mod 7) ]
#        = (prod_j w(n_j)) * sum_T (-1)^|T| prod_j g_T(n_j mod 7)
#        = R(n) * C7(n mod 7).
# We must extract g_T(r) for r in F_7^*.  From chat_T(m) = -sum_{j in T} cbar(j,m), and
#   cbar(j,m) = (e^{-2pi i m j/7} - e^{-2pi i m (j+1)/7})/(2 pi i m)
#             = e^{-2pi i m j/7}(1 - e^{-2pi i m/7})/(2 pi i m).
# So chat_T(m) = -[(1-e^{-2pi i m/7})/(2 pi i m)] * sum_{j in T} e^{-2pi i m j/7}.
# The factor (1-e^{-2pi i m/7})/(2 pi i m): magnitude = |sin(pi m/7)|/(pi|m|)=|w(m)|, phase depends on m mod 7.
# Define for r=m mod 7 (r in 1..6):
#   A(r) = (1 - exp(-2pi i r/7)) / (2 pi i)   [the "/m" pulled out as 1/m]
# then chat_T(m) = -(A(r)/m) * sum_{j in T} exp(-2pi i r j /7).
# So chat_T(m) = (1/m) * [ -A(r) sum_{j in T} zeta^{ -r j} ]  =: (1/m) * h_T(r).
# Thus prod_j chat_T(n_j) = (1/ prod n_j) prod_j h_T(n_j mod 7).
# K(n) = (1/prod n_j) * sum_T (-1)^|T| prod_j h_T(n_j mod 7)  =: (1/prod n_j) * D7(n mod 7).
# This is the EXACT factorization: a pure 1/prod(n_j) reciprocal times a finite mod-7 coefficient D7.
import cmath
def zeta(p):  # exp(-2 pi i / 7)^p
    return cmath.exp(-2j*math.pi*p/7.0)
def A(r):
    return (1 - cmath.exp(-2j*math.pi*r/7.0))/(2j*math.pi)
def h_T(T, r):
    return -A(r)*sum(zeta(r*j) for j in T)

def D7(residues):
    """finite mod-7 coefficient D7(c) = sum_T (-1)^|T| prod_j h_T(c_j),  c_j in 1..6."""
    total = 0+0j
    for cnt in range(7):
        for T in itertools.combinations(range(1,7), cnt):
            prod = 1+0j
            for r in residues:
                prod *= h_T(T, r)
            total += ((-1)**cnt)*prod
    return total

# Direct chat for verification
def chat_T_direct(T, m):
    if m == 0: return 1.0-len(T)/7.0
    if m % 7 == 0: return 0.0
    s = 0+0j
    for j in T:
        a=j/7.0;b=(j+1)/7.0
        s += (cmath.exp(-2j*math.pi*m*a)-cmath.exp(-2j*math.pi*m*b))/(2j*math.pi*m)
    return -s

def K_direct(vals):
    total=0+0j
    for cnt in range(7):
        for T in itertools.combinations(range(1,7),cnt):
            prod=1+0j
            for v in vals: prod*=chat_T_direct(T,v)
            total+=((-1)**cnt)*prod
    return total

def K_factored(vals):
    res = tuple(v%7 for v in vals)
    invprod = 1.0
    for v in vals: invprod /= v
    return invprod * D7(res)

# ---- relation lattice support-6 enumeration (solve dependent coord) ----
def relations_support6(E, L):
    k=len(E); out=[]
    for idxs in itertools.combinations(range(k),6):
        es=[E[i] for i in idxs]
        dep=max(range(6),key=lambda i:abs(es[i])); e_dep=es[dep]
        free=[i for i in range(6) if i!=dep]; efree=[es[i] for i in free]
        ranges=[range(-L,L+1) for _ in range(5)]
        if e_dep==0:
            for vf in itertools.product(*ranges):
                if any(c==0 or c%7==0 for c in vf): continue
                if sum(c*e for c,e in zip(vf,efree))!=0: continue
                for vd in range(-L,L+1):
                    if vd==0 or vd%7==0: continue
                    combo=[0]*6
                    for i,c in zip(free,vf): combo[i]=c
                    combo[dep]=vd
                    out.append((tuple(combo),max(abs(x) for x in combo)))
            continue
        for vf in itertools.product(*ranges):
            if any(c==0 or c%7==0 for c in vf): continue
            s=sum(c*e for c,e in zip(vf,efree))
            if s%e_dep!=0: continue
            vd=-s//e_dep
            if vd==0 or vd%7==0 or abs(vd)>L: continue
            combo=[0]*6
            for i,c in zip(free,vf): combo[i]=c
            combo[dep]=vd
            out.append((tuple(combo),max(abs(x) for x in combo)))
    return out

if __name__=="__main__":
    print("LRC(14) COSET-QUOTIENT DECAY (kps-S12 part 2)")
    print("Step A: VERIFY exact factorization  K(n) = (1/prod n_j) * D7(n mod 7)\n")
    import random
    random.seed(1)
    maxerr=0.0
    for _ in range(2000):
        vals=[]
        while len(vals)<6:
            v=random.randint(-20,20)
            if v!=0 and v%7!=0: vals.append(v)
        kd=K_direct(vals); kf=K_factored(vals)
        maxerr=max(maxerr,abs(kd-kf))
    print(f"  max |K_direct - K_factored| over 2000 random support-6 n = {maxerr:.3e}")
    print("  => factorization EXACT iff this is ~1e-13.\n")

    print("Step B: per-COSET signed reciprocal sums.  correction = sum_c D7(c) * S_c(E,L),")
    print("  S_c(E,L) = sum_{n in Lambda(E), n=c mod7, |n|inf<=L} 1/prod n_j.")
    print("  We track how S_c (and the total) converge per coset.\n")

    from collections import defaultdict
    for E in ([0,1,2,3,4,5,6,7],[0,1,3,5,7,9,11,13]):
        p0=dist_p0(E)
        M7=sum(((-1)**t)*math.comb(6,t)*((1-t/7.0)**(len(E)-1)) for t in range(7))
        target=float(p0)-M7
        print(f"--- E={E}  target p0-M7 = {target:.6f} ---")
        for L in [5,7,9]:
            rels=relations_support6(E,L)
            coset_recip=defaultdict(float)   # c (residue tuple) -> sum 1/prod n
            for combo,ninf in rels:
                c=tuple(v%7 for v in combo)
                ip=1.0
                for v in combo: ip/=v
                coset_recip[c]+=ip
            corr=0+0j
            for c,S in coset_recip.items():
                corr += D7(c)*S
            print(f"   L={L}: #rels={len(rels)}  reconstructed corr = {corr.real:+.6f} (Im={corr.imag:.1e})  frac={corr.real/target:.4f}")
        print()
