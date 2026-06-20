#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) REFRAMING F: C_3 (cube-root, order-3) PARTIAL SYMMETRIZATION of the
relation-lattice correction sum_{n in Lambda(E)} K(n), K(n)=w(n) D7(c), c=n mod 7.

CENTRAL DISTINCTION (the task's premise needs correcting):
  There are TWO different "C_3" structures, and they are NOT the same ring:

  (G) GALOIS C_3 = the order-3 subgroup <sigma_2> = {sigma_1,sigma_2,sigma_4} of
      Gal(Q(zeta_7)/Q) (cyclic order 6).  sigma_a: zeta_7 -> zeta_7^a.
      Its FIXED FIELD is the QUADRATIC subfield Q(sqrt(-7)) (since {1,2,4}=QR).
      => the partial trace  T3(c) := sum_{a in {1,2,4}} sigma_a(D7(c))
         lies in Z[(1+sqrt(-7))/2] = O_{Q(sqrt-7)}, NOT in Z[omega] (Eisenstein).
      The FULL trace (all 6 a) lands in Z (rational): that is HYP-2657.
      The C_3 partial trace is the HALFWAY object: rational over Q(sqrt-7).

  (S) SECTOR C_3 = mult-by-2 permuting the 7 sectors as (0)(1 2 4)(3 5 6).
      This is the dilation of the *coordinate residues* c -> 2c, which by Galois
      equivariance D7(2c)=sigma_2(D7(c)) is the SAME map as (G) applied to D7.
      So (S) and (G) coincide on D7.  Good: only ONE C_3 acts.

  The task asks whether C_3-partial-trace + Z/2-reality (n<->-n) gives a SIGNED
  bound BETTER than absolute (closer to the true ~5x cancellation in |corr|).

We compute, on the actual relation lattices of AP8 and WIDE8:
  (1) Confirm T3(c) in Q(sqrt-7): test sigma_2-fixed.
  (2) The C_3 partial-trace magnitudes |T3(c)| vs absolute sum sum|D7| over the orbit.
  (3) Re-bundle the correction by C_3 dilation orbits and measure the partial-sum
      envelope vs the per-residue absolute bound.  Honest gain factor.
  (4) Compare three bounds on |correction|:
        B_abs  = sum_n |w(n)| |D7(c_n)|                 (absolute, the 5x-lossy one)
        B_re   = sum_n |w(n)| |Re D7(c_n)|              (Z/2 reality only)
        B_c3   = sum_{orbits O} |sum_{n: c in O} w(n) D7(c_n)|  (C_3 bundling, signed)
      and the EXACT |correction| as ground truth.
"""
import sys, itertools, cmath, math
from collections import defaultdict
from fractions import Fraction
sys.path.insert(0, "04-computation")

# ---------- exact Z[zeta_7] arithmetic (basis 1,z,...,z^5; z^6 = -1-...-z^5) ----------
def zexp(e):
    e %= 7
    if e <= 5:
        v=[Fraction(0)]*6; v[e]=Fraction(1); return tuple(v)
    return (Fraction(-1),)*6
def zadd(a,b): return tuple(x+y for x,y in zip(a,b))
def zsub(a,b): return tuple(x-y for x,y in zip(a,b))
def zscale(a,s): return tuple(x*s for x in a)
def zmul(a,b):
    out=[Fraction(0)]*7
    for i in range(6):
        if a[i]==0: continue
        for j in range(6):
            if b[j]==0: continue
            out[(i+j)%7]+=a[i]*b[j]
    c6=out[6]
    return tuple(out[i]-c6 for i in range(6))
ZERO=(Fraction(0),)*6
Z=cmath.exp(2j*math.pi/7)
def znum(a): return sum(complex(a[k])*(Z**k) for k in range(6))
def galois(a,x):
    r=ZERO
    for k in range(6):
        if x[k]!=0: r=zadd(r, zscale(zexp((a*k)%7), x[k]))
    return r

Tlist=[T for r in range(7) for T in itertools.combinations(range(1,7),r)]
SGN={T:Fraction((-1)**len(T)) for T in Tlist}
def sigma_T(T,m):
    r=ZERO
    for t in T: r=zadd(r, zexp((-m*t)%7))
    return r
SIG={(T,m):sigma_T(T,m) for T in Tlist for m in range(1,7)}
PREF={m:zsub(zexp(0), zexp((-m)%7)) for m in range(1,7)}
_Dc={}
def D7(c):
    v=_Dc.get(c)
    if v is not None: return v
    pref=zexp(0)
    for cj in c: pref=zmul(pref,PREF[cj])
    acc=ZERO
    for T in Tlist:
        p=zexp(0)
        for cj in c: p=zmul(p,SIG[(T,cj)])
        acc=zadd(acc, zscale(p,SGN[T]))
    v=zmul(pref,acc); _Dc[c]=v; return v

C3 = [1,2,4]    # the order-3 Galois subgroup = QR
FULL = [1,2,3,4,5,6]

def T3(c):
    """C_3 partial trace: sum_{a in {1,2,4}} sigma_a(D7(c)) = sum_a D7(a*c)."""
    r=ZERO
    for a in C3:
        r=zadd(r, D7(tuple((a*cj)%7 for cj in c)))
    return r

def Tfull(c):
    r=ZERO
    for a in FULL:
        r=zadd(r, D7(tuple((a*cj)%7 for cj in c)))
    return r

def in_quad_field(x):
    """x in Q(sqrt-7) iff fixed by sigma_2 (generator of <2>=QR)."""
    return galois(2,x)==x
def is_rational(x):
    return all(x[k]==0 for k in range(1,6))

# ---------- relation lattice (reused from lrc14_qr_lattice_apply_kps) ----------
def relations_support6(E,L):
    k=len(E); out=[]
    for idxs in itertools.combinations(range(k),6):
        es=[E[i] for i in idxs]
        dep=max(range(6),key=lambda i:abs(es[i])); e_dep=es[dep]
        if e_dep==0: continue
        free=[i for i in range(6) if i!=dep]; efree=[es[i] for i in free]
        for vfree in itertools.product(range(-L,L+1),repeat=5):
            if any(c==0 or c%7==0 for c in vfree): continue
            s=sum(c*e for c,e in zip(vfree,efree))
            if s%e_dep!=0: continue
            vd=-s//e_dep
            if vd==0 or vd%7==0 or abs(vd)>L: continue
            combo=[0]*6
            for i,c in zip(free,vfree): combo[i]=c
            combo[dep]=vd
            out.append(tuple(combo))
    return out

def w_real(vals):
    p=1.0
    for v in vals: p*=1.0/v
    return -p/((2*math.pi)**6)   # i^6=-1, support 6

if __name__=="__main__":
    print("="*72)
    print("PART 1: structural — where does the C_3 partial trace land?")
    print("="*72)
    import random
    random.seed(1)
    samp=[tuple(random.choice([1,2,3,4,5,6]) for _ in range(6)) for _ in range(60)]
    samp+=[(1,1,1,1,1,1),(1,2,4,3,5,6),(3,3,3,3,3,3),(1,2,3,4,5,6)]
    nq=sum(in_quad_field(T3(c)) for c in samp)
    nr_part=sum(is_rational(T3(c)) for c in samp)
    nr_full=sum(is_rational(Tfull(c)) for c in samp)
    print(f"  T3 (partial, a in {{1,2,4}}) in Q(sqrt-7)  [fixed by sigma_2]: {nq}/{len(samp)}")
    print(f"  T3 (partial) actually RATIONAL                              : {nr_part}/{len(samp)}")
    print(f"  Tfull (all 6 a) RATIONAL  (=HYP-2657 full trace collapse)   : {nr_full}/{len(samp)}")
    print("  => partial C_3 trace lives in Q(sqrt-7), NOT generically rational, NOT Z[omega].")

    # magnitudes: |T3(c)| vs absolute sum over the orbit
    print("\n  per-orbit magnitude: |T3| vs sum|D7| over the C_3 orbit (cancellation factor):")
    seen=set(); rows=[]
    for c in samp:
        orb=tuple(sorted(tuple((a*cj)%7 for cj in c) for a in C3))
        key=orb[0]
        if key in seen: continue
        seen.add(key)
        members=set(tuple((a*cj)%7 for cj in c) for a in C3)
        absum=sum(abs(znum(D7(m))) for m in members)
        t3mag=abs(znum(T3(c)))
        if absum>1e-12:
            rows.append((t3mag/absum, t3mag, absum, len(members)))
    rows.sort()
    import statistics
    ratios=[r[0] for r in rows]
    print(f"    orbits sampled={len(rows)}  |T3|/sum|D7| : min={min(ratios):.3f} "
          f"mean={statistics.mean(ratios):.3f} max={max(ratios):.3f}")
    print(f"    (ratio<1 => signed C_3 cancellation; ratio≈1 => no gain at orbit level)")

    print("\n"+"="*72)
    print("PART 2: bounds on the actual relation-lattice correction")
    print("="*72)
    for E,lab in [(list(range(8)),"AP8 consec"),([0,1,3,5,7,9,11,13],"WIDE8 odd-AP")]:
        print(f"\n=== {lab}  E={E} ===")
        for L in [3,4]:
            rels=relations_support6(E,L)
            # exact correction in Z[zeta_7]
            corr=ZERO
            byres=defaultdict(float)    # W_c real weight per residue
            B_abs=0.0; B_re=0.0
            for n in rels:
                c=tuple(v%7 for v in n)
                wv=w_real(n)
                corr=zadd(corr, zscale(D7(c), Fraction(wv).limit_denominator(10**12)))
                byres[c]+=wv
                d=znum(D7(c))
                B_abs+=abs(wv)*abs(d)
                B_re +=abs(wv)*abs(d.real)
            corr_c=znum(corr)
            exact=abs(corr_c)
            # C_3 bundling: group residues into C_3 dilation orbits, sum signed inside orbit
            orbits=defaultdict(lambda: ZERO)
            for c,W in byres.items():
                key=min(tuple((a*cj)%7 for cj in c) for a in C3)
                orbits[key]=zadd(orbits[key], zscale(D7(c), Fraction(W).limit_denominator(10**12)))
            B_c3=sum(abs(znum(v)) for v in orbits.values())
            # also full-F7* bundling (for reference)
            orbits6=defaultdict(lambda: ZERO)
            for c,W in byres.items():
                key=min(tuple((a*cj)%7 for cj in c) for a in FULL)
                orbits6[key]=zadd(orbits6[key], zscale(D7(c), Fraction(W).limit_denominator(10**12)))
            B_full=sum(abs(znum(v)) for v in orbits6.values())
            print(f"  L={L}: #rel={len(rels):6d}  #res={len(byres):4d}  exact|corr|={exact:.4e}")
            print(f"        B_abs ={B_abs:.4e} (x{B_abs/exact:6.2f})   "
                  f"B_re ={B_re:.4e} (x{B_re/exact:6.2f})")
            print(f"        B_c3  ={B_c3:.4e} (x{B_c3/exact:6.2f})   "
                  f"B_full={B_full:.4e} (x{B_full/exact:6.2f})")
            print(f"        C_3 gain over absolute: {B_abs/B_c3:.2f}x   "
                  f"full-trace gain over absolute: {B_abs/B_full:.2f}x")
