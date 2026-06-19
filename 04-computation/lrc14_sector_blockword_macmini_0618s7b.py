#!/usr/bin/env python3
"""
lrc14_sector_blockword_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B continued)

Per-T extremality HOLDS (s7b part B). Now CHARACTERIZE meas(B_T) and prove AP-extremality
of the IE-signed term, per block T.

KEY OBJECT.  For a sector-block T subset Z/7 with |T|=t, define the "keep" indicator on the
circle: K_T(s) = 1 if s notin T (s in Z/7). For offset e>=1, sigma_e(x)=floor(7ex) mod 7.
  B_T(E) = {x: for all e, sigma_e(x) notin T},  meas = integral over [0,1) of prod_e K_T(sigma_e(x)).

FOURIER ON Z/7.  K_T(s) = sum_{n=0}^{6} c_n omega^{n s}, omega=e^{2pi i/7}, c_n = (1/7) sum_{s notin T} omega^{-n s}.
c_0 = (7-t)/7.  And sigma_e(x)=floor(7ex) mod 7, so omega^{n sigma_e(x)} is NOT simply e^{2pi i n e x}
(floor vs the linear phase). The cutting-sequence Fourier needs the SAWTOOTH expansion:
  1_{sigma_e(x)=s} = meas-cell indicator; its Fourier in x has freq multiples of 7e.
Actually the clean statement: as x runs [0,1), sigma_e(x) takes each value 0..6 on e disjoint
intervals each of total length 1/7, ARRANGED as the periodic word (0123456)^e on 7e equal cells.
So prod_e K_T(sigma_e(x)) is a 0/1 word on the common refinement grid; meas(B_T)= (#keep cells)/(lcm grid).

COMBINATORIAL FORM. On the grid of N=lcm_e(7e) cells (cell j = [j/N,(j+1)/N)),
  sigma_e on cell j = floor(7 e (j+0.5)/N) mod 7 = floor( (7e/N)(j+0.5) ) mod 7.
Let N = 7*L where L=lcm(E_nz). Then 7e/N = e/L, and cell j keeps iff for all e: floor(e(j+0.5)/L)*? ...
Cleaner: parametrize y=7x in [0,7). sigma_e(x)=floor(e*7x) mod7 = floor(e*y) mod 7. With y in[0,7).
meas_x(B_T)=(1/7) meas_y{y in[0,7): for all e, floor(e y) mod 7 notin T}.

So define g_T(E) = meas{y in [0,7): for all e in E, (floor(e y) mod 7) notin T}, meas(B_T)=g_T/7.
This is a UNION of arithmetic progressions of intervals problem -- pure combinatorics on words.

THIS SCRIPT:
 (1) Confirm meas(B_T) depends on T only through its CIRCULAR GAP STRUCTURE (rotation-invariant).
 (2) For the per-T extremal, characterize: for fixed t=|T|, which arrangement of T (consecutive
     block vs spread) gives the LARGEST meas(B_T)? And does AP-of-E maximize over E?
 (3) Decompose the t=1 case fully (single missed residue) -- the cleanest. meas(B_{single})
     and prove AP extremal there as a warmup (it is the (-1)^1 term, want MIN).
 (4) Express corr's leading triple term as a relation-lattice sum and confirm 7|n vanishing.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def lcm(a,b): return a*b//gcd(a,b)
def lcml(xs): return reduce(lcm, xs, 1)

def measB_T(E, T):
    E = sorted(set(e for e in E if e!=0)); Tset=set(T)
    if not E: return F(1) if 0 not in Tset else F(0)
    bps=set([F(0),F(1)])
    for e in E:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; ok=True
        for e in E:
            if (int(7*e*xm)%7) in Tset: ok=False; break
        if ok: total+=x1-x0
    return total

def Tclass(T):
    T=sorted(T); n=len(T)
    if n==0: return (0,())
    best=None
    for sh in range(7):
        rot=tuple(sorted((t+sh)%7 for t in T))
        if best is None or rot<best: best=rot
    # represent by gap structure (cyclic gaps between consecutive elements)
    el=sorted(best); gaps=[]
    for i in range(n):
        gaps.append((el[(i+1)%n]-el[i])%7)
    # canonical rotation of gap cyclic word
    gw=min(tuple(gaps[i:]+gaps[:i]) for i in range(n))
    return (n, gw)

# ---------- (1) meas(B_T) is rotation-invariant; depends only on circular gap structure ----------
print("="*92)
print("(1) meas(B_T) depends on T ONLY through its circular gap structure (test on a few E)")
print("="*92)
Etests=[(0,1,2,3,4,5,6,7),(0,1,3,7,12,20,30,44),(0,2,3,5,7,11,13,17)]
for E in Etests:
    bycls={}
    ok=True
    for r in range(1,7):
        for T in itertools.combinations(range(7),r):
            c=Tclass(T); v=measB_T(E,T)
            if c in bycls and bycls[c]!=v: ok=False; print(f"   MISMATCH E={E} T={T} class={c}")
            bycls[c]=v
    print(f"  E={E}: rotation-class-invariance {'HOLDS' if ok else 'FAILS'}, #distinct classes={len(bycls)}")

# enumerate the T-classes for Z/7
print("\n  T-classes on Z/7 (by size and gap word):")
allcls=set()
for r in range(1,7):
    for T in itertools.combinations(range(7),r):
        allcls.add(Tclass(T))
for c in sorted(allcls):
    print(f"     size {c[0]}: gap word {c[1]}")

# ---------- (2) for fixed t=|T|, consecutive-block T vs spread-T: which maximizes meas(B_T)? ----------
print()
print("="*92)
print("(2) For fixed t, which T-shape maximizes meas(B_T)?  (on AP_8 and a random E)")
print("="*92)
for E in [(0,1,2,3,4,5,6,7),(0,1,3,7,12,20,30,44)]:
    print(f"  E={E}")
    for t in range(1,6):
        rows=[]
        seen=set()
        for T in itertools.combinations(range(7),t):
            c=Tclass(T)
            if c in seen: continue
            seen.add(c)
            rows.append((measB_T(E,T),c,T))
        rows.sort(reverse=True)
        cons=tuple(range(t))  # consecutive block
        consv=measB_T(E,cons)
        mx=rows[0]; mn=rows[-1]
        print(f"    t={t}: MAX meas(B_T)={float(mx[0]):.5f} at gapword {mx[1][1]}; "
              f"consec-block={float(consv):.5f}; MIN={float(mn[0]):.5f} at {mn[1][1]}")

# ---------- (3) the single-missed-residue (t=1) term, the cleanest. -----------
print()
print("="*92)
print("(3) t=1 single missed residue: meas(B_{r}) = meas{x: NO e has sigma_e(x)=r}.")
print("    IE term is -sum_r meas(B_r) (each r equivalent by rotation: 7 * meas(B_0class)).")
print("    Want this MINIMIZED (to maximize -term).  Is AP the minimizer over E?")
print("="*92)
def gen_shapes(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out
for k in [8]:
    AP=tuple(range(k)); shapes=gen_shapes(k,k+4)
    T=(1,)  # single residue (rotation rep)
    vals=[(measB_T(E,T),E) for E in shapes]
    mn=min(vals); ap=measB_T(AP,T)
    print(f"  k={k}: t=1 meas(B_T): AP={float(ap):.6f}  global MIN={float(mn[0]):.6f} at {mn[1]}  "
          f"{'AP IS MIN' if ap==mn[0] else 'AP NOT MIN'}")

# ---------- (4) corr leading triple term & 7|n vanishing (relation lattice) -----------
print()
print("="*92)
print("(4) Relation-lattice check: corr(E)=meas(S7)-M7(k). For dissociated E, corr->0.")
print("="*92)
def measS7(E):
    s=F(0)
    for r in range(0,8):
        for M in itertools.combinations(range(7),r):
            s+=F((-1)**r)*measB_T(E,M)
    return s
def M7(k):
    s=F(0)
    for t in range(0,7): s+=F((-1)**t*comb(6,t))*F(7-t,7)**(k-1)
    return s
for E in [(0,1,2,3,4,5,6,7),(0,1,2,5,6,7,10,11),(0,1,3,7,15,31,63,127),(0,1,3,7,12,20,30,44)]:
    k=len(E); s7=measS7(E); corr=s7-M7(k)
    print(f"  E={E}: meas(S7)={float(s7):.6f} corr={float(corr):.6f}")
print("\nDONE.")
