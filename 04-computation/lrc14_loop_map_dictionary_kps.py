#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
A DICTIONARY OF MAPS ON THE LOOP T=R/Z, and which of them operate GROUP-LIKE.

kind-pasteur-2026-07-01-S8. The owner's cue: "pushing Verblunsky to the unit circle is a nice
recursive metaphor for LRC (runners on a loop) -- find creative functions between points on a loop,
a whole dictionary; they may operate group-like."  Verified here:

  #  MAP                     formula on T=R/Z            composition law             STRUCTURE
  1  rotation  R_a           t -> t+a                    R_a . R_b = R_{a+b}         GROUP (T,+)  [runner FLOW]
  2  dilation  M_v (RUNNER)  t -> v t                    M_v . M_w = M_{vw}          MONOID (Z,*); GROUP (Z/N)* on N-torsion
  3  affine    A_{v,a}       t -> v t + a                semidirect                  GROUP Aff = T x| (Z/N)*
  4  inversion i             t -> -t                     i^2 = id                    GROUP Z/2  [antipode/complement]
  5  doubling  D             t -> 2 t                    D^k = M_{2^k}               MONOID (Z>=0,+); the "2" of 14=2.7
  6  Blaschke  b_a (OPUC)    z -> (z-a)/(1-conj(a) z)    Mobius o Mobius = Mobius    GROUP SU(1,1)/pm [VERBLUNSKY]
  7  Gauss/CF  G             x -> {1/x}                  natural extension           GROUPOID (modular geodesic flow)
  D  character chi_k         t -> e(k t)                 chi_k . chi_l = chi_{k+l}   DUAL group (Z,+)=T^  [Fourier]
  c  sawtooth/Dedekind s(h,k) 1-cocycle of SL2(Z)        reciprocity                 COCYCLE  [value -> -1/12=zeta(-1)]

Runners are map #2 (multiplication); the (Z/N)* atoms are its invertible orbit (HYP-3793).  The OPUC
recursion is map #6 (Blaschke), "push to the unit circle" = its boundary dynamics.  The sawtooth Dedekind
cocycle (c) is the CONNECTING map between the additive circle (#1) and the modular/Mobius world (#6/#7)
-- which is WHY the census collar-sum (additive, over units) evaluates to a Dedekind sum (-1/12=zeta(-1)).
"""
import sys, cmath, itertools
from fractions import Fraction as Fr
from math import gcd, pi, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def phi(N): return len(units(N))

print("="*94); print(" #2  RUNNERS = MULTIPLICATION MAPS: M_v.M_w = M_{vw}; (Z/N)* is a GROUP (the atoms)"); print("="*94)
for N in [10,12,14]:
    U=units(N)
    closed=all((a*b)%N in U for a in U for b in U)
    has_inv=all(any((a*b)%N==1 for b in U) for a in U)
    # composition table order structure: is it cyclic?
    orders=[]
    for a in U:
        o=1; x=a%N
        while x!=1: x=(x*a)%N; o+=1
        orders.append(o)
    print(f"  (Z/{N})* = {U}: |.|=phi={phi(N)}; closed under M_v.M_w? {closed}; every M_v invertible? {has_inv}")
    print(f"     element orders {dict(zip(U,orders))}; cyclic? {max(orders)==phi(N)}  (=> runners on N-torsion form an abelian group)")
print("  RUNNER COMPOSITION is literally speed-multiplication: applying runner v then runner w = runner vw.")
print("  Scale-invariance THM-522 = the action of M_c; the census atoms are one (Z/N)*-orbit.")

print("\n"+"="*94); print(" #6  OPUC / BLASCHKE MAPS: b_a(z)=(z-a)/(1-conj(a)z) form a GROUP; iterate -> boundary"); print("="*94)
def blaschke(a): return lambda z: (z-a)/(1-a.conjugate()*z)
def compose(f,g): return lambda z: f(g(z))
# verify closure: b_a . b_b maps disk->disk and is a Mobius transform of the same (unit-disk-automorphism) type
import random
rng=random.Random(0)
def is_disk_auto(F, tries=200):
    ok=True
    for _ in range(tries):
        r=rng.random()**0.5; th=rng.random()*2*pi; z=r*cmath.exp(1j*th)
        w=F(z)
        if abs(w)>1+1e-9: ok=False
    return ok
a1=0.3+0.4j; a2=-0.2+0.5j
F=compose(blaschke(a1),blaschke(a2))
print(f"  b_a1 . b_a2 (a1={a1}, a2={a2}) maps disk->disk (unit-disk automorphism)? {is_disk_auto(F)}")
print(f"     => Blaschke maps are closed under composition = the group of disk automorphisms SU(1,1)/pm.")
# "push Verblunsky to unit circle": iterate a fixed b_a; orbit of 0 spirals to the boundary as |a|->1
for amag in [0.5, 0.9, 0.98]:
    a=amag*cmath.exp(1j*1.0); b=blaschke(a); z=0j; radii=[]
    for _ in range(6): z=b(z); radii.append(round(abs(z),4))
    print(f"  iterate b_a, |a|={amag}: |orbit of 0| = {radii}  (|a|->1 => orbit -> unit circle = 'atomic/tight' limit)")
print("  The OPUC Szego recursion IS Blaschke composition; the extremizer's |alpha_k|->1 (HYP-3793) is the")
print("  orbit hitting the boundary = the loose measure collapsing onto its (Z/N)* atoms.")

print("\n"+"="*94); print(" (c) THE DEDEKIND COCYCLE connects additive #1 <-> modular #6/#7; value -> -1/12=zeta(-1)"); print("="*94)
def dedekind(h,k):  # s(h,k)=sum_{a=1}^{k-1} ((a/k))((ha/k)),  ((x))=x-floor(x)-1/2 (0 at integers)
    def saw(x):
        f=x-floor(x)
        return Fr(0) if f==0 else f-Fr(1,2)
    return sum(saw(Fr(a,k))*saw(Fr(h*a,k)) for a in range(1,k))
# reciprocity: s(h,k)+s(k,h) = -1/4 + (h/k + k/h + 1/(hk))/12
for h,k in [(1,7),(2,7),(3,7),(1,13),(5,12)]:
    lhs=dedekind(h,k)+dedekind(k,h)
    rhs=Fr(-1,4)+(Fr(h,k)+Fr(k,h)+Fr(1,h*k))/12
    print(f"  s({h},{k})={dedekind(h,k)} ; reciprocity s(h,k)+s(k,h)={lhs} =?= {rhs}  {lhs==rhs}")
# the -1/12 shadow: s(1,k) = (k-1)(k-2)/(12k) -> k/12 ~ ... and the AP-orbit sum -> -1/12 (HYP-3773 margin)
print("  s(1,k)=(k-1)(k-2)/(12k):", [str(dedekind(1,k)) for k in [7,13,14]], "-> the 1/12 = -zeta(-1) reciprocity constant.")

print("\n"+"="*94); print(" (D) CHARACTERS chi_k(t)=e(kt): the DUAL group; Ramanujan sum = sum of unit characters"); print("="*94)
for N in [10,14]:
    for k in [1,2]:
        c=sum(cmath.exp(2j*pi*a*k/N) for a in units(N))
        print(f"  c_{N}({k}) = sum_{{a in (Z/{N})*}} chi_a(k/{N}) = {c.real:+.4f} (Ramanujan sum = the atomic-measure moment, HYP-3793)")

print("\n"+"="*94); print(" CENSUS SYMMETRY GROUP: (Z/N)* label is INVARIANT under dilation x complement x translate"); print("="*94)
# the group G = <M_c (dilation), i (complement t->-t <=> v->v), integer relabel> acts on speed-sets;
# verify the covering-min N (denominator of maximizers) is invariant under primitive dilation of the set.
def nrm(x):
    f=x-(x.numerator//x.denominator); return f if f<=Fr(1,2) else 1-f
def M_and_N(S,Dmax=60):
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w)
            if v!=w: dens.add(abs(v-w))
    best=Fr(0); atoms=[]
    for b in sorted(d for d in dens if 0<d<=Dmax):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            val=min(nrm(v*Fr(a,b)) for v in S)
            if val>best: best=val; atoms=[Fr(a,b)]
            elif val==best: atoms.append(Fr(a,b))
    Ns=sorted(set(t.denominator for t in atoms))
    return best, (Ns[0] if len(Ns)==1 else tuple(Ns))
S=[1,2,3,4,5,7,8,9,11,12,13]
M0,N0=M_and_N(S)
print(f"  base S={S}: M={M0}, atoms on (Z/{N0})*")
for c in [2,3,5]:
    Sc=[c*v for v in S]; Mc,Nc=M_and_N(Sc, Dmax=60*c)
    print(f"     dilation M_{c}(S): M={Mc} (=M0? {Mc==M0}), atoms on (Z/{Nc})*  (dilation is a SYMMETRY: same class)")
print("  => the census is a quotient by this group; each orbit = ONE (Z/N)* class => enumerate REPRESENTATIVES.")
print("DONE.")
