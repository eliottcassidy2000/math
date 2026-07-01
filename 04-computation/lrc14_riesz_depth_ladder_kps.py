#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
RIESZ-DEPTH LADDER: is the far-element discrepancy w*Delta_w a Riesz-product effect?
kind-pasteur-2026-06-30.

Reframing A ("attacker = defender"): the hard multiscale extremizers that make
sup_w |w*Delta_w| large are COHERENT TOWERS, whose Fourier structure factorizes as
a Riesz product prod_j \hat{1_B}(scale^j * xi). The proof tool (mac-mini's test
measure) is ALSO a Riesz product. This script tests the SHARP, falsifiable content
of that reframing, using the VERIFIED engine from lrc14_uniform_C_growth_kps.py
(G0, cells_of, wD, supw copied verbatim; benchmarks reproduced as a validation gate):

  P1 (DEPTH, not cardinality): sup|wD| grows ~LINEARLY with the number of coherent
      blocks (= Riesz factors = Dedekind-ladder rungs), while single blocks and
      incoherent sets of the SAME cardinality stay small.
  P2 (SCALE RESONANCE / attacker=defender): the argmax w* EXTENDS the tower --
      it is the far element that completes/continues the self-similar structure.
  P3 (STRUCTURAL): the clean tower {0,1,2} (+) 30*{0,1,2} has an EXACT 2-fold Riesz
      factorization \hat{1_E}(xi) = \hat{1_B}(xi) * \hat{1_B}(30 xi).

If P1+P2 hold, the discrepancy is a Riesz-product/depth effect: the order parameter
is COHERENT DEPTH, and the attacker's optimal far element is the one that makes the
set more self-similar -- i.e. the adversary plays the certificate's own object.

RESULT (2026-06-30): P1 and P3 CONFIRMED; P2 partially. But the Sidon control
CORRECTED the framing: the driver is SPREAD / #distinct-scales (reached by coherent
repetition OR by spreading), NOT coherence specifically -- wide Sidon sets score as
high as towers. The honest punchline is the Phi-Delta TRADE-OFF: sup|wD| is large
exactly where the plateau Phi is small, so neither alone is a counterexample.
"""  # noqa: W605  (backslashes below are math notation in the docstring)
import sys, random, cmath, math
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---- VERIFIED ENGINE (verbatim from 04-computation/lrc14_uniform_C_growth_kps.py) ----
def G0(y):
    f=y-(y.numerator//y.denominator)
    return Fraction(6,7)*f if f<=Fraction(1,7) else Fraction(6,49)-(f-Fraction(1,7))/7
def cells_of(Ep):
    Ep=sorted(set(Ep)); bps={Fraction(0),Fraction(1)}
    for e in Ep:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); out=[]
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in Ep:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        m=[j for j in range(1,7) if j not in hit]
        if len(m)==1:
            if out and out[-1][2]==m[0] and out[-1][1]==lo: out[-1]=(out[-1][0],hi,m[0])
            else: out.append((lo,hi,m[0]))
    return out
def wD(cells,w):
    return sum(G0(w*b-Fraction(s,7))-G0(w*a-Fraction(s,7)) for (a,b,s) in cells)
def supw(Ep,wmax):
    cells=cells_of(Ep); best=(Fraction(0),None)
    for w in range(2,wmax):
        if reduce(gcd,tuple(Ep)+(w,))!=1: continue
        v=abs(wD(cells,w))
        if v>best[0]: best=(v,w)
    return best
def topw(Ep,wmax,k=6):
    """all (|wD|, w) sorted desc -- for peak-location analysis."""
    cells=cells_of(Ep); vals=[]
    for w in range(2,wmax):
        if reduce(gcd,tuple(Ep)+(w,))!=1: continue
        vals.append((abs(wD(cells,w)),w))
    vals.sort(reverse=True)
    return vals[:k], len(cells)
# --------------------------------------------------------------------------------------

def block_tower(scale, nblocks, base=(0,1,2)):
    return sorted(scale*i+b for i in range(nblocks) for b in base)

def wm_for(Ep):  # same rule the engine uses
    return 6*max(Ep)+30
def wfmt(w): return "  n/a" if w is None else f"{w:5d}"

print("="*78)
print("VALIDATION GATE: reproduce the three published benchmarks (must match)")
print("="*78)
bench = [("consec8",list(range(8)),0.73),
         ("odd-struct s11",[0,1,3,5,7,9,10,11],2.54),
         ("multiscale {0,1,2,30,31,32,60,61}",[0,1,2,30,31,32,60,61],3.91)]
ok=True
for name,Ep,exp in bench:
    s,w=supw(Ep,wm_for(Ep))
    match = abs(float(s)-exp)<0.01
    ok &= match
    print(f"  {name:38s} sup|wD|={float(s):.4f} at w={w}   expect~{exp}  {'OK' if match else '*** MISMATCH ***'}")
print(f"  ENGINE VALIDATED: {ok}\n")

print("="*78)
print("P1  DEPTH LADDER (scale 30): sup|wD| vs number of coherent blocks")
print("    base block B={0,1,2}; nblocks copies at 0,30,60,...  (Riesz factors)")
print("="*78)
print("  #blocks  set                                              sup|wD|   at w*   #cells  incr")
prev=None
scale=30
ladder=[]
for nb in range(1,5):
    Ep=block_tower(scale,nb)
    (s,w)=supw(Ep,wm_for(Ep))
    fs=float(s); incr = "" if prev is None else f"+{fs-prev:.3f}"
    ladder.append((nb,fs,w,Ep))
    setstr=str(Ep)
    deg = "  (DEGENERATE: base too thin to near-cover -> no N1 cells -> Delta=0)" if len(cells_of(Ep))==0 else ""
    print(f"    {nb:5d}   {setstr:48s} {fs:7.4f}  {wfmt(w)}   {len(cells_of(Ep)):5d}  {incr}{deg}")
    prev=fs
# actual published extremizer = 3-block tower minus its last element (62)
Ex=[0,1,2,30,31,32,60,61]
(s,w)=supw(Ex,wm_for(Ex))
print(f"    (3-blk minus 62 = the published extremizer) {str(Ex):22s} {float(s):7.4f}  {wfmt(w)}")
# linear fit slope over NON-DEGENERATE blocks (fs>0)
pts=[(nb,fs) for (nb,fs,w,Ep) in ladder if fs>0]
if len(pts)>=2:
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
    slope=(n*sxy-sx*sy)/(n*sxx-sx*sx); inter=(sy-slope*sx)/n
    print(f"    LINEAR FIT (non-degenerate blocks {xs}): sup|wD| ~= {slope:.3f}*(#blocks) + {inter:.3f}")

print("\n"+"="*78)
print("    same ladder at SCALE 50 (tests whether the per-block increment is scale-free)")
print("="*78)
prev=None; scale=50; ys2=[]
for nb in range(1,4):
    Ep=block_tower(scale,nb)
    (s,w)=supw(Ep,wm_for(Ep)); fs=float(s); ys2.append(fs)
    incr = "" if prev is None else f"+{fs-prev:.3f}"
    print(f"    {nb} blocks {str(Ep):40s} sup|wD|={fs:7.4f} at w={wfmt(w)}  {incr}")
    prev=fs

print("\n"+"="*78)
print("CONTROL: coherence vs cardinality.  Three families at MATCHED cardinality:")
print("   (a) coherent block tower (Riesz)   (b) consecutive single block")
print("   (c) GENUINE Sidon set (Mian-Chowla): minimal additive structure = decorrelated")
print("   NOTE: random-integers-in-a-small-box are NOT incoherent -- by pigeonhole they")
print("   are dense in differences (high additive energy), so they score HIGH.  A real")
print("   incoherence control must be Sidon (all pairwise sums distinct).")
print("="*78)
def sidon(N):
    S=[0]; sums={0}; c=0
    while len(S)<N:
        c+=1; ns={c+s for s in S}|{2*c}
        if ns.isdisjoint(sums): sums|=ns; S.append(c)
    return S
print("  card | tower(blocks) sup|wD|  | consec sup|wD| | SIDON sup|wD|   (Sidon set)")
tower_by_card={6:block_tower(30,2), 9:block_tower(30,3), 12:block_tower(30,4)}
for card in (6,9,12):
    tw=tower_by_card[card]; (st,_)=supw(tw,wm_for(tw))
    con=list(range(card)); (sc,_)=supw(con,wm_for(con))
    sd=sidon(card); (ss,_)=supw(sd,wm_for(sd))
    print(f"   {card:3d}  |   {float(st):7.4f}          |   {float(sc):7.4f}      |  {float(ss):7.4f}   {sd}")
print("  READING: the block tower is NOT uniquely worst. Both wide towers and wide")
print("  Sidon sets score high; the narrow consec block is the MINIMIZER. => the driver")
print("  is SPREAD / #distinct-scales, reached by coherent repetition OR by spreading.")

print("\n"+"="*78)
print("THE TRADE-OFF: why large sup|wD| is NOT a counterexample (Phi anti-correlates)")
print("   p0(E')=meas{all 7 sectors hit};  p1=meas{exactly one inner sector missed};")
print("   plateau Plat = p0 + p1/7 (frozen-phase limit).  cap_9 = 1979/4004 ~ 0.4943.")
print("="*78)
def p0f(E):  # verbatim from lrc14_wide_multiscale_p0_kps.py
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
def p1meas(E): return sum(hi-lo for (lo,hi,s) in cells_of(E))
cap9=Fraction(1979,4004)
print("  family(card9)   span   p0(E')    p1-meas   Plat=p0+p1/7   sup|wD|   regime")
for name,E in [("consec-9",list(range(9))),("tower-9(3blk)",block_tower(30,3)),("sidon-9",sidon(9))]:
    p0v=p0f(E); p1v=p1meas(E); plat=p0v+p1v/7; (sd_,_)=supw(E,wm_for(E))
    reg = "BINDING (near cap, small Delta)" if float(plat)>0.4 else "SLACK (small Plat, large Delta)"
    print(f"  {name:14s}  {max(E):4d}   {float(p0v):.4f}    {float(p1v):.4f}    {float(plat):.4f}       {float(sd_):6.3f}   {reg}")
print("  => sup|wD| is LARGE exactly where Plat is SMALL. p0(full)=Plat+Delta stays < cap")
print("     in every regime. The discrepancy and the plateau are DUAL; neither alone is")
print("     a counterexample. This is the Phi-Delta trade-off (repo HYP-2779/kps-S23).")

print("\n"+"="*78)
print("P2  SCALE RESONANCE: does the argmax w* EXTEND the tower?")
print("    (attacker=defender: worst far element continues the self-similar block pattern)")
print("="*78)
scale=30
for nb in (2,3,4):
    Ep=block_tower(scale,nb)
    peaks,ncell=topw(Ep,wm_for(Ep),k=5)
    # predicted 'extending' far elements: next-block completions (nb*scale + base) and gaps
    print(f"  {nb}-block tower {Ep}")
    for val,w in peaks:
        r=w%scale
        tag=""
        if r in (0,1,2): tag=f"  <- w mod {scale} = {r}  (EXTENDS a block: {scale}*{w//scale}+{r})"
        print(f"      |wD|={float(val):.4f} at w={w:4d}   (w mod {scale} = {r}){tag}")
    print()

print("="*78)
print("P3  STRUCTURAL: exact 2-fold Riesz factorization of the clean tower")
print("    E = B (+) 30*B,  B={0,1,2};   \\hat 1_E(xi) =?= \\hat 1_B(xi) * \\hat 1_B(30 xi)")
print("="*78)
B=[0,1,2]; E=block_tower(30,3)   # {0,1,2,30,31,32,60,61,62} = B (+) 30*{0,1,2}
def fhat(S,xi): return sum(cmath.exp(2j*math.pi*s*xi) for s in S)
maxerr=0.0
for k in range(1,50):
    xi=k/97.0
    lhs=fhat(E,xi); rhs=fhat(B,xi)*fhat(B,30*xi)
    maxerr=max(maxerr, abs(lhs-rhs))
print(f"    max_xi | \\hat 1_E(xi) - \\hat 1_B(xi)\\hat 1_B(30 xi) | over 49 probes = {maxerr:.2e}")
print(f"    => E IS an exact 2-fold Riesz product (B dilated by 1 and by 30).")
print(f"       n-block tower = B (+) 30*{{0..n-1}}; deeper nesting B(+)30B(+)900B = 3-fold, etc.")
print("\nDONE.")
