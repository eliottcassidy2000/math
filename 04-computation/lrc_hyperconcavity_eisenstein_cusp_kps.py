#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYPERCONCAVITY: a QUANTITATIVE order parameter for the regularizable/Eisenstein (ordered/AP)
vs un-regularizable/cusp (disordered/generic) dichotomy of HYP-3768/3774/3775.

kind-pasteur-2026-07-01. opus HYP-3775 gave the dichotomy a BINARY geometric criterion:
margin = classical Dedekind sum (=> -1/12=zeta(-1), Eisenstein) EXACTLY when the witness
residues are a scaled INTERVAL (AP); beaters (generic non-AP) are non-modular/cusp. The OPEN
part (mac-mini HYP-3774, opus HYP-3775 sec 5, OPEN-Q-108) is to QUANTIFY the residual -- "how
far below the ordered locus" a disordered set reaches = "how big is the cusp form".

This session's extension uses the owner's "hyperconcavity" hint LITERALLY. Two facts:
  (H) HYPERCONCAVITY: the autocorrelation A_R(t)=#{(x,y) in R^2: x-y=t} of an INTERVAL is the
      FEJER TRIANGLE (m-|t|)_+ -- LOG-CONCAVE, and its Fourier transform is the Fejer kernel
      |hat 1_R|^2 >= 0 (the level-k E2/"cyclotomic magic function" of last session). Order =>
      log-concave autocorrelation => positive-definite => REGULARIZABLE.
  (E) The additive energy E(R)=sum_t A_R(t)^2 = ||hat 1_R||_4^4 (the L4 Fourier / hypercontractive
      moment, last session's best predictor of p0) is MAXIMIZED by the interval/AP and MINIMIZED by
      Sidon/generic. So E is the CONTINUOUS order parameter: max = Eisenstein, drop = cusp.

CHECKABLE CLAIM: interval/AP  <=>  log-concave autocorrelation (0 defect)  <=>  additive energy
maximal  <=>  Dedekind sum closes (regularizable, -1/12). The DEFICIT from the interval = a
computable proxy for the cusp/residual, extending opus's binary criterion to a scalar.
"""
import sys, math, random
from fractions import Fraction as Fr
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def autocorr(R):
    R=list(R); A={}
    for x in R:
        for y in R:
            A[x-y]=A.get(x-y,0)+1
    return A
def add_energy(R):            # E(R) = sum_t A(t)^2 = #{(a,b,c,d): a-b=c-d} = ||hat 1_R||_4^4
    return sum(v*v for v in autocorr(R).values())
def lc_defect(R):
    """hyperconcavity defect: internal zeros + ratio-violations of the autocorrelation on its hull.
    0 <=> log-concave (Fejer triangle) <=> interval/AP-like."""
    A=autocorr(R); ts=sorted(A); lo,hi=ts[0],ts[-1]
    full=[A.get(t,0) for t in range(lo,hi+1)]
    internal_zeros=0; started=False
    for i,v in enumerate(full):
        if v>0: started=True
        if started and v==0 and any(full[j]>0 for j in range(i+1,len(full))): internal_zeros+=1
    ratio_viol=0
    for i in range(1,len(full)-1):
        a,b,c=full[i-1],full[i],full[i+1]
        if a>0 and b>0 and c>0 and b*b < a*c: ratio_viol+=1
    return internal_zeros+ratio_viol, internal_zeros, ratio_viol
def is_interval_upto_scale(R):
    R=sorted(set(R)); d=R[1]-R[0] if len(R)>1 else 1
    return all(R[i+1]-R[i]==d for i in range(len(R)-1))

# ---- Dedekind (for the regularizability / -1/12 side) ----
def saw(x):
    x=Fr(x)
    if x.denominator==1: return Fr(0)
    return x-(x.numerator//x.denominator)-Fr(1,2)
def dedekind(a,D):
    a%=D
    return sum(saw(Fr(i,D))*saw(Fr(a*i,D)) for i in range(1,D))

print("="*94)
print("PART 1  THE AUTOCORRELATION DICHOTOMY: ordered/interval => log-concave (Fejer) => max energy")
print("  E_norm = E(R)/E(interval same size) in (0,1]; defect = hyperconcavity (log-concavity) defect")
print("="*94)
def report(name, R):
    R=sorted(set(R)); m=len(R)
    interval=list(range(m)); Emax=add_energy(interval)
    E=add_energy(R); En=E/Emax; d,iz,rv=lc_defect(R)
    ordered = "ORDERED (interval/AP)" if is_interval_upto_scale(R) else ("~ordered" if d==0 else "DISORDERED")
    print(f"  {name:34s} |R|={m:2d}  E={E:5d}  E_norm={En:.3f}  lc_defect={d:2d} (zeros {iz}, viol {rv})  {ordered}")
    return En,d
print("  -- reference families --")
report("interval {0..9}", range(10))
report("AP step 3 {0,3,..,27}", [3*i for i in range(10)])
report("interval-with-hole (drop 3,4)", [0,1,2,5,6,7,8,9])
report("Sidon (Mian-Chowla) size 8", [0,1,3,7,12,20,30,44])
rng=random.Random(701); report("random size 8 in [0,40]", rng.sample(range(41),8))
print("  -- the actual LRC BEATERS (generic covering-min, opus HYP-3769) --")
beaters={7:[1,2,5,6,7,8], 8:[1,4,5,6,7,11,16], 9:[1,3,4,5,7,11,18,32], 10:[1,2,3,5,6,7,8,9,30]}
Mbeat ={7:Fr(2,13),8:Fr(2,15),9:Fr(4,33),10:Fr(4,37)}
beat_En={}
for n,S in beaters.items():
    En,d=report(f"beater n={n}  {S}", S); beat_En[n]=En
print("  -- the CONSTRUCTION witness residues (interval-core {1..n-2} + antipode -1, opus HYP-3775) --")
con_En={}
for n in range(7,15):
    D=n*n-n+1; Rres=list(range(1,n-1))+[D-1]   # {1,..,n-2} U {-1 mod Phi6}
    En,d=report(f"construction n={n} residues", Rres); con_En[n]=En

print("\n"+"="*94)
print("PART 2  THE REGULARIZABLE SIDE: the construction's Dedekind sum CLOSES (-> -1/12=zeta(-1))")
print("  s(n,Phi6) = -T/(12T+6),  T=n(n-1)/2;  margin=-12 s/n^2;  n^2*margin -> -12 zeta(-1)=1")
print("="*94)
for n in [6,7,10,14,50,200]:
    D=n*n-n+1; T=Fr(n*(n-1),2); closed=Fr(-T,12*T+6)
    if n<=60:
        s=dedekind(n,D); match=(s==closed)
    else:
        s=None; match="(skip enum)"
    margin=Fr(n,D)-Fr(1,n); n2m=n*n*margin
    print(f"  n={n:3d}: s(n,Phi6) closed={float(closed):+.6f}  {'enum '+str(match) if s is not None else match}   "
          f"n^2*margin={float(n2m):.5f} (->1);  Phi6 mod n = {D%n} (==1 => 1-step reciprocity = AP telescoping)")
print("  => construction is REGULARIZABLE: Phi6==1 mod n makes the AP/interval sawtooth sum telescope")
print("     in one reciprocity step -> closed form -> the -1/12 Eisenstein anomaly.")

print("\n"+"="*94)
print("PART 3  THE ORDER PARAMETER QUANTIFIES THE RESIDUAL (the open part)")
print("="*94)
print("  E_norm (order): construction (residues) vs beaters (both at their n):")
for n in sorted(beaters):
    print(f"    n={n}: construction E_norm={con_En[n]:.3f} (log-concave, Eisenstein) "
          f"vs beater E_norm={beat_En[n]:.3f} (disordered, cusp)  M_beater={float(Mbeat[n]):.4f}")
# does beater disorder (1 - E_norm) track its depth below the construction floor n/Phi6?
print("  correlation of beater DISORDER (1-E_norm) with DEPTH (n/Phi6 - M_beater), the residual size:")
xs=[]; ys=[]
for n in sorted(beaters):
    D=n*n-n+1; depth=float(Fr(n,D)-Mbeat[n]); dis=1-beat_En[n]; xs.append(dis); ys.append(depth)
    print(f"    n={n}: disorder(1-E_norm)={dis:.3f}   residual depth(n/Phi6 - M)={depth:.4f}")
# simple Pearson
mx=sum(xs)/len(xs); my=sum(ys)/len(ys)
cov=sum((x-mx)*(y-my) for x,y in zip(xs,ys)); sx=math.sqrt(sum((x-mx)**2 for x in xs)); sy=math.sqrt(sum((y-my)**2 for y in ys))
r=cov/(sx*sy) if sx>0 and sy>0 else float('nan')
print(f"    Pearson(disorder, residual depth) = {r:+.3f}")

print("\n"+"="*94)
print("PART 4  WHY ORDER = REGULARIZABLE: interval autocorr = Fejer triangle => Fourier >= 0 (pos-def)")
print("="*94)
for m in [5,8]:
    R=list(range(m)); A=autocorr(R)
    tri=all(A.get(t,0)==m-abs(t) for t in range(-(m-1),m))
    print(f"  interval size {m}: autocorrelation == triangle (m-|t|)_+ ? {tri}  -> Fourier = |hat 1_R|^2 = Fejer kernel >= 0")
print("  => the ordered (interval/AP) set has a LOG-CONCAVE (Fejer) autocorrelation and a NON-NEGATIVE")
print("     power spectrum -- positive-definite, hence its Dedekind/sawtooth sum telescopes and")
print("     REGULARIZES (-1/12 = zeta(-1), the E2/Eisenstein anomaly, = the level-k Fejer magic function).")
print("     Disorder breaks log-concavity => no telescoping => cusp/un-regularizable residual.")
print("DONE.")
