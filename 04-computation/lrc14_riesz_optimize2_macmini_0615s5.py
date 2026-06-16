#!/usr/bin/env python3
"""
lrc14_riesz_optimize2_macmini_0615s5.py  (mac-mini-2026-06-15-S5; OPEN-Q-104)

DIRECT-GRID (exact) Riesz-product optimization. R(τ)=∏_v (1 + a_v cos 2π v τ) ≥ 0
(need |a_v|≤1). Certificate: ratio = ∫MR/∫R < 1 ⟹ that config is LOOSE.
Simple ∏(1-cos) (a_v=-1) certifies the dilated core (0.955) but NOT the harder
extremizer {1..13}\{6}∪{56} (1.064). Push:
 (1) UNIFORMITY: ∏(1-cos) ratio over many loose primitive mult-of-14 configs — which fail?
 (2) OPTIMIZE the per-speed amplitudes a_v (coordinate descent, direct grid) for the
     extremizer — can we get ratio < 1 (a certificate for the hardest core)?
 (3) try richer R: add 2nd-harmonic factors (1-cos 4πvτ), and a dissociated subset.
"""
import sys, math, itertools, random
from math import cos, pi, gcd

sys.stdout.reconfigure(line_buffering=True)
random.seed(615)

def gcd_all(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

def ratio_direct(S, factors, Q=30030):
    """factors = list of (freq, amp). R=∏(1+amp cos 2π freq τ). ratio=∫MR/∫R."""
    rad=Q//14; sMR=0.0; sR=0.0
    for k in range(Q):
        R=1.0
        for (f,amp) in factors:
            R*=(1+amp*cos(2*pi*f*k/Q))
        if R<=0.0:
            if R<0: R=0.0  # clamp tiny negatives from rounding (factors with |amp|<=1 keep R>=0)
            continue
        M=sum(1 for v in S if ((v*k)%Q)<=rad or ((v*k)%Q)>=Q-rad)
        sMR+=M*R; sR+=R
    return (sMR/sR) if sR>1e-15 else float('inf'), sR/Q

print("="*78)
print("(1) UNIFORMITY of R=∏_v(1-cos 2πvτ): ratio over loose primitive mult-of-14 configs")
print("="*78)
loose_configs = {
  "ext {1..13}\\{6}∪56":  sorted(set([x for x in range(1,14) if x!=6]+[56])),
  "ext {1..13}\\{12}∪84": sorted(set([x for x in range(1,14) if x!=12]+[84])),
  "dilated {1..12}∪14":   sorted(set(list(range(1,13))+[14])),
  "evader 7*{1..12}∪13":  sorted(set([7*i for i in range(1,13)]+[13])),
  "d3 {3,6..36}∪182":     sorted(set([3*i for i in range(1,13)]+[182])),
}
fails=[]
for name,S in loose_configs.items():
    if len(S)!=13: continue
    r,intR = ratio_direct(S, [(v,-1.0) for v in S])
    ok = "CERT<1" if r<1 else "FAILS(≥1)"
    if r>=1: fails.append(name)
    print(f"   {name:22s}: ratio={r:.4f}  {ok}   (∫R={intR:.2e})")
print(f"   --> ∏(1-cos) certifies {len(loose_configs)-len(fails)}/{len(loose_configs)}; FAILS on: {fails}")

print("\n" + "="*78)
print("(2) OPTIMIZE per-speed amplitudes a_v (coord descent) for the EXTREMIZER {..\\6..∪56}")
print("="*78)
Sx = sorted(set([x for x in range(1,14) if x!=6]+[56]))
amps = {v:-1.0 for v in Sx}
def rat(amps_dict):
    return ratio_direct(Sx, [(v,amps_dict[v]) for v in Sx])[0]
cur = rat(amps)
print(f"   start (all a=-1): ratio={cur:.4f}")
for it in range(4):
    improved=False
    for v in Sx:
        best_a=amps[v]; best_r=cur
        for a in [-1.0,-0.9,-0.7,-0.5,-0.3,0.0,0.3,0.6]:
            amps[v]=a; r=rat(amps)
            if r<best_r-1e-5: best_r=r; best_a=a
        amps[v]=best_a; cur=best_r
        if best_a!=-1.0 or True: pass
    print(f"   pass {it+1}: ratio={cur:.4f}  {'CERT<1!' if cur<1 else ''}")
    # stop if converged
print(f"   optimized amplitudes (nonzero-from-default): "
      f"{[(v,round(amps[v],1)) for v in Sx if abs(amps[v]+1.0)>0.05]}")
print(f"   --> best extremizer ratio = {cur:.4f}  {'CERTIFIED LOOSE' if cur<1 else 'still ≥1'}")

print("\n" + "="*78)
print("(3) RICHER R: add 2nd-harmonic factors (1 - b cos 4πvτ) for small speeds")
print("="*78)
for extra in [[(2*v,-0.5) for v in (1,2,3)], [(2*v,-0.4) for v in (1,2,3,4,5)]]:
    facs=[(v,-1.0) for v in Sx]+extra
    r,iR = ratio_direct(Sx, facs)
    print(f"   ∏(1-cos 2πvτ)·∏(1-0.x cos 4π·{[f for f,_ in extra][:3]}..τ): ratio={r:.4f}  {'CERT<1!' if r<1 else ''}")
# also: drop the stranger 56 from the product (it is the dominant/dodge runner)
r2,_ = ratio_direct(Sx, [(v,-1.0) for v in Sx if v!=56])
print(f"   ∏(1-cos) over core WITHOUT the stranger 56: ratio={r2:.4f}  {'CERT<1!' if r2<1 else ''}")
print("\nDONE.")
