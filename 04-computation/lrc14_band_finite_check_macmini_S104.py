#!/usr/bin/env python3
"""mac-mini-S104: the covering-case CLOSURE band. THM-726 Step 2 = finite check outliers<=220; opus
floor = W>W0(~339-475). The gap = multi-killer with largest outlier in (220, ~500]. Verify this band
is all-lonely (M>=1/13), so the finite check PASSES; also check how many are reducible (THM-753 =>
LRC(<=13)) vs need the exact-Q certificate. This pins whether the covering case closes on the band."""
from fractions import Fraction as F
import numpy as np, random
def M_num(S,dens=2_000_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def is_safe_peel_exists(S):
    for v in sorted(S,reverse=True):
        core=[u for u in S if u!=v]
        if len(core)<2: return True
        mu0,t0=M_num(core,dens=400_000)
        if min((v*t0)%1.0,1-(v*t0)%1.0)>=mu0-2e-4: return True
    return False
def covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
one13=1/13
rng=random.Random(104)
# multi-killer families with largest outlier in the BAND (220, 500]: interval cores {1..k} + outliers
band=[]; below13=[]; tot=0; reducible=0
for _ in range(40000):
    k=rng.choice([9,10,11]); core=list(range(1,k+1))
    n_out=13-k
    outs=set()
    # outliers >=15 carrying missing moduli, with the LARGEST in (220,500]
    while len(outs)<n_out: outs.add(rng.choice([q*rng.randint(15,36) for q in [13,14,12,11,84,182,26,28]]))
    outs=sorted(outs)[:n_out]
    if not outs or max(outs)<=220 or max(outs)>500: continue
    S=sorted(set(core+outs))
    if len(S)!=13 or not covering(S): continue
    tot+=1
    Mm,_=M_num(S)
    if Mm < one13-1e-4: below13.append((round(Mm,4),S))
    if is_safe_peel_exists(S): reducible+=1
    band.append(Mm)
    if tot>=400: break
print(f"BAND multi-killer families (largest outlier in (220,500]): {tot}")
if band:
    print(f"  M: min={min(band):.4f} median={sorted(band)[len(band)//2]:.4f} max={max(band):.4f}  (1/13={one13:.4f}, 1/14={1/14:.4f})")
    print(f"  families with M < 1/13 (would need care): {len(below13)}")
    for M,S in sorted(below13)[:6]: print(f"     M={M}: {S}")
    print(f"  REDUCIBLE (safe-peel => LRC(<=13), THM-753): {reducible}/{tot} ({100*reducible/tot:.1f}%)")
print()
print("=> if the band is all M>=1/13 (or reducible => LRC(<=13)), the finite check PASSES and the")
print("   covering case closes: single(THM-724) + multi[Step2<=220 + THM-751/753 monotonicity + this band")
print("   + opus floor >W0]. Fully closing = executing this bounded band exactly (kps/opus exact-Q).")
