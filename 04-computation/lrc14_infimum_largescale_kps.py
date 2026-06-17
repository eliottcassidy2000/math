#!/usr/bin/env python3
"""
lrc14_infimum_largescale_kps  (kind-pasteur, both sides)

DISPROVE-SIDE idea fuel + PROVE-side bound test: can large-element families push
the lonely measure BELOW 1/1260, or even toward 0?

Three probes:
 [A] single perturbation drop-12 -> w for ALL w up to 5000 (float), find the
     L(w) curve and its minimum; show it does NOT dip below 1/1260 and -> 6/7*?
     Actually as w->inf, dropping 12 and adding huge w: the added speed's danger
     arcs become a fine uniform comb of measure 1/7, decoupling. So L(drop12,w)
     -> meas(G_12) * (1 - 1/7)?? test the limit.
 [B] 2-element perturbation: replace {11,12} or {12,13} etc by two large coordinated
     speeds (e.g. w, 2w or w, w+1) and see how low L goes. Scan w up to 300.
 [C] the "resonance" families: {1..11,13, m} for m a large multiple structure,
     and {1..11,13,24}-type sporadics scaled. Find smallest L vs element size.

Goal: locate exactly WHERE (which family, which direction) any L<1/1260 or L->0
would have to live. Save the L(w) curves.
"""
from fractions import Fraction as F
import math

def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort(); tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    arcs.sort(); t=F(0); cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else: t+=ch-cl; cl,ch=lo,hi
    return F(1)-(t+ch-cl)

base=list(range(1,14))

# [A] L(drop12 -> w) curve, w up to 5000
print("[A] single perturbation drop-12 -> w, L as w grows",flush=True)
rest=[x for x in base if x!=12]
minL=(2.0,None); samples=[]
for w in range(14,5001):
    if w in rest: continue
    lf=L_float(rest+[w])
    if lf<minL[0]: minL=(lf,w)
    if w in (14,20,24,36,48,60,72,100,200,500,1000,2000,5000):
        samples.append((w,lf))
print("   sample L(w):")
for w,lf in samples:
    print(f"     w={w:5d}: L={lf:.6e}")
print(f"   min over w<=5000: L={minL[0]:.6e} at w={minL[1]}")
# limit: G_12 measure times fraction not covered by a 'random' fine comb
G12meas=float(L_exact([x for x in base if x!=12])) # this is meas of complement = e-gap? careful
# Actually L_exact(rest) = lonely measure of the 12-element set rest = meas(G_12)
print(f"   meas(G_12)=L_exact(rest)= {L_exact(rest)} = {float(L_exact(rest)):.6e}")
print(f"   6/7 * meas(G_12) = {F(6,7)*L_exact(rest)} = {float(F(6,7)*L_exact(rest)):.6e}  (decoupling limit)")

# [B] coordinated 2-element large perturbations
print("\n[B] 2-element perturbations replacing {a,b} by large coordinated speeds",flush=True)
print("    drop {11,12}, add {w, 2w}:")
bestB=(2.0,None)
for w in range(14,301):
    S=[x for x in base if x not in (11,12)]+[w,2*w]
    if len(set(S))!=13: continue
    lf=L_float(S)
    if lf<bestB[0]: bestB=(lf,(w,2*w))
print(f"      best L={bestB[0]:.6e} at (w,2w)={bestB[1]}")
print("    drop {12,13}, add {w, w+1}:")
bestB2=(2.0,None)
for w in range(14,301):
    S=[x for x in base if x not in (12,13)]+[w,w+1]
    if len(set(S))!=13: continue
    lf=L_float(S)
    if lf<bestB2[0]: bestB2=(lf,(w,w+1))
print(f"      best L={bestB2[0]:.6e} at (w,w+1)={bestB2[1]}")
print("    drop {12,13}, add {w, 2w}:")
bestB3=(2.0,None)
for w in range(14,301):
    S=[x for x in base if x not in (12,13)]+[w,2*w]
    if len(set(S))!=13: continue
    lf=L_float(S)
    if lf<bestB3[0]: bestB3=(lf,(w,2*w))
print(f"      best L={bestB3[0]:.6e} at (w,2w)={bestB3[1]}")

# [C] scale a tight sporadic: {1..11,13,24}*nothing -- instead try {1..11,13,m}
print("\n[C] family {1..11,13,m}: L vs m (single perturbation of AP at 12)",flush=True)
fam=[x for x in base if x!=12]
rows=[]
for m in range(14,400):
    if m in fam: continue
    lf=L_float(fam+[m])
    rows.append((lf,m))
rows.sort()
print("    smallest L in family {1..11,13,m}, m<=399:")
for lf,m in rows[:12]:
    Le=L_exact(fam+[m])
    print(f"      m={m:4d}: L={Le}={float(Le):.6e}")

# crucial: does ANY of [B] beat 1/1260? exact-confirm the winners
print("\n[D] exact-confirm any 2-element winner below 1/1260:",flush=True)
TARGET=F(1,1260)
def check(label, S):
    Le=L_exact(S)
    rel="<" if Le<TARGET and Le>0 else ("=" if Le==TARGET else (">" if Le>TARGET else "=0"))
    print(f"    {label}: L={Le}={float(Le):.6e}  {rel} 1/1260")
if bestB[1]: check(f"drop{{11,12}} add {bestB[1]}", [x for x in base if x not in (11,12)]+list(bestB[1]))
if bestB2[1]: check(f"drop{{12,13}} add {bestB2[1]}", [x for x in base if x not in (12,13)]+list(bestB2[1]))
if bestB3[1]: check(f"drop{{12,13}} add {bestB3[1]}", [x for x in base if x not in (12,13)]+list(bestB3[1]))
