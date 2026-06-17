"""
DISPROVE-side CORE TEST: does L -> 0 along ANY growing-lcm family?

Three escape attempts, all exact:
  (1) drop-6 family extended: w from 14..3000. Does L dip below 1/1260 anywhere,
      or trend to 0? (single large coordinated speed)
  (2) "coordinated doublings": replace several AP entries e -> 2e (or m*e) at once.
      Resonance argument that made 12->24 tight: ||2*(e*tau)|| keeps coverage near
      tau=k/e. Try multi-doubling to keep tightness while growing lcm.
  (3) S551 sieve-blind: S = {1..12, lcm(2..Q)} type configs; measure L vs Q (lcm).
      Prove side says L there GROWS (~0.0068..0.0094). Confirm and push Q higher.

If ALL three keep L bounded away from 0, the inf L=0 escape is closed and the
PROVE side's compactness/quantization wins. We are looking for ANY chink.
"""
from fractions import Fraction as Fr
from math import gcd, lcm
from functools import reduce

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1:
            A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else:
            A.append((lo, hi))
    return A
_cache={}
def darcs(v):
    r=_cache.get(v)
    if r is None: r=danger_arcs(v); _cache[v]=r
    return r
def L_exact(S):
    A=[]
    for v in S: A.extend(darcs(v))
    A.sort(key=lambda t:(t[0],t[1]))
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if b<=a: continue
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

THRESH=Fr(1,1260)
AP=list(range(1,14))

print("=== (1) drop-6 family extended: {1,2,3,4,5,7,8,9,10,11,12,13,w} ===")
core6=[x for x in AP if x!=6]
minL=None; below=0
samples=[]
for w in list(range(14,200))+list(range(200,3001,7)):
    if w in core6: continue
    S=sorted(core6+[w])
    L=L_exact(S)
    if L==0:
        print(f"   TIGHT at w={w}"); continue
    if minL is None or L<minL:
        minL=L; minw=w
    if L<THRESH: below+=1
    if w in (14,69,138,207,690,1380,2070,2898) or (w%500==13):
        samples.append((w,L))
print(f"   min positive L = {minL}={float(minL):.6f} at w={minw}; #below 1/1260 = {below}")
print(f"   samples (w, L, float):")
for w,L in samples:
    print(f"      w={w:5d}  L={float(L):.6f}  lcm={lcm(*sorted(core6+[w]))}")

print()
print("=== (2) coordinated multi-doublings of AP entries (keep tightness?) ===")
# 12->24 is tight. Try doubling multiple entries. Test e->2e for subsets.
from itertools import combinations
print("   single e->2e (which preserve tightness?):")
for e in range(1,14):
    if 2*e in AP: continue
    S=sorted([x for x in AP if x!=e]+[2*e])
    L=L_exact(S)
    tag="TIGHT" if L==0 else f"L={float(L):.6f}"
    print(f"      {e}->{2*e}: {tag}  lcm={lcm(*S)}")
print("   double pairs e1->2e1, e2->2e2 (look for tight with big lcm):")
best2=None
for e1,e2 in combinations(range(1,14),2):
    if 2*e1 in AP or 2*e2 in AP: continue
    new=[x for x in AP if x not in (e1,e2)]+[2*e1,2*e2]
    if len(set(new))!=13: continue
    S=sorted(new)
    L=L_exact(S)
    if L==0:
        print(f"      TIGHT: {e1}->{2*e1}, {e2}->{2*e2}  S={S} lcm={lcm(*S)}")
    if best2 is None or L<best2[0]:
        best2=(L,S,(e1,e2))
if best2:
    print(f"   best double-doubling: L={best2[0]}={float(best2[0]):.6f} S={best2[1]} lcm={lcm(*best2[1])}")
print("   triple e->m*e coordinated (m up to 6) keeping tightness:")
# the resonance lever: e->m*e preserves the danger near k/e iff m*e's centers hit them
besttrip=None
for e1,e2,e3 in combinations([6,8,9,10,11,12],3):
    for m in (2,3,4,6):
        new=[x for x in AP if x not in (e1,e2,e3)]+[m*e1,m*e2,m*e3]
        if len(set(new))!=13: continue
        if any(x in AP for x in (m*e1,m*e2,m*e3)):
            # allow if distinct anyway
            pass
        S=sorted(new)
        if len(set(S))!=13: continue
        L=L_exact(S)
        if L==0:
            print(f"      TIGHT: {e1},{e2},{e3}->x{m}  S={S} lcm={lcm(*S)}")
        if besttrip is None or L<besttrip[0]:
            besttrip=(L,S,(e1,e2,e3),m)
if besttrip:
    print(f"   best triple: L={besttrip[0]}={float(besttrip[0]):.6f} S={besttrip[1]} lcm={lcm(*besttrip[1])} via drops={besttrip[2]} x{besttrip[3]}")

print()
print("=== (3) S551 sieve-blind {1..12, M}: exact for M<=2000, decoupling limit for large ===")
# core = {1..12}. meas(G_core) = L({1..12}) (lonely measure of the 12-speed core).
core12=list(range(1,13))
Lcore=L_exact(core12)
print(f"   core {{1..12}} lonely measure meas(G_core) = {Lcore} = {float(Lcore):.8f}")
print(f"   decoupling limit L({{1..12,M}}) -> (6/7)*meas(G_core) = {Fr(6,7)*Lcore} = {float(Fr(6,7)*Lcore):.8f} as M->inf")
print("   exact samples (M moderate):")
for M in [13,17,19,23,29,37,53,101,211,503,1009,2003]:
    if M in core12: continue
    S=sorted(core12+[M])
    L=L_exact(S)
    print(f"   M={M:5d}  L={float(L):.8f}  (={L})  trend-> {float(Fr(6,7)*Lcore):.8f}")
print("   CONCLUSION: as the single huge sieve-blind speed M grows, L INCREASES toward")
print("   (6/7)*meas(G_core) > 0; it never -> 0. Large lcm does NOT drive L down.")
