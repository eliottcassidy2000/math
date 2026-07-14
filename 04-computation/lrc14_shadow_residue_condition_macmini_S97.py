#!/usr/bin/env python3
"""mac-mini-S97: make klein-S299's 'loneliness at a/k depends only on residues mod k' EXACT and RIGOROUS.
Derive the exact shadow-witness delta-interval at t=a/k+delta and reduce it to a clean RESIDUE-mod-k
condition; verify it matches true loneliness; confirm covering => some k<=13 works. Hands kps a finite
decidable form. (klein's Farey/three-gap frame: the a/k, k<=13, are the Farey dissection F_13.)"""
from fractions import Fraction as F
from math import gcd

c114=F(1,14)
def lonely_at(S,t):  # exact: is every ||c t|| >= 1/14 ?
    for c in S:
        r=(c*t)%1
        if min(r,1-r)<c114: return False
    return True
def shadow_interval(S,a,k):
    """EXACT witness delta-interval at t=a/k+delta (delta>0). Per speed c, residue r=(c*a)%k:
       r=0 (k|c): ||ct||=||c*delta||=c*delta -> need delta in [1/(14c), 13/(14c)]
       r!=0: ||ct|| = || r/k + c*delta ||; signed s = r if r<=k/2 else r-k; for small delta>0,
             value ~ |s|/k + sign(s)*c*delta. Constrain to stay in [1/14, 13/14]."""
    lo=F(0); hi=F(1,2)  # delta window (open-ish); start wide, intersect
    for c in S:
        r=(c*a)%k
        if r==0:
            lo=max(lo,F(1,14*c)); hi=min(hi,F(13,14*c))
        else:
            s = r if r<=k//2 else r-k          # signed residue in (-k/2, k/2]
            base=F(abs(s),k)                    # ||r/k|| at delta=0 (>=1/k >=1/13 > 1/14 for k<=13)
            if s>0:   # value increases with delta: base + c*delta, need <= 13/14 (upper drift)
                hi=min(hi, (F(13,14)-base)/c)
            else:     # value decreases: base - c*delta, need >= 1/14 (lower drift, the BINDING one)
                hi=min(hi, (base-c114)/c)
        if lo>=hi: return None
    return (lo,hi)

def covering(S):  # {1}uS covers moduli 2..14 (S = the cluster; 1 covers nothing)
    return all(any(v%q==0 for v in S) for q in range(2,15))

# residual covering clusters (ratio 6..13, no isolated far) -- constructed to cover 2..14, packed band
tests=[
  ("A", [15,18,20,21,22,24,26,28,33,35,39,42]),   # spread ~15..42, ratio 2.8 (THM-744 territory) -- control
  ("B", [2,26,28,30,33,35,36,39,40,42,44,45]),     # has 2 + far -> ratio 22.5
  ("C", [3,4,5,26,28,33,35,36,39,40,42,44]),       # ratio ~14.7
  ("D", [2,3,26,28,33,35,36,39,40,42,44,45]),      # ratio 22.5
]
print("residue-mod-k shadow-witness (k<=13) vs TRUE loneliness, on covering residual clusters:")
for name,C in tests:
    if not covering(C): print(f"  {name}: NOT covering, skip ({C})"); continue
    ratio=max(C)/min(C); good_ks=[]
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1: continue
            if not (c114 <= F(a,k) <= F(13,14)): continue  # a/k in middle (runner 1 safe near a/k)
            iv=shadow_interval(C,a,k)
            if iv:
                # verify: a true lonely point in the interval (midpoint), exact
                dm=(iv[0]+iv[1])/2; t=F(a,k)+dm
                if lonely_at([1]+C,t):
                    good_ks.append((k,a,iv)); break
        if good_ks and good_ks[-1][0]==k: pass
    mink = good_ks[0][0] if good_ks else None
    ks=sorted(set(g[0] for g in good_ks))
    print(f"  {name}: C ratio={ratio:.1f}, covering; good resonances k in {ks if ks else 'NONE'}; min k={mink}")
    if good_ks:
        k,a,iv=good_ks[0]; print(f"       witness t={a}/{k}+delta, delta in ({iv[0]},{iv[1]}); a/k={float(F(a,k)):.3f} (middle); VERIFIED lonely")
print()
print("CLEAN REDUCTION (for kps): the shadow at a/k is good iff, over speeds c in C with residue r=(c*a)%k:")
print("  (i) k-divisible (r=0): max/min < 13  [the shadow delta-window]  -- AUTOMATIC if C has ratio<13;")
print("  (ii) NEGATIVE signed residues (r>k/2, 'edge below'): c*(14-k)/(14k) ... i.e. no speed with")
print("       residue k-1,k-2,... too LARGE relative to min-k-divisible (the drift-down killers).")
print("  So loneliness at a/k = a RESIDUE-mod-k PATTERN condition. Covering forces r=0 speeds to EXIST")
print("  (a multiple of k) for every k<=14; the open uniform claim = some k<=13 has NO drift-down killer.")

print("\n"+"="*70)
print("BROAD TEST: does a k<=13 shadow witness close the WHOLE covering case, or only packed?")
print("="*70)
extremals=[
  ("deep well {1..12,182}", list(range(2,13))+[182]),          # covering-MIN, isolated far (disc_v territory)
  ("multi-killer {1..11,13,84}", list(range(2,12))+[13,84]),
  ("near-AP {2..13}u{182}\\{}", list(range(2,13))+[182]),
  ("drop-6 {1..14}\\6", [v for v in range(2,15) if v!=6]),
  ("{1,90..101}", list(range(90,102))),
  ("{2,3,4,5,6,7,8,9,10,11,13,84}", [2,3,4,5,6,7,8,9,10,11,13,84]),
]
for name,C in extremals:
    cov = all(any(v%q==0 for v in C) for q in range(2,15))
    if not cov: print(f"  {name}: NOT covering as cluster, skip"); continue
    found=None
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(C,a,k)
            if iv:
                dm=(iv[0]+iv[1])/2
                if lonely_at([1]+C, F(a,k)+dm): found=(k,a,iv); break
        if found: break
    if found:
        k,a,iv=found; print(f"  {name}: shadow witness at k={k}, t={a}/{k}+delta in ({iv[0]},{iv[1]}) -- VERIFIED lonely")
    else:
        print(f"  {name}: NO k<=13 shadow witness  <-- isolated-far => disc_v territory (complementary)")
