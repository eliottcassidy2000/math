#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE CROSS-PRIME APEX FAMILY: how the core gap-spectrum varies with the apex prime p (klein-S21).

For odd primes p: the binding DOUBLET gap = 2-2cos(pi/p) = 4 sin^2(pi/2p) (lambda_min of the signless
Laplacian of C_p, positive iff p odd). The FLAT/difference-set cores: the QR is a Paley (p,(p-1)/2,*)
difference set IFF p = 3 mod 4 (so the apex 7=3mod4 is the Paley case). Map: #distinct gap values, the
min (doublet), the QR spectrum, and the p mod 4 split.
"""
import math, cmath, itertools
from collections import Counter

def spectrum(O,p,W):
    return [abs(sum(W**((k*x)%p) for x in O))**2 for k in range(p)]
def gap(O,p,W):
    sp=spectrum(O,p,W); return min(sp[k] for k in range(1,p))
def is_flat(O,p,W):
    sp=[round(x,6) for x in spectrum(O,p,W)]; return len(set(sp[1:]))==1

print("="*88)
print(" CROSS-PRIME APEX FAMILY: doublet gap, QR structure, #gap values, by apex prime p")
print("="*88)
print(f"\n {'p':>3} {'p mod4':>6} {'doublet gap 4sin^2(pi/2p)':>26} {'#distinct gaps':>15} {'QR flat? (Paley diff-set)':>26}")
for p in [3,5,7,11,13,17,19]:
    W=cmath.exp(2j*math.pi/p)
    db=2-2*math.cos(math.pi/p)
    QR=set(pow(a,2,p) for a in range(1,p))   # quadratic residues
    qr_flat = is_flat(QR,p,W)
    # distinct gap values: sample over all subsets is 2^p (too big for p>=13); sample sizes 1..min(p,7) reps
    gaps=set()
    if p<=13:
        for r in range(1,p):
            # sample: all r-subsets is C(p,r); cap to avoid blowup
            combs=itertools.combinations(range(p),r)
            cnt=0
            for O in combs:
                gaps.add(round(gap(set(O),p,W),5)); cnt+=1
                if cnt>20000: break
    ndist = len(gaps) if gaps else None
    print(f" {p:>3} {p%4:>6} {db:>26.6f} {str(ndist):>15} {str(qr_flat):>26}")

print("\n"+"="*88)
print(" THE DOUBLET-GAP SEQUENCE 4sin^2(pi/2p) (the apex obstruction by prime) and its asymptotics:")
print("="*88)
for p in [3,5,7,11,13,17,19,23,29]:
    db=2-2*math.cos(math.pi/p); approx=math.pi**2/p**2
    print(f"   p={p:>2}: 4sin^2(pi/2p) = {db:.6f}   (~ pi^2/p^2 = {approx:.6f}, ratio {db/approx:.4f})")
print("   => the apex obstruction shrinks like pi^2/p^2; SMALLER odd apex prime = LARGER floor.")
print("      p=7 (the LRC(14) apex) gives 0.198 -- the 3rd-largest (after p=3,5).")

print("\n"+"="*88)
print(" p mod 4 SPLIT (the Paley/difference-set theme, CLAUDE.md nu_2=0 <=> p=3 mod4):")
print("="*88)
for p in [3,5,7,11,13,17,19]:
    W=cmath.exp(2j*math.pi/p); QR=set(pow(a,2,p) for a in range(1,p))
    sp=[round(x,4) for x in spectrum(QR,p,W)]
    flat=is_flat(QR,p,W)
    note = "PALEY DIFFERENCE SET (flat, |QR-hat(k)|^2=(p+1)/4 const)" if flat else "QR symmetric (-1 in QR), NOT a diff-set; 2 spectral values"
    print(f"   p={p:>2} (={p%4} mod4): QR={sorted(QR)} -> nonzero spectrum {{{','.join(str(v) for v in sorted(set(sp[1:])))}}}  {note}")
print("   => p=3 mod4 (7,11,19): QR is a Paley difference set (FLAT). p=1 mod4 (5,13,17): QR NOT flat.")
print("      the apex 7=3mod4 is in the Paley/difference-set/flat-QR class (the gap-2 cores exist).")
