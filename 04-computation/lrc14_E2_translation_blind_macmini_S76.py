#!/usr/bin/env python3
"""mac-mini-S76: the SHARP reason E2 is blind (opus-S182 redirect). E2 = additive energy
#{a+b=c+d} is TRANSLATION-INVARIANT; loneliness L is dilation-invariant but NOT translation-
invariant. So E2 CANNOT separate the AP from a translate, but E3/Schur (#{a+b=c}) and L can.
=> the fine-scale (1/14) invariant is E3/Schur, NOT E2/Freiman. Confirms my pairwise-blind
finding is the E3/Schur redirect."""
from itertools import product
c=1.0/14
def E2(A):  # additive energy #{a+b=c+d}
    from collections import Counter
    s=Counter(a+b for a in A for b in A); return sum(v*v for v in s.values())
def Schur(A):  # #{(a,b,c) in A^3 : a+b=c}  (3-term, translation-SENSITIVE)
    Aset=set(A); return sum(1 for a in A for b in A if a+b in Aset)
def L(A,res=200000):
    n=0
    for j in range(res):
        t=(j+0.5)/res
        if all(min((v*t)%1,1-((v*t)%1))>=c for v in A): n+=1
    return n/res
print("set                    | E2(add.energy) | Schur(a+b=c) | L(lonely)")
print("-"*66)
for nm,A in [("AP {1..13}",list(range(1,14))),
             ("translate {6..18}",list(range(6,19))),
             ("translate {11..23}",list(range(11,24))),
             ("dilate 2*{1..13}={2..26 even}",[2*k for k in range(1,14)]),
             ("{1..11,13,84} (cov)",[*range(1,12),13,84])]:
    print(f"{nm:22s} | {E2(A):14d} | {Schur(A):12d} | {L(A):.5f}")
print("\n=> E2 IDENTICAL for AP and its translates (translation-invariant) but L DIFFERS")
print("   (dilation-inv only). So E2/additive-energy/Freiman is TRANSLATION-BLIND at the 1/14")
print("   scale. Schur(a+b=c) and L both drop off the AP => the correct invariant is E3/SCHUR.")
print("   (Dilate 2*{1..13}: same E2 as AP up to scaling, L=0 too -- dilation-invariance holds.)")
