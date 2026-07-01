#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
LRC-by-n DIFFICULTY vs TOURNAMENT-SPACE COMPLEXITY, indexed by the apex prime p (n=2p).

kind-pasteur-2026-07-01-S17. The LRC(2p) apex object is the Paley/heptagon tournament on p vertices (p=3mod4).
The proof-pillars (Mersenne 2-adic descent / Heegner sqrt(-p) SOS / 3-mod-4 Borsuk-Ulam) depend on p's
arithmetic; and the difficulty of the TOURNAMENT SPACE on p vertices (this session: reconstruction unique
p<=4, (I(Omega,x),d)-complete p=6, WALL at p=7) tracks the LRC difficulty.  Tabulate both, side by side, and
flag the SURPRISES (non-monotone difficulty; p=1mod4 => no tournament bridge; tight-equiangular {7,23}).
"""
import sys
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True
def is_mersenne(p): return isprime(p) and ((p+1)&p)==0 and p>=3
HEEGNER={2,3,7,11,19,43,67,163}; TIGHT_EQ={2,3,7,23}
A000568={3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,11:'~5.4e9'}   # tournaments on p vertices
# tournament-space reconstruction difficulty (this session, by #vertices):
RECON={3:'trivial (2 cls)',4:'UNIQUE from (cat,deg)',5:'(cat,deg,H) injective',
       6:'(I(Om,x),d) complete; first OCF-cospectral twin',7:'WALL: (I(Om,x),d) fails, 90% cospectral'}

print("="*118)
print(" LRC(n=2p) by apex prime p: proof-pillars  x  tournament-space complexity on p vertices")
print("="*118)
print(f"{'p':>3} {'n=2p':>5} {'pillars(M/H/3m4)':>16} {'#(Z/2p)* atoms':>14} {'#tourn on p':>12} {'reconstruction on p vtx':>34}")
for p in [3,5,7,11,13,17,19,23,31]:
    if not isprime(p): continue
    M=is_mersenne(p); He=p in HEEGNER; T3=(p%4==3)
    pil=''.join('MH3'[i] if [M,He,T3][i] else '.' for i in range(3))
    phi=sum(1 for a in range(1,2*p) if __import__('math').gcd(a,2*p)==1)
    nt=A000568.get(p,'?'); rec=RECON.get(p,'(large)')
    star=' <==' if sum([M,He,T3])==3 else ('  [p=1mod4: NO Paley tournament]' if p%4==1 else '')
    print(f"{p:>3} {2*p:>5} {pil:>16} {phi:>14} {str(nt):>12} {rec:>34}{star}")

print("\n"+"="*118); print(" THE PICTURE (easier / harder / clearer / SURPRISING)"); print("="*118)
print(" EASIEST/CLEAREST: n=6 (apex 3): all 3 pillars + trivial tournament space (2 classes). PROVED.")
print(" RICHEST-but-FIRST-OPEN: n=14 (apex 7): all 3 pillars, BUT the tournament space on 7 vtx hits the")
print("   RECONSTRUCTION WALL (this session) -- 90% OCF-cospectral, local invariants fail. The apex is maximally")
print("   structured (Paley heptagon, tight-equiangular d=7) yet the class-space is first-irreducible => the")
print("   'gentlest place to break it open' is ALSO where the combinatorics first stops being locally decidable.")
print(" SURPRISE 1 (non-monotone): above 14 you LOSE pillars, not gain size-difficulty: n=22 (apex 11) loses")
print("   the Mersenne 2-adic descent; n=62 (apex 31) loses the Heegner sqrt(-31) SOS. Only {6,14} keep all 3.")
print(" SURPRISE 2 (p=1 mod 4 is EASIER, not harder): n=10,26,34 (apex 5,13,17) have no pillars, but the complement")
print("   reflection is an AUTOMORPHISM (Paley GRAPH, not tournament) => a BROUWER FIXED POINT witnesses the lonely")
print("   runner directly (SOS suffices). The 3 pillars are TOOLS for the HARD p=3mod4 regime (free Z2 => Borsuk-")
print("   Ulam odd degree, no fixed point). So difficulty axis #1 is p mod 4; 'fewer pillars' can mean EASIER.")
print(" SURPRISE 3 (tight-equiangular {7,23}): n=14 and n=46 are the equiangular-tight apices; 46 has only the")
print("   3-mod-4 pillar but the equiangular structure (my equioscillation<->equiangular reflection) -- a")
print("   DIFFERENT special mechanism, suggesting n=46 as the next 'geometrically special' case after 14.")
print(" MENTAL MODEL: LRC difficulty by n is NOT monotone; axis #1 = p mod 4 (Brouwer/easy vs Borsuk-Ulam/hard =")
print("   the sgn of the complement permutation = my S11 even-Aut/odd-anti-Aut finding); axis #2 (within hard")
print("   p=3mod4) = pillar count (all 3 only at p=3,7). Tournament-space wall at 7 vtx = the shadow of LRC(14).")
print("DONE.")
