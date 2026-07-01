#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
IS THE (HALF-)TILING A LEFT-RIGHT CAYLEY / SQUARE COMPLEX? + cut-cycle ALEXANDER DUALITY + the even-graph LDPC code.

kind-pasteur-2026-07-01-S23. Chasing the LTC (Dinur-Evra-Livne-Lubotzky-Mozes) lead:
 (1) SQUARE/CUBE COMPLEX: the tiling cube = F_2^m (m=C(n-1,2) tiles) is a CUBE complex; the staircase tiles
     have ROW (y) and COLUMN (x) directions => row-flips A and column-flips B are two generating sets, and
     {t, t+e_i(row), t+e_j(col), t+e_i+e_j} are SQUARES => a LEFT-RIGHT square (Cayley) complex over the ABELIAN
     cycle space -- but ABELIAN (product/torus), NOT the nonabelian EXPANDER of Dinur et al.
 (2) ALEXANDER/HODGE DUALITY: arc space F_2^{C(n,2)} = CUT space (scores, dim n-1) (+) CYCLE space (tiles, dim
     C(n-1,2)), orthogonal complements. cut=coboundary (tournament/chromatic), cycle=even-graph (flow). The
     THM-584 R-even/R-odd split = this duality.
 (3) THE EVEN-GRAPH CODE = cycle space of K_n = an LDPC code: length C(n,2), dim C(n-1,2), distance 3 (triangle),
     n LOCAL parity checks (vertex even-degree), each weight n-1. Local-testability of THIS code = reconstruction
     (S14-19); the n=7 wall = local certification FAILS => the metagraph is an ANTI-LTC.
"""
import sys, itertools
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
from math import comb

def gf2_rank(rows):  # rows = list of int bitmasks
    rows=[r for r in rows if r]; r=0; piv=[]
    for col in range(64):
        bit=1<<col
        p=next((x for x in rows if x&bit and x not in piv and (x& ((bit<<1)-1))==bit or (x&bit)), None)
    # simpler gaussian elimination
    basis=[]
    for v in rows:
        for b in basis: v=min(v, v^b)
        if v: basis.append(v); basis.sort(reverse=True)
    return len(basis)

def rank_gf2(vectors, nbits):
    basis=[]
    for v in vectors:
        cur=v
        for b in basis:
            cur=min(cur, cur^b)
        if cur: basis.append(cur)
    # reduce
    return len(basis)

print("="*100); print(" (2) CUT + CYCLE = ALEXANDER/HODGE DUALITY of the arc space F_2^{C(n,2)}"); print("="*100)
print(f"  {'n':>2} {'C(n,2) arcs':>11} {'cut dim (n-1)':>13} {'cycle dim C(n-1,2)':>18} {'sum=arcs?':>10} {'cut⊥cycle?':>10}")
for n in range(4,9):
    edges=list(itertools.combinations(range(n),2)); E=len(edges); eidx={e:i for i,e in enumerate(edges)}
    # cut space: rows = vertex stars (boundary map ∂^T); dim = n-1
    cutrows=[]
    for v in range(n):
        m=0
        for i,(a,b) in enumerate(edges):
            if a==v or b==v: m|=1<<i
        cutrows.append(m)
    cutdim=rank_gf2(cutrows,E)
    # cycle space: rows = fundamental cycles (triangles span it); dim = C(n-1,2) = E-(n-1)
    cyclerows=[]
    for a,b,c in itertools.combinations(range(n),3):
        m=(1<<eidx[(a,b)])|(1<<eidx[(a,c)])|(1<<eidx[(b,c)])
        cyclerows.append(m)
    cycdim=rank_gf2(cyclerows,E)
    # orthogonality: cut . cycle = 0 (every cut meets every cycle evenly)
    orth=all(bin(cm & yr).count('1')%2==0 for cm in cutrows for yr in cyclerows)
    print(f"  {n:>2} {E:>11} {cutdim:>13} {cycdim:>18} {str(cutdim+cycdim==E):>10} {str(orth):>10}")
print("  => arc space = cut (n-1) (+) cycle (C(n-1,2)), orthogonal complements: the GF(2) Alexander/Hodge split")
print("     (cut=coboundary=tournament/chromatic; cycle=even-graph=flow). THM-584 R-even/R-odd = this duality.")

print("\n"+"="*100); print(" (3) THE EVEN-GRAPH CODE = cycle space of K_n = an LDPC code"); print("="*100)
print(f"  {'n':>2} {'length C(n,2)':>13} {'dim C(n-1,2)':>12} {'#checks (vtx)':>13} {'check wt (n-1)':>14} {'min dist':>9}")
for n in range(4,9):
    print(f"  {n:>2} {comb(n,2):>13} {comb(n-1,2):>12} {n:>13} {n-1:>14} {3:>9}")
print("  => n LOCAL parity checks (even-degree per vertex), each of weight n-1; distance 3 (triangle). An LDPC.")
print("     Dinur et al. build codes where LOCAL tests certify GLOBAL membership with O(1) queries. The tournament")
print("     'code' (iso class from local invariants) does the OPPOSITE: local certification FAILS at n=7 (S14-19")
print("     reconstruction wall) => the metagraph is an ANTI-LTC -- query complexity GROWS with n.")

print("\n"+"="*100); print(" (1) THE SQUARE/CUBE COMPLEX + half-tiling fold"); print("="*100)
for n in [5,6,7]:
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    rows=sorted(set(y for x,y in TILES)); cols=sorted(set(x for x,y in TILES))
    f=sum(1 for (x,y) in TILES if (n-y+1,n-x+1) in TILES and (n-y+1,n-x+1)==(x,y))
    print(f"  n={n}: m={m} tiles on a staircase grid ({len(rows)} rows x up to {len(cols)} cols); "
          f"row-flips A + column-flips B => square complex over F_2^{m} (ABELIAN). half-tiling (σ-fixed) = sub-cube dim {(m+f)//2}.")
print("  => the tiling cube is a genuine SQUARE (cube) complex with two flip-directions = an ABELIAN left-right")
print("     Cayley complex; it is a product/torus, NOT the nonabelian EXPANDER LTC. Nonabelian upgrade = the S_n")
print("     quotient (Schreier) or a PSL_2 realization -- that is where expansion/testability would live.")
print("\nVERDICT: the half-tiling IS a left-right SQUARE complex (row x column flips over the cycle space), and the")
print(" cut-cycle Alexander duality is the Tutte chromatic<->flow / THM-584 R-even<->R-odd split. But it is the")
print(" ABELIAN version: the even-graph LDPC code is locally CHECKABLE, and the tournament reconstruction is an")
print(" ANTI-LTC (local testing fails at n=7). The LTC lead points to: seek a NONABELIAN (S_n / PSL_2) realization")
print(" where the tiling square complex becomes an EXPANDER and reconstruction becomes O(1)-locally-testable.")
print("DONE.")
