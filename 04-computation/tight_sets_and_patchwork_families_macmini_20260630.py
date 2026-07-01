#!/usr/bin/env python3
"""
mac-mini-2026-06-30-S62
=======================
CLASSIFY (A) the NON-DIFFERENCE-CLOSED tight sets and (B) the PATCHWORK COVERING families,
and UNIFY the LRC tight locus with the covering-min via the "skip-and-patch" operation.

Tight set: M(S) = 1/n exactly (M(S) = max_t min_{v in S} ||v t||), the extremal LRC configs.
Difference-closed tight sets = dilated APs (THM-560 / opus HYP-3749). This script classifies
the NON-diff-closed residual and connects it to the covering-min construction {1..n-2, n(n-1)}.

Findings (verified below):
  A1. Single-swap tight census (n=5..16): non-diff-closed single-swap tight sets exist only at
      n=5 (skip2,add7), n=6 (skip2,add9), n=8 (skip6,add12), n=14 (skip12,add24).
  A2. GW-DOUBLING FAMILY: skip (n-2), add 2(n-2) is TIGHT  <=>  n = 2 (mod 6)  (6 | n-2).
      Verified n=5..26: tight exactly at n=8,14,20,26. The patch 2(n-2) is the smallest
      multiple of (n-2) outside the AP; it kills the resonance hole the skip opens at D=n-2,
      and tightness survives iff the Jacobsthal window 6|(n-2) holds.
  A3. Multi-swap sporadics exist (n=8 two-swap {1,4,5,6,7,11,13}); n=5,6 are small exceptions.
  B.  PATCHWORK UNIFICATION: both canonical LRC14 sets are single-element patchworks of AP{1..13}
      with the patch a MULTIPLE of the skipped speed (killing its resonance hole):
        - skip 12 (=n-2), patch 24 =2(n-2)     -> M=1/14   (LRC floor, tight, NON-covering)  [GW]
        - skip 13 (=n-1), patch 182=lcm(13,14) -> M=14/183 (covering-min, covering)  [construction]
      Among ALL skip+multiple-patch single-swaps, ONLY these two hit 1/14 and 14/183 respectively.
"""
from fractions import Fraction as F
from math import gcd, lcm
from functools import reduce
from itertools import combinations

def Mexact(S):
    Sg = sorted(set(S)); cand = set()
    for i in range(len(Sg)):
        for j in range(len(Sg)):
            for d in (Sg[i]-Sg[j], Sg[i]+Sg[j]):
                if d > 0:
                    for k in range(1, d): cand.add(F(k, d))
    best = F(0); at = None
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in Sg)
        if g > best: best = g; at = t
    return best, at

def covers(S, n): return all(any(v % q == 0 for v in S) for q in range(2, n+1))
def diff_closed(S):
    Sset = set(S); return all(abs(a-b) in Sset for a in S for b in S if a != b)
def prim(S): return reduce(gcd, S) == 1

# --------------------------------------------------------------------------------------
print("="*80)
print("A1. SINGLE-SWAP TIGHT CENSUS (AP {1..n-1}, remove one k, add one r<=4n): the sporadics")
print("="*80)
for n in range(5, 17):
    AP = list(range(1, n)); tgt = F(1, n); found = []
    for k in AP:
        base = [x for x in AP if x != k]
        for r in range(2, 4*n+1):
            if r in AP: continue
            S = tuple(sorted(base+[r]))
            if not prim(S): continue
            if Mexact(S)[0] == tgt: found.append((k, r))
    tag = "  AP-only (apex-clean)" if not found else ""
    desc = "; ".join(f"skip{k}+add{r}(r%n={r%n},{'=2k' if r==2*k else ''})" for k, r in found)
    print(f"  n={n:2d}: {len(found)} sporadic single-swap(s){tag}   {desc}")

print()
print("="*80)
print("A2. GW-DOUBLING FAMILY: skip (n-2), add 2(n-2).  TIGHT  <=>  n = 2 (mod 6)")
print("="*80)
for n in range(5, 27):
    k = n-2; r = 2*k
    S = sorted([x for x in range(1, n) if x != k]+[r])
    M, _ = Mexact(S); t = (M == F(1, n))
    print(f"  n={n:2d}: 6|(n-2)? {str((n-2)%6==0):5}  skip {k:2d} add {r:2d} -> M={str(M):8}"
          f"  {'TIGHT, ndc=%s, cover=%s'%(not diff_closed(S), covers(S,n)) if t else ''}")

print()
print("="*80)
print("A3. MULTI-SWAP sporadic at n=8 (two-swap), + n=5,6 small exceptions")
print("="*80)
n = 8; AP = list(range(1, n)); tgt = F(1, n); two = set()
for rem in combinations(AP, 2):
    base = [x for x in AP if x not in rem]
    for a, b in combinations(range(2, 30), 2):
        if a in AP or b in AP: continue
        S = tuple(sorted(base+[a, b]))
        if len(set(S)) != n-1 or not prim(S): continue
        if Mexact(S)[0] == tgt: two.add(S)
print(f"  n=8 two-swap tight sets: {sorted(two)}  (all ndc={all(not diff_closed(s) for s in two)})")
for n, S in ((5, (1, 3, 4, 7)), (6, (1, 3, 4, 5, 9))):
    print(f"  n={n} small sporadic {S}: M={Mexact(S)[0]} ndc={not diff_closed(S)} "
          f"(residues mod {n}: {sorted(x%n for x in S)})")

print()
print("="*80)
print("B. PATCHWORK UNIFICATION (n=14): skip s of AP{1..13}, patch = a MULTIPLE m*s of s.")
print("   The patch kills the resonance hole at D=s. Only two hit the critical M-values.")
print("="*80)
n = 14; AP = list(range(1, n)); floor = F(1, n); cmin = F(14, 183)
print(f"  floor 1/14={float(floor):.5f}   covering-min 14/183={float(cmin):.5f}")
print(f"  {'skip s':>6} | multiples m*s giving notable M")
for s in AP:
    base = [x for x in AP if x != s]; notes = []
    for m in range(2, 16):
        r = m*s
        if r in base or r > 500: continue
        S = tuple(sorted(base+[r]))
        if not prim(S): continue
        M, _ = Mexact(S); cv = covers(S, n)
        tag = 'TIGHT' if M == floor else ('COVMIN' if (cv and M == cmin) else ('cov' if cv else ''))
        if tag in ('TIGHT', 'COVMIN'): notes.append(f"{r}:{M} <<{tag}>>")
    if notes: print(f"  {s:>6} | {'  '.join(notes)}")
gw = sorted([x for x in AP if x != 12]+[24]); con = sorted([x for x in AP if x != 13]+[182])
print(f"\n  GW  = AP skip 12 (=n-2), patch 24 = 2(n-2)     : M={Mexact(gw)[0]} cover={covers(gw,n)}  [tight LRC floor]")
print(f"  CON = AP skip 13 (=n-1), patch 182=lcm(13,14)  : M={Mexact(con)[0]} cover={covers(con,n)}  [covering-min]")
print(f"  Both patches are multiples of the skipped speed: 24=2*12, 182=14*13.")
print("\nDONE.")
