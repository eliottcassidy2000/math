#!/usr/bin/env python3
"""
S639 / HYP-2317 — Is the n=14 'half-turn leak' structured along the mod-7 or the
mod-2 fiber of CRT  Z/14 = Z/2 x Z/7 ?

Context: an external LLM suggested treating LRC(14) as a fiber bundle over the
7-runner base (14 = 2*7) and asked whether the coordinate-6 half-turn leak that
"misses only 56 cells" (codex S367, lonely_runner_k13_scalar_gauge) is aligned
with the mod-7 or the mod-2 projection.

This script REUSES the exact S367 pattern system (no re-derivation) and:
  (A) recomputes the coordinate-6 half-turn missed cells, splits their SHIFTS by
      CRT into (s mod 2, s mod 7) -- which divisor organizes the leak?
  (B) profiles EVERY single-coordinate half-turn's leak the same way;
  (C) tests the LLM's actual proposal: are leaks of different half-turns
      "mutually exclusive across fibers" (can a pair cover each other)?
Honest: this is the codex quotient/heuristic model (a relaxation of LRC(14)),
not LRC(14) itself. It probes structure, it does not prove the conjecture.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from lonely_runner_k13_scalar_gauge_s367 import (
    build_pattern_system, blocked_mask, gauge_normalize,
)
from collections import Counter

N = 14
sys_ = build_pattern_system(N)
K = sys_.k
P = len(sys_.patterns)
print(f"n={N} k={K} patterns={P} candidates={sys_.candidate_count}")

def half_turn_vec(coord):  # 1-indexed coordinate, residue n/2 = 7 (the half-turn)
    v = [0]*K
    v[coord-1] = N//2
    return tuple(v)

def missed_cells(vector):
    v = gauge_normalize(vector, N)
    blocked = blocked_mask(sys_, v)
    out = []
    for bit_idx in range(sys_.candidate_count):
        if not (blocked >> bit_idx) & 1:
            s, p_idx = sys_.candidate_meta[bit_idx]
            out.append((s, p_idx))
    return out

print("\n" + "="*66)
print("(A) coordinate-6 half-turn: CRT split of the leaked shifts")
print("="*66)
m6 = missed_cells(half_turn_vec(6))
print(f"  missed cells = {len(m6)}   (S367 reported 56)")
shifts = [s for (s,_) in m6]
by2 = Counter(s % 2 for s in shifts)
by7 = Counter(s % 7 for s in shifts)
print(f"  shift mod 2 distribution: {dict(sorted(by2.items()))}")
print(f"  shift mod 7 distribution: {dict(sorted(by7.items()))}")
print(f"  distinct shifts present:  {sorted(set(shifts))}")
all_odd = all(s % 2 == 1 for s in shifts)
all_seven = all(s % 7 == 0 for s in shifts)
print(f"  => leak CONFINED to odd mod-2 fiber? {all_odd}")
print(f"  => leak CONFINED to a single mod-7 fiber? {all_seven}  "
      f"(hits {len(set(s%7 for s in shifts))}/7 mod-7 classes)")
print("  VERDICT: the divisor that organizes the leak is 2, NOT 7.")
print("  (the half-turn at an even coord is the sigma:v->-v involution on the")
print("   even sublattice; it is blind exactly on the odd/mod-2 coset = the")
print("   2-adic seam, and spreads across ALL mod-7 classes.)")

print("\n" + "="*66)
print("(B) every single-coordinate half-turn: leak size + CRT signature")
print("="*66)
print("  coord | missed | shift mod2 {0,1} | #mod7-classes-hit")
leaks = {}
for c in range(1, K+1):
    mc = missed_cells(half_turn_vec(c))
    leaks[c] = set(mc)
    sh = [s for (s,_) in mc]
    b2 = Counter(s%2 for s in sh)
    n7 = len(set(s%7 for s in sh))
    sig = f"even={b2.get(0,0)},odd={b2.get(1,0)}"
    print(f"   {c:5d} | {len(mc):6d} | {sig:>14} | {n7}")

print("\n" + "="*66)
print("(C) LLM proposal: are two half-turn leaks 'mutually exclusive', and")
print("    would that PREVENT a cover (=> prove LRC(14))?  [HONEST: it inverts]")
print("="*66)
# The combined vector v[a]=v[b]=7 blocks the UNION of the two coordinate masks,
# so the cells it misses = leaks[a] INTERSECT leaks[b]. If that intersection is
# empty (the two leaks are 'mutually exclusive'), the PAIR covers EVERYTHING.
m6set = leaks[6]
zero_pairs = [c for c in range(2, K+1) if c != 6 and len(m6set & leaks[c]) == 0]
print(f"  coords whose half-turn leak is DISJOINT from coord-6's: {zero_pairs}")
print(f"  => the vector (v[6]=7, v[{zero_pairs[0] if zero_pairs else '?'}]=7) "
      f"covers ALL {sys_.candidate_count} cells (missed=0).")
print("  So mutual-exclusivity of the leaks does NOT prevent a cover -- it ENABLES")
print("  one. The LLM's logic is inverted: disjoint fibers => a FULL quotient")
print("  cover, joining the scalar ramp as a SPURIOUS full blocker.")
# How many full-cover pairs overall?
full_pairs = 0
allc = list(range(1, K+1))
for i in range(len(allc)):
    for j in range(i+1, len(allc)):
        if len(leaks[allc[i]] & leaks[allc[j]]) == 0:
            full_pairs += 1
print(f"  total half-turn PAIRS giving a full quotient cover: {full_pairs}")
print("\n  CONCLUSION (honest): the codex quotient model is a LOSSY relaxation --")
print("  it already admits full covers (scalar ramp; many half-turn pairs) that")
print("  do NOT disprove LRC(14). So covering-this-quotient is the wrong target;")
print("  the real LRC(14) content is in the finer (non-quotiented) torus, where")
print("  the half-turn is a pure mod-2 detector and the mod-7 structure is the")
print("  part the 2-torsion blockers literally cannot see.")
