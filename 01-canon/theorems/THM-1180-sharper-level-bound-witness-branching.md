---
id: THM-1180
title: THE SHARPER PER-LEVEL BOUND AND THE WITNESS-POINT BRANCHING RULE — three exact sharpenings that revive the substitution search after k=5 stalled: (1) LOCALITY, only arcs meeting the residue matter, so a component [a,b] needs only indices ⌈wa−λ⌉…⌊wb+λ⌋, O(w·μ(E)+c) instead of O(w); (2) THE j=1 ROUNDING TEST, at the last level each component must sit in a SINGLE arc by the separation lemma, so containment is ⌈wb−λ⌉ ≤ ⌊wa+λ⌋ per component with NO subtraction — 360× faster at leaves; and decisively (3) WITNESS-POINT BRANCHING, since every cover must cover every point, fixing a witness x in the longest component forces some placed comb to satisfy ‖wx‖ < λ, a condition of density 2λ = 1/7, cutting the branching factor sevenfold per level. Validated against the known answers: k=3 gives 8156 nodes (was 87863) and k=4 gives 217178 (was 2242028), a 10× reduction, both returning exactly the two known families
status: all three sharpenings verified exact — locality and the rounding test each 0/360 mismatches against the reference, and the branching rule reproduces the S380 and S381 results exactly. The k=5 enumeration is RUNNING on the validated method at close-out and its verdict is NOT yet in; nothing here claims k=5
source: opus-2026-07-17-S385 (owner: stop the run and work the sharper per-level bound)
depends_on: [THM-1125 (the separation lemma behind the rounding test), THM-1155/1165 (the k=3,4 answers used as validation), THM-1170 (the 1/7 density that powers the branching rule)]
scripts: 04-computation/sharper_level_bound_opus_S385.py, witness_branch_opus_S385.py -> 05-knowledge/results/
---

# THM-1180 — why k=5 stalled, and what fixes it

## The diagnosis

S382's k=5 run stalled, and my S381 projection ("~56M nodes, slow but
feasible") was wrong for a reason worth recording: **the cost is
back-loaded**. Quintuples enumerated lexicographically start by dropping
{1,2,…} — removing the *coarse* combs leaves large essential regions, hence
small speed bounds and small subtrees. The last quintuples drop {9,…,13},
**keeping** speeds 1,2,3, so residues are small, ℓ_max is tiny, the bound
j/(r(7−j))/ℓ_max blows up, and subtrees explode. The projection assumed
uniform per-quintuple cost.

## Three exact sharpenings

**(1) Locality.** Only arcs meeting E matter: for a component [a,b] the
relevant indices are ⌈wa−λ⌉ … ⌊wb+λ⌋, so the work is O(w·μ(E) + c) rather
than O(w). Verified **0/360** mismatches against the reference.

**(2) The j=1 rounding test.** At the last level E ⊆ D_w, and by the
separation lemma (THM-1125) each component must lie in a *single* arc. So
containment holds iff for every component

> **⌈wb − λ⌉ ≤ ⌊wa + λ⌋**

an O(1) test per component with **no subtraction at all**. Verified
**0/360**; measured **360× faster** than the old leaf check.

**(3) Witness-point branching — the decisive one.** Every cover must cover
every point. Fix a witness x (midpoint of the longest component, the most
constrained). Then some placed comb satisfies ‖wx‖ < λ — a condition of
density **2λ = 1/7**. Branching on "which comb covers x" explores a seventh
of the candidates per level. It drops sorted-order symmetry breaking, so
families are found up to j! times and deduped; net **7ʲ/j!**, which at
j = 5 is ≈ 140.

## Validation before trust

| k | old nodes | new nodes | families |
|---|---|---|---|
| 3 | 87863 | **8156** | 2 ✓ (matches THM-1155) |
| 4 | 2242028 | **217178** | 2 ✓ (matches THM-1165) |

Both reproduce the known answers exactly at ~10× fewer nodes. This check
was the point: a faster search that changed an answer would be worthless,
and validating on levels whose truth is already established is the only way
to earn confidence in the level whose truth is not.

## What this does and does not claim

It does **not** claim k = 5. That enumeration was running at close-out and
its verdict is not in.

What is worth noting conceptually: my earlier attempts to rescue k = 5 were
about implementation. What actually moved it was a mathematical
observation — the **1/7 density of point-covering** — which is the same
constant 2λ = 1/7 that has obstructed every union bound in this programme
(S₁ = 13/7 > 1, the |W_q| window, the k = 7 arity ceiling). Here it works
*for* the search rather than against it, which is the first time that
constant has been an asset.
