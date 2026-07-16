---
id: THM-921
title: THE IL_n SHELF CENSUS — exhaustive at n = 6 (2^15 codes): ground states (Q = 3 = Guy, 384 codes) and ONE shelf (Q = 4, 768 codes); the frustration profiles are EXACT integers by quadruple gap-class, and THE PALINDROMIC (REFLECTION-EVEN) CLASS (1,1,2,2) IS FRUSTRATED ONLY ON THE SHELF (3 per minimum) AND NEVER ON THE GROUND — the shelf/ground separation of the crossing-energy landscape is precisely "does the reflection-even quadruple class carry frustration": the owner's stratum conjecture confirmed in the sharpest form, tying the crossing landscape to THM-869's shelves and the residue-six fixed locus
status: exhaustive machine-exact at n = 6; profiles integral (ground: 2×(1,1,1,3) + 1×(1,2,1,2) per min; shelf: 3×(1,1,2,2) + 1×(1,2,1,2)); the alternating class (1,2,1,2) is the unavoidable core (1 per minimum at every energy). ε-TABLE NOTE (honest negative): the crude slice-at-N certificate is too lossy for small triples (values ≫ margins); the referee residual for residue six sharpens to: feed THM-920's lemma the TRUE per-triple λ₂ (computable from the lattice basis) instead of the generic floor — named, mechanical
source: mac-mini-2026-07-16-S122 (owner: run the IL_n shelf census against the reflection-even stratum; aim to complete LRC(14))
depends_on: [THM-920 (landscape formalization + the λ₂ lemma), THM-869 (the shelf template), boxeph T1545]
script: 04-computation/iln_shelf_census_epsilon_tables_macmini_S122.py -> 05-knowledge/results/iln_shelf_census_epsilon_tables_macmini_S122.out
---

# THM-921 — the shelf census

n = 6, exhaustive over 2^15 page codes: local minima at Q = 3 (= Guy Z(6); 384 codes)
and Q = 4 (768 codes) only. Frustration by cyclic gap-class of the quadruple:

| energy | (1,1,1,3) | (1,2,1,2) alternating | (1,1,2,2) PALINDROMIC |
|---|---|---|---|
| ground Q = 3 | 2 per min | 1 per min | **0** |
| shelf Q = 4 | 0 | 1 per min | **3 per min** |

The alternating class carries exactly one frustration at every minimum (the unavoidable
core); the ground spends its remaining budget on the asymmetric class (1,1,1,3); **the
shelf is the configuration that dumps its excess onto the reflection-even palindromic
class** — the crossing landscape's shelves are reflection-even-frustrated, exactly the
pattern of THM-869 (upset-saturated shelves) and the residue-six fixed locus (mass on
{s, 6−s}). Three landscapes, one law: **excess above the ground concentrates on the
involution-fixed stratum.** Next: n = 7 census (2^21, feasible), the shelf count sequence,
and whether the ground's (1,1,1,3)-profile predicts Guy's floor-product combinatorially.
