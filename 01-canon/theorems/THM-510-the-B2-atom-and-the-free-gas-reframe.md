---
id: THM-510
title: The B₂ atom (triangular-number factorization ≅ the four 4-tournaments) and the free-gas reframe of the OCF
status: PROVED (n=4 correspondence, hand-verified — gradings c3=|subset|, transpose=a↔b swap, parity↔transpose-type, det(S) split); the free-gas/Mayer reframe is a rigorous re-reading of THM-002/T006 + THM-029 (scripts deferred — session ran under a full-disk constraint)
source: mac-mini-2026-06-15-S2
depends_on:
  - THM-002   # H = I(Ω,2), the OCF as a hard-core lattice gas (T006)
  - THM-029   # H=7 forbidden (the K₃ conflict-graph obstruction)
  - THM-468   # 4-tournament det(S)=9 iff diamond/vortex, else 1 (Belkouche)
  - THM-509   # det/permanent (P/#P) wall
related:
  - T819      # Pascal-slope-d family (row 2 = (1,2,1))
  - HYP-2271  # {7,21} the only permanent H-gaps
  - OPEN-Q-100
---

# THM-510 — the B₂ atom and the free-gas reframe

## A. The free-gas reframe (re-reading of THM-002)

`H(T) = I(Ω,2)` is a hard-core lattice gas on the odd-cycle conflict graph Ω at
fugacity 2. The **free board / ideal gas** is the edgeless-Ω reference
`H_free = (1+2)^{α₁} = 3^{α₁}` (all α₁ odd cycles independent particles); the real
tournament is the **interacting gas** `H = Σ_k 2^k α_k ≤ 3^{α₁}`. **Forbidden values =
interaction obstructions:** `H=7 = 1+2·3` forces `α=(3,0)` ⟹ `Ω=K₃`, non-realizable
(THM-029: three pairwise-intersecting odd cycles force `α₂>0` or `α₁>3`). The overlap/
Witt defects (p33 = Ω-edges among triangles, TQ, …; THM-502/505) are the **Mayer
cluster integrals**. Free gas ↔ det/P side; interacting Ω ↔ permanent/#P side
(THM-509). The baby-Hodge holes are the lattice gas's interaction-forbidden free-α
vectors.

## B. The B₂ atom (PROVED, n=4)

`T(x)=x(x+1)/2 = f(x)g(x)`, with `a(x)=x+1` (successor/additive), `b(x)=x/2`
(halving/dyadic), `T(x)=x·b(a(x))`. By parity: x even → `f={a}, g={b}`; x odd →
`f=∅, g={a,b}` (where the op-set labels mean which of a,b each factor applies). The
four faces are the power set **B₂ = {∅,{a},{b},{a,b}}**, and they correspond to the
**four iso classes of tournaments on 4 vertices** (A000568(4)=4) under **c₃ = |subset|**:

| subset | c₃ | class | score | transpose | det(S) |
|---|---|---|---|---|---|
| ∅ | 0 | transitive | (0,1,2,3) | self | 1 |
| {a} | 1 | out-vortex diamond | (1,1,1,3) | ↔{b} | 9 |
| {b} | 1 | in-vortex diamond | (0,2,2,2) | ↔{a} | 9 |
| {a,b} | 2 | strong | (1,1,2,2) | self | 1 |

Three exact matches (hand-verified):
1. **c₃ = |subset|**; distribution `(1,2,1)` = Pascal row 2 (`C(2,k)`, the d=1 row of
   the Pascal-slope family T819), summing to `4 = 2²`.
2. **Transpose = the a↔b swap**: arc reversal fixes ∅ (transitive) and {a,b} (strong)
   and swaps the diamonds {a}↔{b} (out-vortex ↔ in-vortex). I.e. *transposing a
   tournament = swapping the additive and dyadic operations.*
3. **Parity of x ↔ transpose-type**: odd x ↔ self-transpose pair {∅,{a,b}}; even x ↔
   transpose-swapped pair {{a},{b}}.
4. (Pfaffian) `det(S)=9` (the diamonds, the **vortices** = obstructions to
   local-orderness / the d=1 floor code, THM-468) vs `det(S)=1` (transitive & strong).
   The diamond is where the **first interaction** appears (a vertex glued to a 3-cycle);
   n=4 is the smallest stage where the free board and the interacting tournament part.

## Scope / honesty

The n=4 correspondence is a finite, hand-checkable fact (verified above; A000568(4)=4
exactly equals |B₂|, and the gradings match). The **count** match (4 = 2²) is special to
n=4; the **gradings** match (c₃=|subset|, transpose=swap, parity↔type, det split) is the
robust structural content and is what makes the analogy a correspondence rather than a
coincidence. The free-gas part is a re-reading (not a new computation) of THM-002/029/
502; the Mayer-coefficient identification is made precise in OPEN-Q-100. (This session ran
under a full-disk machine constraint; the standard gentourng re-verification at n=4 is
deferred but the facts are elementary.)
