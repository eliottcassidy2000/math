---
id: THM-4472
title: "The four-vertex reading of 3n+-1: the two converse diamonds (1,1,1,3)/(0,2,2,2) are the forward map and the inverse tree of one AM-fair pair (time reversal), not the two sheets; converse = negation o time reversal"
status: >
  PROVED + INDEPENDENTLY AUDITED.
  For the sheet 3n+b (b = +-1), take the AM-fair pair {s_o, s_e = s_o + b}
  (THM-4470) and its images T(s_o), T(s_e). The tournament on these four
  vertices, with the two map arcs x -> T(x) and every other pair oriented
  smaller -> larger, is (0,2,2,2), a 3-cycle over a sink with H = 3, iff
  b s_o > 0. Otherwise it is transitive with H = 1. So H = 2 + sgn(b s_o).
  The inverse-tree tournament (map arcs reversed) is the CONVERSE class,
  (1,1,1,3) resp. transitive. The reflection rho(x) = s_o + s_e - x, which
  preserves the vertex set because the pair sum is preserved, carries the
  forward tournament onto the converse of the backward one.
  Negation x -> -x fixes the class. Reversing all six arcs equals negation
  composed with time reversal. No local rule (map arcs plus any mix of order
  and fixed-label orientations; 512 rules) gives 3n+1 and 3n-1 the converse
  pair on a common side of 0.
  So the owner's "3 and +-1 = the two diamonds that swap when all arcs are
  inverted" is REAL as time reversal. 3n+b inverts to (y-b)/3, so the
  offset flips sign under inversion: forward is the sink, the inverse tree
  the source. It is REFUTED as the sheet swap. "3 = H(C3)" is an analogy:
  every odd multiplier q >= 3 gives the same diamond.
source: collatz-procgen-20260922 session, tournament lane (2026-09-24), from the owner's claim "the '3' and '+-1' is an underlying tournament on four vertices ... the 2 isomorphism classes which swap when all 6 arcs are inverted"; the coordinator's probe; audited and promoted by the session orchestrator 2026-09-24
depends_on:
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
related:
  - 01-canon/theorems/THM-001-redei
  - 01-canon/theorems/THM-002-ocf
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (transport theorem)
script: 04-computation/experiments/procgen_tourn_20260924_four_vertex.py
script_audit: 04-computation/experiments/procgen_tourn_20260924_orchestrator_check.py
output: 05-knowledge/results/procgen_tourn_20260924.out
output_audit: 05-knowledge/results/procgen_tourn_20260924_orchestrator_check.out
script_sha256: bb042ed26b3c3fc1d567395ba9532423c4c2cc87a662a5ac676e32fafa8186c9
script_audit_sha256: 1f3e9827777a85cb8e9b29bad4f0ee8913707d39707419a19a5cc2eeb0950b50
output_sha256: ebe55081f3f18391b34e089387d3fc84b794a0d6f8abff50bc5774a3d851987d
output_audit_sha256: 0d44d2ff44195d6b63b7867ff4c4845d22bc2657b668ffdc0f5fae4fb48c5960
hash_basis: raw LF bytes
audit: >
  The orchestrator's probe predated the lane and agrees with it.
  Independent orchestrator code (written without reading the lane's
  scripts) checks 232 quadruples over both sheets and both sides of 0: the
  forward class, H, the converse relation to the backward class, and the
  reflection identity rho(F) = converse(B). The lane's full pipeline was
  re-run, and its output is byte-identical. The 512-rule census
  (Proposition B) was not re-implemented; its conclusion follows for the
  intrinsic constructions from the exact classification above.
---

# THM-4472 -- the four-vertex reading of `3n+-1`

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_tourn_20260924_four_vertex_sheets_redei](../../05-knowledge/results/procgen_tourn_20260924_four_vertex_sheets_redei.md).

## 1. Statement

**The four tournaments on four vertices** (score sequence, number of Hamiltonian paths `H`):

| class | description | `H` | under converse |
|---|---|---|---|
| TT `(0,1,2,3)` | transitive | 1 | fixed |
| `(0,2,2,2)` | 3-cycle over a sink | 3 | swapped with `(1,1,1,3)` |
| `(1,1,1,3)` | source over a 3-cycle | 3 | swapped with `(0,2,2,2)` |
| `(1,1,2,2)` | strong | 5 | fixed |

**Theorem.** Let `T_b(n) = n/2` (n even) and `(3n+b)/2` (n odd), `b = +-1`, and take the pair `{s_o, s_e}` with `s_o` odd and `s_e = s_o + b`. This is the pair whose sum `T_b` preserves (THM-4470). Put `V = {s_o, s_e, T_b(s_o), T_b(s_e)}` (four distinct points once `|s_o| >= 5`).
1. **Forward.** Take the map arcs `s_o -> T(s_o)` and `s_e -> T(s_e)`, and orient every other pair smaller -> larger.
   * If `b s_o > 0`: the class is `(0,2,2,2)` and `H = 3`.
   * If `b s_o < 0`: the class is transitive and `H = 1`.
   * Hence `H = 2 + sgn(b s_o)`.
2. **Backward (inverse tree).** Reversing the map arcs gives the converse class: `(1,1,1,3)` or transitive.
3. **Reflection.** `rho(x) = s_o + s_e - x` preserves `V`, because `T(s_o) + T(s_e) = s_o + s_e`, and `rho(F) = converse(B)`.
4. **Negation.** `x -> -x` sends `(3n+b, s_o)` to `(3n-b, -s_o)`. It keeps `b s_o`, hence the class. Reversing all six arcs is exactly negation composed with time reversal.

## 2. Proof sketch

Take `3n+1`, `i >= 2`, `V = {i < 2i-1 < 2i < 3i-1}`.
* **Forward.** The arcs are `2i-1 -> 3i-1` and `2i -> i`, plus the order arcs. So `i -> 2i-1 -> 2i -> i` is a 3-cycle and `3i-1` is a sink.
* **Backward.** The arcs become `3i-1 -> 2i-1` and `i -> 2i`, plus the order arcs. So `i` is a source over the 3-cycle `2i-1 -> 2i -> 3i-1 -> 2i-1`.
* **The other cases.** The negative pairs and the `3n-1` pairs follow by negation (item 4), since negation reverses exactly the order arcs.
* **Reflection.** `rho` reverses the order and swaps `s_o <-> s_e` and `T(s_o) <-> T(s_e)`, so it keeps the map arcs and flips the order arcs. Composing with the backward construction, which flips only the map arcs, gives the full converse. ∎

## 3. Reading

* **The owner's two converse classes are the two directions of time.**
  * `3n+b` inverts to `(y-b)/3`: the offset changes sign under inversion.
  * The forward map is "3-cycle over a sink" and the inverse tree is "source over a 3-cycle".
  * They coincide with the owner's "hidden choice" objects of the tree programme: `(N-1)/3` is the backward `-1`.
* **The sheet swap does not exchange them.**
  * `3n+1` and `3n-1` on a common side of 0 give tournaments differing in at most one arc: one diamond and one transitive.
  * The diamond appears exactly on the side where `b s_o > 0`, i.e. the side carrying the contracting 2-cycle, where the sign law makes positive cycles contracting.
* **"3 = H(C3)" is an analogy.** Every odd multiplier `q >= 3` gives the same diamond.
* **Redei's parity has no Collatz counterpart.**
  * Redei's `H` is odd because the OCF involution has exactly one fixed point.
  * The 2-adic periodic points of `T` (`2^p` per period, the Banach fixed points of THM-4471 §4) are paired by the **free** involution `x_w <-> x_(w-bar)`, which commutes with `T`.
  * The two parity mechanisms are opposite.
