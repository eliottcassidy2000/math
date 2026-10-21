---
id: HYP-9351
title: "After THM-4739: does second-order port information (scores plus curl profile of the partner) determine the merged metagraph for every n? NO already at n = 6 (L2 >= 200,000 candidates, L3 exactly 48,000); the information ladder L1 < L2 < L3 < ...; the remaining consecutive-port patterns of the (n-5,4) layer ({2,n-3,n-2}, {n-4,n-3,n-2}) and the epsilon = 0 ledger for k >= 4; electrical invariants of the metagraph (resistance from the transitive node, Kirchhoff index) and their laws"
status: OPEN (pre-audit); finite evidence stated per item
source: mac-mini-2026-10-21-electric
depends_on:
  - 01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md
related:
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
  - 05-knowledge/hypotheses/HYP-9347-three-move-law-layer-theorems-fundamental-black-graph-and-the-klein-dictionary-beyond-order-4.md
---

# HYP-9351 -- the information ladder of the metagraph

## (a) Second-order determination

THM-4739 (D): first-order port data (scores) determine the class-level metagraph at `n = 4`
but leave 4 candidates at `n = 5`; second-order port data (scores + curl profile of the
partner) force every port class at `n = 5`; at `n = 6` L2 leaves 153 free port classes (>= 200,000
solutions) and L3 (second scores added) leaves 37 free port classes and exactly 48,000
class-level solutions, 3 of which project to the real merged metagraph. **Claim (a1).** For each
`n` some finite level of the ladder (L4 = the orbit type of the full curl 2-form, L5 = the
switching class, ...) determines the metagraph; **(a2)** the level needed grows with `n`; the
first test is L4 at `n = 6` (does the full cyclic-triangle hypergraph of the partner force
all 403 port classes?). The 48,000 third-order impostors at `n = 6` are the explicit measure
of the information the lines carry beyond third order.
At the class level `(scores, curl profile)` is not a complete invariant (one collision at
`n = 5`, eleven groups at `n = 6`, 107 at `n = 7`), so any determination at L2 rests on the
port structure, not on class invariants: the question is how much of the fibre's internal
information the lines carry.

## (b) The epsilon = 0 ledger for k >= 4

THM-4739 (L)(i) is exhaustive for `k = 4, 5` at `n <= 12` and proved for `k = 3` by the packet.
**Claim.** The permutation argument (the product of the positional reversals of `Q1` and `P`
would be an automorphism with an orbit through `a` too small for `|I| = k - 1` in-neighbours)
proves the first word black for all `n >= k + 5` and all `k`; the second and third words follow
the packet's section 2 (`n > 4k - 2`) and section 3 witnesses.

## (c) The other consecutive-port patterns

`residual_cover.out`: for `I = {2, n-3, n-2}` (closed) the inherited `C` ledger covers the half
with `b -> y`; the residual half needs `5` (n = 9) or `7` (n = 10) certificate templates under
the greedy cover, the largest a block transport `(0,1,...,6,9,7,8), (4,5,6,9,7,3,0,1,2,8),
(4,5,6,7,3,0,1,2,8,9)`; for `I = {n-4, n-3, n-2}` (closed and open) the inherited `A`-cuts and
the mirrored-TT3 ledger cover `1920/2048` and the residual `128` need `3-7` templates. **Claim.**
Each residual splits by at most two further middle arcs into uniformly repaired sub-cases; the
first test is the separating-arc analysis restricted to the residual halves.

## (d) Electrical invariants

Resistance distance from the transitive node and the Kirchhoff index of `M_n` (`3, 34.8, 157.2,
2610.6` for `n = 4..7`). **Questions.** Is the resistance from the transitive node a function of
`(t_u, c3)` asymptotically? Does the Kirchhoff index have a closed form in `T(n), SC(n), 2^m`?
Which classes maximise the resistance (the "electrically remote" tournaments: at `n = 7` the
maximum is `2.08`)? Is the effective conductance between the two pure-blue nodes of THM-4730's
alphabet (`(1,1,1)` and the `T5`-centred stacks) a metagraph invariant with a law?

## (e) The determination problem of HYP-9345 (g), restated

The real metagraph is one of four L1-consistent multigraphs at `n = 5`. The three impostors
are honest electrical networks (same degrees, same closed/open bookkeeping, same c3 and V
potentials and currents) that no first-order law can exclude; the second-order law does.
Listing the impostors explicitly and the first `n` where even L3 leaves impostors is the
concrete next step.
