---
id: THM-4733
title: "The divisor plane of p^2 q r versus the marked-path tiling cube: the 3-4-3 and 4-2-4 pictures and their two reflections. (A) On the tiling cube {0,1}^m of tournaments with a marked Hamiltonian path (equivalently graphs on the ordered set [n-1], trijection_reset_20261020) the three natural involutions are the complement line kappa (graph complementation, the metagraph line K_P), the mirror sigma (order reflection = converse), and their product sigma kappa = path reversal R_P (reverse the path arcs only, read along the reversed path). (B) The Klein action (C, Q, CQ) = (divisor complement, prime exchange q <-> r, their product) on the ten proper divisors of p^2 q r, with rank levels 3-4-3 and sign levels 4-2-4, is EXACTLY the action (kappa, sigma, sigma kappa) on the four mirror-moved tiles of the order-five cube folded along one mirror orbit: D(p^2 q r) = {0,1}^4/(p_1 <-> p_2) with a = number of backward tiles in the folded orbit and (b, c) the bits of the other orbit; rank = Hamming weight, 4-2-4 = sign of (folded weight) - (unfolded weight). The same Klein set also embeds equivariantly as a ten-element invariant subset of that four-tile cube. (C) Obstruction: no equivariant injection of D(p^2 q r) into the full tiling cube of any order n >= 3 (nor into graphs on an ordered set with at least two points) exists, in either assignment of the two reflections, because a mirror-fixed tile exists: kappa is fixed-point free while Q fixes four divisors, and sigma kappa = path reversal is fixed-point free on genuine presentations while CQ fixes pq and pr. The divisor picture is therefore a quotient and a sub-object of the fixed-tile-free part of the cube, never of actual marked tournaments; the two central divisors pq, pr are the shadow of the forbidden path-reversal-symmetric presentations"
status: >
  PROVED (pre-audit; elementary): Proposition A (the three involutions of the cube and their tournament meanings),
  Theorem B (the equivariant fold and the equivariant embedding), Theorem C (the fixed-point obstruction).
  FINITE-EXACT: divisor_plane_check.out (the fold and embedding verified on the 16 patterns, the Klein orbit types,
  the obstruction count on the full n = 5 cube).
  CITED: trijection_reset_20261020 (README, arithmetic.md: the (f,g) = (a-b-c, b-c) chart, the Klein action with
  fixed counts 0, 4, 2, the tetrahedral obstruction), endpoint_three_repairs_20261020h (the trijection and K_P as graph
  complement), THM-643 (sigma = converse, f = floor((n-1)/2) fixed tiles), THM-4729 Proposition K (R_P = (K_P)^op),
  THM-4729 L1 (no antisymmetric tiling), procgen_tcpc_20261001 section 4.3 and collatz_finite_address_lattice_20261015c
  (the earlier 3-4-3 readings).
  NOT PROVED / not claimed: any isomorphism-preserving meaning of the fold (identified presentations are in general
  non-isomorphic tournaments); any detour, minor or colouring consequence; any arithmetic consequence for F = U + S.
source: mac-mini-2026-10-20-debtcones (owner's question: can the two reflections of the 3-4-3 / 4-2-4 picture connect explicitly to graph complementation and tournament path reversal; find a map preserving real structure or prove the obstruction)
depends_on:
  - 01-canon/theorems/THM-643-blue-black-line-parity-structure.md
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 05-knowledge/results/trijection_reset_20261020/arithmetic.md
related:
  - 05-knowledge/results/trijection_reset_20261020/README.md
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
  - 05-knowledge/results/procgen_tcpc_20261001_tournament_clock_prime_collatz.md
scripts:
  - 04-computation/experiments/endpoint_layers_20261020/divisor_plane_check.py
outputs:
  - 04-computation/experiments/endpoint_layers_20261020/divisor_plane_check.out
---

# THM-4733 -- the divisor plane is the fixed-tile-free shadow of the tiling cube

**PROVED (pre-audit) + FINITE-EXACT** as itemized in the status.

## 0. The two pictures

*Divisor plane* (trijection_reset_20261020/arithmetic.md). The ten proper nontrivial
divisors `d = p^a q^b r^c` (`0 <= a <= 2`, `b, c in {0,1}`, `(a,b,c)` not `(0,0,0)` or
`(2,1,1)`) carry the chart `(f, g) = (a-b-c, b-c)`. Rank `a+b+c = 1, 2, 3` gives the
levels `3-4-3`; the sign of `f` gives `4-2-4`. Three commuting involutions:
`C(a,b,c) = (2-a, 1-b, 1-c)` (divisor complement, `f -> -f`, `g -> -g`, no fixed point),
`Q(a,b,c) = (a, c, b)` (exchange of `q, r`, `g -> -g`, fixed `p, p^2, qr, pqr`),
`CQ(a,b,c) = (2-a, 1-c, 1-b)` (`f -> -f`, fixed `pq, pr`). Klein orbits:
`{p, pqr}`, `{p^2, qr}`, `{pq, pr}`, `{q, r, p^2 q, p^2 r}`.

*Tiling cube* (01-canon/definitions.md, THM-643). Tournaments on `[n]` containing the base
path `n -> ... -> 1` are the tilings `t in {0,1}^m`, `m = C(n-1,2)`; the trijection of
trijection_reset_20261020 identifies them with simple graphs on the ordered set `[n-1]`
(tile `(a,b)` = pair `{b, a-1}`, bit 1 = edge). Involutions: the complement `kappa`
(flip every tile; the metagraph line `K_P(T) = T^op xor E(P)`; graph complementation), the
mirror `sigma` (tile `(a,b) -> (n+1-b, n+1-a)`; order reflection of `[n-1]`), which fixes
`f = floor((n-1)/2) >= 1` tiles for `n >= 3`.

## 1. Proposition A: the three involutions and their tournament meanings

(i) `T_{sigma t}` is the converse `T_t^op` relabelled along the reversed path (THM-643,
MISTAKE-065). (ii) `T_{kappa t} = K_P(T_t)`, the complement line. (iii) `T_{sigma kappa t}`
is the **path reversal** `R_P(T_t) = T_t xor E(P)` (reverse the path arcs only), relabelled
along the reversed path. *Proof of (iii).* Reversing the path arcs of `T_t` and relabelling
by `i -> n+1-i` restores the base path; a non-path arc `a -> b` (`a > b`, bit 0 at tile
`(a,b)`) becomes `n+1-a -> n+1-b`, which runs from the smaller to the larger label of the
mirrored tile `(n+1-b, n+1-a)`, i.e. bit 1 there. So the new tiling is `sigma(t)^c =
sigma(kappa t)`. (THM-4729 Proposition K states the same as `R_P = (K_P)^op`.) `#`

So the Klein four-group `<kappa, sigma>` acting on presentations is exactly
`{identity, complement line, converse, path reversal}`. Its fixed points: `kappa` has none;
`sigma` fixes the blue (grid-symmetric) tilings; `sigma kappa` fixes nothing when `f >= 1`
(a mirror-fixed tile would have to equal its own complement; THM-4729 L1).

## 2. Theorem B: the divisor plane is the folded four-tile cube

Take `n = 5`: six tiles, `sigma` fixes `(4,2)` and `(5,1)` and moves
`(3,1) <-> (5,3)`, `(4,1) <-> (5,2)`. On the sub-cube `{0,1}^4` of the four moved tiles
write a pattern as `(p_1, p_2, q, r) = (t(3,1), t(5,3), t(4,1), t(5,2))`; there
`sigma(p_1,p_2,q,r) = (p_2,p_1,r,q)` and `kappa` flips all four bits. Define the **fold**

```text
Phi(p_1, p_2, q, r) = p^{p_1 + p_2} q^q r^r     (a = p_1 + p_2 in {0,1,2}, b = q, c = r).
```

`Phi` is a surjection of the 16 patterns onto the 12 divisors of `p^2 q r` (fibres of size
1 or 2, the size-2 fibres over `a = 1`), and it is equivariant:

```text
Phi(kappa x) = C(Phi x),      Phi(sigma x) = Q(Phi x),      Phi(sigma kappa x) = CQ(Phi x).
```

Indeed `kappa` sends `(p_1,p_2,q,r)` to `(1-p_1,1-p_2,1-q,1-r)`, so `a -> 2-a`, `b -> 1-b`,
`c -> 1-c`; `sigma` preserves `p_1 + p_2` and exchanges `q, r`. Removing the bottom
`(0,0,0,0)` and top `(1,1,1,1)` (weights 0 and 4) leaves the ten proper divisors with rank
= Hamming weight (`3-4-3` = weights `1, 2, 3` of the folded cube) and
`f = a - b - c` = (weight on the folded orbit) `-` (weight on the unfolded orbit), whose
sign gives `4-2-4`. The `CQ`-fixed divisors `pq, pr` are the images of the four
path-reversal-fixed patterns `p_1 + p_2 = 1, q + r = 1`.

**Embedding.** The same Klein set also sits inside the four-tile cube: the ten patterns
`0000, 1111, 1100, 0011, 1010, 0101, 1000, 0100, 0111, 1011` are a `<kappa,sigma>`-invariant
subset isomorphic to `D(p^2 q r)` as a Klein set (`p -> 0000`, `pqr -> 1111`,
`p^2 -> 1100`, `qr -> 0011`, `pq -> 1010`, `pr -> 0101`, `q -> 1000`, `r -> 0100`,
`p^2 r -> 0111`, `p^2 q -> 1011`). This is an abstract isomorphism of Klein sets (the
sigma-fixed patterns `0000, 1111, 1100, 0011` carry the four `Q`-fixed divisors, the
`sigma kappa`-fixed pair `1010, 0101` carries `pq, pr`, and one free orbit carries the rest);
it is NOT a section of `Phi` (`Phi(0000)` is the unit divisor), so the fold and the embedding
are two different equivariant comparisons. Both are checked pattern by pattern in
divisor_plane_check.out. `#`

Hence: **divisor complement = the complement line (graph complementation), prime exchange
= the converse (order reflection), and their product = tournament path reversal**, exactly,
on the mirror-moved part of the order-five cube.

## 3. Theorem C: the obstruction on genuine presentations

Let `X_n = {0,1}^m` be the full tiling cube of order `n >= 3` (or the graphs on an ordered
set with at least two points) with `(kappa, sigma)`. There is **no** injective map
`D(p^2 q r) -> X_n` intertwining `(C, Q)` with `(kappa, sigma)`, nor with `(sigma, kappa)`.

*Proof.* Equivariant maps send fixed points of an involution to fixed points of its image.
(1) If `C -> kappa`, `Q -> sigma`: then `CQ -> sigma kappa`, and the `CQ`-fixed divisors
`pq, pr` must land on `sigma kappa`-fixed tilings; none exist, since the mirror-fixed tile
`(n,1)` (or any of the `f >= 1` fixed tiles) would have to equal its own complement
(THM-4729 L1). (2) If `C -> sigma`, `Q -> kappa`: the four `Q`-fixed divisors must land on
`kappa`-fixed tilings, and `kappa` has none. In graph language: no graph on an ordered set
`[m]`, `m >= 2`, equals the complement of its order reflection, because the pairs
`{i, m+1-i}` are reflection-fixed. `#`

So the divisor picture can be matched with real structure only after the mirror-fixed
tiles are forgotten: it is the quotient (Theorem B) and a sub-object (the embedding) of the
**fixed-tile-free** part of the order-five cube. On actual marked tournaments, path reversal
has no fixed presentation, which is exactly why `pq, pr` have no counterpart: they are the
shadow of the forbidden path-reversal-symmetric presentations. The tetrahedral obstruction of
trijection_reset_20261020 (every tetrahedral involution fixes two pair letters, `C` fixes
nothing) is the same phenomenon seen from the geometric side.

## 4. What the fold does and does not preserve

The fold identifies the two patterns `(1,0,q,r)` and `(0,1,q,r)`, i.e. the presentations in
which exactly one of the mirror pair `(3,1), (5,3)` is backward. These are in general
non-isomorphic tournaments (the identification is by the count of backward tiles in one
mirror orbit, not by an automorphism), so the fold preserves the Klein action, the rank and
sign gradings, and nothing finer: no `H`, no score sequence, no detour ledger. This agrees
with the lane's verdict that "no divisor-to-tournament detour map is asserted"; the present
theorem states precisely which structure IS transported (the four-group with its gradings)
and why nothing on genuine presentations can be (Theorem C).

## 5. Reading the other 3-4-3

The lift note of endpoint_three_repairs_20261020h uses `3,4,3` for the ordered substitution
`C3[TT3, V, TT3]` (`n = 10`) and `4,2,4` for a stable-tree mark allocation. Those are
different objects; in the present theorem the integers are Hamming-weight levels of a folded
four-cube, and the common integer ten is `|D(p^2 q r)| = 12 - 2 = 2^4 - 4 - 2`
(sixteen patterns, four folded away, two extremes removed).

## 6. Files

`04-computation/experiments/endpoint_layers_20261020/divisor_plane_check.py` / `.out`:
the fold, its equivariance, the embedding, the Klein orbit types, and the obstruction counts
(`sigma kappa`-fixed tilings on the full `n = 5` cube: 0; on the four-tile cube: 4).
