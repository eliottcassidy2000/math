# THM-455 — Transitive Subtournaments under Skew-Sylvester Doubling: the interleaving obstruction, and the Mersenne tower as an explicit Erdős–Moser family (STUB → filling this session)

**Status:** Parts (1)–(2) DERIVED (proofs below, computational verification running);
parts (3)–(4) PENDING branch-I computation. Honest stub per protocol.
**Source:** kind-pasteur-2026-06-09-S2 (T769, HYP-2356/2357).
**Related:** THM-447/448 (doubling, tower, link recursion), Erdős–Moser problem (1964),
mac-mini THM-453 Part A (592's red target = order-faithful transitive subtournament).

## (1) Lower bound: trans(D(T)) ≥ trans(T) + 1 (PROVED)

Let S = {v₁ → ⋯ → v_t} be transitive in T with source v₁. In D(T), append the source's twin:
v₁ → v₁' (twin arc), and v₁' → v_j for all j ≥ 2 (cross arc (i',j) follows T: v₁ beats v_j).
So v₁, v₁', v₂, …, v_t is transitive of size t+1. ∎

## (2) The interleaving obstruction (DERIVED — verification pending)

Let S ⊆ V(D(T)) be transitive, S₁ = unprimed part, S₂ = primed part. S₁ is transitive in T,
S₂ is transitive in T^op (so |S₁|, |S₂| ≤ trans(T), giving the trivial trans(D(T)) ≤ 2·trans(T)).
The cross arcs (both directions follow T, twins i → i') impose: in the transitive order of S,
each primed v_j' must sit exactly in the slot after v_j (relative to the unprimed chain), while
the primed elements among themselves must appear in REVERSED T-order. Consequence (derived for
chain-comparable configurations): if v_i', v_j' ∈ S with i ≺ j in the relevant transitive order,
then NO unprimed index in (i, j] can lie in S — primed elements must cluster outside the
unprimed span, except for a single twin splice at the boundary.

**For T transitive this forces trans(D(T)) = trans(T) + 1 exactly** (computation of the
max over (U, P) splits: |S| = n + 1 achieved by U = [1..k], P = {k..m}', and no more).

**Conjecture (testing now):** trans(D(T)) = trans(T) + 1 for ALL T — i.e. the doubling adds
exactly one transitive unit, the slowest possible growth. If true, iterated doubling from any
seed is an explicit family with trans = trans(seed) + #doublings ~ log₂ n.

## (3) The tower recursion (PENDING — branch I)

trans(T) = 1 + max_v trans(T[N⁺(v)]), and the tower's links include its own previous level
(B₀(T_{2m−1}) = T_{m−1}, THM-448(e)). Per tower level the floors give: trans(T7) = 3 (Paley T₇
is the unique largest TT₄-free tournament), trans(T15) ≥ 5, trans(T31) ≥ 6, trans(T63) ≥ 6.
Question: does the tower step (double + border) add exactly +2 (rate 2log₂ n, Paley-like) or
sometimes +1 (rate log₂ n, which would be remarkable)? Exact values being computed.

## (4) Comparison table (PENDING — branch I, with literature verification)

trans of: tower T7/T15/T31/T63, Paley_31 control, random controls; known extremal values for
TT_k-free tournaments (Sanchez-Flores etc., web-verified citations to be inserted).
