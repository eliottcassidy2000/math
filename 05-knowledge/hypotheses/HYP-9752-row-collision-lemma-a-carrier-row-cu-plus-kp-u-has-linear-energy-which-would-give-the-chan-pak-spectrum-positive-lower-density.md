---
id: HYP-9752
title: "Row collision lemma: for a carrier row G(U) = cU + k p(U) + phi [U odd], the collisions with top split level l and common prefix w number at most C 2^l/(1 + w); hence its collision energy is O(2^L) and, by Cauchy-Schwarz, the Hamiltonian-path spectrum has positive lower density"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-growthbound, wave 29 of chanpak_waves_20261113.
  PROVED: the top-split reduction (each collision is determined by its suffix pair and one parity bit) and the implication
  lemma ⇒ linear energy ⇒ positive lower density.
  NUMERICAL (exact counts to L = 24): the energy per 2^L converges, to 0.2825 for 18U + 3p + 2[odd], 0.0053 for
  189U + 11p and 1.267 for 2U + p. At L = 20 the maximum of count·(1 + w)/2^l is 14 for the 3-cycle row and 6 for 2U + p.
  Concurrent (credited): the 7-21 session (commit 0cf4e0a, notes/06) measured one-hub energy ≈ 10.08·2^m and noted
  that E = O(2^m) gives positive density. chanpak_dyadic_collision_ports_20261115 proved the split identity and
  short-displacement bounds.
session: mac-mini-2026-11-15-growthbound
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 29)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w29_row_collisions.py, w29_row_collision_lemma_test.py
related:
  - HYP-9722 (carrier forms, density one), THM-5192 (A) (carrier form), chanpak_dyadic_collision_ports_20261115
---

# HYP-9752 - the row collision lemma

**Rows are counts.** Setting Y = 1 in a carrier form (THM-5192 (A)) gives a row

    G(U) = c U + k p(U) + φ·[U odd] + const,

and every value is the Hamiltonian-path count of an explicit tournament. For example, the 3-cycle core gives
18U + 3p(U) + 2[U odd] + 9.

**Top-split reduction (PROVED).** Let U ≠ U' < 2^L collide, with highest differing bit at position l − 1. Write
U = 2^l w + b and U' = 2^l w + b', with b ≥ 2^(l−1) > b'. The split identity
p(2^l w + b) = p(w) + p(b) + (2^l − 1 − rev_l b) w + [b < 2^(l−1)]·(w mod 2) (chanpak_dyadic_collision_ports) turns
G(U) = G(U') into

    k d w + k (w mod 2) = c (b − b') + k (p(b) − p(b')) + φ (b₀ − b'₀),   d = rev_l(b) − rev_l(b').

Here d is odd, hence nonzero. So each suffix pair (b, b') and parity bit yields at most one prefix w.

**Lemma (CONJECTURE).** Let N(l, w) be the number of ordered collision pairs with split level l and prefix w. Then

    N(l, w) ≤ C·2^l/(1 + w)

for a constant C depending on (c, k, φ).

**Consequences (PROVED from the lemma).**
- *Linear energy.* E(L) = Σ_l Σ_{w < 2^(L−l)} N(l, w) ≤ C Σ_l 2^l (1 + (L − l) ln 2) = O(2^L).
- *Positive density.* Cauchy–Schwarz gives #distinct values ≥ 4^L/(2^L + E) ≥ 2^L/(1 + O(1)). The values lie below
  (c + k)·2^L, so the row has positive lower density. These are tournament counts, so the Hamiltonian-path spectrum
  has positive lower density among the odd numbers. The previous unconditional lower bounds are of order X^(3/4).

**Why the lemma is plausible.** For large w the term −k w rev_l(b) dominates the suffix equation. A collision then
needs |rev_l b − rev_l b'| ≲ (c + 2k) 2^l/(k w), and among the about 4^l/w such pairs exact equality has density
about 2^(−l). The measured profile has exactly this shape: per-level totals grow linearly in the prefix length, and
count·(1 + w)/2^l stays bounded.

**Prediction P102.** For the 3-cycle row, max over (l, w) of N(l, w)(1 + w)/2^l stays ≤ 16 at L = 26, and
E(L)/2^L stays in [0.27, 0.29].

## Update (wave 30): the lemma as a twisted collision

Peel the top suffix bit: β = 2^(l−1) + x and β' = y, with x, y < 2^(l−1). The identities
rev_l(β) = 1 + 2 rev_(l−1)(x) and p(2^(l−1) + x) = p(x) + 2^(l−1) − 1 − rev_(l−1)(x) + [x < 2^(l−2)] turn a collision at
(l, w) into an equation between two twisted digital functions:

    A_w(x) − B_w(y) = κw + κ w₀ + κ − (c + κ)2^(l−1) − κ[x < 2^(l−2)] − φ(x₀ − y₀),
    A_w(x) = c x + κ p(x) − κ(2w + 1) rev(x),    B_w(y) = c y + κ p(y) − 2κw rev(y).

So N(l, w) ≤ C 2^l/(1 + w) says that the value sets of A_w and B_w, each of size 2^(l−1), meet a bounded shift in about
2^l/w places. That is the density a random pair of sets would have, since both are spread over a range of length
about κw·2^l.
- *Trivial bound.* Bounding rev(y) by the window forced by B_w gives 2^l/w solutions *per x*: a factor 2^l too weak.
- *What a proof needs.* Some equidistribution of the digital function c·rev(ρ) + κp(rev(ρ)) under the twist
  ρ ↦ −κwρ. That is the domain of the dyadic mixing and voltage-graph methods of the concurrent packets.

The lemma remains OPEN.

**Partial follow-up, 2026-11-15.** The [reverse-prefix theorem](../results/chanpak_local_rows_span3_20261115/ROWS.md) proves that low split levels and small reversed displacements contribute sublinear collision energy in the stated ranges. The [mesoscopic row limit](../results/chanpak_local_rows_span3_20261115/ROW_LIMIT.md) gives the continuous mass but does not control the distant exact collisions. The full hypothesis remains open.


**Further partial follow-up, 2026-11-15.** The [two-ended lattice reduction](../results/chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md) excludes every negative reversed displacement for the physical rows and proves that spans below $(1/(\log_2 6-1/3)-\epsilon)L$ contribute $o(2^L)$ by a slope-adapted digit split. Infinite positive-slope collision rays persist; the same-direction long-span estimate and linear energy remain open.

## Update (wave 31): long prefixes never collide (PROVED)

In the twisted form, write δ = rev(y) − rev(x). Then

    2κwδ = R₀ + κ rev(x) − [c(x − y) + κ(p(x) − p(y))],

with |R₀| ≤ κw + (c + κ)2^(l−1) + O(κ + φ) and the bracket bounded by (c + κ)2^(l−1) + O(1).
- *δ vanishes.* For w ≥ W₀(l) = ((2c + 3κ)/(2κ))·2^l + O(1) the right side is less than 2κw in absolute value, so
  δ = 0 and y = x.
- *No solution.* The equation then reduces to −κ rev(x) = κw + κw₀ + κ − (c + κ)2^(l−1) − κ[x < 2^(l−2)]. Its right
  side is positive once w > (c + κ)2^(l−1)/κ + O(1), so there is no solution.

Hence **N(l, w) = 0 for w ≥ W₀(l)**. Every collision is "balanced": U and U' share at most about half their leading
bits, so the row is injective on every block [2^l w, 2^l(w+1)) with w ≥ W₀(l), roughly on intervals of length √U.

Exact check (3-cycle row, L = 20): no collision has w ≥ 8·2^l. The lemma's remaining content is the balanced regime
w < W₀(l), where the measured collisions live (prefix lengths 4–7 bits at L = 20). There the claim is still
N(l, w) ≤ C 2^l/(1 + w).

**Sharpness, 2026-11-15 (klein-2026-11-15-rc, [THM-5232](../../01-canon/theorems/THM-5232-reverse-complement-collisions-every-odd-core-with-its-reverse-complement-collides-on-an-anti-diagonal-so-the-row-collision-bound-is-sharp-and-one-explicit-family-carries-linear-energy.md)).**
The bound shape is attained, so the lemma cannot be strengthened:
- p(rc x) = p(x) for the reverse-complement rc_a x = 2^a − 1 − rev_a x;
- every odd core collides with its reverse-complement along an anti-diagonal when k | c 2^(j−1);
- this gives N(l, w) ≥ (c/2k)·2^l/w on a positive-log-proportion set of prefixes in every dyadic block, so C ≥ 3
  for the 3-cycle row.

The same explicit family has linear size: 44% of the 3-cycle row's collisions at L = 18. A proof of the upper
bound should subtract such lattice-ray families and bound only the generic single-point lines. The lemma itself
remains OPEN.
