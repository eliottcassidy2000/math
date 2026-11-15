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
