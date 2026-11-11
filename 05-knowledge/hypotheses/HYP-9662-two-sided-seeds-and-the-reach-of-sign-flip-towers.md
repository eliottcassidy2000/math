---
id: HYP-9662
title: "Two-sided seeds and the reach of the sign-flip towers. (A) Every odd c has a seed and an anti-seed: pairs (B, O) with offset +1, -1 and slope c 2^j, so that both lines 1 + c 2^s and c 2^s - 1 are eventually path counts. (B) Density form: the odd integers M that are the slope of a seed (or of an anti-seed) have positive lower density; this would give a positive lower density of Hamiltonian path counts, since 1 + 2M and 2M - 1 are path counts. (C) The towers over transitive pairs reach every odd c other than a finite set, for each sign. (D) A gadget (X, P) with kappa_X(P) = 0 has |P| = 2 and r_X(P) = -H(X); with H(X) = 1 it is a transitive pair at distance 1 or 2"
status: >
  OPEN. (A): the +1 half is HYP-9658 (A). FINITE-EXACT evidence: both hold for every odd c <= 4097 (seeds: THM-5126
  and odd_c_seed_symmetry_20261110, re-certified in THM-5130 (E); anti-seeds: THM-5130 (E), 2,049 certificates). (B): OPEN; the certificates have odd slope for 1,225 (seeds) and
  1,187 (anti-seeds) of the 2,049 odd c <= 4097, and at least 3,722 of the 4,095 odd N in [3, 8191] lie on a
  certified ray. (C): OPEN and uncertain; pure towers of depth <= 4 in the box of the pre-registration reach 1,612
  (with the two extra transitive starts: 1,649) of the 2,049 odd c <= 4097 with offset +1; at a fixed depth the
  reachable set is thin. (D): FINITE-EXACT for orders <= 7 (THM-5130 (G)); OPEN in general.
  None of these implies the Chan-Pak conjecture: rays of offset +1 or -1 do not exhaust the spectrum (15 lies on
  none, THM-5130).
session: kind-pasteur-2026-11-10-S22
source: 05-knowledge/results/kps_s22_sign_flip_towers_20261110/README.md
related:
  - 01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md
  - 05-knowledge/hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md
  - 05-knowledge/results/sequence_gf_20261110/SYNTHESIS.md (Session249: the Mahler source R_* of the conjectured spectrum)
---

# HYP-9662 - two-sided seeds and the reach of the sign-flip towers

Notation of THM-5130: a pair `(B, O)` has offset `r = H(B) - W_B(O)` and slope `M = U_B(O)`; the closed tournament
`B => TT_t => v`, `v -> O`, has `r + 2^t M` Hamiltonian paths. A seed has `r = 1`, an anti-seed `r = -1`.

## Statements

**(A) Two-sided seeds.** For every odd `c` there are a seed and an anti-seed whose slopes have odd part `c`.

**(B) Density form.** Let `Sigma_+` (`Sigma_-`) be the set of odd integers that are the slope of a seed (an
anti-seed). Then `Sigma_+` and `Sigma_-` have positive lower density.

Consequence if true: `1 + 2M` is a path count for `M` in `Sigma_+` and `2M - 1` for `M` in `Sigma_-`, so the
spectrum would have positive lower density. (The best proved lower bound is of order `X^(2/3)`, THM-4887.) The
conjectured spectrum, in the Mahler form of Session249, needs `Sigma_+` to contain every odd number except 3 and 5
whenever the point `1 + 2c` is to come from a ray; that is not required for Chan-Pak, because a point may be realised
off the rays.

**(C) Reach of pure towers.** For each sign, all but finitely many odd `c` occur as the odd part of the slope of a
tower over a transitive pair (THM-5130 (C)), at some depth.

**(D) Neutral gadgets.** If `kappa_X(P) = Pi_2(X) - 2 W_X(P) - 2 W3_X(P) = 0` then `|P| = 2` and
`W_X(P) = 2 H(X)`. If moreover `H(X) = 1`, then `P` is a pair at distance 1 or 2 of a transitive tournament.

## Evidence and scope

- (A): THM-5130 (E), odd `c <= 4097`, both signs; 1,969 seeds and 1,959 anti-seeds come from pairs of order at most
  9, the others from one or two tower levels over such pairs.
- (B): among the least-order certificates, 1,225 seeds and 1,187 anti-seeds have odd slope (`c <= 4097`); the
  certificates were chosen for least order, so these are lower bounds.
- (C): pre-registered run E1 (depth <= 4, `x <= 10`, `m <= 4`): `N_2 = 805`, `N_4 = 1612` of 2,049. Without a seed in
  that run: 5, 11, 23, 29, 73, 85, 91, 95, 127, ... (29 is reached from the transitive pair start that E1 omitted).
  Run E3 (depth <= 6, `x <= 5`, `m <= 2`): 1,592 of 2,049; 5, 73 and 85 appear at depth 6 (5 with slope `5 * 2^32`),
  11, 23, 91, 95, 125, 127 do not.
  Heuristic against a quick proof: at depth `k` the slope is a sum of boundedly many terms `2^x 3^y` with
  coefficients that are products of 3, 5 and 13, so each fixed depth gives a set of density zero; (C) can only hold
  through unbounded depth.
- (D): all gadgets with `X` of order at most 7 (neutral.py).

## Pre-registered tests for a later session

1. Odd `c` in `(4097, 16385]`: (A) predicts a seed and an anti-seed for each, found among pairs of order at most 10
   and at most two tower levels. A value of `c` with neither after pairs of order 11 are included would count
   against (A) in the form "seeds are plentiful", not against (A) itself.
2. (B): the proportion of odd `M <= X` that are slopes of seeds of order at most 12 does not fall below 0.5 for
   `X = 2^13, 2^14, 2^15`. A proportion falling like a power of `X` would refute the density form at these orders.
3. (D): order 8.

## What would refute

- (A): an odd `c` and a proof that no pair has offset `+1` (or `-1`) and slope `c 2^j`. The forbidden values 7 and
  21 already exclude the slopes 3; 5 and 10 (seeds) and 1, 2, 4; 11 (anti-seeds), so such obstructions exist for
  single slopes; none is known for a whole line.
- (D): a gadget with `kappa = 0` and `|P| != 2`.
