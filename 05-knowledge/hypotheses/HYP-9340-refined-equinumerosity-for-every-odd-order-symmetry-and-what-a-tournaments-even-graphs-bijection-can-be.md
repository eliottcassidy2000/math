---
id: HYP-9340
title: "(a) For every odd-order subgroup H <= S_n, the tournaments with an H-symmetry and the even graphs with an H-symmetry are equinumerous up to isomorphism (proved for Sylow H by THM-4721; exhaustive for every odd-order class, n <= 8, and fifteen classes at n = 9). (b) What an exact bijection between tournaments and even graphs can be: the Walsh change of basis is the canonical linear one; equivariant set bijections exist and are canonical exactly for odd-order symmetry groups and exist (non-canonically) for AGL(1,p) on prime circulants; no S_n-equivariant one (THM-4531), conjecturally for every n >= 5 (HYP-9172); no uniform multiplier rule for prime circulants. (c) Sequence conjectures: odd graphs 0, 1, 2, 7, 22, 100, 588, 5466 (b-file to n = 60), labelled even graphs 1, 1, 4, 20, 396, 16432, 1312056, 200582824, connected odd graphs 0, 1, 1, 4, 13, 69, 455, 4704 are not in OEIS (2026-10-19); A334335 = connected even graphs extended to n = 128"
status: >
  (a) OPEN (NUMERICAL): true for every odd-order subgroup class of S_n for n <= 8 and for fifteen classes at n = 9,
  including Z3 wr Z3 of order 81, Z3^3, Z9 and four actions of Z3xZ3 (2954 = 2954 for the 3-cycle type) (THM-4721, subgroups_fix*.out:
  cyclic of every odd cycle type, Z3xZ3, Z7:Z3, Z15), PROVED when H is a Sylow p-subgroup (THM-4721 Theorem S) and
  at the level of Aut-classes of conjugates for every odd H (THM-4721 Proposition Q). The Sylow proof needs
  multiplicity one and the Frattini step, both of which fail separately for non-Sylow H (the tournament with
  Aut = Z3xZ3 on 6 vertices has two classes of 3-cycle modules; K_{3,3} is odd but N(Z3)-even) and cancel in total.
  (b) PROVED facts assembled (THM-4721 Theorems W, G', Corollary C; THM-4531 O3); OPEN: a uniform natural rule for
  the AGL(1,p)-equivariant circulant bijection; HYP-9172 for all n.
  (c) FINITE-EXACT (census n <= 8, Burnside formula for odd graphs to n = 60, inverse Euler transform of the repo's
  A000568 b-file to n = 128); OEIS lookups 2026-10-19 (JSON API, generic UA): odd graphs, labelled even graphs,
  connected odd graphs: no match; A334335 found with 19 terms and the comment "it appears ... connected even graphs
  (verified n <= 10)", now PROVED by THM-4721 Theorem A.
source: mac-mini-2026-10-19-evengraphs
related:
  - 01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md
  - 01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md
  - 05-knowledge/hypotheses/HYP-9172-no-natural-bijection-for-twisted-mallows-sloane-at-any-n-ge-5.md
  - 05-knowledge/hypotheses/HYP-9339-the-two-three-axis-numbers-fixed-points-and-relations.md
note: 05-knowledge/results/even_graphs_20261019_sylow_refinement_components_walsh_duality.md
---

# HYP-9340 — refined equinumerosity for every odd-order symmetry, and what an exact bijection can be

## (a) The refined conjecture

For `H <= S_n` of odd order let `tau_T(H)` (resp. `tau_E(H)`) be the number of isomorphism types of tournaments (resp. even graphs) on `n` vertices whose automorphism group contains a conjugate of `H`.

**Conjecture.** `tau_T(H) = tau_E(H)` for every `n` and every odd-order `H`.

Evidence: every odd-order subgroup class of `S_n`, `n <= 8`, and fifteen classes at `n = 9` (THM-4721; at `n = 9`: `Z_3` of type `(3,1^6)`: 2954/2954, `Z_3 x Z_3` with two fixed points: 86/86, regular `Z_3 x Z_3`: 1/1, `Z_3^3`: 2/2, `Z_3 wr Z_3`: 1/1, `Z_9`: 3/3). E.g. at `n = 8`: `Z_7`: 4/4, `Z_15`: 2/2, `Z_5`: 12/12, `Z_3` of type `(3,3,1,1)`: 48/48, `Z_3` of type `(3,1^5)`: 280/280, `Z_3 x Z_3`: 16/16, `Z_7 x| Z_3`: 2/2. Proved for Sylow `p`-subgroups (THM-4721 Theorem S), hence for all odd `H` of prime order `p` with `n < 2p`, and for prime circulants (`n = p`).

What is known for general odd `H` is the multiplicity-weighted form (THM-4721 Proposition Q): `sum_T m_T(H) = sum_{X even} m_X(H) + sum_{X odd, N(H)-even} m_X(H)`, `m_x(H)` = number of `Aut(x)`-classes of conjugates of `H` in `Aut(x)`. The conjecture says the two defects — tournaments with several classes of `H`-subgroups, and odd graphs all of whose odd automorphisms lie outside `N(H)` — always balance. At `n = 6`, `H = Z_3` of type `(3,1,1,1)` both defects are 1 (the two-module tournament; `K_{3,3}`). A proof would presumably pair these structures, or run Möbius inversion over the odd-order subgroup lattice; the reduction to a subgroup-lattice identity is recorded in the note (section 6).

Hostile control: the per-permutation version is false (THM-4721 section 6), so the conjecture is genuinely about orbits, not fixed points.

## (b) What an exact bijection can be

Royle et al. ask for "a reasonably natural bijection between the sets of unlabelled even graphs and tournaments". The repo's answer, assembled:

1. **The canonical object is linear.** The signed even-subgraph counts `Phi_X` (`X` even) and the type indicators are two `Z[1/2]`-bases of the same space (THM-4721 Theorem W), related by an invertible integer matrix with `det = +-2^k` and two orthogonality relations. Every isomorphism invariant of tournaments is a unique combination of signed even-subgraph counts.
2. **Equivariant set bijections exist exactly as the counting predicts.** For a symmetry group `Gamma` of odd order there is a canonical `Gamma`-equivariant bijection (reverse a `Gamma`-invariant tournament along the graph). For `AGL(1,p)` acting on prime circulants one exists but is not canonical (Hall on the lattice of odd multiplier subgroups). For `S_n` none exists for `5 <= n <= 15` (THM-4531), conjecturally for all `n >= 5` (HYP-9172).
3. **Open.** A uniform rule for the circulant case (the multiplier derivatives `D_c` fail for each fixed `c`; a rule would have to send the `K_max`-invariant tournaments to the empty graph for every `p`). A bijection respecting the refined counts of (a), i.e. mapping `H`-symmetric tournaments to `H`-symmetric even graphs for every odd `H` simultaneously, is the right target: by (a) the counts allow it, and any such map is automatically compatible with the Sylow structure.

## (c) Sequences (candidates for OEIS submission, with the THM-4720 and S48 lists)

| sequence | terms (n = 1, 2, ...) | status |
|---|---|---|
| odd graphs = A000088 − A000568 | 0, 1, 2, 7, 22, 100, 588, 5466, 83132, 2272112, 115244616, 10982861424 | not in OEIS; b-file to n = 60 (`odd_graphs_bfile.txt`) |
| labelled even graphs | 1, 1, 4, 20, 396, 16432, 1312056, 200582824 | not in OEIS; census n <= 8 |
| connected even graphs = A334335 | 1, 0, 1, 2, 8, 43, 398, 6413, 184596, 9540998, ... | interpretation now PROVED; b-file to n = 128 (`connected_even_graphs_bfile.txt`) |
| connected odd graphs = A001349 − A334335 | 0, 1, 1, 4, 13, 69, 455, 4704 | not in OEIS; census n <= 8 |
| strong tournaments | 1, 0, 1, 1, 6, 35, 353, 6008 | A051337; `EULER(A334335) = INVERT(A051337) = A000568` |
| even 4-, 6-coloured complete graphs | = THM-4720's `T_4`, `T_6` | THM-4721 Theorem K |
