---
id: THM-4818
title: "Anti-automorphism cycles have length 2 (mod 4); symmetry-first enumeration of the order-11 symmetric primes with two exact validations (cross-check for the reserved THM-4816/4817 ledger); a census of the dipole law's potentials on every class of orders 4-8: arc connectivity is exact, the singleton distance exceeds delta^0 exactly when the condensation leaves no eligible endpoint (a rigid non-self-converse failure exists at order 8), the score potential L has generic local minima for the transitive target and no linear negation-invariant potential is exact there; the class distance from P_7 to the transitive class is 6 = L(P_7); L is the A_{n-1} root-word metric on sorted scores; every regular 11-tournament is exactly five all-black moves from a singleton plus a strong core"
status: >
  PROVED (elementary; klein-2026-10-27, not independently audited): (A) the 2 mod 4 law; (C2) the condensation
  criterion in the direction "no eligible vertex => distance > delta^0"; (C4) the LP dominates every linear
  negation-invariant potential bound; (D) the root-word reading and the two-sided sorted-score bound for an arbitrary
  target. FINITE-EXACT: (B) the 270,999 symmetric prime classes of order 11 with their twelve cycle-index types, passing
  two exact identities (self-converse primes 259,307 = the value forced by OEIS A002785(11); sum over non-rigid primes
  of (1 - 1/|Aut|) = 1337728/165 = P(11) - P_lab(11)/11!); the order-12 ledger of the packet given four OEIS inputs
  (CITED: A000568(12), A002785(12), A051337(12), A054946(12)); (C) the class move graphs of orders 4-8 (12, 56, 456,
  6,880 classes; every directed Hamiltonian path of every class) with exact distances to the non-strong, singleton
  and transitive loci; (C2)'s converse at orders <= 8; (C4) the LP at orders 5 and 6; (D) the labelled BFS over all
  2^21 tournaments on 7 vertices; (E) all 1,223 regular 11-tournaments. CITED: THM-4811 (A)-(D) (the dipole law and the
  lambda / delta^0 exactness for strong sources), THM-4810 Theorem L, THM-4814 (A) (the automorphism lemma and the
  species identities), THM-4803, Camion 1959, OEIS as listed. RESERVED ELSEWHERE: the order-11 type table and the
  order-12 ledger as theorems are mac-mini-2026-10-27-dipole's THM-4816 / THM-4817; this entry records the method,
  the validations and the data as an independent cross-check and takes no credit for the ledger.
  CORRECTION (2026-10-27, mac-mini-2026-10-27-dipole, MISTAKE-683): in (B) the exact values are
  P_lab(12) = 68,487,766,269,146,864,640 and l(12) = 35,524,962,710,023,558,019 (the printed values are double-precision
  roundings; THM-4817); in (C2) the order-8 class S_4 => S_4 (scores 6,6,5,5,1,2,2,1) is rigid but self-converse.
session: klein-2026-10-27 (codex/session-overnight-20261027)
source: 05-knowledge/results/order11_prime_types_n12_ledger_dipole_census_20261027/README.md
depends_on:
  - 01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md
  - 01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md
  - 01-canon/theorems/THM-4814-order-ten-prime-cycle-index-types-order-eleven-ledger-skew-hadamard-nuclear-norm-extremality-and-black-core-distance-of-regular-tournaments.md
related:
  - 05-knowledge/hypotheses/HYP-9424-singleton-distance-is-delta0-iff-an-eligible-endpoint-exists-linear-potentials-never-exact-for-the-transitive-target-and-the-order-12-types.md
  - 05-knowledge/hypotheses/HYP-9422-regular-tournaments-exactly-m-black-moves-from-a-singleton-plus-strong-core-nuclear-norm-maximizers-off-the-hadamard-orders-and-the-order-11-types.md
  - 05-knowledge/results/external_preprints_mub6_local_smoothing_eth_merge_20261027.md
  - 07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md
scripts:
  - 05-knowledge/results/order11_prime_types_n12_ledger_dipole_census_20261027/scripts/ (gen_sym11.py, anti_types_check.py, tourn_types_dedup.c, species11_fast.py, species12b.py, movegraph.c, analyze_mg.py, singleton_failures.py, lp_dipole.py, lbfs7.c, blackcore.py, run.sh)
---

# THM-4818 — the 2 mod 4 law, the symmetric primes of order 11, and what the dipole law's potentials do on every class of orders 4–8

Conventions as in THM-4811: `T` a tournament on `n` vertices, `s(v)` out-degrees, `d = s − (n−1)/2`, `K` the skew
matrix, a legal word `Q = (a, …, b)` a directed Hamiltonian path, the native move `τ_Q T` keeps the arcs of `Q` and
reverses every other pair, so `K' = −K + 2A_Q` and `d' = −d + e_a − e_b` (the dipole law). Distances are on isomorphism
classes unless "labelled" is said.

## (A) Anti-automorphism cycles have length 2 (mod 4)

Let `ψ` be an anti-automorphism (`u → v ⟺ ψ(v) → ψ(u)`). Applying `ψ` `k` times to an arc `u → v` gives
`ψ^k(u) → ψ^k(v)` for even `k` and `ψ^k(v) → ψ^k(u)` for odd `k`. If `ψ^k` swaps `u` and `v` with `k` even, the arc
`u → v` forces `v → u`: impossible. In a cycle of length `2m` the pair `{x, ψ^m x}` is swapped by `ψ^m`, so `m` is
odd. Two fixed points `u, v` give `u → v ⟹ v → u`; a cycle of odd length `≥ 3` makes an odd power of `ψ` fix two
points of the cycle while acting as an anti-automorphism. Hence: **all cycles of `ψ` have length `≡ 2 (mod 4)`, with
at most one fixed point.** (THM-4814 stated only "all cycles even". The self-complementary-graph analogue, lengths
`≡ 0 (mod 4)`, differs because a complement-reversing step returns the pair after an odd number of steps.) At order
11 the admissible types are `2^5 1`, `6 2^2 1`, `10 1`; `scripts/anti_types_check.py` verifies the law by orbit
propagation for every partition of every `n ≤ 12`.

## (B) The symmetric primes of order 11 by invariant enumeration (cross-check)

By THM-4814's lemma a prime tournament of order 11 with nontrivial `Aut` has an automorphism of type `3^2 1^5`,
`3^3 1^2`, `5^2 1` or `11`; by (A) a self-converse one has an anti-automorphism of type `2^5 1`, `6 2^2 1` or
`10 1`. Enumerating the invariant tournaments of these seven permutations (one bit per orbit on pairs; the free
five-point block from the 12 classes of 5-tournaments; 2-cycles of `ψ` oriented `p → q`) gives 1,316,960 candidates,
which the deduplicating nauty census reduces in 7 s to **270,999 distinct symmetric prime classes** (259,307
self-converse; |Aut| = 3: 12,128; 5: 23; 11: 3; 55: 1), in twelve cycle-index types (packet §1). Two exact
validations hold: the self-converse prime count equals the value forced by A002785(11) through THM-4814's
decomposition identity, and `Σ_{|Aut|>1} (1 − 1/|Aut|) = 1337728/165 = P(11) − P_lab(11)/11!` with
`P_lab(11)` from the labelled species identity. The four classes with an 11-cycle are the circulants on `Z/11`
(`(32 + 8)/10 = 4` by Burnside), all self-converse. The resulting order-12 ledger (packet §2) gives, with the four
OEIS inputs, `P(12) = 142,980,435,791` (A000568(12) and A051337(12) agree on the non-strong count 1,798,161,433),
1,170,949 self-converse primes, `c_s(12) = 76,154,750,809`, `c(12) = 76,601,960,166`,
`P_lab(12) = 68,487,766,269,146,865,664`, `ℓ(12) = 35,524,962,710,023,558,531`. The theorem-level statements of the
table and ledger are THM-4816 / THM-4817 (reserved); the method and the two validations are this entry's.

## (C) The potentials on every class of orders 4–8

The class move graph (every legal word of every class, canonicalised) is undirected because the move is an
involution. Exact BFS distances versus the potentials:

1. **`λ` is exact** for the non-strong locus at every class of orders 4–8, both colours and black-only (THM-4811
   (B), (C) on classes; `R5` is the one black-unreachable class).
2. **Condensation criterion for the singleton locus.** A legal word starts in the initial strong component and
   ends in the terminal one. By the dipole law `|d(w)|` grows only when `w` is the start and `s(w) ≤ (n−1)/2`, or the
   end and `s(w) ≥ (n−1)/2`. Call `w` *eligible* if `s(w) = δ⁰` and `w` lies in the initial component, or
   `n−1−s(w) = δ⁰` and `w` lies in the terminal one. **If no vertex is eligible the singleton distance exceeds
   `δ⁰`** (the first move cannot increase any `|d|` from its maximum `(n−1)/2 − δ⁰`, and a vertex below the
   maximum needs more than `δ⁰` moves). At orders ≤ 8 the converse holds: the 0, 1, 3, 16 failures at `n = 5, 6, 7,
   8` are exactly the classes without an eligible vertex, all non-strong, all at distance `δ⁰ + 1 = 2`; among them
   the rigid, non-self-converse class with scores `6,6,5,5,1,2,2,1` (components 4 + 4). Exactness fails through
   forced structure, not forced symmetry; for strong sources every vertex is eligible and THM-4811 (D) applies.
3. **The transitive target.** `L = ½Σ|s_(i) − i|` is exact on 6/12, 23/56, 274/456, 5,449/6,880 classes; its
   local minima (classes with no `L`-decreasing move) number 22, 110, 904 at `n = 6, 7, 8`, mostly rigid
   non-self-converse (748 of 904 at `n = 8`); the distance exceeds `L` by up to 3 at `n = 7`. The transitive class
   has one legal word, blue, so it has exactly one neighbour and infinite black distance from everything else.
4. **No linear negation-invariant potential is exact for the transitive target.** A word of length `ℓ` gives
   `(K' − (−1)^ℓ K)/2 = Σ_i (−1)^{ℓ−i} A_{Q_i}`; dropping integrality and the requirement that `Q_i` be a path of
   the current tournament leaves the LP `min Σ_Q c_Q` subject to `Σ_Q c_Q A_Q = D`, `c ≥ 0`, over the `n!` directed
   paths of `K_n`. Its dual feasible points are the skew weight matrices `W` with `|⟨W, A_Q⟩| ≤ 1` for all `Q`, i.e.
   exactly the linear potentials with dipole-Lipschitz constant 1 (scores, cuts, sorted scores, and the linear
   functionals whose suprema give `λ`, `L` and the trace-norm tariff), so the LP bound dominates every bound of that
   kind. Minimised over labellings and rounded to the parity of `ℓ`, it is exact at `n = 5` except on the three
   classes with `|Aut| = 3` and has gaps on rigid non-self-converse classes at `n = 6`. The binding constraint is
   that the kept path is directed in the current tournament.
5. **Blue moves are not loops on classes** (32/32, 88/96, 780/780, 4,500/4,512 blue moves change the class at
   `n = 5, 6, 7, 8`).

## (D) P_7 is six moves from the transitive class; `L` is the root-word metric

The three regular 7-tournaments are at class distance 6 = `L` from `TT_7`. The labelled BFS over all 2^21
tournaments from the standard `P_7` puts 2,520 of the 5,040 transitive labellings at distance 6 and 2,520 at
distance 7 (the `i < j` labelling among the latter); so THM-4810's "P_7 at distance 7" is a labelled statement and
`L` is exact at `P_7` (MISTAKE-682 for the reflection's inference). Reading: `d` is a sum-zero vector, the dipole
`e_a − e_b` a root of `A_{n−1}`, a word of length `ℓ` writes `d_ℓ − (−1)^ℓ d_0` as a signed sum of `ℓ` roots, and
the least number of roots summing to a sum-zero integer vector `v` is `‖v‖_1/2`; minimising over labellings is the
sorted matching. Hence for any target `T'`,
`d(T, T') ≥ min( ceil_even ½Σ|s_(i)(T) − s_(i)(T')|, ceil_odd ½Σ|s_(i)(T^op) − s_(i)(T')| )`,
which is THM-4810's `L` when `T' = TT_n`.

## (E) Regular 11-tournaments

Every one of the 1,223 regular classes of order 11 has an all-black word of exactly `m = 5` native moves to a
singleton plus a strong core (THM-4814's search on `gentourng -d5D5 11`); HYP-9422 (1) is FINITE-EXACT through
order 11 and is the regular case of THM-4811 (D2) one order beyond its census.

Proofs, tables and reproduction: the packet README; the external-preprint hooks (MUB(6) and the Eisenstein unit
`1 + 2ω`; local smoothing/Kakeya; ETH) are in the companion note and change no status here.
