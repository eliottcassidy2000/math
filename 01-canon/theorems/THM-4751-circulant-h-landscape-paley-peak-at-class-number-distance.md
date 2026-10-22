---
id: THM-4751
title: "The circulant H-landscape: a Paley peak at class-number distance from the rotational peak. (1) For a prime q = 3 mod 4, q > 3, the Paley tournament P_q and the rotational tournament R_q = Cay(Z_q, {1,...,(q-1)/2}) are B_q = ((q-1)/2 - (2 - (2/q)) h(-q))/2 difference-class flips apart in the circulant class graph (minimum over multipliers), and agree on a fraction 1/2 + (2 - (2/q)) h(-q)/(q-1) of their arcs (Dirichlet; the half-row imbalance of the Ogg note, the circulant companion of THM-4666). (2) B_q <= 1 iff q in {3, 7, 11}; B_q >= 2 for every q >= 19 (elementary). (3) FINITE-EXACT, odd n <= 23: the H-maximal circulant is R_n or a one-flip neighbour of R_n (P_7, Z_9{1,2,3,5}, P_11); the Paley class is a strict local maximum of H at q = 7, 11, 19, 23 and 31 (all its flip-neighbours are isomorphic). (4) The one-flip conjecture HYP-9362(a) would imply that P_q is never the circulant H-maximum for q >= 19, the circulant 'only if' of HYP-9357. (5) DICTIONARY: among primes q = 3 mod 4 the spectral ring Z[(1+sqrt(-q))/2] of P_q is norm-Euclidean (covering radius^2 = (q+1)^2/(16q) < 1) exactly for q = 3, 7, 11, and a PID exactly for the Heegner primes 3, 7, 11, 19, 43, 67, 163."
status: >
  PRE-AUDIT (mac-mini-2026-10-21-fanopsp, 2026-10-21).
  PROVED: statements 1 (from Dirichlet's class number formula; the overlap A_q - B_q = (2 - (2/q))h(-q) is classical
  and already in the repo as the "half-row imbalance" of P_q in 05-knowledge/results/ogg_triangular_chowla_20260927.md;
  new here is the multiplier-minimised flip distance), 2 (elementary, least quadratic non-residue bound), 4 (formal
  consequence), 5 (classical facts on each side).
  FINITE-EXACT: statement 3 (complete multiplier-class tables of circ_counts.c for odd n <= 23, one H computation at
  q = 31; scripts below), and the check of statement 1 for all 1135 primes 7 <= q < 20000.
  PRIOR ART IN THE REPO:
  - THM-4666 (the natural order of P_p has ranking excess p h(-p));
  - the Ogg note (h(-p) as the half-row imbalance);
  - THM-481 (Paley tournaments generate the extended QR codes and the Gleason ring; Kim-Sole 2008). That theorem is the
    proved form of the Paley -> E_8 / Golay bridge used in the results note.
  CAUTION: statement 5 is a DICTIONARY. Every condition that is a threshold in (11, 19) on primes q = 3 mod 4 yields
  {3, 7, 11}, because no such prime lies between 11 and 19.
session: mac-mini-2026-10-21-fanopsp (owner prompt: 1729, 24 dimensions, Heegner numbers, Niemeier/Leech, {7,21} vs {8,24})
source: 05-knowledge/results/lattices_codes_tournaments_20261021.md
scripts:
  - 04-computation/experiments/lattice_bridges_20261021/class_number_distance.gp (+ .out)
  - 04-computation/experiments/lattice_bridges_20261021/circ_flip_distance.py, circ_landscape.py, paley_flip_margin.py (+ .out; data ../tournament_pseudoprimes_20261021/circ_*.out)
  - 04-computation/experiments/lattice_bridges_20261021/P31_oneflip.out (ie_ham.c)
depends_on:
  - Dirichlet's class number formula for Q(sqrt(-q)), q = 3 mod 4, q > 3: (2 - (2/q)) h(-q) = sum_(a < q/2) (a/q) (CITED)
  - the classification of norm-Euclidean imaginary quadratic fields (d = 1, 2, 3, 7, 11) and of class number one (Heegner-Baker-Stark) (CITED)
related:
  - 05-knowledge/hypotheses/HYP-9362-one-flip-conjecture-paley-peak-and-the-7-21-to-8-24-dictionary.md
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
  - 01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md
  - 01-canon/theorems/THM-481-paley-gauge-generates-gleason.md
---

# THM-4751 — The circulant H-landscape: Paley peak and rotational peak

**One line.** Among circulant tournaments on `Z_q`, the number of Hamiltonian paths has a rotational peak and a separate Paley peak. The two peaks lie `B_q` difference-class flips apart, and `B_q` is the class number of `Q(√−q)` up to an affine change. The Paley peak is the global one only when it is adjacent to the rotational peak, which happens exactly for `q = 3, 7, 11`.

## Notation

* `q` is a prime with `q ≡ 3 (mod 4)`, and `m = (q − 1)/2`.
* `QR`, `NQR` are the nonzero squares and non-squares mod `q`.
* `P_q = Cay(Z_q, QR)` is the Paley tournament and `R_q = Cay(Z_q, {1, …, m})` the rotational tournament.
* A circulant tournament `Cay(Z_n, S)` contains exactly one of `k`, `n − k` for each `k ≤ m`. The **circulant class graph** has as vertices the classes of `S` up to multipliers `a ∈ Z_n^×`. Two vertices are adjacent when the sets differ in one difference class (one "flip").
* `A_q = |QR ∩ [1, m]|`, `B_q = m − A_q`.

## 1. Overlap and flip distance (PROVED)

By Dirichlet's class number formula, `A_q − B_q = (2 − (2/q)) h(−q)` for `q > 3`. This half-row imbalance is already in the repo (Ogg note). Hence

    A_q = (m + (2 − (2/q))h(−q))/2,   B_q = (m − (2 − (2/q))h(−q))/2.

**Arcs.** `P_q` and `R_q` share the arc `i → i + k` (`k ≤ m`) iff `k ∈ QR`. So they agree on a fraction `1/2 + (2 − (2/q))h(−q)/(q − 1)` of their arcs.

**Flip distance.** The minimum over multipliers `a` of `#{k ≤ m : −k ∈ aQR}` is `B_q`.
* `a ∈ QR` fixes `QR`, and since `−1 ∈ NQR` the count is then `B_q`.
* `a ∈ NQR` gives `A_q > B_q`, by `h > 0`.

**Checked.** All 1135 primes `7 ≤ q < 20000` (PARI `qfbclassno`; `class_number_distance.out`).

**Heegner primes** (`h = 1`). `B_q = 1, 1, 3, 9, 15, 39` at `q = 7, 11, 19, 43, 67, 163`. These are the Paley tournaments least correlated with `R_q`.

**Companion.** THM-4666 measures `h(−p)` with the linear order on `P_p` (ranking excess `p·h(−p)`). Here the class number measures the distance to the cyclic order's tournament `R_q`.

## 2. Adjacency only for q = 3, 7, 11 (PROVED, elementary)

`B_3 = 0` and `B_7 = B_11 = 1`. `B_q ≥ 2` for every `q ≥ 19`:

* **`q ≡ 3 (mod 8)`.** 2 and 8 are non-residues, and both are `≤ m` once `q ≥ 17`.
* **`q ≡ 7 (mod 8)`.** With `n_0` the least non-residue, `n_0` and `2n_0` are non-residues.
  * The elementary bound is `n_0 < 1/2 + √(q + 1/4)`. With `m_0 = ⌈q/n_0⌉`, the number `m_0n_0 − q ∈ (0, n_0)` is a residue, so `m_0` is a non-residue and `m_0 ≥ n_0`.
  * This gives `2n_0 ≤ m` for `q ≥ 31`.
  * At `q = 23`, `n_0 = 5` and 10 is a non-residue.

## 3. The landscape (FINITE-EXACT)

Source: the complete multiplier-class tables (`circ_counts.c`) for odd `n ≤ 23`.

**Global maximum.**
* It is `R_n` for `n = 13, …, 23`.
* It is a one-flip neighbour of `R_n` for `n = 7` (`P_7`), `n = 9` (`Z_9{1,2,3,5}`, the class of `Z_9{2,3,4,8}`) and `n = 11` (`P_11`).

**Strict local maxima.**
* `n = 7, 9, 11`: the global maximum only.
* `n = 13, 15`: `R_n` only.
* `n = 17`: `R_17` and one class at distance 2.
* `n = 19`: `R_19`, `P_19` (at `B_19 = 3`) and one class at distance 2.
* `n = 21`: 5 classes.
* `n = 23`: 6 classes, among them `P_23` (at `B_23 = 4`).

**The Paley peak.** Residue multipliers are automorphisms of `P_q` and permute its `m` difference classes simply transitively, so all flip-neighbours of `P_q` are isomorphic. The margins `H(neighbour)/H(P_q)` (`paley_flip_margin.out`) are:

| `q` | margin |
|---|---|
| 7 | 0.92593 |
| 11 | 0.97825 |
| 19 | 0.99246 |
| 23 | 0.99542 |
| 31 | 0.99822 |

The `q = 31` value comes from `H = 19841263005616787646549771` against `H(P_31) = 19876729650849384218015865`. All margins are below 1, so `P_q` is a strict local maximum.

## 4. Consequence for HYP-9357 (PROVED, conditional)

If the one-flip conjecture HYP-9362 (a) holds, then by §2, `P_q` is not the circulant `H`-maximum for any prime `q ≥ 19`.

## 5. Euclidean and Heegner primes (DICTIONARY)

The nontrivial eigenvalues of `P_q`, `(−1 ± √−q)/2`, generate `O_K` with `K = Q(√−q)`.

* **Euclidean.** `O_K` is norm-Euclidean exactly for `q = 3, 7, 11`. Its deep hole in `C` has squared distance `(q + 1)²/(16q)`, and this is `< 1` iff `q < 7 + √48`.
* **Heegner.** `O_K` is a PID exactly for `3, 7, 11, 19, 43, 67, 163`.

So `{3, 7, 11}` is cut out by covering radius, while the Heegner primes 19, 43, 67, 163 are PIDs without the Euclidean property, and their Paley tournaments are far from rotational (§1). This joins the list of characterizations of `{3, 7, 11}` in HYP-9357 as a DICTIONARY entry (see the caution in the status).
