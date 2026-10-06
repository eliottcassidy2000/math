---
id: THM-4553
title: "The point-stripping ladder: PSL(2,7) on P^1(F_7) (2-transitive, no invariant tournament) > Borel = Aut(P_7) (its only invariant tournaments are P_7 and its reverse) > split torus <x -> 2x> = Aut(P_7 - 0), the automorphism group of THM-4524's unique all-odd 6-tournament and the Collatz Frobenius of the S15 nineteenth note; the octonionic J at e_0 is the Fano matching q -> 3q (QR_7 -> NQR_7); beta(P_7) = 6 with 63 minimum sets (56 Hall obstructions of three types and 7 exotic stars)"
status: "PROVED (i)-(iv); FINITE-EXACT + LEAN-CHECKED (v) and the orders in (i)-(ii); independent audit: see the results note, section 10"
session: mac-mini-2026-10-06-sixseven
source: 05-knowledge/results/sixes_and_sevens_20261006.md
scripts:
  - 04-computation/experiments/sixseven_20261006_structure.py (sections E, F, G; + .out)
  - 04-computation/lean/standalone/sixseven_20261006_paley_seven_certificate.lean (native_decide; no sorry)
related:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
  - 01-canon/theorems/THM-4552-exceptional-knight-tori-six-is-a-two-place-product-seven-is-a-paley-torus.md
  - 05-knowledge/results/collatz_paley_bridge_20261001.md
  - 05-knowledge/hypotheses/HYP-9167-paley-minus-a-vertex-is-redei-rigid.md
  - 05-knowledge/hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md
---

# THM-4553 — the point-stripping ladder

`P_7` is the Paley tournament on `F_7` (`x -> y` iff `y - x ∈ QR_7 = {1, 2, 4}`), and `P_7 - 0` is its
subtournament on `F_7^*`.

**(i) Eight to seven.**

* `PSL(2,7)` (order 168) acts 2-transitively on the 8 points of `P^1(F_7)`, so it preserves no tournament.
* The stabiliser of `∞` is the Borel subgroup `B = {x -> ax + b : a ∈ QR_7, b ∈ F_7}`, of order 21.
* `B = Aut(P_7)`.
* `B` has exactly two orbitals on ordered pairs of distinct points of `F_7`, so the `B`-invariant tournaments on `F_7` are exactly `P_7` and its reverse.

**(ii) Seven to six.**

* The stabiliser of `∞` and `0` is the image of the split torus, `H = {x -> ax : a ∈ QR_7} = <x -> 2x>`, of order 3. (We write `H`, not `T`, to keep `T` for the Collatz map.)
* `H = Aut(P_7 - 0)`. By THM-4524, `P_7 - 0` is the unique tournament on 6 vertices with every arc on an odd number of Hamiltonian paths.
* `H` is the group `<x2>` of the nineteenth S15 Collatz note (Proposition 5): the stabiliser of the connection set `QR_7` in `Aut(P_7)`. It acts on the code `QR_7 = {1, 2, 4}` of the trivial cycle as the Collatz map `T(n) = n/2, 3n+1` acts on `{1, 4, 2}`.
* The 6 non-trivial translations `x -> x + b` (with the identity, the unipotent radical of `B`) are the symmetries the note's carry kills.

**(iii) The octonionic `J` (PROVED).** Take the octonions with imaginary units `e_0, ..., e_6` and `e_x e_{x+1} = e_{x+3}` (indices mod 7; this table is alternative). Then:

* the three Fano lines through 0 are `{0, q, 3q}` for `q ∈ QR_7`, namely `{0,1,3}, {0,2,6}, {0,4,5}`;
* left multiplication `J = L_{e_0}` on `T_{e_0} S^6 = span(e_1, ..., e_6)` satisfies `J e_q = e_{3q}` and `J e_{3q} = -e_q` for `q ∈ QR_7`;
* `H` permutes these three complex lines cyclically.

So the almost-complex structure of `S^6` at a basis point is the matching of the six remaining Fano points that sends the out-neighbours `QR_7` of 0 in `P_7` to the in-neighbours `NQR_7`.

**(iv) The shape of `P_7 - 0` (PROVED).** `QR_7` and `NQR_7` each induce a cyclic triangle. Every arc between them goes `QR_7 -> NQR_7`, except the three antipodal arcs `u -> -u` (`u ∈ NQR_7`), which form the antipodal orbit proved odd in THM-4524 (recorded in HYP-9167).

**(v) Blocking (FINITE-EXACT + LEAN-CHECKED).**

* Deleting any 5 arcs of `P_7` leaves a Hamiltonian path.
* Exactly 63 six-arc sets leave none, so `β(P_7) = 6 = N - 1 = hall(P_7)`.
* 56 of the 63 are Hall obstructions in HYP-9168's sense (they kill every 1-path-cycle factor):
  * 21 two-source sets;
  * 21 two-sink sets;
  * 14 sets making three vertices share a single out-neighbour (or a single in-neighbour).
* The other 7 are the stars (strip a vertex). They are exotic minima: the isolated vertex plus the two cyclic triangles of `P_7 - v` is a surviving 1-path-cycle factor, and Hamiltonian paths die because the isolated vertex has no arcs.

## Proofs

(i)

* `PSL(2, q)` is 2-transitive on `P^1(F_q)` (CITED, standard). A tournament invariant under a 2-transitive group would contain both `(x, y)` and `(y, x)`.
* The stabiliser of `∞` is the image of `[[a, b], [0, a^{-1}]]`, acting by `x -> a^2 x + ab`. Its elements are `x -> cx + d`, `c ∈ QR_7`.
* These preserve the square class of `y - x`. So the orbitals are `{y - x ∈ QR}` and `{y - x ∈ NQR}`, and `B ⊆ Aut(P_7)`.
* `|Aut(P_7)| = 21` by exhaustion over `S_7` (Python and Lean).

(ii)

* The stabiliser of 0 in `B` is `H`: `x -> cx`, `c ∈ QR_7`, and it fixes `P_7 - 0`.
* `|Aut(P_7 - 0)| = 3` by exhaustion over `S_6` (Python and Lean), with the three maps `x -> x, 2x, 4x`.
* The Collatz identification is Proposition 5 of `collatz_paley_bridge_20261001.md`.

(iii) The lines are the translates `{x, x+1, x+3}`; those through 0 are `x = 0, 6, 4`. The products are read off the
table, and `e_0 (e_0 e_q) = -e_q` by alternativity. Multiplication by 2 maps `{1,3} -> {2,6} -> {4,5} -> {1,3}`. ∎

(iv) A direct check of the 15 arcs: `1 -> 2 -> 4 -> 1`, `3 -> 5 -> 6 -> 3`, and the `NQR -> QR` arcs are `3 -> 4`, `5 -> 2`, `6 -> 1`. ∎

(v) The Lean file checks the following (the classification into types is in the Python mirror only):

* all `C(21,5) = 20349` five-sets by a subset dynamic program for Hamiltonian paths;
* the star;
* the count 63 among the `C(21,6) = 54264` six-sets.

The Python mirror classifies the 63 by type, and checks that the stars leave a 1-path-cycle factor.

## Remarks

* Stripping one point twice takes `PSL(2,7)` (no tournament) to `Aut(P_7)` (the Paley tournament) to `Aut(P_7 - 0)` (the all-odd tournament, with the Collatz Frobenius as its whole symmetry). Note 19's contrast of maximal symmetry against no symmetry becomes Borel against torus in this language.
* The queen/knight 4+4 split of the 8 slopes of the 7x7 torus (THM-4552) is **not** `PSL(2,7)`-invariant: its stabiliser has order 4 (8 in `PGL(2,7)`).
* **Modular-curve reading (CITED, standard; typed DICTIONARY).**
  * The Borel `B` is the image mod 7 of `Γ_0(7)` (index `8 = |P^1(F_7)|` in `PSL_2(Z)`).
  * The torus `H` is the image of `Γ_0(7) ∩ Γ^0(7)`, which is conjugate to `Γ_0(49)` by `τ -> 7τ`.
  * So the ladder is the tower `X(1) <- X_0(7) <- X_0(49)` (degrees 8, then 7). `X_0(49)` has genus 1: it is the CM elliptic curve 49a1, with complex multiplication by `Q(sqrt(-7))`.
  * `QR_7 = {1, 2, 4}` (the trivial cycle's code) is the set of residues mod 7 of the primes that split in `Q(sqrt(-7))`, and the Paley eigenvalues `(-1 ± sqrt(-7))/2` generate its ring of integers.
  * No Collatz statement follows. The identification only says which Langlands object (a GL(1) Hecke character of `Q(sqrt(-7))`, induced to the weight-2 form of level 49) carries the symmetry the trivial cycle keeps.
* The knight torus `7 x 7` uses `μ_8`, the norm-one subgroup of the non-split torus `F_49^*` (THM-4552 (ii)). The Collatz trivial cycle uses `H`, the image of `SL_2`'s split torus, inside the split torus `F_7^* × F_7^*` of `GL_2`. The cuspidal and principal-series representations of `GL_2(F_7)` are indexed by characters of these two maximal tori (Deligne–Lusztig, CITED). This is a DICTIONARY only, with no predicted statement.
* Compared with the knight torus (THM-4552 (v)): there the star *is* the Hall (2-factor) obstruction and the census finds nothing else. In `P_7` the star is exotic, and the Hall obstructions of sizes 2 and 3 also cost `2k = 6`.
