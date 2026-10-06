---
id: THM-4552
title: "Exceptional knight tori: the 6x6 knight torus is the two-place product C_4 x Paley(9) (twin squares x, x+(3,3)); the 7x7 knight torus lives in F_49, with knight moves = the norm-5 coset of mu_8 (non-squares), nightrider = complement of Paley(49), queen = Paley(49), and an order-8 rotation; direction classes are transversal iff gcd(n,30)=1"
status: "PROVED (i)-(iii), (vi); FINITE-EXACT (iv)-(v); CITED in (vi): Weil bound for Kloosterman sums, finite Euclidean graphs; independent audit: see the results note, section 10"
session: mac-mini-2026-10-06-sixseven
source: 05-knowledge/results/sixes_and_sevens_20261006.md
scripts:
  - 04-computation/experiments/sixseven_20261006_structure.py (+ .out, ALL CHECKS PASSED)
  - 04-computation/experiments/sixseven_20261006_knight_torusblock_n.c (+ sixseven_20261006_knight_torusblock_n.out)
  - 04-computation/experiments/sixseven_20261006_knight_torus_block.py (CP-SAT lazy-cut cross-check)
related:
  - 01-canon/theorems/THM-4550-chessboard-lonely-runner-dictionary-and-the-knight-weave-law.md
  - 01-canon/theorems/THM-4553-the-point-stripping-ladder-psl27-borel-torus-and-the-collatz-frobenius.md
  - 05-knowledge/hypotheses/HYP-9211-knight-torus-blocking-number-is-seven-by-stars-only.md
  - 05-knowledge/hypotheses/HYP-9168-hp-blocking-number-equals-hall-bound.md
  - 05-knowledge/results/collatz_paley_bridge_20261001.md
  - 05-knowledge/results/collatz_connectivity_from_rigidity_20261001.md
---

# THM-4552 — the exceptional knight tori

Let `K = {(±1, ±2), (±2, ±1)}` and `G_n = Cay(Z_n^2, K)`, the knight graph on the `n x n` torus. For
`n >= 5` the 8 moves are distinct mod `n`, so `G_n` is 8-regular. It is edge-transitive under translations and
the point group `D4`, which acts regularly on `K`.

**(i) Six: a two-place product (PROVED).** Under the CRT map `Z_6^2 -> Z_2^2 x Z_3^2`, the knight set
maps bijectively onto `S_2 x S_3`, where:

* `S_2 = {(1,0), (0,1)}` is the rook step mod 2;
* `S_3 = {±(1,1), ±(1,-1)}` is the bishop step mod 3.

Hence

`G_6 ≅ Cay(Z_2^2, S_2) x Cay(Z_3^2, S_3) = C_4 x Paley(9)` (tensor product; `Paley(9) = K_3 □ K_3`).

Because `C_4 = K_{2,2}` has twins, every square `x` has the same 8 neighbours as `x + (3,3)`.

**(ii) Seven: a Paley torus (PROVED).** Identify `Z_7^2` with `F_49 = F_7[i]`, with norm `N(a,b) = a^2 + b^2`.

1. `K = {z : N(z) = 5} = (1 + 2i) μ_8`, where `μ_8 = ker(N: F_49^* -> F_7^*)`. Every knight move is a non-square of `F_49^*`.
2. The nightrider step set `F_7^* K` equals `{z : N(z) ∈ NQR_7}`. The queen step set (slopes `0, ∞, ±1`) equals `{z ≠ 0 : N(z) ∈ QR_7} = (F_49^*)^2`. So the queen torus graph is `Paley(49)` and the nightrider torus graph is its complement, both `srg(49, 24, 11, 12)`. The 8 slopes of `P^1(F_7)` split as `{0, ∞, 1, 6}` (queen) and `{2, 3, 4, 5}` (knight).
3. Multiplication by `ω = 5 + 5i` (order 8, `ω^2 = i`) is an automorphism of `G_7` fixing 0, a "45-degree rotation".

**(iii) Transversality (PROVED).** The four direction classes of `G_n` (cosets of `<v>` for `v = (1,2), (2,1), (1,-2), (2,-1)`, each an `n`-cycle) are pairwise transversal iff `gcd(n, 30) = 1`. Transversal means `Z_n^2 = <v> ⊕ <w>`, so each cycle of one class meets each cycle of the other exactly once. The least such `n > 1` is 7.

**(iv) Symmetry census (FINITE-EXACT).**

* `|Aut(G_n)|` (nauty) is `8 n^2` for `n = 9, 11, 12, 13, 14`. The exceptions are:
  * `28800` at `n = 5` (`G_5 = K_5 □ K_5`);
  * `2^18 · 144` at `n = 6`;
  * `16 · 49` at `n = 7`;
  * `512 · 64` at `n = 8`;
  * `16 · 100` at `n = 10`.
* For `n <= 26`, the linear maps of `Z_n^2` preserving `K` form `D4` (order 8) except at `n ∈ {5, 6, 7, 8, 10}`. An element of order 8 occurs only at `n = 5` and `n = 7`.

**(v) Blocking (FINITE-EXACT).** For `n = 5, 6, 7, 8` the fewest deletions leaving `G_n` without a Hamiltonian cycle is `β(n) = 7`. The minimum sets are exactly the `8 n^2` stars (all but one move at a square). At `n = 8` this covers all 8,637,487,551 six-sets and all 359,895,314,625 seven-sets through a fixed move.

**(vi) Seven is a Kloosterman graph; the Ramanujan knight tori (PROVED; CITED Weil).**

* `G_7 = Cay(F_49, {z : N(z) = 5})` is the finite Euclidean graph `E_7(5)`: squares at squared distance `5 = 1^2 + 2^2` (Medrano–Myers–Stark–Terras, *Finite analogues of Euclidean space*, CITED). Among primes `p`, the knight set is the whole sphere `{x^2 + y^2 = 5}` only at `p = 7` (`p + 1 = 8`).
* The eigenvalue of `G_7` at the character `(a_1, a_2)` is `-Kl_7(3 N(a))`, where `Kl_7(c) = Σ_{x ∈ F_7^*} e((x + c/x)/7)` and `N(a) = a_1^2 + a_2^2`. So the nontrivial spectrum is the six Kloosterman sums over `F_7`, each with multiplicity 8.
* By Weil's bound `|Kl_q(c)| <= 2 sqrt q` (CITED), every nontrivial eigenvalue is at most `2 sqrt 7 = 2 sqrt(d - 1)`. So `G_7` is a Ramanujan graph: the Ramanujan bound for degree `d = q + 1 = 8` is Weil's bound for `q = 7`.
* **Census:** the knight tori that are Ramanujan are exactly `n ∈ {5, 6, 7, 8, 10}`, the same five `n` as in (iv). For `n >= 12`, the character `(1, 0)` gives `4 cos(2π/n) + 4 cos(4π/n) >= 2 + 2 sqrt 3 > 2 sqrt 7`; `n = 9, 11` are checked directly.

## Proofs

(i) Reduce the 8 moves mod 2 and mod 3. The pairs `((a,b) mod 2, (a,b) mod 3)` are the 8 elements of
`S_2 x S_3`; this is a direct listing. A Cayley graph on `A x B` with connection set `S_A x S_B` is the tensor product of
`Cay(A, S_A)` and `Cay(B, S_B)`.

* `Cay(Z_2^2, S_2)` is the 4-cycle `00-10-11-01`.
* `Cay(Z_3^2, S_3)` is the union of the two classes of slope-`±1` lines, each a `K_3`. They are transversal, so the graph is `K_3 □ K_3 = Paley(9)`.

In `C_4`, `g` and `g + (1,1)` have equal neighbourhoods. The CRT image of `(3,3)` is `((1,1), (0,0))`. ∎

(ii)

1. `N` is surjective onto `F_7^*` with kernel of order `48/6 = 8`, so `{N = 5}` has 8 elements. The 8 knight moves have `N = 5`.
2. `z ∈ (F_49^*)^2` iff `z^24 = 1` iff `N(z)^3 = 1` (because `N(z) = z^8`) iff `N(z) ∈ QR_7 = {1, 2, 4}`. Since `5 ∈ NQR_7`, every knight move is a non-square.
3. A nonzero `z` lies on one of 8 lines through 0, each containing 6 nonzero elements. `N` is constant on a line up to the square class (`N(kz) = k^2 N(z)`).
   * The knight slopes `2, 4, 5, 3` (from `(1,2), (2,1), (1,-2), (2,-1)`) have norm class `NQR`.
   * The queen slopes `0, ∞, 1, 6` have norms `1, 1, 2, 2`, in `QR`.

   The two families are 4 + 4 of the 8 lines, which gives the two step sets.
4. `Paley(49) = Cay(F_49, (F_49^*)^2)` is self-complementary, with parameters `(49, 24, 11, 12)`.
5. `ω^2 = (5 + 5i)^2 = 50i = i`, so `ω` has order 8 and `N(ω) = 50 = 1`. Hence `ω K = K`. ∎

(iii) Two classes are transversal iff `v, w` form a basis of `Z_n^2` iff `det(v, w)` is a unit mod `n`. The six
determinants are `-3, -4, -5, -5, -4, 3`, all units iff `gcd(n, 2·3·5) = 1`. ∎

(vi) Eigenvalues of a Cayley graph on `Z_7^2` are the character sums `λ(a) = Σ_{s ∈ K} e((a_1 s_1 + a_2 s_2)/7)`. Write `a_1 s_1 + a_2 s_2 = Re(ā s) = Tr(ā s / 2)` in `F_49 = F_7[i]` (`Tr(z) = z + z^7 = 2 Re z`). Substituting `w = ā s / 2` turns the sum over `N(s) = 5` into a sum over `N(w) = 5 N(a)/4 = 3 N(a)`. The classical norm-torus identity `Σ_{N(w) = c} e(Tr(w)/q) = -Kl_q(c)` (`q` odd, `w ∈ F_{q^2}`; CITED, e.g. Iwaniec–Kowalski, chapter 11) finishes the proof. The census uses the closed formula `λ(a, b) = 2 Σ cos(2π(a x + b y)/n)` over `(x, y) ∈ {(1,2), (2,1), (1,-2), (2,-1)}`. ∎

(iv), (v): see the scripts. The blocking census:

* For each `k = 6, 7` and each `k`-set of moves through a fixed move `e0`, nested loops keep the sub-pool of a pool of about 1000 Hamiltonian cycles avoiding the prefix.
* The last move is tested against the AND of the sub-pool masks.
* Uncertified sets go to an exact DFS, and every cycle found is re-verified and pooled.

The CP-SAT lazy-cut model independently gives `β(5) = 7` with stars only.

## Remarks

* **Typed bridges.**
  * (i) is the "two places" of the eighteenth S15 Collatz note (`n < 6^L` is fixed by its residues mod `2^L`, `3^L`) at `L = 1`. This is an ANALOGY: the parts are a rook and a bishop, not halving and tripling.
  * (ii) puts the knight on the `NQR_7` side and the queen on the `QR_7` side. That is the same quadratic character as the nineteenth note's real/2-adic codes of the trivial cycle, an EXPLAINED COINCIDENCE with no map between cycles and leaps.
* **Langlands-adjacent reading of (vi) (typed).** `Kl_7(c)` is the trace of Frobenius at `c` on Deligne's rank-2 Kloosterman sheaf over `G_m/F_7`, whose Riemann hypothesis is the Weil bound (CITED). Kloosterman sheaves are geometric-Langlands eigensheaves (Heinloth–Ngô–Yun, CITED). So the 7x7 knight torus is a finite graph whose spectrum is a Frobenius-trace table, and whose expansion is a Riemann hypothesis. This is an exact identification, with no consequence for tours or for Collatz.
* Six fails every field-like structure that seven has: transversality (iii), toroidal queens (`gcd(n, 6) = 1`, Pólya, CITED), and orthogonal Latin squares / `AG(2, 6)` (Tarry, CITED).
