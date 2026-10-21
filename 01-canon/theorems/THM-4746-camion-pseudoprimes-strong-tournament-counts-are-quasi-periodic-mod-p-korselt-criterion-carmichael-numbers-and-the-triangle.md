---
id: THM-4746
title: "Camion pseudoprimes. Let s(m) be the number of strong (equivalently, by Camion's theorem, Hamiltonian) labelled tournaments on m vertices (A054946) and c(m) the number of circulant tournaments on Z_m (2^((m-1)/2) for odd m, 0 for even m). (1) For every odd prime p and every m >= 1, s(m+p-1) = (2/p) s(m) (mod p), and the same holds for D = s - c. (2) Hence s(kp) = (2/p)^k s(k) (mod p) and s(p) = (2/p) (mod p): the strong count obeys Euler's criterion. (3) Modulo 3 the strong exponential generating function is sin x, so s(m) = 2^((m-1)/2) (mod 3) for odd m and s(m) = 0 (mod 3) for even m >= 2; modulo 5 the odd case also holds for every odd m. (4) Korselt-type criterion: if n is composite with squarefree odd part, then s(n) = c(n) (mod n) iff p divides D(r_p) for every odd prime p | n, where r_p = n (mod p-1) taken in 1..p-1. (5) D(1) = D(2) = D(3) = 0 (K1, no strong K2, the triangle) and D(r) > 0 for r >= 4, so every Carmichael number, every 2p and every 3p is a Camion pseudoprime, and so is every 3ab with b >= 5 and a = 3b-2 both prime (all residues in the triangle class 3), in particular 4^q - 1 whenever the Mersenne number 2^q - 1 and the Wagstaff number (2^q+1)/3 are both prime; 21 = 3*7 is one and the base-2 pseudoprime 341 = 11*31 is not (its residue at 31 is the Wagstaff factor 11, and 31 does not divide D(11)), while 3*341 = 1023 = 2^10 - 1 is. (6) The tournament-count test 2^C(n,2) = c(n) (mod n) is the Fermat test to the base c(n) = 2^((n-1)/2); its only universal residue is r = 1, so the triangle (s(3) = c(3) = 2) is exactly what separates the two tests. (7) The same quasi-periodicity, with (3/p), holds for strong semicomplete digraphs; connected graphs (unordered components) have none."
status: >
  PROVED: statements 1-7 (elementary; statement 1 by a group-ring argument in the Hurwitz ring of binomial
  convolution over F_p or F_(p^2); statement 2 also by the rotation action; statement 3 by an integer identity;
  statement 4 by the Chinese remainder theorem, statement 1 and the free action of 2-subgroups of rotations).
  FINITE-EXACT: periodicity for every odd prime p <= 199 (m <= 3p+5); the zero-class table Z_p for p <= 199; the
  criterion against direct computation for all 472 odd squarefree composites n <= 2500 (prime factors <= 199) and
  all 142 even n <= 400 with squarefree odd part; the direct list of all odd composite Camion pseudoprimes n <= 1200;
  the squarefree list to 20000 (2145 members; all nine Carmichael numbers <= 20000 pass); the Camion-Wieferich primes
  (D(p) = 0 mod p^2) for p <= 3000, which are exactly 3 and 13; every 2^k - 1 with k <= 64 except k = 59 (undecided),
  the passes being exactly k = 4, 10, 14, 26, 34, 38, 62; the controls of statement 7 for 5 <= p <= 23. NOT CHECKED: novelty against the literature (OEIS A054946 records no congruence;
  the sequences of statements 5-6 are not in the OEIS as of 2026-10-21). OPEN: HYP-9358.
session: mac-mini-2026-10-21-fanopsp (owner prompt: tournament structure, the 7/21 gaps, converse symmetry, 341 = 11*31 and {2,3,11}; "think pseudoprimes")
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_psp.c (+ camion_psp_1200.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_split.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_period.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_verify.c (+ camion_verify.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_even_sq.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_D_factor.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_list.c (+ camion_list_20000.out, camion_squarefree_list.txt)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_controls.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/mersenne_camion.py (+ .out)
depends_on:
  - Camion 1959 (a tournament is strong iff it has a Hamiltonian cycle) — CITED, used only for the name and the reading s(m) = Hamiltonian tournaments
  - Wright's recurrence 2^C(m,2) = sum_k C(m,k) s(k) 2^C(m-k,2) (OEIS A054946) — standard; reproved in one line below
related:
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (the same ordered strong decomposition makes H multiplicative)
  - 01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md (path-rigid strong atoms are K1 and C3)
  - 01-canon/theorems/THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md
  - 01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md
  - 01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md
  - 05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-wieferich-type-primes-sporadic-zero-classes-prime-squares.md
  - 05-knowledge/results/trijection_four_nodes_three_edges_20261020.md (kind-pasteur: |Tour(m)/Z_m|, the Burnside count read in statement 6)
---

# THM-4746 — Camion pseudoprimes

**One line.** At a prime order p, the rotations make the number of strong labelled tournaments congruent to the number of circulants, which is `(2/p)` by Euler's criterion. The composite orders that still pass are governed by a Korselt-type criterion. Its three free residues `1, 2, 3` are K1, the absent strong K2, and the triangle. So every Carmichael number, every `2p` and every `3p` passes. That includes `21`, but not the classical base-2 pseudoprime `341`.

## Notation

* `s(m)`: strong labelled tournaments on `[m]`: `1, 0, 2, 24, 544, 22320, 1677488, …` (A054946). By Camion's theorem these are exactly the labelled tournaments with a Hamiltonian cycle.
* `T(m) = 2^C(m,2)`: all labelled tournaments.
* `c(m)`: circulant tournaments on `Z_m` with its fixed labelling. `c(1) = 1`, `c(m) = 2^((m−1)/2)` for odd `m`, and `c(m) = 0` for even `m ≥ 2` (the pair `{x, x+m/2}` cannot be oriented rotation-invariantly).
* `D(m) = s(m) − c(m)`. `D(1) = D(2) = D(3) = 0`, `D(4) = 24`, `D(5) = 540`, `D(6) = 22320`, `D(7) = 1677480`. For `r ≥ 4`, `D(r) > 0`: a strong labelled tournament that is not regular is not a circulant, and such tournaments exist for every `r ≥ 4`.
* `ε = (2/p)` for an odd prime `p`.
* **Camion pseudoprime:** a composite `n` with `s(n) ≡ c(n) (mod n)`.

**Wright's recurrence.** Every tournament is the ordered sum of its strong components. Splitting off the top component gives

    T(m) = Σ_(k=1..m) C(m,k) s(k) T(m−k)   (m ≥ 1).                          (W)

## 1. Quasi-periodicity modulo p (PROVED)

**Theorem 1.** For every odd prime `p` and `m ≥ 1`: `s(m+p−1) ≡ ε s(m)` and `D(m+p−1) ≡ ε D(m)` (mod `p`).

*Proof.*

* **The ring.** Let `K = F_p` if `ε = 1`, and `K = F_(p²)` otherwise. Fix `g₀ ∈ K` with `g₀^(p−1) = ε`.
  * Let `H(K)` be the ring of sequences `a = (a_m)_(m≥0)` over `K` with binomial convolution `(a∗b)_m = Σ_k C(m,k) a_k b_(m−k)`. Its unit is `e = (1, 0, 0, …)`.
  * For `x ∈ K` put `E_x = (x^m)_m`, with `0⁰ = 1`. The binomial theorem gives `E_x ∗ E_y = E_(x+y)`.
  * So for the additive subgroup `G = g₀F_p` (order `p`), `ι : K[G] → H(K)`, `[x] ↦ E_x`, is a ring homomorphism.
* **`T` is a combination of exponentials.** `T(m+p−1)/T(m) = 2^(m(p−1)) · 2^((p−1)(p−2)/2) ≡ ε^(p−2) = ε` for every `m ≥ 0`.
  * The sequences with `a_(m+p−1) = ε a_m` form a `(p−1)`-dimensional space: they are determined by `a_0, …, a_(p−2)`.
  * The `p−1` sequences `E_(g₀a)`, `a ∈ F_p^*`, lie in it. They are independent: their values at `m = 0..p−2` form a Vandermonde matrix with distinct nodes.
  * So `τ := (T(m) mod p)_m = ι(t)` with `t = Σ_(a≠0) c_a [g₀a]`.
* **Inversion.** By (W), `(e − σ) ∗ τ = e`, where `σ = (s(m) mod p)_m` and `σ_0 = 0`.
  * The augmentation of `t` is `τ_0 = 1`. The augmentation ideal of `K[G]` is nilpotent (`([x] − 1)^p = [px] − 1 = 0`), so `t` is a unit of `K[G]`.
  * `ι(t^(−1)) ∗ τ = e`. Since `τ_0 = 1`, `τ` has a unique inverse in `H(K)`, so `e − σ = ι(t^(−1)) = Σ_(x∈G) d_x E_x`.
* **Conclusion.** For `m ≥ 1` the term `x = 0` vanishes, so `s(m) ≡ −Σ_(a≠0) d_(g₀a) (g₀a)^m`. Multiplying each term by `(g₀a)^(p−1) = ε` gives `s(m+p−1) ≡ ε s(m)`.
* **The circulant count.** `c(m+p−1) ≡ ε c(m)`: for odd `m`, `2^((p−1)/2) ≡ ε`; for even `m`, both sides are `0`. ∎

**The mechanism is the linear order of the strong components** (control, `camion_controls.out`).

* Strong semicomplete digraphs have the same ordered decomposition, with `3^C(m,2)`. They obey `s₃(m+p−1) ≡ (3/p) s₃(m)` for `5 ≤ p ≤ 23`, with the same proof.
* Connected labelled graphs have unordered components: `exp`/`log` instead of `1/(1−S)`. They have no period or antiperiod `L < 2p` modulo `p` for any `p ≤ 23`.
* This is the same ordered decomposition that makes `H` multiplicative over strong components in THM-1370 (1), and thereby produces the gaps `7` and `21`.

## 2. Euler's criterion and the rotation lemma (PROVED)

**Corollary 2.** `s(kp) ≡ ε^k s(k) (mod p)` for `k ≥ 1`. In particular `s(p) ≡ (2/p) (mod p)`.

* *First proof.* Apply Theorem 1 `k` times to `kp = k + k(p−1)`.
* *Second proof (rotations).*
  * The rotations of order `p` act on labelled strong tournaments on `Z_(kp)`. Each orbit has size 1 or `p`.
  * A tournament fixed by them has every strong component a union of cosets, because the condensation is a linear order with no nontrivial automorphism.
  * Running (W) over cosets, with `2^(k(kp−1)/2) ≡ ε^k 2^C(k,2)` replacing `T(k)`, shows that the fixed count is `≡ ε^k s(k)`.
* `s(p) ≡ 2^((p−1)/2)` is the statement "strong tournaments ≡ circulants" at prime order. Theorem 1 reads it as Euler's criterion.

## 3. The moduli 3 and 5 test nothing (PROVED)

**Mod 3.** `T(m) ≡ (−1)^C(m,2)` is the coefficient sequence of `cos x + sin x`.

* The claim is `σ ≡ (coefficients of sin x)`, i.e. `s(m) ≡ (−1)^((m−1)/2)` for odd `m` and `0` for even `m ≥ 2`.
* The binomial convolution of `sin` with `cos + sin` is `sin x cos x + sin² x`. Its coefficients for `m ≥ 1` are the integers `2^(m−1)(s_m − c_m)`.
* (W) therefore requires `c_m + s_m ≡ 2^(m−1)(s_m − c_m) (mod 3)`, where `c_m, s_m` are the coefficients of `cos` and `sin`.
  * For even `m`: `s_m = 0` and `1 + 2^(m−1) ≡ 0`.
  * For odd `m`: `c_m = 0` and `2^(m−1) ≡ 1`.
* Hence `s(m) ≡ 2^((m−1)/2) (mod 3)` for every odd `m`.

**Mod 5.** The odd residues modulo `p − 1 = 4` are `1` and `3`, and `D(1) = D(3) = 0`. So `s(m) ≡ 2^((m−1)/2) (mod 5)` for every odd `m`. For even `m`: `s(m) ≡ 0` if `m ≡ 2 (mod 4)`, and `s(m) ≡ ±24` if `4 | m`.

**From 7 on.** `p = 7` is the first prime with an odd residue other than `1, 3`, namely `r = 5`. There `D(5) = 540 ≡ 1 (mod 7)`: the prime 7 is the first to "see" the 544 strong labelled 5-tournaments.

## 4. The Korselt-type criterion (PROVED)

**Theorem 4.** Let `n ≥ 4` be composite with squarefree odd part. Then `s(n) ≡ c(n) (mod n)` iff `p | D(r_p)` for every odd prime `p | n`, where `r_p ∈ {1, …, p−1}` and `r_p ≡ n (mod p−1)`.

*Proof.*

* **The 2-part.** If `2^a ∥ n` with `a ≥ 1`, then `c(n) = 0`. The rotation subgroup of order `2^a` acts freely on labelled tournaments on `Z_n`: every nontrivial element has a power equal to the rotation by `n/2`, which reverses each pair `{x, x+n/2}`. So `2^a | s(n)`.
* **Odd `p`.** By Theorem 1, `D(n) ≡ ε^((n − r_p)/(p−1)) D(r_p) (mod p)`.
* **CRT** combines the parts. ∎

**Zero classes.** `Z_p = {r : p | D(r)}` always contains `{1, 2, 3}`.

* Sporadic extra classes occur at a rate consistent with chance. Examples (odd `r`): `p = 43: 39`; `47: 17`; `59: 21`; `67: 65`; `71: 51`; `79: 73`; `101: 19`; `103: 37, 43`; `113: 109`; …. The table is in `camion_verify.out`.
* Even sporadic classes come from the factors of `D(r)`. For example `D(6) = 2⁴·3²·5·31`, so `r = 6` is a zero class mod 31, and `186 = 6·31` passes.
* `D(r)` for `r ≤ 16` is factored in `camion_D_factor.out`.
* The criterion agrees with direct computation in every tested case: 472 odd squarefree composites `≤ 2500` and 142 even `n ≤ 400`.

## 5. Consequences (PROVED)

* If `n ≡ 1, 2` or `3 (mod p−1)` for every odd prime `p | n` (odd part squarefree), then `n` is a Camion pseudoprime. In particular:
  * **every Carmichael number** (squarefree, `p − 1 | n − 1`): `561, 1105, 1729, 2465, 2821, 6601, 8911, …` all pass;
  * **every `2p` and every `3p`** (`p` an odd prime): `r_p = 2`, resp. `r_p = 3`, and the factor 3 is free (§3).
* **A Chernick-type family (PROVED).** If `b ≥ 5` and `a = 3b − 2` are both prime, then `n = 3ab` is a Camion pseudoprime.
  * `a ≡ 1 (mod b−1)` gives `n ≡ 3 (mod b−1)`.
  * `3b = a + 2` gives `n ≡ 3 (mod a−1)`.
  * The factor 3 contributes the residue `1`.
  * Both large primes therefore sit in the triangle class `r = 3`.
  * Members: `195 = 3·5·13`, `399 = 3·7·19`, `1023 = 3·11·31`, `1443 = 3·13·37`, `4623 = 3·23·67`, `12099 = 3·37·109`, `16383 = 3·43·127`, `19599 = 3·47·139`, ….
* **Mersenne–Wagstaff case.** With `a = 2^q − 1` and `b = (2^q + 1)/3` (so `a = 3b − 2`), `4^q − 1 = 3ab` is a Camion pseudoprime whenever both are prime: `q = 5, 7, 13, 17, 19, 31, 61, 127, …`.
  * These are the exponents with both a Mersenne and a Wagstaff prime, which is the New Mersenne Conjecture's double case.
  * `mersenne_camion.out` decides every `2^k − 1`, `k ≤ 64`, except `k = 59` (undecided: a residue is too large to evaluate). The passes are exactly `k = 4, 10, 14, 26, 34, 38, 62`. The composite `2^22 − 1` fails at 89 and `2^46 − 1` fails at 178481.
* **The two Cipolla numbers split.**
  * *Why.* At its Mersenne factor `a`, the Cipolla number `ab = (4^q−1)/3` has residue `r = b`, the Wagstaff factor itself.
    * It passes there only if `a | D(b)`. For `q = 3`, `b = 3` is the triangle class.
    * The extra factor 3 in `4^q − 1` moves the residue to `3`.
  * `21 = 3·7 = (4³−1)/3` is a Camion pseudoprime.
  * `341 = 11·31 = (4⁵−1)/3` is not: `r_31 = 11`, and `31 ∤ D(11) = 2⁵·3²·5·11·2226142198097`.
  * `341` is the smallest base-2 Fermat pseudoprime, and `21` fails that test. So the strong test and the Fermat test disagree in both directions on the first two Cipolla numbers.
* **Prime squares** need the residue modulo `p²`. `D` is not periodic modulo `p²`. Direct computation: `9` and `169` pass; `25, 27, 49, 81, 121, 125` fail.
  * FINITE-EXACT (odd `m ≤ 700`): `D(m) ≢ 0 (mod 9)` exactly for `m ≡ 1 (mod 6)`, `m ≥ 7`. These classes contain no multiple of 3, so every odd multiple of 9 up to 700 passes modulo 9. No proof is claimed.
  * The "Camion–Wieferich" primes, with `D(p) ≡ 0 (mod p²)` (`p²` passes), are `3` (since `D(3) = 0`) and `13` (`13² | D(13)`), and no others up to 3000 (`camion_list_20000.out`).
* **Data.**
  * Odd composite Camion pseudoprimes `≤ 1200` (direct): `9, 15, 21, 33, 39, 45, 51, 57, 63, 69, 87, 93, 105, 111, 123, 129, 141, 159, 169, 177, …, 1105, …`.
  * The ones not of the form `3p` up to 20000 are listed in the results note: `105, 195, 231, 399, 561, 595, 1023, 1105, 1239, 1443, 1729, …`.

## 6. The tournament-count test and the triangle (PROVED)

**Statement.** For odd composite `n`, `T(n) ≡ c(n) (mod n)` iff `ord_n(2) | (n−1)²/2`, i.e. iff `n` is a Fermat pseudoprime to the base `c(n) = 2^((n−1)/2)`.

* It contains every odd base-2 pseudoprime (`341, 561, 645, 1105, …`). It also contains composites that are not base-2 pseudoprimes: `85 = (4⁴−1)/3`, `133`, `481`, `511 = 2⁹ − 1`, `949`, ….
* Its Korselt form has only the universal residue `r = 1`: `T(r) ≡ c(r) (mod p)` iff `2^((r−1)²/2) ≡ 1`, and `r = 3` gives `4 ≡ 1`, i.e. `p = 3`.
* The strong test has the extra universal residue `3` because `s(3) = c(3) = 2`: both strong labelled triangles are circulants. `T(3) = 8 ≠ c(3)`. **So the triangle is exactly what separates the two tests.**
  * `21` passes the strong test only.
  * `85` and `341` pass the count test only.
  * Carmichael numbers pass both.

**Burnside reading.** For prime `m`, the number of rotation classes of all, resp. strong, tournaments on `Z_m` is `(T(m) + (m−1)c(m))/m`, resp. `(s(m) + (m−1)c(m))/m`.

* The count-test pseudoprimes, resp. Camion pseudoprimes, are exactly the composite odd `m` at which these prime-order formulas still return integers.
* The first formula is kind-pasteur's `|Tour(m)/Z_m|` (trijection note, 2026-10-20) evaluated as if `m` were prime.

## 7. What this does and does not say

* **Says.** The rotation symmetry of tournaments is an honest Fermat test with a Korselt criterion and Carmichael-number behaviour. Its free residues are the orders at which every strong tournament is a circulant (`1, 3`, and vacuously `2`).
* **The same atoms recur.** These are the strong atoms singled out elsewhere in the canon:
  * the path-rigid atoms `{K1, C3}` of THM-2477 / THM-4741;
  * the first members of the pure-blue alphabet `{K1, C3, T5}` of THM-4730;
  * the substituted towers `C_3^[k]` of the converse-Gluck classification (THM-4743 / Session24).
* **Does not say.** Nothing here explains why `7` and `21` are not Hamiltonian-path counts (THM-1370 does), and nothing touches LRC.
  * The shared mechanism is the ordered strong decomposition; the shared atoms are `K1` and `C3`.
  * The remaining coincidences (`7` is the first prime that sees strong 5-tournaments; `21` is a Camion pseudoprime) are recorded as ANALOGY in the note.
