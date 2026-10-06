# Ramanujan's constant and the trivial Collatz cycle: the code's Gauss period is the Heegner point of discriminant −7 (`j = −15^3`), the only class-number-one CM point a binary code can carry, and under `T` the trivial cycle is the only rational cycle that carries it; Ramanujan's `Q(sqrt(−163))` is invisible to every binary code because 2 is inert there; Δ mod 7 (Ramanujan's own split of the primes into `{1,2,4}` and `{3,5,6}` mod 7), `G_7 = 2^(1/4)` and the 64-series, the Klein quartic over `F_8`, Chudnovsky's 7 and 127, typed

2026-10-06, session opus-2026-10-06-S17 (worktree `codex/session-ramanujan-heegner7-20261006`).

Owner's prompt (verbatim): "consider how ramanujan's constant relates with the split primes {1,2,4} mod 7 and the trivial collatz cycle's code and any other surprisingly relevant ideas you can dig up and merge in"

**Status.**

- **PROVED (elementary, on classical inputs):**
  - Proposition 1: the code's Gauss period is the Heegner point `τ_7 − 1`.
  - Lemma 2: the norm bound for split primes.
  - Proposition 3: visibility, i.e. which imaginary quadratic fields lie in a decomposition field of 2.
  - Proposition 4: the Gauss period of a doubling orbit is a CM point only when `m` is squarefree and `<2>` has index 2, and then it is `(μ(m) ± sqrt D)/2` with `D` fundamental. Statement and proof outline supplied by the audit; proof written out here.
  - Proposition 5: under `T`, the trivial cycle is the only rational cycle whose code is a class-number-one CM point. Proof supplied by the audit.
  - Proposition 6: Frobenius of ordinary curves over `F_2`.
  - Proposition 7: `7 | #E_7(F_(2^k))` iff `3 | k`, and `E_7[7]` is rational iff `21 | k`.
  - The root-selection derivation of `f(sqrt(−7))^24 = 2^12` (section 5.2).
  - The exact degrees 6 and 176 (squarefree denominators, by Proposition 4's argument).
- **FINITE-EXACT:** every check in [the script](../../04-computation/experiments/ramanujan_heegner7_trivial_cycle_20261006.py) (output [`.out`](ramanujan_heegner7_trivial_cycle_20261006.out), ending `ALL CHECKS PASSED`):
  - the Heegner table, including the four non-maximal orders;
  - the Gross–Zagier factorizations and norms;
  - the Deuring pattern for `5 <= p < 72`;
  - all 16 smooth Weierstrass equations over `F_2`;
  - `#E_7(F_(2^k))` for `k <= 8` (direct) and `k <= 36` (recurrence), and the order of `π` in `(O_K/7)^*`;
  - the Klein quartic counts for `p < 110` and over `F_(2^n)`, `n <= 6`;
  - `τ(n)` mod 7 and mod 49 for `n <= 6000`;
  - the Chudnovsky integer identity;
  - the class-number scan;
  - the census of all doubling orbits mod odd `m < 3001` (Propositions 4–5);
  - the exact Gauss-period degrees, by reduction modulo a prime `l ≡ 1 (mod m)`.
- **NUMERICAL (60 to 200 digits):**
  - `j(g)`, the `γ_2` values, `f(sqrt(−7))^24` and `k_7`;
  - the near-integers and Ramanujan's series;
  - the Landau–Ramanujan constant `C_7`;
  - the Weber-unit census.
- **CITED:**
  - Heegner–Baker–Stark;
  - Gross–Zagier 1985;
  - Deuring;
  - Schneider 1937;
  - Weber's `j = (f^24 − 16)^3/f^24` (Cox §12);
  - Ramanujan 1914;
  - Ramanujan's τ manuscript (Berndt–Ono 1999), Wilton 1930, Kolberg 1962, Swinnerton-Dyer 1973;
  - Klein 1879, Elkies 1999, Fité–Lorenzo García–Sutherland 2018, and Serre's bound;
  - the Chudnovskys 1988;
  - Landau 1908 and Selberg–Delange;
  - Rabinowitsch 1913;
  - Construction A.
- **AUDITED:** an independent adversarial audit (2026-10-06) recomputed every claim with its own code and found six errors and several typing and citation issues, all corrected in place (section 9; MISTAKE-568). The audit also contributed Propositions 4 and 5.
- **Collatz: OPEN.** The merge is an **EXPLAINED COINCIDENCE, extended**. The Collatz-relevant facts are faces of two:
  - 2 is a nontrivial square mod 7, equivalently `ord_7(2) = 3`, equivalently 2 splits in `Q(sqrt(−7))`;
  - `h(−7) = 1`.

  Sections 5.3–5.4 and 6.2 are context, not faces of these. Nothing constrains any other integer cycle.

Inherits:

- [collatz_paley_bridge_20261001.md](collatz_paley_bridge_20261001.md), Propositions 1–4:
  - the trivial cycle `1 → 4 → 2` of `T(n) = n/2, 3n + 1` has word `100`;
  - its real code is `{1,2,4}/7 = QR_7` and its 2-adic code is `−{1,2,4}/7`, with numerators `NQR_7`;
  - under the shortcut map `T_1` the period-3 role passes to `{−5, −7, −10}`;
  - only length 3 makes a single cycle code a full Paley set; Gersonides' `2^2 − 3 = 1` makes this one integral.
- The S650 Euler/Heegner horizon (`heegner_prime_horizon_tournament_s650.out`, HYP-2225): the lucky primes `{2, 3, 5, 11, 17, 41}` are the primes `p` with `x^2 + x + p` prime through `x = p − 2`, and `d = 4p − 1` maps them to `{7, 11, 19, 43, 67, 163}`.
- The S16 Langlands ladder ([square11_octic_langlands_collatz_20261006.md](square11_octic_langlands_collatz_20261006.md), section 3).

## 0. Answers in brief

1. **The link, in one line.**
   - The trivial cycle's code orbit `{1, 2, 4}/7` has Gauss period `g = ζ + ζ^2 + ζ^4 = (−1 + sqrt(−7))/2`, with `ζ = e^(2πi/7)`.
   - `g + 1 = (1 + sqrt(−7))/2` is the Heegner point of discriminant −7, so `j(g) = −3375 = −15^3`.
   - Ramanujan's constant comes from the same function at the Heegner point of discriminant −163, a root of Euler's `x^2 − x + 41`: `e^(π sqrt 163) = −j(τ_163) + 744 − 196884 e^(−π sqrt 163) + ... = 640320^3 + 744 − 7.4993×10^(−13)`.
   - So the trivial cycle's code and Ramanujan's constant come from one function, `j`, at two of the nine Heegner points. These are the maximal orders; there are thirteen class-number-one discriminants counting the orders −12, −16, −27, −28.
   - The two points sit at the two ends of Euler's list of lucky primes: `p = 2` (`x^2 + x + 2` is the minimal polynomial of `g`) and `p = 41`.
2. **Why `{1, 2, 4}`.** The split residues mod 7 form the subgroup generated by the split prime 2 itself, because `ord_7(2) = 3 = (7 − 1)/2`. The trivial cycle's code is a doubling orbit. Hence it is the decomposition group of 2 in `Q(ζ_7)`, i.e. the set of residues of the primes that split in `Q(sqrt(−7))`.
3. **Why Ramanujan's field is on the other side.**
   - A prime that splits in a class-number-one field `Q(sqrt(−d))` (`d` squarefree, `d ≡ 3 (mod 4)`) is at least `(d + 1)/4` (Lemma 2). So 2 splits only for `d = 7`, and every prime below 41 is inert in `Q(sqrt(−163))`.
   - Consequently `Q(sqrt(−7))` is the only class-number-one imaginary quadratic field contained in any decomposition field of 2 (Proposition 3).
   - At the level of points: every CM Gauss period of a doubling orbit is `(μ(m) ± sqrt D)/2` with `D` fundamental (Proposition 4). So the only class-number-one CM point any binary code (a rational periodic point of the doubling map) carries is `τ_7`, up to translation and conjugation.
   - `Q(sqrt(−163))` is seen by no binary code.
   - Faces of the same fact, all verified:
     - among all thirteen class-number-one discriminants, `j` is odd exactly for −7 and −28, the two orders of `Q(sqrt(−7))`;
     - equivalently (Deuring), among the nine Heegner CM curves only the discriminant −7 curve `49a1` has ordinary reduction at 2;
     - Gross–Zagier forces `j(τ_7) = −3^3 5^3`, excluding the prime 2 exactly because 2 splits.
4. **The 2-adic face.**
   - Every ordinary elliptic curve over `F_2` has Frobenius `(±1 ± sqrt(−7))/2` (all 16 smooth Weierstrass equations checked).
   - The trace-one curve has Frobenius roots `−g, −ḡ`, the negatives of the Paley `P_7` eigenvalues. Its canonical lift is the CM curve with `j = −15^3` (Deuring).
   - Over `F_(2^k)` it has a rational point of order 7 iff `3 | k`. Frobenius acts on the kernel of `sqrt(−7)` as `×4`, the map by which `T` moves along the trivial cycle. This is exact but tautological: both are the order-3 element of `<2>` acting on `Z/7`.
5. **Ramanujan's other objects meet `Q(sqrt(−7))`.**
   - **The τ function.**
     - `τ(p) ≢ 0 (mod 7)` iff `p ≡ 1, 2, 4 (mod 7)`, and then `τ(p) ≡ 2p`, because `τ(p) ≡ p(1 + p^3)` and `p^3 ≡ (p/7)`. Verified for `p <= 6000`.
     - Ramanujan's unpublished manuscript states `τ(n) ≡ nσ_3(n) (mod 7)` and itself sorts the primes into `p ≡ 1, 2, 4` and `p ≡ 3, 5, 6 (mod 7)` (Berndt–Ono 1999, eqs. (6.2)–(6.4)).
   - **The class invariant and the 64-series.** `G_7 = 2^(1/4)`, i.e. `f(sqrt(−7))^24 = 2^12`. It gives:
     - `e^(π sqrt 7) = 2^12 − 24 − 276 e^(−π sqrt 7) − ... = 4072 − 0.0679`, the d = 7 near-integer in Weber's normalization; the `j`-normalization gives `4119 − 47.07`, no near-integer;
     - Ramanujan's series `16/π = Σ (42k + 5) ((1/2)_k/k!)^3 / 64^k`, his 1914 eq. (29), annotated by him with `q = e^(−π sqrt 7)`, with `64 = G_7^24 = 1/(4 k_7^2 k_7'^2)`; verified to 50 digits;
     - `e^(π sqrt 28) = 255^3 − 744 − 0.01187`, from `j(sqrt(−7)) = 255^3`.
   - **The Landau–Ramanujan constant** of the values of `x^2 + xy + 2y^2` (built from 7, the split primes, and squares of the inert ones): `C_7 = 0.72472`.
   - **Chudnovsky's linear coefficient** `545140134 = 2 · 3^2 · 7 · 11 · 19 · 127 · 163` satisfies `12 × 545140134 = sqrt(163 (1728 − j(τ_163)))`.
     - Gross–Zagier admits only primes dividing some `163 − y^2` that are non-split in `Q(sqrt(−163))` and `Q(i)`.
     - 7 qualifies because `163 ≡ 2 ≡ 3^2` is a square mod 7, i.e. `163 mod 7` lies in `{1, 2, 4}`. 127 qualifies as `127 = 163 − 6^2`.
     - That 7 and 127 are Mersenne numbers is NUMEROLOGY.
6. **The Klein quartic** `x^3 y + y^3 z + z^3 x = 0`.
   - `#X(F_p) = p + 1 − a_p(E_7)(1 + χ(p) + χ̄(p))` for `p < 110`, with `χ` the cubic character mod 7.
   - Over `F_2`, Frobenius cycles the three CM factors because `ord_7(2) = 3`, so `L(T) = 1 + 5T^3 + 8T^6`. This model has 24 points over `F_8`, Serre's bound for genus 3.
   - Optimality over `F_8` is a property of this model. Elkies' `F_2`-model is already optimal over `F_4`.
7. **E8 and 744.**
   - The chain: Hamming code with zeros `α, α^2, α^4` → E8 → `θ_E8 = E_4` → `γ_2 = E_4/η^8` → `γ_2((3 + sqrt(−163))/2) = −640320`.
   - `744 = 3(240 + 8)`.
   - Exact but cosmetic for Collatz: E8 is unique.
8. **Verdict.**
   - **Under `T`.** The trivial cycle is the only rational cycle whose code (real or 2-adic) is a class-number-one CM point (Proposition 5). Below code denominator 3001, the only CM codes of any class number are the trivial cycle's (`D = −7`) and the `1/5` cycle's (`D = −15`).
   - **Under `T_1`.** The point `τ_7` recurs: on the −5 cycle, and on rational cycles with code denominator `m = 7p`, of which there are 21 below 3001, e.g. the cycle of `1/11`, word `110000`. So the uniqueness depends on the coding.
   - **Collatz content.** The Collatz-specific input is still only Gersonides' integrality. This complements, rather than sharpens, the Paley note's Proposition 3. Collatz is untouched.
   - **Langlands ladder.** Everything stays abelian over `Q(sqrt(−7))`: the CM forms are induced from Hecke characters, and Δ mod 7 is reducible. This is consistent with S16's search-based observation that the exact arithmetic in the Collatz notes is GL(1).

## 1. The two codes, recalled

From the Paley note (Proposition 1):

- The parity map `Φ(x) = Σ (T^i(x) mod 2) 2^i` sends the trivial cycle to `−{1, 2, 4}/7` in `Z_2`. For example `Φ(1) = 1 + 8 + 64 + ... = −1/7`, rechecked mod `2^40` in `.out` §3.
- Read as a binary fraction, the word `100` gives the real code `{4, 1, 2}/7`.

So the real numerators are `{1, 2, 4} = QR_7` and the 2-adic numerators are `{3, 5, 6} = NQR_7`. A word's 2-adic numerator orbit is minus the real numerator orbit of the reversed word, so its Gauss period is the complex conjugate.

## 2. The code's Gauss period is a Heegner point (PROVED + NUMERICAL)

**Proposition 1.** Let `ζ = e^(2πi/7)` and `g = ζ + ζ^2 + ζ^4`, the Gauss period of the real code orbit.

1. `g = (−1 + sqrt(−7))/2`. Hence `g^2 + g + 2 = 0`, `g ḡ = 2`, and `g − ḡ = 1 + 2g = sqrt(−7)`.
2. `Im g = sqrt(7)/2 > 0`, and `g + 1 = τ_7 := (1 + sqrt(−7))/2`.
3. `j(g) = j(τ_7) = −3375 = −15^3`.
4. The 2-adic code orbit has Gauss period `ḡ`, and `−ḡ = g + 1 = τ_7`.

*Proof.*

1. `g + ḡ = −1`, and Gauss's quadratic sum gives `g − ḡ = Σ (a/7) ζ^a = i sqrt 7`.
2. Immediate from 1.
3. `j(τ + 1) = j(τ)`, and `h(−7) = 1` makes `j(τ_7)` a rational integer. The value was recomputed to 60 digits.
4. Conjugation. ∎

**The comparison with Ramanujan's constant.** For squarefree `d ≡ 3 (mod 4)`, `τ_d = (1 + sqrt(−d))/2` is a root of `x^2 − x + (d + 1)/4`, i.e. Euler's polynomial `x^2 + x + p` shifted, with `p = (d + 1)/4`. At `q = e^(2πiτ_d) = −e^(−π sqrt d)` the expansion `j = q^(−1) + 744 + 196884 q + ...` gives

`e^(π sqrt d) = −j(τ_d) + 744 − 196884 e^(−π sqrt d) + O(e^(−2π sqrt d))` (`.out` §1):

| `d` | `j(τ_d)` | `e^(π sqrt d) − (−j + 744)` | 2 in `Q(sqrt(−d))` |
|---|---|---|---|
| 7 | `−15^3` | `−47.07` | split |
| 11 | `−32^3` | `−5.86` | inert |
| 19 | `−96^3` | `−0.222` | inert |
| 43 | `−960^3` | `−2.23×10^(−4)` | inert |
| 67 | `−5280^3` | `−1.34×10^(−6)` | inert |
| 163 | `−640320^3` | `−7.4993×10^(−13)` | inert |

So the trivial cycle's code sits at the start of the Heegner list, where the near-integer is worst, and Ramanujan's constant at the end. Both are values of `j` at the root of an Euler polynomial: `x^2 − x + 2` and `x^2 − x + 41`, the extreme lucky primes of S650.

(Ramanujan did not mention `e^(π sqrt 163)` (MathWorld; section 10). The list is the Heegner list.)

## 3. Visibility, and which codes carry CM points (PROVED + FINITE-EXACT)

**Lemma 2 (norm bound).** Let `d > 3` be squarefree, `d ≡ 3 (mod 4)`, `h(−d) = 1`. If a prime `p` splits in `K = Q(sqrt(−d))`, then `p >= (d + 1)/4`.

*Proof.*

- `O_K = Z[τ_d]` has norm form `N(x + yτ_d) = x^2 + xy + ((d + 1)/4) y^2`.
- A split prime is the norm of a principal prime `α ∉ Z`, so `y ≠ 0`.
- Then `p = (x + y/2)^2 + (d/4) y^2 >= d/4`. As `d/4` is not an integer, `p >= (d + 1)/4`. ∎

Equality holds for every Heegner `d >= 7`: the smallest split primes are `2, 3, 5, 11, 17, 41` (`.out` §4). This is the S650 horizon read from the field side. The audit notes the conclusion also holds for the non-maximal `d = 27`, where the smallest split prime is `7 = 28/4`.

**Corollaries.**

- (i) In a class-number-one imaginary quadratic field, 2 splits only for `d = 7`.
  - For `d ≡ 3 (mod 4)`: `(d + 1)/4 <= 2` forces `d <= 7`, and `d = 3` has 2 inert.
  - `d = 1, 2` ramify 2.
  - Directly: for `d ≡ 7 (mod 8)`, `d >= 15`, the form `2x^2 + xy + ((d + 1)/8) y^2` is reduced and not principal, so `h(−d) >= 2`.
  - The scan of squarefree `d ≡ 7 (mod 8)` up to 20000 is a redundant check.
- (ii) Every prime below 41 is inert in `Q(sqrt(−163))`.

**Proposition 3 (visibility).** Let `m` be odd and `H = <2> ⊂ U = (Z/m)^*`, the decomposition group of 2 in `Q(ζ_m)/Q`. An imaginary quadratic field `K` of discriminant `D` lies in `Q(ζ_m)^H` iff `|D|` divides `m` and 2 splits in `K`. Hence:

- among class-number-one fields only `Q(sqrt(−7))` occurs, iff `7 | m`;
- `Q(sqrt(−163))`, like every other Heegner field, lies in no decomposition field of 2.

*Proof.*

- `K ⊂ Q(ζ_m)` iff `|D|` divides `m`.
- `K` is fixed by `H` iff its character has `χ_D(2) = +1`. As `D` is odd, 2 is unramified in `K`, so this says 2 splits.
- Then apply Corollary (i). ∎

**Proposition 4 (CM Gauss periods; statement and outline from the audit).** Let `m > 1` be odd, `C = rH` a coset of `H = <2>` in `U = (Z/m)^*`, and `η_C = Σ_(c ∈ C) ζ_m^c`. Then `η_C` is imaginary quadratic iff:

- `m` is squarefree,
- `[U : H] = 2`, and
- the quadratic character `ψ` with kernel `H` is odd.

In that case `η_C = (μ(m) ± sqrt D)/2`, where `D ≡ 1 (mod 8)` is the fundamental discriminant of `ψ`. So `η_C` generates the maximal order of `Q(sqrt D)`, and `j(η_C)` is a singular modulus of class number `h(D)`.

*Proof.*

1. **Expansion.** Expanding the indicator of `C` in characters of `U/H`, with `k = [U : H]`:
   `η_C = (1/k) Σ_χ χ̄(r) G(χ)`, where `G(χ) = Σ_(c ∈ U) χ(c) ζ_m^c`.
2. **Gauss sums.** If `χ` is induced from a primitive `χ*` of conductor `f`, then `G(χ) = μ(m/f) χ*(m/f) G(χ*)`. This is nonzero iff `m/f` is squarefree and coprime to `f`.
3. **The degree.**
   - The conjugates of `η_C` are the `η_(aC)`, and `a ↦ η_(aC)` has Fourier coefficients `χ̄(r) G(χ)/k`.
   - `σ_a` fixes `η_C` iff `χ(a) = 1` for every `χ` in `S = {χ : G(χ) ≠ 0}`.
   - So the degree of `η_C` is the order of the group generated by `S`.
4. **Squarefree.** Degree 2 forces `S ⊆ {χ_0, ψ}` with `ψ ∈ S` quadratic. Then `ψ ∈ S` forces `m/f_ψ` squarefree and coprime to the squarefree `f_ψ`, so `m` is squarefree.
5. **Index 2.** For squarefree `m`, every `G(χ) ≠ 0`, so `S` is the whole dual of `U/H` and degree 2 means `k = 2`.
6. **The value.** Then `η_C = (μ(m) + ψ(r) G(ψ))/2`, with `G(ψ) = ±sqrt D` (Gauss). This is imaginary iff `ψ(−1) = −1`. Finally `2 ∈ H = ker ψ` gives `D ≡ 1 (mod 8)`.
7. **Converse.** It is the same computation read backwards. ∎

*Census* (`.out` §5, all odd `m < 3001`):

- 269 denominators have an imaginary quadratic Gauss period, every one squarefree, of index 2, with `D` fundamental `≡ 1 (mod 8)`.
- Class number one occurs only for `D = −7`, at `m = 7` and at the 21 values `m = 7p` (`p` with 2 a primitive root and `p ≢ 1 (mod 3)`): 21, 35, 77, 203, …, 2933. This is exactly the list Proposition 4 predicts.
- The prime case `m = p ≡ 7 (mod 8)` with `ord_p(2) = (p − 1)/2` gives `(−1 + sqrt(−p))/2` with `h(−p) = 1, 3, 5, 7, 5, 5, 11, …` for `p = 7, 23, 47, 71, 79, 103, 167, …`.

**Proposition 5 (under `T`, class number one isolates the trivial cycle; proof from the audit).** Among all rational cycles of `T(n) = n/2, 3n + 1`, the trivial cycle `{1, 4, 2}` is the only one whose real or 2-adic code has a class-number-one CM Gauss period.

*Proof.*

1. **Which words.** `T`-cycles are coded by periodic words without `11`, since `3n + 1` is even. Every such word codes a rational cycle (parity conjugacy). A word contains `11` iff some rotation's binary value lies in `(3/4, 1)`, i.e. iff its coset (in lowest terms) has an element `c` with `3m/4 < c < m`.
2. **Which denominators.** By Proposition 4 and Lemma 2, `D = −7` and `m ∈ {7} ∪ {7p}`. The two cosets of `H` are the classes `χ_(−7) = ±1` of `U`, since `H ⊆ ker χ_(−7)` and both have index 2.
3. **`m = 7`.** `{3, 5, 6}` contains 6 (word `110`); `{1, 2, 4}` avoids `(21/4, 7)`. That is word `100`, the trivial cycle.
4. **`m = 7p`.** Both cosets meet `(3m/4, m)`.
   - `p = 3`: 16 (`+1`), 17 (`−1`).
   - `p = 5`: 29 (`+1`), 27 (`−1`).
   - `p >= 11`: the interval has length `7p/4 > 14`, so it contains two full blocks of 7 consecutive integers, i.e. 6 integers of each `χ_(−7)`-sign. At most 2 of them are multiples of `p`, since the length is below `2p`.
5. **2-adic codes.** The 2-adic code of a word is the conjugate of the real code of its reversal, which also avoids `11`. ∎

The census confirms this below 3001: the only CM cosets that are `T`-codes are `{1, 2, 4}/7` (`D = −7`) and `{1, 2, 4, 8}/15` (`D = −15`, `h = 2`, the cycle `1/5 → 8/5 → 4/5 → 2/5`, word `1000`).

Under `T_1`, which realizes every word, Proposition 5 fails:

- `τ_7` is carried by the −5 cycle (word `110`), by `{1/5, 4/5, 2/5}` (word `100`), and by `1/11 → 7/11 → 16/11 → 8/11 → 4/11 → 2/11` (word `110000`, code orbit `{1, 2, 4, 8, 16, 11}/21`, Gauss period exactly `τ_7`);
- the orbits mod 23, 47, 71, … are `T_1`-codes carrying class numbers 3, 5, 7, ….

## 4. The 2-adic face: ordinary curves over `F_2` (PROVED + FINITE-EXACT + CITED)

**Proposition 6.**

1. An elliptic curve over `F_2` is ordinary iff its trace `a` is odd. Then `a = ±1`, and the Frobenius satisfies `x^2 ∓ x + 2 = 0`, of discriminant −7.
2. For `a = +1` the roots are `1 + g = −ḡ` and `1 + ḡ = −g`, the negatives of the nontrivial Paley `P_7` eigenvalues.
3. **Canonical lift.**
   - For an ordinary curve every endomorphism commutes with Frobenius, and `Z[(1 + sqrt(−7))/2]` is maximal, so `End(E_0) = O_K`.
   - Deuring's lifting theorem gives a canonical lift with CM by `O_K`, so `j = j(τ_7) = −3375 ≡ 1 (mod 2)`.

FINITE-EXACT (`.out` §6): there are 16 smooth Weierstrass tuples over `F_2`.

- 8 are ordinary, all with `a_1 = 1`, traces `±1` (4 each), and discriminant −7.
- 8 are supersingular, with traces `−2, 0, 2` in counts 2, 4, 2.

**Proposition 7.** For `E_7 = 49a1` (mod 2: `y^2 + xy = x^3 + x^2 + 1`, trace 1), `7 | #E_7(F_(2^k))` iff `3 | k`. The full `E_7[7]` is `F_(2^k)`-rational iff `21 | k`.

*Proof.*

1. **The kernel of `sqrt(−7)`.**
   - `#E_7(F_(2^k)) = N(1 − π^k)` with `π = 1 + g`.
   - `sqrt(−7) = 1 + 2g` is the prime above 7, so `7 | N(1 − π^k)` iff `π^k ≡ 1 (mod sqrt(−7))`.
   - From `1 + 2g ≡ 0` we get `g ≡ 3 (mod sqrt(−7))`, so `π ≡ 4`, of order 3 in `F_7^*`. Also `π − π̄ = sqrt(−7)`, so `π̄ ≡ π`.
2. **The full `E[7]`.**
   - It is rational iff `π^k ≡ 1 (mod 7)`.
   - `π^3 − 1 = sqrt(−7) · g` has `sqrt(−7)`-valuation exactly 1, so `π` has order 21 in `(O_K/7)^*`. ∎

Checked (`.out` §6):

- the counts `2, 8, 14, 16, 22, 56, 142, 288` for `k = 1..8`;
- the divisibility for `k <= 36`;
- the order 21;
- `v_7 #E_7(F_(2^21)) = 3`.

**Reading.**

- The kernel of `sqrt(−7)` is the kernel of the rational 7-isogeny of `49a1`. Frobenius at 2 acts on it as `×4 = 2^(−1)`.
- Its nonzero orbits `{P, 2P, 4P}` and `{3P, 5P, 6P}` are the two code cosets, permuted as `T` permutes `{1, 4, 2}` (Paley note, Proposition 2.4). This is exact but tautological.
- The trivial cycle's period is the degree `[F_8 : F_2]` over which `E[sqrt(−7)]`, not the whole `E[7]`, becomes pointwise rational.

## 5. Ramanujan's other objects (FINITE-EXACT + NUMERICAL + CITED)

### 5.1 Δ mod 7 sees the split primes

**Statement.** `τ(n) ≡ n σ_3(n) ≡ n σ_9(n) (mod 7)` for all `n`. This is in Ramanujan's unpublished manuscript (Berndt–Ono 1999, eq. (6.2)); Wilton (1930) gave the first published proof of the mod-7 vanishing. For a prime `p ≠ 7`, `τ(p) ≡ p(1 + p^3)`, and `p^3 ≡ (p/7)` by Euler's criterion, with exponent `3 = (7 − 1)/2 = ord_7(2)`. Hence:

- `τ(p) ≡ 0 (mod 7)` iff `p ≡ 3, 5, 6 (mod 7)` (or `p = 7`);
- `τ(p) ≡ 2p (mod 7)` iff `p ≡ 1, 2, 4 (mod 7)`.

**In Ramanujan's own words.** The manuscript's eq. (6.4), which decides when `τ(n) ≢ 0 (mod 7)` from the factorization of `n`, distinguishes exactly the classes `p ≡ 1, 2, 4` and `p ≡ 3, 5, 6 (mod 7)`. So Ramanujan himself wrote down the split classes the owner asked about.

**Checks** (`.out` §12, `n <= 6000`):

- the congruence holds for every `n`;
- Kolberg's refinement mod `7^2` holds for every `n ≡ 3, 5, 6 (mod 7)` (Kolberg 1962);
- the zero set of `τ(p)` mod 7 is the residues `{0, 3, 5, 6}`.

**Galois reading.** The semisimplified mod-7 representation of Δ is `ω ⊕ ω^4 = ω ⊗ (1 ⊕ ω^3)`, and `ω^3(Frob_p) ≡ (p/7) = (−7/p)` is the character of `Q(sqrt(−7))`. 7 is one of Δ's exceptional primes `2, 3, 5, 7, 23, 691` (Swinnerton-Dyer 1973).

**Typing.** Δ is also the denominator of `j = E_4^3/Δ`. That one form carries both this congruence and Ramanujan's constant is ANALOGY, not mechanism.

### 5.2 `G_7 = 2^(1/4)`, Ramanujan's 64-series, and the d = 7 near-integers

**Derivation of `f(sqrt(−7))^24 = 2^12` (PROVED from cited inputs).**

- `h(−28) = 1` gives `j(sqrt(−7)) = 255^3`.
- Weber's identity is `j = (f^24 − 16)^3/f^24`.
- `(X − 16)^3 = 255^3 X` has roots `4096` and `−2024 ± 765 sqrt 7`, both negative.
- `f(sqrt(−7))^24 = e^(π sqrt 7) Π (1 + e^(−(2n−1)π sqrt 7))^24 > 0`, so it equals `2^12`. Equivalently `f(sqrt(−7)) = sqrt 2` (Weber) and `G_7 = 2^(−1/4) f(sqrt(−7)) = 2^(1/4)` in Ramanujan's notation.

**Consequences (NUMERICAL).**

- **The near-integer.** `x f^24 = Π (1 + x^(2n−1))^24 = 1 + 24x + 276x^2 + 2048x^3 + 11202x^4 + ...` with `x = e^(−π sqrt 7)`. Hence
  `e^(π sqrt 7) = 2^12 − 24 − 276x − 2048x^2 − ... = 4072 − 0.06790`.
  The `j`-normalization fails because `196884 x = 48.35`.
- **The 64-series.** `k_7 = θ_2^2/θ_3^2 = (3 − sqrt 7)/(4 sqrt 2)` and `4 k_7^2 k_7'^2 = 1/64 = G_7^(−24)`. Ramanujan's 1914 eq. (29), annotated by him with `q = e^(−π sqrt 7)` and `2kk' = 1/8`, is
  `16/π = Σ_(k>=0) (42k + 5) ((1/2)_k/k!)^3 / 64^k`,
  verified to 50 digits. The d = 7 member of his `1/π` family is powered by `2^6`, as Chudnovsky's d = 163 series is by `640320^3`.
- **The order of conductor 2.** `e^(π sqrt 28) = 255^3 − 744 − 0.011874`. Here `γ_2 = 255 = 2^8 − 1` at `sqrt(−7)` and `γ_2 = −15` at `(3 + sqrt(−7))/2`. The Mersenne shape of 255 comes from Weber's identity, `255 = (2^12 − 2^4)/2^4`.
- **Weber-unit census (NUMERICAL).** For `m ≡ 7 (mod 8)` (2 split), `f(sqrt(−m))/sqrt 2` is an algebraic unit for all eight cases tested (`.out` §11):
  - `m = 7` (value 1);
  - 15 (`w^3 = φ`);
  - 23 (`w^3 = w + 1`, the plastic number);
  - 31, 47, 55, 71;
  - 39 (minimal polynomial `X^12 − 3X^9 − 4X^6 − 2X^3 − 1`, found by the audit and by a search in `w^3`).

  We expect this to be classical (Weber) but did not locate the source, so it is typed NUMERICAL. At `m = 7` the unit is 1 because `h(−28) = 1`.

### 5.3 The Landau–Ramanujan constant of the split primes (CITED method, constant computed)

**Setup.** Let `B_7(x)` count the `n <= x` of the form `a^2 + ab + 2b^2`, i.e. the norms from `Z[g]`: the integers whose prime factors `≡ 3, 5, 6 (mod 7)` occur to even powers.

**The constant.** Selberg–Delange gives `B_7(x) ~ C_7 x/sqrt(log x)` with

`C_7 = π^(−1/2) (L(1, χ_(−7)) · (1 − 1/7)^(−1) · Π_(q ≡ 3,5,6 (7)) (1 − q^(−2))^(−1))^(1/2) = (sqrt(7)/6 · Π (1 − q^(−2))^(−1))^(1/2) = 0.7247195`.

**Accuracy.** The Euler product is truncated at `2·10^6`; the same truncation reproduces `K = 0.7642236` to `6×10^(−9)`. The audit, using accelerated products, gets `C_7 = 0.72471952146864625690` and `K` to 20 digits.

**Counts.** `B_7(x) sqrt(log x)/x = 0.787848, 0.771482, 0.761244, 0.754799` at `x = 10^4, …, 10^7`, decreasing slowly toward `C_7`.

### 5.4 Chudnovsky's 7 and 127 (FINITE-EXACT + CITED)

**The identity.** `1728 − j(τ_163) = 640320^3 + 1728 = 163 m^2`, with `m = 40133016 = 2^3 3^3 · 7 · 11 · 19 · 127`. Moreover `12 × 545140134 = 163 m`, i.e. `B = sqrt(163(1728 − j(τ_163)))/12`. Checked in exact integer arithmetic, with the series gaining 14.2 digits per term.

**Gross–Zagier reading.** For `j(τ_163) − j(i)` (discriminants −163 and −4), the admissible primes:

- divide some `(163 · 4 − x^2)/4 = 163 − y^2`, and
- are non-split in both `Q(sqrt(−163))` and `Q(i)`.

All primes of `m` qualify (`.out` §13).

**Where 7 and 127 come from.**

- `127 = 163 − 6^2`.
- 7 divides `163 − y^2` for `y = 3, 4, 10, 11`, because `163 ≡ 2 ≡ 3^2 (mod 7)` is a square mod 7, i.e. `163 mod 7` lies in `{1, 2, 4}`. For primes `≡ 3 (mod 4)` this is equivalent to 7 being inert in `Q(sqrt(−163))`. Deuring then gives `j(τ_163) ≡ 1728 (mod 7)`, the unique supersingular value; in fact mod 49.

**Typing.** That 7 and 127 are Mersenne primes is NUMEROLOGY.

## 6. The Klein quartic and E8 (FINITE-EXACT + NUMERICAL + CITED)

### 6.1 The Klein quartic `x^3 y + y^3 z + z^3 x = 0`

**Symmetry.**

- Its automorphism group over C is `PSL(2, 7) ≅ GL(3, 2)`, of order 168 (Klein 1879). Over Q only the cyclic coordinate shift survives, consistent with the audit's count of exactly 3 projective automorphisms over `F_5`.
- `GL(3, 2)` is the collineation group of the Fano plane, whose lines are the translates `{1, 2, 4} + b`, and the automorphism group of the Hamming code.
- The two 3-dimensional characters of `PSL(2, 7)` take the values `g, ḡ` on the two classes of elements of order 7 (Elkies 1999, §1.3).

**Point counts** (`.out` §8):

- For every prime `p < 110`, `p ≠ 7`: `#X(F_p) = p + 1 − a_p(E_7)(1 + χ(p) + χ̄(p))`, with `χ` a cubic character mod 7.
  - So `#X(F_p) = p + 1` unless `p ≡ 1 (mod 7)`.
  - The untwisted `p + 1 − 3a_p` fails exactly at the split primes `≡ 2, 4 (mod 7)`.
  - This is what one expects if `Jac(X)` is Q-isogenous to the restriction of scalars of `E_7` from `Q(ζ_7)^+`, i.e. `L(X) = L(E_7) L(E_7 ⊗ χ) L(E_7 ⊗ χ̄)`. That would follow from the identity at all `p` by Faltings; we checked `p < 110` only (NUMERICAL).
  - Geometrically `Jac(X) ~ E^3` (Elkies 1999, citing Ekedahl–Serre). Over Q the standard model's Jacobian is not isogenous to a cube (Fité–Lorenzo García–Sutherland 2018, §3), whereas Elkies' model (1.22) is.
- Over `F_(2^n)`, for this model: 2 is inert in `Q(ζ_7)^+` (`ord_7(2) = 3`), so Frobenius cyclically permutes the three CM factors.
  - Hence `L(T) = 1 + 5T^3 + 8T^6` (using `π^3 + π̄^3 = −5`), i.e. counts `3, 5, 24, 17, 33, 38` for `n = 1..6`, all checked directly.
  - Over `F_8`, `24 = 9 + 3⌊2 sqrt 8⌋` is Serre's bound, so this model is maximal there.
  - Elkies' `F_2`-model (a twist, with Frobenius polynomial `(T^2 − T + 2)^3`) has `0, 14, 24, …` points and is already maximal over `F_4`. So "optimal at degree 3" is a property of this model only.

### 6.2 E8, `γ_2` and 744

**The construction chain (FINITE-EXACT).**

- The cyclic Hamming code with generator `1 + x + x^3` has zeros `α, α^2, α^4`, i.e. zero set `{1, 2, 4}`.
- Its weights are `1, 7, 7, 1`; the extended `[8, 4, 4]` code has 14 words of weight 4.
- Construction A gives E8 with 240 roots, so `θ_E8 = E_4`.
- The integer `q`-series give `q^(1/3) γ_2 = 1 + 248q + 4124q^2 + 34752q^3 + ...` with `248 = 240 + 8`. Cubing gives `q j = 1 + 744q + 196884q^2 + ...`, so `744 = 3 × 248`.

**At the Heegner point (NUMERICAL).** With the branch `η = e^(2πiτ/24) Π (1 − q^n)`:

- `γ_2((3 + sqrt(−163))/2) = −640320` to 55 digits;
- at `(1 + sqrt(−163))/2` it is `640320 e^(−iπ/3)`, because `γ_2(τ + 1) = e^(−2πi/3) γ_2(τ)`.

**Typing.** The dependence on `{1, 2, 4}` is cosmetic, because E8 is unique.

## 7. Hostile controls and numerology

1. **Coding dependence.**
   - Under `T_1`, the trivial cycle `{1, 2}` (word `10`) has code `{1, 2}/3` with Gauss period −1: no CM point.
   - `τ_7` is carried instead by the −5 cycle and by rational cycles with `m = 7p` (section 3).
   - Under `T`, Proposition 5 holds and the trivial cycle is the unique carrier.
2. **No other known integer cycle carries a CM point** (`.out` §5). The Gauss periods are:
   - −1 for the `−1` cycle of `T` (rational);
   - degree 6 for the `−5` cycle of `T` (`m = 31`), PROVED since `m` is squarefree;
   - degree 176 for the `−17` cycle of `T_1` (`m = 2047 = 23·89`), PROVED likewise;
   - degree 7776 for the `−17` cycle of `T` (`m = 262143`), FINITE-EXACT by reduction mod a prime `l ≡ 1 (mod m)`;
   - and the fixed points 0 and `−1` of `T_1` are trivially rational.

   By Schneider (1937), `j` at these periods, conjugated into the upper half-plane, is transcendental.
3. **Higher class numbers.** Under `T`, below 3001 only the `1/5` cycle (`D = −15`, `h = 2`) joins the trivial cycle. Under `T_1`, every index-2 doubling orbit is a code (e.g. mod 23, 47, 71, with `h = 3, 5, 7`).
4. **Numerology, typed:**
   - **Heegner numbers as `|2^A − 3^l|`.** `1, 2, 3, 7, 11, 19` have this form; `43, 67, 163` do not. 11 of the integers `1..20` have it.
   - **`163 ≡ 2 (mod 7)`.** This is not empty: it is what puts 7 into `j(τ_163) − 1728` and into Chudnovsky's `B` (section 5.4). But a residue landing in `{1, 2, 4}` happens half the time a priori, and it carries no Collatz content.
   - **`640320 ≡ 2 (mod 7)`.** No content.
   - **Mersenne shapes** `15, 255, 7, 127`. Each is explained by Weber or Gross–Zagier, not by codes.
   - **`2^18 − 1 = 3^3 · 7 · 19 · 73`.** The generic factorization of a Mersenne number.
   - **The 24s.**
     - `τ(2) = −24` and the 24 in `2^12 − 24` are the same 24: the exponent of `η` in `Δ = η^24`, since `f(τ)^24 = −Δ((τ + 1)/2)/Δ(τ)`.
     - Only the 24 points of `X(F_8)` (`= 9 + 3·5`) is independent of them.
5. **What would have counted as Collatz leverage:** a mechanism tying the existence of an integer cycle to CM data, or a CM invariant separating the trivial cycle from hypothetical long cycles by something other than its period word. Nothing here does that.
   - Proposition 5 is a statement about all rational `T`-cycles at once, decided by the word.
   - Integrality is still Gersonides'.

## 8. Open threads (typed)

1. **All CM codes under `T` (FINITE + PROVABLE).** Are the trivial cycle (`D = −7`) and the `1/5` cycle (`D = −15`) the only rational `T`-cycles whose codes are CM points of any class number?
   - Verified below 3001.
   - A proof for all `m` needs both `ψ`-classes of units to meet `(3m/4, m)`. Pólya–Vinogradov plus a finite check should settle it.
2. **`Q(sqrt(−23))` (OPEN, low prior).**
   - It is the smallest imaginary quadratic field in which both 2 and 3 split (`d ≡ 23 (mod 24)`; `h = 3`).
   - Its Weber invariant is `f(sqrt(−23))/sqrt 2 = ρ`, the plastic number.
   - Δ mod 23 sees its Hilbert class field (Wilton).
   - Whether it deserves a "Collatz-adapted CM field" reading is open.
3. **The canonical-lift coordinate (OPEN, low prior).**
   - In the embedding with `v_2(g) = 1`, the Frobenius unit root is `1 + g ∈ Z_2^*`, and `Φ(1) = (1 + 2g)^(−2)`.
   - Whether the Serre–Tate coordinate has any meaning for the Bernstein–Lagarias conjugacy is open; no mechanism is known.
4. **A certificate for the 64-series (FINITE, cheap, low value).** A Lean bracket check of its first terms against `16/π`, a companion to S16's certificate.

## 9. Audit

An independent adversarial audit (2026-10-06, blind subagent, its own code) confirmed every computation.

- **Confirmed with independent methods:**
  - `j` via theta functions, `γ_2` via the pentagonal `η`, `k_7` via `K'/K`;
  - `τ` via Jacobi's `η^3`;
  - Klein counts to `n = 9`;
  - `C_7` and `K` to 20 digits.
- **It found six errors, corrected above:**
  - (E1) §7.3 mixed the codings. The orbits mod 23, 47, 71 contain `11`, so they are codes of `T_1`-cycles, not `T`-cycles, and "class number one isolates the trivial cycle" was unproved and false under `T_1`. It is now Proposition 5 (true under `T`, proof from the audit) with the `T_1` counterexamples.
  - (E2) "nine class-number-one CM points" omitted the four non-maximal orders; the odd-`j` statement is now over all 13.
  - (E3) The 7-torsion over `F_8` is `E[sqrt(−7)]`; the full `E[7]` needs `F_(2^21)`.
  - (E4) Two of the "unrelated" 24s are the same 24, the `η`-exponent.
  - (E5) `g ≡ 3` holds mod `sqrt(−7)`, not mod 7.
  - (E6) The rounding `0.7879` should be `0.7878`.
- **Typing and citations:**
  - the point-level claim now rests on Proposition 4 (statement and outline from the audit);
  - the ×4 sentence is labeled tautological;
  - "faces of three" became two (`ord_7(2) = 3` and `χ_(−7)(2) = 1` are equivalent);
  - §6.2's "EXACT" was retyped;
  - the Gauss-period degrees were upgraded;
  - the restriction-of-scalars sentence was retyped and the Fité–Lorenzo García–Sutherland reference added;
  - the Chudnovsky mechanism (Gross–Zagier norms) replaced "small inert primes";
  - Lemma 2 gained "squarefree";
  - "optimal over `F_8`" is now stated for this model only;
  - "Ramanujan's list" became the Heegner list;
  - Plouffe was credited for the name; Kolberg and Wilton were credited;
  - "sharpens" became "complements".
- **Logged as MISTAKE-568.**

## 10. Reproduction and references

`python3 04-computation/experiments/ramanujan_heegner7_trivial_cycle_20261006.py` (a few minutes; mpmath, sympy, numpy). Output: [ramanujan_heegner7_trivial_cycle_20261006.out](ramanujan_heegner7_trivial_cycle_20261006.out), 15 sections, ending `ALL CHECKS PASSED`.

References:

- **Ramanujan's constant and class number one:**
  - Hermite (1859) noticed the near-integer `e^(π sqrt 163)`;
  - Simon Plouffe coined the name "Ramanujan's constant", after Martin Gardner's April 1975 hoax column (MathWorld, "Ramanujan Constant"); Ramanujan himself did not mention `e^(π sqrt 163)`;
  - Heegner 1952; Baker 1966; Stark 1967;
  - G. Rabinowitsch (1913).
- **Ramanujan:**
  - *Modular equations and approximations to π*, Quart. J. Math. 45 (1914): class invariants, the `1/π` series, eq. (29);
  - B. Berndt and K. Ono, *Ramanujan's unpublished manuscript on the partition and tau functions* (1999), eqs. (6.2)–(6.4).
- **τ congruences:**
  - J. R. Wilton (1930);
  - O. Kolberg, Math. Scand. 10 (1962);
  - H. P. F. Swinnerton-Dyer, *On ℓ-adic representations and congruences for coefficients of modular forms*, LNM 350 (1973).
- **Singular moduli and CM:**
  - B. Gross and D. Zagier, *On singular moduli*, J. reine angew. Math. 355 (1985);
  - M. Deuring (1941);
  - T. Schneider, Math. Ann. 113 (1937);
  - D. Cox, *Primes of the form x^2 + ny^2*, §12.
- **The Klein quartic:**
  - F. Klein, Math. Ann. 14 (1879);
  - N. Elkies, *The Klein quartic in number theory*, in *The Eightfold Way*, MSRI Publ. 35 (1999);
  - F. Fité, E. Lorenzo García and A. Sutherland, arXiv:1712.07105 (2018), §3;
  - J.-P. Serre (1983), the bound `q + 1 + g⌊2 sqrt q⌋`.
- **1/π series:** D. and G. Chudnovsky, *Approximations and complex multiplication according to Ramanujan* (1988).
- **Lattices:** Conway–Sloane, *Sphere Packings, Lattices and Groups* (Construction A).
- **Landau–Ramanujan:** E. Landau (1908).
