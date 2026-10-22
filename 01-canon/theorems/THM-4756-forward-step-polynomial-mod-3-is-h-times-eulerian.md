---
id: THM-4756
title: "Every tournament is H copies of the transitive tournament modulo 3: the forward-step polynomial F(T,x) = sum over permutations of x^(number of T-forward steps) satisfies F(T,x) = H(T) A_n(x) (mod 3), with A_n the Eulerian polynomial; A_n = (x-1)^(n-1) (n odd), (x-1)^(n-2)(x+1) (n even) mod 3. Proof: F(T,x) = sum_m c_m(T)(x-1)^(n-m) with c_m the ordered path-cover numbers, c_m = 0 mod 3 for m >= 3 (factor m!), and c_2 = ((-1)^n + 1) H mod 3 from the exact doubled Redei-Berge identity sum_(X subset V) H(T[X])H(T[V\\X]) = sum_R 4^|R| 2^(n-|V(R)|) (R over sets of disjoint odd cycles). Also F(T,x) = (1+x)^(n-1) mod 2. This completes the proof of THM-086 (and its Eulerian corollary); at the orders n = 2*3^j + 1 the mod-3 pattern collapses to three terms."
status: >
  PROVED (pre-audit, mac-mini-2026-10-22-sierp). Exhaustively verified for every tournament with n <= 9 (198,946 classes,
  fwd_poly_mod3.out), and the exact identity of step 2 for every tournament with n <= 6 (uT11_identity.out).
  It upgrades THM-086 (kind-pasteur S37, "PROOF SKETCH + VERIFIED n <= 10") to a complete proof of its Theorem A
  (3 | c_j for j < val(n)) and of its Eulerian corollary (3 | A(n,k) implies 3 | F_k(T)). Step 2 is the principal
  specialisation U_T(1,1) of the Redei-Berge / Grinberg-Stanley function in power sums (GS-OCF bridge), read modulo 3.
  The failure for p = 5 noted in THM-086 is explained: modulo 5 the weights 2^(-|C|) depend on the cycle length.
session: mac-mini-2026-10-22-sierp
source: 05-knowledge/results/sierpinski_steiner_collatz_tournaments_20261022.md
scripts:
  - 04-computation/experiments/sierpinski_ocf_20261022/fwd_poly_mod3.c (+ fwd_poly_mod3.out)
  - 04-computation/experiments/sierpinski_ocf_20261022/uT11_identity.py (+ .out)
  - 04-computation/experiments/sierpinski_ocf_20261022/fwd_mod3.py (+ .out)
depends_on:
  - the OCF H(T) = I(Omega(T), 2) (Grinberg-Stanley) applied to every induced subtournament
related:
  - 01-canon/theorems/THM-086-universal-taylor-zeros-mod3.md
  - 01-canon/theorems/THM-085-f-omega-mod9-universal.md
  - 01-canon/theorems/THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md
---

# THM-4756 — modulo 3, every tournament is H copies of the transitive one

Let `T` be a tournament on a vertex set `V` with `|V| = n`. Let `F(T, x) = Σ_σ x^(fwd σ)`, the sum over all `n!` orderings `σ` of `V`, where `fwd σ` counts consecutive pairs `σ_i → σ_(i+1)`. Let `c_m(T)` be the number of ordered `m`-tuples of vertex-disjoint directed paths covering `V`, so that `c_1 = H(T)`.

**Step 1 (Taylor coefficients are path covers).**
* Since `x^f = Σ_j C(f, j)(x − 1)^j`, the coefficient of `(x − 1)^j` is `Σ_σ C(fwd σ, j)`.
* A pair (ordering, chosen `j` of its forward steps) is the same as an ordered list of `n − j` directed paths: cut at the unchosen steps.
* Hence `F(T, x) = Σ_(m=1)^n c_m(T)(x − 1)^(n−m)`.

**Step 2 (doubled Rédei–Berge identity, exact).** Put `H(∅) = 1` and apply the OCF to `T[X]` and to `T[V \ X]`.
* A pair of independent sets of odd cycles inside `X` and inside `V \ X` is the same as an independent set `R` of `Ω(T)` split into two parts. The vertices outside `V(R)` may lie on either side. So

      Σ_(X ⊆ V) H(T[X]) H(T[V\X]) = Σ_R 2^|R| · 2^|R| · 2^(n − |V(R)|) = Σ_R 4^|R| 2^(n − |V(R)|).

* More generally, the sum over ordered `k`-tuples of disjoint possibly empty sets covering `V` is `Σ_R (2k)^|R| k^(n − |V(R)|)`. This is `U_T(1^k)`.

**Step 3 (modulo 3).** Modulo 3, `4 ≡ 1`, `2 ≡ −1`, and every cycle has odd length, so `2^(n − |V(R)|) ≡ (−1)^n (−1)^|R|`.
* Step 2 is therefore `≡ (−1)^n I(Ω, −1) ≡ (−1)^n H(T) (mod 3)`.
* Removing the two terms `X = ∅, V` leaves `c_2(T) ≡ ((−1)^n − 2) H ≡ ((−1)^n + 1) H`.
* For `m ≥ 3`, `c_m = m!·(unordered covers) ≡ 0`.

So

    F(T, x) ≡ H(T)(x − 1)^(n−1)            (n odd),
    F(T, x) ≡ H(T)(x − 1)^(n−2)(x + 1)     (n even)      (mod 3).

For the transitive tournament (`H = 1`), `F` is the Eulerian polynomial `A_n`. Hence `F(T, x) ≡ H(T) A_n(x)` for every `T`. Checked independently against the Eulerian numbers for `n ≤ 40`.

**Modulo 2.** Every `c_m` with `m ≥ 2` is even and `H` is odd (Rédei). So `F(T, x) ≡ (1 + x)^(n−1)`: every tournament carries the binary Sierpinski row.

**Consequences.**
* For odd `n`, `F_k(T) ≡ H·(−1)^(n−1−k) C(n−1, k) (mod 3)`. The ternary Sierpinski (Lucas) pattern of Pascal's triangle, scaled by `H`, is the forward-step distribution of every tournament.
* At `n = 2·3^j + 1` it collapses to `F ≡ H(1 + x^(3^j) + x^(2·3^j))`.
* The pasted endpoint collapse at lengths `3^j + 1` does not occur: those `n` are even, and the even formula carries the factor `x + 1`.
