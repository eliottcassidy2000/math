# Sierpinski, Steiner, Collatz and Hamiltonian paths (2026-10-22, mac-mini-2026-10-22-sierp)

**Owner prompt.**
* How `3^k + 1` and `k^3 + 1` relate to Collatz `3n + 1`.
* A pasted note on the Steiner quasigroup `x ⋆ y = −x − y (mod 3)`: lengths `3^k + 1` collapse to their endpoints; a three-state refinement `S_a, S_b, S_c` of `H`; `H ≡ I(Ω, −1) (mod 3)`; profiles for 7 and 21; a "minimal genus four" analogy.
* `{25, 1}` against `{7, 21}`, the prominence of 27 in Collatz, and earlier "path of least resistance / natural recursion" work.
* Be bold.

**Lanes.**
* Session29 (`05-knowledge/results/session29_framed_profiles_20261022/`) took the Collatz side: powers of 3 and the odd parts of `3^k + 1` as cylinders, a test at 27, the profile hierarchy, and genus-four separators. That work is cited, not duplicated.
* This note is the tournament side.
* Canon:
  * THM-4754 (3-adic Sierpinski tournament and Collatz);
  * THM-4755 (blow-up gains through the skew spectrum, the lexicographic law, the ternary tower product);
  * HYP-9364 (the self-blow-up principle).
* Scripts: `04-computation/experiments/sierpinski_ocf_20261022/`.

**Status words** follow the repo: PROVED, CITED, FINITE-EXACT, DERIVED (a computation through a standard limit theorem, not a formal proof), NUMEROLOGY.

## 0. Summary

| item | status | where |
|---|---|---|
| the pasted reduction is the binomial average in `F_3` (`⋆` = midpoint), and collapses to the endpoints at depth `3^k` | PROVED (classical) | §1 |
| in tournaments the pattern lives in the forward-step polynomial: `F(T,x) ≡ H(T)·A_n(x) (mod 3)` for every tournament (`A_n` Eulerian), i.e. `H` times a Pascal row mod 3, collapsing to three terms at `n = 2·3^k + 1`; this completes THM-086 | PROVED (THM-4756), exhaustive `n ≤ 9` | §1, THM-4756 |
| the pasted list of profiles for `H = 21` contains two impossible profiles | corrected; the full lists are prior art | §2 |
| `x ↦ 3x + b` are the self-similarities of the 3-adic tournament `C_3^[∞]`, halving reverses every arc, and Collatz trajectories flip the arc between them once per halving | PROVED, new as a Collatz dictionary | §3, THM-4754 |
| the OCF weight `2^k = 3J_k + (−1)^k` is a one-step Collatz relation (Jacobsthal split) | PROVED | §3 |
| exact `H` of the ternary towers through 729 vertices; `v_3` law | FINITE-EXACT; values to 243 vertices are Session24's | §4 |
| a tower's blow-up gain is the level product `Π ((3^(2i+1)+1)/(3^(2i+1)−1))^(3^i)`, whose second factor is `((27+1)/(27−1))^3` | PROVED, given the determinant formula | §4, THM-4755 |
| Wormald's gain is `det(dI − S)/|det(idI − S)|`, and his constant is `cosh(1)/cos(1)` | Bhattacharjee 2026; rederived independently here | §5, THM-4755 |
| the tower's own normalised `H` converges to the same limit `κ_S`, to 7 digits | FINITE-EXACT plus extrapolation | §4, HYP-9364 |
| the same principle predicts `R_n → cosh 1/cos 1`, `P_q → e`, and `H(R_q)/H(P_q) → 1.0506` | data agree to 1e-3, 5e-5, 1e-3 | §6, HYP-9364 |
| number board for `{25, 1}`, `{7, 21}` and 27 | NUMEROLOGY, labelled | §7 |

## 1. The pasted Steiner reduction, and where it lives in tournaments

**The reduction.** In `F_3`, `1/2 = −1`, so `x ⋆ y = −x − y = (x + y)/2`: the third point of the line through `x, y` (midpoint).
* One row of the triangle multiplies by `−(1 + shift)`. A string of length `n` reduces to `(−1)^(n−1) Σ C(n−1, i) x_i = 2^(−(n−1)) Σ C(n−1, i) x_i`, the binomial (random-walk) average of the string.
* The pasted formula and the endpoint collapse at `n − 1 = 3^k` (from `(1 + z)^(3^k) = 1 + z^(3^k)`) are correct, and Lucas gives the converse.

**The pasted tournament statistic is weak.** It labels each vertex by out-degree mod 3 and reduces the label string of every Hamiltonian path.
* For a regular tournament all labels are equal, so every path gets the same value.
* At `n = 3^k + 1`, which is always even, it only counts paths by the labels of their two endpoints. It is an invariant, but carries little.

**Where the Sierpinski pattern really lives** (THM-086, kind-pasteur S37, a proof sketch verified for `n ≤ 10`, now proved in full by THM-4756).
* Let `F(T,x) = Σ_k F_k x^k` count permutations by their number of `T`-forward steps. For odd `n`, `F(T,x) ≡ H·(x − 1)^(n−1) (mod 3)`: a row of Pascal's triangle mod 3, scaled by `H`.
* So the forward-step distribution of every tournament carries the Lucas/Sierpinski zero pattern (Désarménien's factorisation of the Eulerian polynomial mod 3).
* At the orders `n = 2·3^k + 1 = 3, 7, 19, 55, 163, …` it collapses to three terms: `F ≡ H(1 + x^(3^k) + x^(2·3^k))`. Checked on all 456 classes at `n = 7` (`fwd_mod3.out`).
* These orders contain the Heegner numbers 3, 7, 19, 163 (NUMEROLOGY).
* At the pasted lengths `n = 3^k + 1` (even), THM-086 has two parameters, and the pattern is flat plus endpoints rather than a collapse.
* **THM-4756 (new) proves it all**, in the form `F(T,x) ≡ H(T)·A_n(x) (mod 3)`: modulo 3, every tournament is `H` copies of the transitive one.
  * The Taylor coefficients of `F` at `x = 1` are the ordered path-cover numbers `c_m`.
  * `c_m ≡ 0` for `m ≥ 3` (factor `m!`).
  * `c_2 ≡ ((−1)^n + 1)H`, from the exact identity `Σ_X H(T[X])H(T[V∖X]) = Σ_R 4^|R| 2^(n−|V(R)|)`.
  * Modulo 2 the same argument gives `F ≡ (1 + x)^(n−1)`, the binary Sierpinski row.
  * Exhaustive check on all 198,946 tournaments with `n ≤ 9`.

**`H` mod 3.** `H ≡ I(Ω, −1) (mod 3)`, and `I(Ω, −1)` is minus the reduced Euler characteristic of the independence complex of `Ω`. Both are prior art: HYP-1028..1030, THM-467 and T274. The pasted deletion recurrence written with `⋆` is the same identity.
* Prior work also rules out any congruence route to the 7/21 gaps. Session25 shows that every odd residue modulo every `M` is a strong `H`-value, and on strong classes of orders 5–8, `I(Ω, −1)` is never 0.

## 2. The profiles for 7 and 21

* Since `2^k ≡ (−1)^k`, any `α_3 ≥ 1` forces `α_2 ≥ 3` (the three pairs inside a disjoint triple) and `α_1 ≥ 3`.
* So the profiles `(6, 0, 1)` and `(4, 1, 1)` in the pasted list are impossible. The arithmetic list for `H = 21` is `(10, 0), (8, 1), (6, 2), (4, 3)` with `α_3 = 0`; for `H = 7` it is `(3, 0)` alone.
* The repo already has these and more:
  * THM-343 (`H = 7`, Lean-formalised);
  * THM-079, THM-115 and THM-1370 (all `n`, via Busch's strong floor `f = 3, 5, 9, 15, 25, 45` and the "crevasses" `7 ∈ (5, 9)`, `21 ∈ (15, 25)`);
  * the graph-level list in `camion_busch_gaps_polyhedra_collatz_20261002.md`: `K_3`; and `K_10`, `K_8 − e`, `K_6 − 2K_2` or `K_6 − P_3`, `P_4` or `K_3 + K_1`.
* Exhaustive data `n ≤ 9` (`analyze1.out`): no realized profile solves either equation.

## 3. The 3-adic Sierpinski tournament and Collatz (THM-4754)

**The tower.** `S_k = Cay(Z/3^k, {x : unit part of x ≡ 1 (mod 3)}) = C_3^[k]`, with `S_∞` on `Z_3`.
* Its three self-similarities are `x ↦ 3x + b`. So `3x + 1` is one of them, and its fixed point `−1/2 = …111` (base 3) is the all-ones address.
* Multiplying by a unit `≡ 2 (mod 3)` reverses every arc, so halving does.

**The arc-flip law.** If `x ≡ y (mod 2^t)`, then `C^t(y) − C^t(x) = 3^o (y − x)/2^e`, so the `S_∞`-arc between the two trajectories flips exactly once per halving.

**Reading.** "2 ≡ −1 (mod 3)" is one statement in several guises:
* halving is the converse operation on the 3-adic tournament;
* the OCF fugacity 2 becomes the Euler-characteristic fugacity −1;
* the Steiner rule is the midpoint.

The repo had no Pascal-mod-3 Collatz object (survey), so this dictionary is new. It is exact but proves nothing about Collatz convergence.

**The Jacobsthal split.** `H = 3 Σ J_k α_k + I(Ω, −1)` with `J_k = (2^k − (−1)^k)/3`.
* The even-index values `J_(2m) = (4^m − 1)/3 = 1, 5, 21, 85, 341` are exactly the odd numbers that `3n + 1` sends to a power of 4. In particular the forbidden `21 = J_6` goes to 64.
* THM-1425 already notes that the OCF on a path conflict graph produces Jacobsthal numbers (`I(P_4, 2) = 21`).

## 4. The ternary towers: exact counts and the product constant (THM-4755)

**Exact values.**
* `H(S_k)`: 3, 3159, 417382500592116023859, then 98- and 404-digit numbers at 81 and 243 vertices. These come from the C_3 substitution recursion (THM-3121; Session24 already had them).
* New here: a 1554-digit value at 729 vertices.
* `v_3(H(S_k)) = (3^k − 1)/2 = v_3|Aut|`, except `k = 2`, where it is one more.

**Normalised counts.**
* `κ(S_k) = H·2^(n−1)/n!`: 2, 2.228571, 2.572358, 2.713530, 2.764891, 2.782564.
* Successive differences shrink by about 1/3. Richardson in `1/n` gives 2.79150957576.

**The predicted limit.** By THM-4755, the blow-up gain of `S_k` is the level product, which tends to

    κ_S = Π_(i≥0) ((3^(2i+1) + 1)/(3^(2i+1) − 1))^(3^i) = 2 · (28/26)^3 · (244/242)^9 ⋯ = 2.79150988819834…,

* This is within 3.1e-7 of the extrapolation.
* The factors are exactly the numbers of the prompt, `3^(odd) ± 1`: `4/2`, then `28/26`, then `244/242`, …. The second factor is `(27 + 1)/(27 − 1) = (7 + 21)/(25 + 1)`, cubed.
* Its appearance is forced by the spectrum of `C_3[C_3]`: the within-block modes have `σ² = 3` against `d² = 81`.

**How the limit was first found.** The heuristic route was a renormalisation.
* The cut density `α` of path covers maps to the run density `(2 + α)/3` of the blocks, so `1 − α` shrinks by 1/3 per level.
* Each level contributes the Gaussian density ratio `(3 + (1−α)²)/(3 − (1−α)²)`.
* The determinant formula of §5 turns this into a theorem about blow-up gains.
* Wormald's Theorem 2 (`E[H(C_3[random blocks])] ∼ 2g`) is the first factor.

## 5. Wormald's constant and Bhattacharjee's formula

**Wormald's setting** (arXiv:2610.16245, posted 2026-10-12). Blow up a regular `T*` on `d` vertices by random tournaments. The gain `ρ_T(T*)` is `√(d/H)` for a Hessian determinant `H`.
* He conjectured that `d/H` is a rational square (Conjecture 6), that `ρ(R_d)` has a limit (Conjecture 2), and that `h(n)/g(n) → c ≈ 2.855958` (Conjecture 5).

**Bhattacharjee's release** (github.com/creelie/Breakthrough-problems-1, v1.0.0, 2026-10-14, unreviewed). It claims:
* `H = d det(d²I + S²)/det(d²I − S²)`, with `S` the skew adjacency matrix, settling Conjecture 6 modulo his Proposition 3.1;
* the cyclic closed form `((d+1)^d + (d−1)^d)/((d+i)^d + (d−i)^d)` and the limit `cosh(1)/cos(1)`, settling Conjecture 2;
* the doubly regular lower bound.

**This session found the same formula independently** before finding the release.
* The weighted block walk has transition matrix `(J + S)/d`, which is `S/d` on the mean-zero space for every regular `T*`. The local limit theorem then gives `ρ = Π (d² + σ²)/(d² − σ²)`.
* It reproduces every regular entry of Wormald's table exactly (`wormald_check*.out`).
* As `d → ∞` on `R_d` it gives `Π_(j odd) (1 + 4/π²j²)/(1 − 4/π²j²) = cosh 1/cos 1 = 2.8559578926`. His fits gave 2.85595787 and 2.85595789.

**New relative to the release:**
* the lexicographic law and the tower product (§4);
* the spectral moment expansion `log ρ = −Σ tr(S^(4m+2))/((2m+1) d^(4m+2))`;
* the observation that Wormald's six other 9-vertex "circulants" are not tournaments as printed and lie below the Jensen bound;
* the self-blow-up principle (§6).

## 6. The self-blow-up principle (HYP-9364)

**Statement.** For regular `T` on `n` vertices, `κ(T) = H(T)·2^(n−1)/n! = ρ(T)(1 + O(1/n))`.

**Evidence.**
* The tower: 7 digits, with deviation about `1.35/n`.
* All circulant classes: `κ/ρ ∈ [0.954, 0.982]` at `n = 25`, with both edges rising like `1/n`.

**It does not fix the ordering of classes** (rank correlation near 0). That ordering is set by the packing corrections of THM-4751 (7).

**Predictions, checked against exact data:**
* `H(R_n)·2^(n−1)/n! → cosh(1)/cos(1)`. Fits in `1/n` over `n ≤ 31` give 2.8536 and 2.8572, bracketing 2.85596. This contradicts the rotational half of HYP-9028, which predicts `e`.
* `H(P_q)·2^(q−1)/q! → e`. The 3-point fit gives 2.718332, against 2.718282 (HYP-9274 (A)).
* `H(R_q)/H(P_q) → cosh(1)/(e cos 1) = 1.050648`. The 3-point fit gives 1.04953 (HYP-9362 (c), now with a conjectured value).

**Picture.**
* `e` is the gain of a flat spectrum: Paley, and random regular tournaments.
* `cosh(1)/cos(1)` is the gain of the most concentrated spectrum (rotational), with `σ_j/n → 2/(πj)` on odd `j`.
* The ternary tower sits between them, at `κ_S`.
* In each case the exponent is a sum of `artanh(σ²/n²)` over the skew spectrum, with leading term exactly 1.

## 7. Number board (NUMEROLOGY unless stated)

* `3³ + 1 = 28 = 7 + 21` and `3³ − 1 = 26 = 25 + 1`. Their ratio `28/26` is the second factor of `κ_S`. This is exact and forced by the spectrum (§4), so the appearance is a theorem; the reading of `28` and `26` as the sums of the prompt's pairs is not.
* `3^k + 1 = k^3 + 1` only at `k = 3`.
  * `3³ + 1 = 28 = C(9)` is the Collatz image of `9 = 3²`.
  * The shortcut image is `U(9) = 7`.
  * `k³ + 1 = (k + 1)Φ_6(k)` and `Φ_6(k) = 1 + 2C(k, 2)`, so `7 = Φ_6(3) = (3³+1)/4` and `21 = Φ_6(5) = (5³+1)/6`.
  * Theorem P of `collatz_one_out_edge_systems_20260930.md` has the oddpart identity for `n³ + 1`.
* `{7, 13, 21} = {111_2, 111_3, 111_4}` are three iterations of `x ↦ bx + 1` from 0. The Collatz member `13 = (3x+1)^3(0)` is attained as `H`, and its neighbours are the two forbidden values. HYP-1043 (`Φ_3(2^k)`, which breaks at 73) and the F_21 note are prior art.
* `21` is a multiple of 3, so it has no Collatz predecessor under `U`; `U(21) = 1`.
* `H(R_7) = 175 = 7·25` and `H(P_7) = 189 = 7·27`. The per-vertex counts of the two vertex-transitive 7-tournaments are 25 and 27.
* **27 in Collatz.** 111 steps: 41 odd and 70 halvings. Prior analysis is THM-4504 and THM-4502 (the 41-tail family).
  * The 70 halvings equal the last coordinate of the Weyl vector `(0, …, 24 | 70)` of `II_(25,1)`, the cannonball `Σ_(i≤24) i² = 70²`. The repo had no such link (survey).
  * Exact tournament form: a regular tournament on `2N + 1` vertices has `c_3 = Σ_(i≤N) i²`, so by Watson `c_3` is a square only for regular tournaments on 1, 3 and 49 vertices, and `c_3(R_49) = 70²`.
  * The link to 27's halvings is NUMEROLOGY.

## 8. "Path of least resistance" and "natural recursion" in the repo

* "Least resistance" means series-parallel electrical reduction (`paley_least_resistance_runs_one_modular_form_20261014.md`).
* "Natural recursion" is the source-and-sink step `n → n + 2` (`07-reflections/the-two-staircases.md`).
* In this note the natural recursion is the C_3 substitution, and the "path of least resistance" is the weighted block walk: forward moves weigh 2, stays 1, backward moves 0. Its spectrum fixes the constant.

## 9. Open

* Prove HYP-9364 (a) from the exact recursion.
* Find the right error term in the principle.
* Find a statistic that orders circulant classes at fixed `n`.
* Push `H(R_33)`, `H(R_35)`.
* Decide THM-086 in full; it is still a proof sketch.

## 10. Later results of the same session

**THM-4756, second level.** The `k`-coloured identity `U_T(1^k) = Σ_R (2k)^|R| k^(n−|V(R)|)` at `k = 3` gives a congruence for the number `u_3` of unordered covers by three directed paths:
* `u_3(T) ≡ H(T) + hc(T) (mod 3)` for odd `n`, where `hc` counts directed Hamiltonian cycles;
* `u_3(T) ≡ 0 (mod 3)` for even `n`.

The reason is that modulo 9 only Hamiltonian cycles survive in `U_T(1,1,1)`. Checked on every tournament with `n ≤ 8` (`u3_congruence.out`). Deeper levels need cycle-length-graded data, so the clean mod-3 picture stops here.

**Wormald's Conjecture 1 for regular bases, d ≤ 13** (FINITE-EXACT; `wormald_conj1.out`, `wormald_conj1_13.out`). Over every regular `d`-tournament, the blow-up gain `ρ = Π (d² + σ²)/(d² − σ²)` has a unique maximum at the rotational `R_d`:

| `d` | regular classes | `ρ(R_d)` | runner-up |
|---|---|---|---|
| 5 | 1 | 2.3157894737 | — |
| 7 | 3 | 2.4630381803 | 2.3993783994 (Wormald's non-circulant example) |
| 9 | 15 | 2.5475143943 | 2.4981158008 (`C_3[C_3]`) |
| 11 | 1223 | 2.6021808258 | 2.5682251025 |
| 13 | 1,495,297 | 2.6404182020 | 2.6170585197 |

The minimum equals the doubly regular bound when a doubly regular tournament exists (`d = 7, 11`). Bhattacharjee's release does not treat Conjecture 1.

**The self-blow-up principle on all regular tournaments** (`principle_regular.out`). At `n = 11` (1223 classes), `κ/ρ ∈ [0.917, 0.980]`, the same band as for circulants. But the rank correlation between `κ` and `ρ` is strongly negative at `n = 9, 11` (−0.72, −0.74). This is the small-`n` regime where flat spectra win `H` (THM-4751 (7)). For circulants it is near 0 from `n = 15` on.

**The Sierpinski interval tournaments `Σ_N`** (the 3-adic tournament on `{0, …, N−1}`).
* **Recursion.** `Σ_N = C_3[Σ_(N_0), Σ_(N_1), Σ_(N_2)]` with `N_b = #{x < N : x ≡ b (mod 3)}`. So the C_3 substitution law computes `H(Σ_N)` for every `N` (`sierpinski_interval.py`, Held–Karp check for `N ≤ 14`).
* **The sequence** `1, 1, 3, 5, 11, 45, 159, 621, 3159, 14877, 79359, 515997, 3210655, 22181917, …` is in `sierpinski_interval_90.out` and is not in OEIS. Neither is the tower sequence `3, 3159, 417382500592116023859, …`.
* **Oscillation.** The normalised value has a base-3 oscillation. Multiples of 3 (balanced blocks) sit higher, the peaks are at `3^k`, and the pasted lengths `3^k + 1` sit in the dip right after each peak.
* **Automorphisms.** By nauty, `Aut(Σ_N)` is a 3-group for every `N ≤ 90`. Its exponent is mirror-symmetric about `3^k`: around 9 for offsets `t ≤ 6`, around 27 over the tested range, and around 81 for `t ≤ 9` (`interval_aut.out`). Since `|Aut|` divides `H`, `v_3(H(Σ_N)) = v_3|Aut(Σ_N)| + δ(N)` with `δ ≤ 5` for `N ≤ 90`; this explains the near-symmetry of `v_3(H)` about 27 and 81.

**The tower constant in Collatz form** (an exact rewriting; reading it as Collatz is NUMEROLOGY):

    κ_S = Π_(n = 1, 9, 81, 729, …) ((3n + 1)/(3n − 1))^(√n).

The numerators `4, 28, 244, 2188` are the Collatz images `3n + 1` of the even powers of 3, i.e. of the vertex counts of `S_0, S_2, S_4, S_6`.

**The automorphism groups of `Σ_N`, exactly** (PROVED, pre-audit; `interval_aut_recursion.out` agrees with nauty for every `N ≤ 90`).

*Recursion.* `|Aut(Σ_N)| = 3^(v(N))` with

    v(1) = v(2) = 0,    v(N) = v(N_0) + v(N_1) + v(N_2) + [3 | N].

* `C_3` is a prime quotient, so by modular decomposition (Gallai) every automorphism of `C_3[A, B, C]` permutes the three blocks by a rotation of `C_3`.
* A rotation is possible exactly when the blocks are isomorphic. Here that means equal sizes, i.e. `3 | N`.

*Mirror symmetry.* `v(3^k + t) = v(3^k − t)` for `0 ≤ t ≤ 2·3^(k−1)`, by induction on `k`.
* The two numbers split into residue-class offsets with the same multiset `{⌈t/3⌉, ⌈(t−1)/3⌉, ⌈(t−2)/3⌉}`.
* Each part is at most `2·3^(k−2)`, inside the previous level's range.
* The first failure is at `t = 2·3^(k−1) + 1 = 7, 19, 55, 163, 487` (checked through `k = 6`).

*Consequence.* Since `|Aut|` divides `H`, `v_3(H(Σ_N)) = v(N) + δ(N)` with `0 ≤ δ ≤ 5` for `N ≤ 90`. This explains the near-symmetry of the valuations of the Hamiltonian counts about 27 and 81.
