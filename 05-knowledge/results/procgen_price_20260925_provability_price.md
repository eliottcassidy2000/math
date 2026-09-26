# The price of provability: explicit L-step-descent trees at flip density at most 2ρ_L ≤ 2^(1−0.05L) (HYP-9136 proved), a 2^(−0.774L) lower bound, and the cycle −5 → −7 → −10 as the second obstruction

**Status.**
- **UPDATE 2026-09-26:** [THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md) proves the matching sharp exponent, `delta_L >= 2^(-(1-h)L-O(sqrt(L log L)))`, by critical growth bands and affine integer capacity. HYP-9137 is PROVED; the open-price language and weaker lower bounds below are historical. Collatz remains OPEN.
- **PROVED** (hand proofs below; the consequences of every lemma are also checked by code):
  - **(A) HYP-9136.** For every `L >= 8` there is an explicit member `G_L` of `P_L`. It is a tree whose only cycle is `{1,2}`, and its upper flip density is at most `2ρ_L <= 2^(1-(1-h)L)`, with `1 - h = 0.050044`. Hence `δ_L <= 2^(1-(1-h)L) -> 0`. This is the upper half of the coordinator's sharp conjecture, with the conjectured exponent.
  - **(B) The undecided density.** `ρ_L = |Bad_L|/2^L` satisfies `2^(-(1-h)L)/poly(L) <= ρ_L <= 2^(-(1-h)L)` (Chernoff bound and cycle lemma), where `h = h(log_3 2)`.
  - **(C) Lower bound.** Every member of `P_L`, with any threshold `n_0`, has upper flip density at least `ρ_L / W_L`. Here `W_L <= sum_(k<L) g_k ~ c λ^L` with `λ = (3+√13)/4`, so `δ_L >= 2^(-(0.7737+o(1))L)`.
  - **(D) SHEET.** The same construction and proof work on the `3n−1` sheet (constants change sign). Provable `3n−1` approximants exist at the same price, and they are trees: all three `3n−1` cycles are destroyed.
  - **(E) Structure.**
    - A rescue whose partner is forced up and then comes down in 2 steps is always available, except on the 2-adic neighbourhood `n ≡ −5 (mod 2^9)` of the bad cycle `−5 → −7 → −10`. There the partner needs 7–8 steps.
    - Every member of `P_L` must flip pairs at every 2-adic depth near every bad 2-adic cycle.
- **FINITE-EXACT.**
  - Exact `ρ_L` for `L <= 40`; `|Bad_20| = 27,328` reproduces the foundry and brackets value.
  - The densities of `G_L` for `L <= 40` on the pairs `<= 5·10^6`. The ratio `δ(G_L)/ρ_L` falls from `0.917` (`L = 8`) to `0.856` (`L = 40`).
  - The lemma checks: 66 runs, both sheets, every `L` in `8..40`.
  - The DRIFT data.
- **OPEN.**
  - The matching lower bound `δ_L >= 2^(-(1-h)L-o(L))`. The best proved exponent is `0.774`.
  - `δ_L` for `4 <= L <= 7`, where there is only finite-range evidence.
  - A positive floor for the price of the `5x+1` family (DRIFT). The evidence is that the construction's price tracks the `5x+1` undecided density, which tends to `0.176`.
- **REFUTED.** Nothing.

Session `collatz-procgen-20260922` (price lane, 2026-09-25). Scripts:
- `04-computation/experiments/procgen_price_20260925_{density,lower,run}.py`;
- the C program `procgen_price_20260925_greedy.c`.

Output: [procgen_price_20260925.out](procgen_price_20260925.out) (18 s, peak 314 MB). Parents:
- [THM-4470](../../01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md);
- [HYP-9136](../hypotheses/HYP-9136-provable-pairing-price-tends-to-zero.md);
- [pairings note §2.4](procgen_brackets_20260924_pairings_transitions.md);
- [choice ladder](collatz_procgen_20260922_choice_ladder.md).

## 0. Setting and results

`T(n) = n/2` (n even), `(3n+1)/2` (n odd). Pairs are `{2i−1, 2i}`, with one bit `ε_i` each. A number `n` in pair `i` goes up by `i` iff (`n` odd) XOR `ε_i`. So flipping pair `i` sends the odd member `2i−1` **down** to `i−1` and the even member `2i` (its **partner**) **up** to `3i`. `P_L` is the set of members in which every `n >= 3` falls below itself within `L` steps. `δ_L` is the infimum of the flip densities over `P_L`. `P_L ⊂ P_(L+1)`, so `δ_L` is non-increasing.

| quantity | value | status |
|---|---|---|
| `ρ_L` (Collatz-undecided density) | `2^(-(1-h)L+O(log L))`, `1-h = 0.050044`; exact values below | PROVED / FINITE-EXACT |
| upper bound | `δ_L <= 2ρ_L <= 2^(1-(1-h)L)`, `L >= 8`, via the explicit tree `G_L` | **PROVED** |
| measured price of `G_L` | `0.917 ρ_L` (L=8) → `0.856 ρ_L` (L=40) | FINITE-EXACT |
| lower bound | `δ_L >= ρ_L/W_L = 2^(-(0.774+o(1))L)` | PROVED |
| sharp lower bound `2^(-(1-h)L-o(L))` | — | OPEN |

**The exponent from the choice-ladder count.**
- `Bad_L` is the set of classes mod `2^L` whose parity word keeps `3^(a_k) > 2^k` for all `k <= L`.
- Its growth `|Bad_L| = 2^(hL+O(log L))` is the finite form of `dim_H Bad_inf = h(log_3 2) = 0.94995`.
- **Upper bound.** `ρ_L <= P(Bin(L,1/2) > L log_3 2) <= 2^(-(1-h)L)` (Chernoff, since `D(p||1/2) = 1 - h(p)` bits).
- **Lower bound.** By the cycle lemma, at least `C(L, ⌈pL⌉+1)/L` words stay strictly above the line.
- **Data.** `ρ_L L^(3/2) 2^((1-h)L)` stays in `[3.3, 11]` for `12 <= L <= 4000` and drifts slowly upward (DP to `L = 4000`), consistent with `ρ_L ≍ L^(-3/2) 2^(-(1-h)L)`.
- `ρ_8 = 0.0742`, `ρ_20 = 0.0261`, `ρ_40 = 0.00582`, `ρ_100 = 2.4·10^-4`.

## 1. Why the naive constructions fail (motivation, FINITE-EXACT)

- **Flip every bad odd `x` at step 0.** This rule is periodic mod `2^(L-1)` in the pair index, so the pair-0 obstruction applies (THM-4470 §5). The partner `x+1` then rides up the rising run.
- **Partner problem.** If the flipped point `v` is not `≡ 2 (mod 3)`, its partner `v+1` is visited by other orbits, which are deflected upward. Early prototypes left 1–2% of all `n` without an `L`-step descent.
- **The three facts that make a clean construction possible.**
  1. **Isolation.** An odd up-image `v = T(p)` is `≡ 2 (mod 3)`, so its partner `v+1 ≡ 0 (mod 3)`. A multiple of 3 is entered only from its double, or from a flipped even member `2j → 3j`. Such partners are visited only by themselves.
  2. **Processing order.** A flip made while certifying `m`, at a point below `2m`, sends every later visitor below `m`.
  3. **Freezing.** Once a certificate is frozen, no later flip can touch it.

## 2. The construction `G_L` and the theorem

**Construction.**
- Process `n = 3, 4, 5, ...` in order. Every pair is FREE or FROZEN (with a bit). The **default path** `D(n)` uses the frozen bits and reads free bits as 0, i.e. as Collatz.
- If `D(n)` descends within `L` steps, freeze the pairs of its points before the descent. This is `n`'s **certificate**.
- Otherwise `n` is **rescued**: set one free pair to 1 and freeze the new certificate. Take the first applicable option, in the order **A, F, B** if `n ≡ 3 (mod 8)` and **F, B, A** otherwise:
  - **A** = `pair(n)`: `n` goes down at once.
  - **F** = `pair(w)` for the first point `w` of `D(n)` that is `T(v)` of an odd `v` going up, with `w ≡ 3 (mod 4)` and `w < 2n`.
  - **B** = `pair(T(n))`: `n → T(n) → (3n−1)/4`.
- The program also has a partner rescue P and a depth-first fallback. The theorem shows that for `L >= 8` neither is ever used, and the program confirms this: `P = dfs = 0` for every `L` in `8..40`.
- The **owner** of a flip is the `n` whose rescue made it.

**Theorem.** Let `L >= 8`.
- Every `n >= 3` is certified, so the construction never gets stuck and `G_L ∈ P_L`.
- `ε_1 = 0`, so `G_L` is a tree with the single cycle `{1,2}`.
- The flipped pairs have upper density at most `2ρ_L`.
- Hence `δ_L <= 2ρ_L <= 2^(1-(1-h)L)`, and **HYP-9136 holds**.

**Proof.** The proof is by induction over the processing order. The inductive hypothesis is that all flips made so far are of types A, F and B.

**(G1) Suffix property.** A frozen certificate never changes. If a certificate from `s` passes `x` before descending, then `x` descends below `s < x` along it.

**(G2) Flip geometry.** A flip with owner `m` sits at an odd point `v`:
- `v = m` for A;
- `m < v < 2m` for F and B.

In all cases `(v−1)/2 < m`, and the pair index `i = (v+1)/2` lies in `(m/2, m]`.

**(G3) A partner is visited before descent only by itself.**
- **F and B.** Here `v = T(p)` with `p` on the owner's certificate going up, so `pair(p)` is frozen at 0.
  - The partner is `u = v+1 ≡ 0 (mod 3)`. Before descending, a path meets no flipped odd point (see G4).
  - So its only possible entries into `u` are from `2u`, or from the flipped even member `2u/3 = p+1`. The latter needs `pair(p)` flipped, which is impossible.
  - Going backwards along any path, the entries into a partner come from strictly larger numbers (`2u`, or a larger partner `2^k(p+1)`). So a path that visits `u` before descending started at `u` or above it.
- **A.** Here `u = m+1`.
  - Its odd preimage `(2m+1)/3 < m` would have been processed first and would have frozen `pair(m)`, so A would not have been available.
  - Later numbers `n' > m` satisfy `n' >= u`.

**(G4) Default paths.** `D(n)` follows `T` until it descends or meets a flipped point.
- At a flipped **odd** point `v` (owner `m < n`) it moves to `(v−1)/2 < m < n`, so it descends.
- By (G3), a flipped **even** point (a partner) is met only if `n` is that partner.
- Hence a rescued `n` is odd, and its Collatz path stays above `n` for `L` steps. So `n ∈ Bad_L` (up to finitely many small `n`), and `n ≡ 3 (mod 4)`.
- `n ≡ 1 (mod 4)` is never rescued: `T(n)` is even, and the odd member `T(n)−1 ≡ 1 (mod 3)` of its pair is never flipped before `n`.

**(G5) Odd rescues always exist.** If an odd `n` needs rescue, `pair(T(n))` is free, so B is available. A certifies `n` in 1 step and B in 2 steps. F is accepted only when its path descends within `L`; otherwise the construction moves on to B.
- *Why `pair(T(n))` is free.*
  - (a) `T(n)+1 = 3(n+1)/2` is entered only from its double (larger numbers) or from `n+1`. The latter would mean `pair(n)` is flipped, and then `n` goes down at once.
  - (b) If an earlier certificate passes `T(n)` at position `j >= 1`, then `D(n)` follows it and descends within `1 + L − j <= L` steps.
  - (c) `T(n)` cannot be an earlier A-point (its owner would be `T(n) > n`) or B-point (its owner would be `n` itself). If it were an F-point, then `n` would lie on that certificate.
- For `n ≡ 7 (mod 8)`, the first F-candidate is `w = T(n)` itself.

**(G6) Partners of A- and F-flips come down in 2 steps.** Here `v ≡ 3 (mod 4)`, and `u = v+1 → 3u/2 = T(v)+1 → 3u/4`.
- This needs `pair(T(v))` unflipped when `u` is processed.
- A flip there would need an owner `T(v) > u` (A), or `v` (B: but a flipped `v` descends at once), or an F-certificate through `v` going up.
- An F-certificate through `v` going up is impossible:
  - before the flip it would have frozen `pair(v)` at 0;
  - after the flip, `v` goes down.

**(G7) B happens only at `n ≡ −5 (mod 2^9)`.** Let `y_k = T^k(n)`. B is used only for `n ≡ 3 (mod 8)`, with A and F failed; `D(n)` is the Collatz path.
- Since `n` does not descend within `8 <= L` steps, the following hold (for `n >= 25`):
  - `y_3` odd, since otherwise `y_4 < n`;
  - `y_4` odd, since otherwise `y_5 = (27n+23)/32 < n`;
  - `y_6` odd, since otherwise `y_7 < n`;
  - `y_7` odd, since otherwise `y_8 = (243n+319)/256 < n`.
- If `y_5` were odd, then `y_4` (`< 2n`, `≡ 3 mod 4`) would be the first F-candidate, and it descends at step 5.
- If `y_8` were odd, then `y_7` would be, descending at step 8.
- Both candidates are free by (G8), so F fails only if `y_5` and `y_8` are even.
- So the parity word starts `110110110`, i.e. `n ≡ 507 (mod 512)`, the class of `−5`.
- The only `n < 25` with `n ≡ 3 (mod 8)` are 3, 11 and 19, and they descend by default.

**(G8) The candidates `y_4` and `y_7` are free.** A frozen pair would require one of the following.
- **`y_k + 1` visited earlier.** Excluded as in (G5a): its flipped entry `y_(k−1)+1` would make `D(n)` descend at `y_(k−1)`.
- **An earlier flip of `y_k`.** A and B are impossible, and F reduces to the next case.
- **A visit to some `y_i` (`i <= k`) by an earlier certificate from `s < n` at position `j`.**
  - If `j >= i`, then `D(n)` descends within `L` along that certificate. This is a contradiction, since `n` needs rescue.
  - If `j < i` (a **late visit**), look at the first common point `z = y_(i0)` of the two paths. There the two paths arrive through different predecessors.
    - `z ∈ {y_3, y_6}` (`≡ 1 mod 3`) has only one preimage below `2n`.
    - `z ∈ {y_1, y_2}` would force `s >= n`.
    - For `z = y_4` or `y_5`, the path of `s` would have to reach `2y_4 ≈ 3.375n` in 2 steps, or `2y_5 ≈ 5.06n` in 3. Each step multiplies by at most `3/2 + 1/(2s)`, so this is impossible.
    - For `z = y_7`, the path of `s` must reach `2y_7 = (243n+319)/64` by 4 or 5 up-steps. The start may be an odd `s`, or a partner `s` with first step `3s/2`. The endpoints are `(81s+c)/16` with `c ∈ {65, 38}`, and `(243s+c')/32` with `c' ∈ {211, 130}`. Equality forces `81(4s−3n) ∈ {59, 167}` or `243(2s−n) ∈ {−103, 59}`, which is impossible.

**(G9) Partners of B-flips come down within 8 steps.** The partner `u = y_1 + 1 = (3n+3)/2` follows

`u → y_2+1 → 3y_3+2 → 3y_4+2 → 3y_5+2 → 3y_6+1 → y_7 → y_8 → y_9`,

where `y_9 = (729n+1085)/512 < u`. These identities follow from the prefix `110110110`.
- For example, `T(y_2+1) = (3y_2+4)/2 = 3y_3+2`.
- The path re-enters `n`'s own orbit at `y_7`, the `−7`-type point of the next lap of the `−5` cycle.
- The even points `3y_5+2`, `3y_6+1` and `y_8` are not partners (G3). Each odd point either moves by `T` or is flipped, and a flipped odd point makes the path descend at once (G4).
- So `u` descends within 8 steps, and within 7 if `y_7` was itself flipped by the rescue of `y_6`.

**Conclusion.**
- **Certification.** Odd `n` is certified by default or by (G5). Even non-partners descend in 1 step. Partners descend in 2 steps (G6) or at most 8 steps (G9). So neither P nor the fallback is ever called, and by induction all flips are of types A, F and B.
- **Tree.** No path from `n >= 3` visits 1 or 2 before descending, and every flipped pair index is `>= 2`. So pair 1 is never touched and `ε_1 = 0`.
- **Density.** Each rescued `n` lies in `Bad_L` (with finitely many exceptions) and flips one pair `i ∈ (n/2, n]`. So `#{flipped i <= X} <= |Bad_L ∩ [1, 2X+1)| + O_L(1)`. ∎

**Remarks.**
- `G_L` is explicit and computable: `ε_i` is final after processing `n <= 2i+1`. Pairs `<= N/2` are final once `n <= N` has been processed.
- The bound `2ρ_L` is crude. The measured density is `0.86–0.92 ρ_L` (table in §4). About 35% of `Bad_L` numbers ride on earlier certificates and need no flip of their own; for example, at `L = 16` there are 211,180 rescues against about 322,570 bad `n <= 10^7`.

## 3. The lower bound (PROVED)

**Proposition.** Let `M ∈ P_L` with any threshold `n_0`. Let `W_L = max_v sum_(k<L) sum_(n: T^k n = v) 3^(a(n,k))/2^k`, where `a(n,k)` is the number of odd steps among the first `k`. Then the upper flip density of `M` is at least `ρ_L/W_L`.

**Proof.**
- Take `n >= n_0` in `Bad_L`. If no pair of `n, T n, ..., T^(L−1) n` is flipped, then `M` follows `T` for `L` steps and `n` does not descend. So the first flipped point `v(n)` exists.
- Let `A` be the set of members of flipped pairs; its upper density equals the flip density `δ̄`. Weight `n` by `1/n = (v/n)/v`, with `v/n <= (1+ε) 3^a/2^k` for large `n`. Then

  `sum_(n in Bad_L, n <= X) 1/n <= (1+ε) W_L · sum_(v in A, v <= (3/2)^L X + O(1)) 1/v`.

- The left side is `(ρ_L+o(1)) ln X`. The right side is `(1+ε) W_L (δ̄+o(1)) (ln X + O(L))`. Hence `δ̄ >= ρ_L/W_L`.
- Logarithmic density absorbs the fact that a flip at height `R` serves `n` at `v/R`.

**The bound on `W_L`.** In the backward tree of `T`:
- a node `≡ 2 (mod 3)` has two children, `2v` (`≡ 1`) and `(2v−1)/3` (any residue);
- a node `≡ 1` has only `2v`;
- a node `≡ 0` has only `2v`.

Hence the weighted count satisfies `g_k <= 1.5 g_(k−1) + 0.25 g_(k−2)`, with `g_0 = 1` and `g_1 = 2`. So `W_L <= sum_(k<L) g_k ~ c λ^L`, where `λ = (3+√13)/4 = 1.65139`.

**Result.** The lower-bound exponent is `(1−h) + log_2 λ = 0.7737`.

| L | `ρ_L` | bound `sum g` | `δ_L >=` | construction `δ(G_L)` |
|---|---|---|---|---|
| 8 | 0.0742 | 99.3 | `7.5·10^-4` | 0.0680 |
| 16 | 0.0323 | 5602 | `5.8·10^-6` | 0.0285 |
| 24 | 0.0171 | `3.1·10^5` | `5.5·10^-8` | 0.0149 |
| 40 | 0.00582 | `9.5·10^8` | `6.1·10^-12` | 0.00499 |

The brute-force maximum of the weighted count over `v <= 20000` is well below the tree bound: 80.6 against 99.3 at `L = 8`.

**Why the bound is weak.** A flip at a point `v` whose backward tree is "rich" could, in principle, serve exponentially many bad orbits, and flips at great heights are cheap. A proof of the sharp exponent must show two things:
- **Hubs are rare.** The multiplicity is `O(1)` on average; a second-moment or 3-adic equidistribution argument would be needed.
- **High flips do not finish the job.** After a flip at height `R >= 2`, the orbit is still at height `R/2 >= 1`.

The coordinator's heuristic ("each flip repairs boundedly many classes") holds for the construction: each flip repairs its owner and the orbits riding on its certificate. It is not proved as a bound over all members.

## 4. Data (FINITE-EXACT, `N = 10^7`, pairs `<= 5·10^6`, no stuck `n`, every `n <= 10^7` descends within `L`)

| L | `δ(G_L)` | `ρ_L` | ratio | 3n−1 sheet | ratio | rescues A / F / B / P |
|---|---|---|---|---|---|---|
| 4 | 0.2134 | 0.1875 | 1.138 | 0.2115 | 1.128 | 123812 / 884442 / 318384 / 212257 |
| 5 (= 6) | 0.1267 | 0.1250 | 1.014 | 0.1208 | 0.966 | 70144 / 712602 / 85273 / 56851 |
| 7 | 0.0961 | 0.1016 | 0.946 | 0.0902 | 0.888 | 50602 / 605813 / 41725 / 6523 |
| 8 | 0.06805 | 0.07422 | 0.917 | 0.06326 | 0.852 | 32988 / 458779 / 10525 / 0 |
| 12 | 0.04921 | 0.05518 | 0.892 | 0.04577 | 0.829 | 21916 / 338605 / 3992 / 0 |
| 16 | 0.02847 | 0.03226 | 0.883 | 0.02648 | 0.821 | 11728 / 198214 / 1238 / 0 |
| 20 | 0.02279 | 0.02606 | 0.874 | 0.02122 | 0.814 | 9136 / 159127 / 898 / 0 |
| 24 | 0.01488 | 0.01708 | 0.871 | 0.01389 | 0.813 | 5745 / 104278 / 516 / 0 |
| 32 | 0.008368 | 0.009627 | 0.869 | 0.007781 | 0.808 | 3162 / 58744 / 280 / 0 |
| 40 | 0.004985 | 0.005823 | 0.856 | 0.004753 | 0.816 | 1865 / 35218 / 170 / 0 |

**The table.**
- The rescue columns are for the Collatz sheet.
- Equal rows reflect `P_5 = P_6`, `P_8 = P_9`, and so on (the Terras lengths).
- For `L <= 7` the partner rescue P is used; the theorem does not cover these `L`, and there is no stuck `n` to `10^7`.

**Lemma checks.** Code in the C program, every `L` in `8..40`, both sheets, `N = 2·10^6`; and `L ∈ {8,12,16,24,32,40}` at `N = 10^7`. Every rescued `n` is odd, lies in the bad residue mod 4, and has a Collatz path that does not descend within `L`. B occurs only at `n ≡ 3 (mod 8)` with prefix `110110110`. Every flipped index lies in `(n/2, n]`.

| partner descent time | 2 | 7 | 8 | other |
|---|---|---|---|---|
| A/F partners, L = 8 | 333,207 | — | — | 0 |
| B partners, L = 8 | — | 108 | 6,908 | 0 |

**Cross-check.** An independent Python re-implementation reproduces the C densities exactly for `L ∈ {4,8,12,16,20}`, `N = 2·10^5` (PART 4).

## 5. Where the price is paid

**The two obstructions, `−1` and `−5`.**
- **`−1`.** Near `−1` (the rising runs, `n ≡ −1 mod 2^t`) the rule "F at `T(n)`" handles the run completely. Each run point `≡ 7 (mod 8)` flips its successor, and the partner comes down in 2 steps.
- **`−5`.** In the neighbourhood `n ≡ −5 (mod 2^9)` of the bad cycle `−5 → −7 → −10` (word `(110)^∞`, factor `(9/8)` per lap):
  - no odd up-image `≡ 3 (mod 4)` exists;
  - `A` is blocked by `(2n+1)/3`.
  So the construction must flip the `−7`-point, whose partner is not 2-step. The partner comes down because its orbit re-merges with the cycle one lap later (G9).
- **Where F never works.** Bad words with no factor `111` at all grow like `2^(0.3L)` (script part 3b). On them no F-rescue ever exists. `B` is needed only when `A` is also blocked, which forces the `−5` class.

**Obstruction lemma (PROVED).**
- **Statement.** Let `c` be a 2-adic periodic point of `T` with more than `p log_3 2` odd steps per period `p` (a bad cycle), rotated so that its word is prefix-bad. For each depth `D >= L`, every `n ≡ c (mod 2^D)` is in `Bad_L`, so some pair `pair(T^k n)` with `k < L` is flipped.
- **Consequence 1.** The flip set of any member of `P_L` accumulates 2-adically at the pair indices of every bad cycle: `−1`, `−5`, the `−17` cycle, `(110)^k 10`, and so on.
- **Consequence 2.** No member of `P_L` has a flip set that is periodic, or locally constant near all of these points. This generalises THM-4470 §5, whose obstruction is the case `c = −1`.

**The coordinator's heuristic.** Flips are needed "essentially only on Bad classes". This is exact for `G_L`: every flip has a `Bad_L` owner, and the partner flips never needed for `L >= 8`.

## 6. Controls

**SHEET (3n−1; offset 1; pairs `{2i, 2i+1}`; `U(n) = (3n−1)/2` for odd `n`).**
- **Transfer.** The negation `T(−n) = −U(n)` carries the construction over to `U`:
  - pairs go to pairs;
  - up-images of `U` are `≡ 1 (mod 3)`;
  - partners `w−1 ≡ 0 (mod 3)` are isolated;
  - F uses `w ≡ 1 (mod 4)`, B is for `n ≡ 5 (mod 8)`, and the prefix becomes `110110110`, i.e. `n ≡ 5 (mod 512)` at the actual `3n−1` cycle `5 → 7 → 10`;
  - the constants `c_w` change sign, and the obstructions `81 ∤ ±59, ±167` and `243 ∤ ±103, ±59` persist.
- **Rescue A is never available** on this sheet: the partner `n−1` of `n` is processed first and freezes `pair(n)`. The proof only uses "A blocked", so it is unaffected. So the theorem holds for `U` with the same bound `2ρ_L`; `ρ_L` is identical by the Terras bijection under negation. The lemma checks pass for every `L` in `8..40`.
- **What this means.**
  - Provable `3n−1` approximants exist at price `<= 2^(1-(1-h)L)`, measured at `0.81–0.85 ρ_L`.
  - They are **trees**: with threshold 2, the cycles `{5,7,10}` and the two cycles with minima 17 are destroyed.
  - `3n−1` has three cycles, so the price of provability cannot see cycles, any more than the pair sums of THM-4470 can see tree-ness. A finite set of flips is density-free.

**DRIFT (5x+1 family: up `v → (5v+[v odd])/2`, down `v → ⌊v/2⌋`; not sum-preserving).**
- **Undecided density.** `β_L` falls to a positive limit: `0.273, 0.222, 0.206, 0.196, 0.181, 0.176` at `L = 8, 20, 30, 40, 100, 400`.
- **Construction.** The same sequential construction, with A and a search of up to 8 flips (N = 3·10^5):

  | L | price | price / `β_L` | stuck `n` |
  |---|---|---|---|
  | 8 | 0.265 | 0.97 | 326 |
  | 12 | 0.196 | 0.80 | 8 |
  | 20 | 0.149 | 0.67 | 1 |
  | 30 | 0.131 | 0.64 | 0 |

- **Reading.** The price tracks the undecided density, which does **not** tend to 0. For `L <= 20` the repertoire is not even complete (there are stuck `n`), so these are not members of `P_L`.
- **No floor is proved.** The weighted-multiplicity bound gives only `β_L/W_L → 0`, so "the price stays bounded away from 0" is **evidence, not a theorem**.
- **The contrast.** The AM–GM fairness of the `3x+1` pairing (THM-4470 §1) is exactly what makes the drift negative, the undecided density tend to 0, and hence the price tend to 0.

## 7. Meaning

- **Where Collatz sits.** Collatz is at flip-density distance **0** from divergent members (THM-4470 §4) **and** from provable trees (this note). Both are now theorems. The conjecture sits on the common boundary of two dense regions of the pairing cube.
- **Is "the price decays like `2^(−0.05L)`" a theorem?**
  - The **upper** half is a theorem: `δ_L <= 2^(1−0.050044 L)`.
  - The exponent `1−h` is also the exact exponent of the undecided density, and the construction pays `0.86–0.92` of it.
  - The **lower** half is not: the proved exponent is `0.774`.
  - "The price decays exactly like `2^(−(1−h)L+o(L))`" remains a conjecture, supported by the construction and by the averaged heuristic `E[multiplicity] = O(1)`.
- **What this does not do.**
  - Nothing here touches Collatz itself. Density-zero changes flip its truth value (DEFECT).
  - `G_L` is provable by construction, because its orbits are forced down on `Bad_L`. It says nothing about orbits Collatz itself leaves undecided.
- **What it does.**
  - It locates the cost of provability exactly: one flip per Collatz-undecided number that no earlier certificate already covers.
  - It shows that the only non-local repair needed is at the second bad cycle `−5`, after the first bad cycle `−1` (the pair-0 obstruction).

## 8. Reproduction

`python3 04-computation/experiments/procgen_price_20260925_run.py > 05-knowledge/results/procgen_price_20260925.out`

This compiles `procgen_price_20260925_greedy.c` into a temporary directory, takes about 18 s and peaks at 314 MB. The parts are:
- (1) `procgen_price_20260925_density.py`: exact `ρ_L`, `L <= 40` exact and `<= 4000` in floating point; the compression `E[1/Rmax]`.
- (2) `G_L` on both sheets for `L <= 40` at `N = 10^7`, the lemma checks, and DRIFT.
- (3) `procgen_price_20260925_lower.py`: `λ`, the `W_L` bound and its brute-force check, the no-111 count, `β_L`.
- (4) The Python cross-check.

Single runs:
- build with `cc -O2 -o greedy 04-computation/experiments/procgen_price_20260925_greedy.c`;
- run `./greedy L N offset q [MAXF]`, with `offset` 0 or 1 and `q` = 3 or 5.

## 9. Hypothesis bookkeeping (no files created)

**HYP-9136 → PROVED** by §2: `δ_L <= 2^(1-(1-h)L)` for `L >= 8`, with an explicit tree. It is a candidate for promotion to a THM after independent audit, especially of (G3), (G7), (G8) and (G9). Its secondary claim `δ_2 = 0.2907…` is untouched.

New candidates, for the coordinator to number:
- **(Q1) Sharp lower bound.** `δ_L >= 2^(-(1-h)L-o(L))`. Equivalently: no member of `P_L` beats `ρ_L` by an exponential factor. Proved exponent `0.774`.
- **(Q2) Constant.** `δ(G_L)/ρ_L` converges, with data `0.917 → 0.856`. Whether `δ_L/ρ_L` stays bounded below is (Q1) with polynomial loss.
- **(Q3) DRIFT floor.** `inf_L δ_L^(5x+1) > 0`. The evidence is that the price is about `0.64 β_L` and `β_L → 0.176`.
- **(Q4) Small L.** `G_L ∈ P_L` for `4 <= L <= 7`. Verified to `10^7`, where P-rescues are used; not proved.
