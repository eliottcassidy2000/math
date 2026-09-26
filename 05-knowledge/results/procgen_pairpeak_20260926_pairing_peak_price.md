# The pairing price at the peak: one high flip re-merges, a reflected barrier does not; the private pairing price is ρ_L·exp(−Θ(L^(1/3))); the consistent price δ_L (HYP-9140) stays OPEN

Lane `pairpeak`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_pairpeak_20260926_{lib,run}.py` and `procgen_pairpeak_20260926_greedy.c`.
Output: [procgen_pairpeak_20260926.out](procgen_pairpeak_20260926.out) (one runner; every `ok:` line is a `check(...)` that aborts on failure; 100 checks, `ALL CHECKS PASSED`).

## Status

**Summary.** HYP-9140 (`δ_L <= poly(L) ρ^peak_L`) is **not** proved, and neither is the weaker `δ_L <= ρ_L exp(−c L^(1/3))`. What is proved:
- a single high flip does not work (it re-merges with probability 1/2);
- many flips at stopping times of the flipped orbit do work, each costing a factor `λ = 0.585`;
- hence the pairing price is peak-discounted **when every undecided source may use its own flips** (the *private* price `π_L`).

The remaining obstruction is consistency between sources. In the tested constructions it takes the form of blocking by earlier frozen certificates.

- **PROVED** (full hand proofs in §§1–6):
  - **Lemma 1 (coupling of a single flip).**
    - Setup: `v` is odd, `x_t` is the Collatz orbit of the flipped image `(v−1)/2`, and `y_t` is the Collatz orbit of `T(v) = (3v+1)/2`.
    - Then `y_t = 3^(a_t) x_t + b_t`, with an explicit four-case transducer started at `(a,b) = (1,2)`.
    - If the Collatz word of `y_1` begins `1^r 0 0`, the two orbits **merge** at time `r+3`. At a stopping time this has Haar probability exactly `1/2`, and `1−c = 0.369` under the tilted measure `Q`.
    - Consequence: post-flip freshness (Lemma F of the peak note) holds only for the flipped continuation on its own; conditioned on the source being undecided it fails. A single flip at a high point of an undecided orbit does not rescue it whenever the continuation begins `1^r00` and the orbit stays above `3n+2` for `r+2` steps.
    - This corrects the "Why it is plausible" paragraph of HYP-9140.
  - **Lemma 2 (bijection).** Suppose the flip decision at time `t` depends only on the modified parity prefix, and may use the slope level `S'_t = U_t − ct`. Then the modified parity word of `n` is a bijective image of `n mod 2^L`; this is iterated Lemma F. With the tilt, each flip multiplies the weight of a surviving word by `λ = (1−c)/c = 0.585`.
  - **Theorem A (reflected barrier).**
    - The process `Π^(m)` flips every odd point of its own orbit at slope level `>= m−1`.
    - It leaves undecided a proportion `N_m(L)/2^L <= 2^(−(1−H)L) ρ_b^floor(L/n_m)` of the classes, where `ρ_b = 1−(1−λ²)/14 = 0.953013` and `n_m = ceil(2(m+1)²/σ²)`.
    - No analysis of the transducer is needed.
  - **Theorem B (private pairing price).**
    - Definition: `π_L` = least average cost `lim sup X^(−1) Σ_(n<=X) min_F Σ_(i∈F) n/i` when each `n` may choose its own set `F` of flipped pairs (a private member), with no consistency between different `n`.
    - `2ρ^peak_L <= π_L <= 2·3^(1−m) L ρ_L + 2N_m(L)/2^L` for every `m >= 2`, and a multi-scale bound (§4).
    - `π_L <= 2^(−(1−H)L) exp(−(γ−o(1)) L^(1/3))` with `γ = 0.3574`.
    - With THM-4480: `π_L = ρ_L exp(−Θ(L^(1/3)))`. **P1 is answered negatively for the private pairing price.**
  - **Theorem C (conditional).** Assume the *Robin inequality* `N_m(L) <= K_0 A_(m+K)(L)` for all `m`, where `A_M(L)` counts undecided words with peak `< 3^M`. Then `π_L <= (K_0 3^K (27L+18) + 4L 2^(−L)) ρ^peak_L`, i.e. HYP-9140 holds for the private price.
  - **Lemma 6 (domination).** Suppose every flipped odd point `v` is an up-image with an unflipped predecessor pair and satisfies one-flip domination `T^(j−1)((v−1)/2) <= T^j(v)` for `1 <= j <= L`. Then every `n` that is not a flipped partner satisfies `F^t(n) <= T^t(n)` until it descends, so no Collatz-good `n` is harmed.
- **FINITE-EXACT.**
  - Lemma 2 and the DP: for `L ∈ {8,10,12,14}` and `m ∈ {2,…,5}`, all `2^L` classes give distinct modified words. The survivors equal the exact DP count `N_m(L)`, and slope survival is equivalent to actual non-descent. The private member reproduces the `Π^(m)` path on all `2^14` classes (`L = 14`, `m = 3`).
  - Theorem A against the exact integer DP: for every `1 <= L <= 180` and every `m`, the maximum ratio is `0.518`.
  - Robin inequality with `K = 1`: `max N_m(L)/A_(m+1)(L) = 1.003255` over `1 <= L <= 180` and all `m`, attained at `(L,m) = (50,7)`. With `K = 2` the maximum is `1.0000889` at `(58,9)`. So `K_0 = 1` is false, but `K_0 = 1.0033` holds on this range.
  - Merge law: exactly `2^14 − 1` of the `2^15` odd classes mod `2^16` begin `1^r00` (`r+2 <= 15`).
  - Consistent construction (MODE 0): valid sections of members of `P_L` for all `3 <= n <= 10^6` (`L = 8, 10, …, 40, 48, 56, 64`) and `n <= 10^7` (`L = 16, 32`). There is no stuck `n`, and every rescued `n` is Collatz-bad (no harm).
  - Simulator check: MODE 2 reproduces THM-4475's `G_L` densities `0.068052, 0.049190, 0.028426, 0.022694` (`L = 8, 12, 16, 20`) exactly.
- **VERIFIED** (float DP, relative error `~1e−13`):
  - Theorem A holds at `L = 256, 512, 1024` for `m <= 40`.
  - The Robin ratio is `<= 1.0005` at `L = 256, 512, 1024` for `m <= 60`.
  - The rigorous multi-scale bound on `π_L` is `1.26 ρ_L` (`L = 256`), `0.143 ρ_L` (`L = 512`) and `0.0065 ρ_L` (`L = 1024`).
- **EMPIRICAL.**
  - Single flip at the highest odd point of an undecided orbit:

    | `L` | rescues | re-merges | other |
    |---|---|---|---|
    | 16 | 0.334 | 0.395 | 0.271 |
    | 32 | 0.250 | 0.548 | 0.202 |
    | 48 | 0.199 | 0.636 | 0.165 |

  - Sampled private multi-scale barrier cost: `6.5, 8.5, 9.1 × 2ρ^peak` at `L = 32, 64, 128`, which is `0.41, 0.15, 0.029 × ρ_L`.
  - The consistent construction's density is `2.0 → 3.4 → 5.3 × 2ρ^peak` (`L = 8, 32, 64`), i.e. `0.35 → 0.11 → 0.046 × 2ρ_L`.
  - Interference: `8–11%` of rescues (`L >= 32`) lose their private best scale, and `>= 95%` of these are blocked by an earlier frozen certificate.
  - One-flip domination fails for `0.20, 0.29, 0.36` of random barrier-type points (`L = 16, 32, 64`).
- **OPEN.**
  - HYP-9140 for the consistent price `δ_L`.
  - Even `δ_L <= ρ_L exp(−c L^(1/3))` for some `c > 0`.
  - The Robin inequality for all `L`.
  - The sharp constant of `π_L`: conjecturally `κ_3 = 2.108`, as for `ρ^peak`.
  - Whether the barrier construction of §7 (or the peak lane's greedy) is never stuck for `n > 10^7`.
- **CORRECTED.** HYP-9140 says "after a flip … the continuation is a fresh Collatz word … it falls … with probability tending to 1". For undecided sources this is false. The flipped and unflipped continuations are coupled by Lemma 1 and re-merge with probability `1/2` per flip. The single-flip rescue rate *decreases* with `L` (table above).
- Nothing here bears on the Collatz conjecture. These are fixed-horizon modification prices. No novelty or priority claim.

## 0. Setting and the three prices

`T(n) = n/2` (n even), `(3n+1)/2` (n odd); `c = log_3 2`, `σ² = c(1−c) = 0.232857`, `λ = (1−c)/c = 0.584963`, `1−H = 1 − H_2(c) = 0.050044`. Words `u ∈ {0,1}^L`, `S_j = e_j − jc`, `w_j = 3^(S_j)`; `Bad_L`, `ρ_L = |Bad_L|/2^L`, `w*(u) = max_(j<L) w_j`, `ρ^peak_L = 2^(−L) Σ_(Bad_L) 1/w*` (peak note §1). `Q` = iid Bernoulli(`c`) letters; tilt identity `2^(−L) = Q(u) 2^(−(1−H)L) λ^(S_L)` (peak note §4).

Pairing family (THM-4470/4475): pairs `{2i−1, 2i}`; flipping pair `i` sends `2i−1 → i−1` and the partner `2i → 3i`. `P_L`: every `n >= 3` falls below itself within `L` steps; `δ_L` = least upper density of flipped pair indices over `P_L`.

| price | definition | known |
|---|---|---|
| arbitrary edits `ε_L` | edits may send a point anywhere | `ρ^peak/M*_L <= ε_L <= ρ^peak_L` (THM-4480) |
| **private pairing** `π_L` | each `n` gets its own pairing member (§4) | `2ρ^peak_L <= π_L <= 2^(−(1−H)L) e^(−(0.3574−o(1))L^(1/3))` (**this note**) |
| consistent pairing `δ_L` | one member for all `n` | `ρ^peak/M*_L <= δ_L <= 2ρ_L` (THM-4475/4480) |

In the arbitrary-edit setting consistency is free, because an edit sends everything to 1. In the pairing family a flip only halves. The flipped orbit then has to be carried down (§§1–3), and it is shared with every other orbit through the same points (§§6–7).

## 1. Why one high flip is not enough: the coupling lemma

**Lemma 1 (coupling).** Let `v` be odd, and put `x_1 = (v−1)/2` (the flipped image) and `y_1 = (3v+1)/2 = T(v)`. Let `x_(t+1) = T(x_t)`, `y_(t+1) = T(y_t)`, and `(a_1, b_1) = (1, 2)`. While `a_t >= 1`, define:

| case | `(a_(t+1), b_(t+1))` |
|---|---|
| (E) `y_t`, `x_t` both even | `(a_t, b_t/2)` |
| (O) `y_t`, `x_t` both odd | `(a_t, (3b_t + 1 − 3^(a_t))/2)` |
| (U) `y_t` odd, `x_t` even | `(a_t + 1, (3b_t+1)/2)` |
| (D) `y_t` even, `x_t` odd | `(a_t − 1, (b_t − 3^(a_t−1))/2)` |

Then, up to and including the first `t` with `a_t = 0`:
- `y_t = 3^(a_t) x_t + b_t`, and
- `x_t ≡ y_t + b_t (mod 2)`.

*Proof.*
- **Start.** `y_1 = 3(v−1)/2 + 2 = 3x_1 + 2`.
- **Parity.** `3^a x ≡ x (mod 2)`.
- **Step.** Write `y = 3^a x + b`.
  - (E): `y/2 = 3^a (x/2) + b/2`.
  - (O): `(3y+1)/2 = 3^a (3x+1)/2 + (3b + 1 − 3^a)/2`.
  - (U): `(3y+1)/2 = 3^(a+1)(x/2) + (3b+1)/2`.
  - (D): `y/2 = 3^(a−1)(3x+1)/2 + (b − 3^(a−1))/2`.
- **Exactness of the divisions.**
  - In (E) and (O), `b ≡ y − x ≡ 0`; in (O) also `3b + 1 − 3^a ≡ 0`.
  - In (U) and (D), `b` is odd. ∎

**Corollary 1.1 (merging).** Suppose the parities of `y_1, y_2, …` begin `1^r 0 0` (`r >= 0`). Then `x_(r+3) = y_(r+3)`, and the two orbits coincide from then on. If they begin `1^r 0 1`, then `(a_(r+3), b_(r+3)) = (2, 2)`.

*Proof.*
- While `y_t` is odd and `(a,b) = (1,2)`: `b` is even, so `x_t` is odd, and (O) returns `(1, (6+1−3)/2) = (1,2)`.
- At the first even `y` (time `r+1`): `x` is even, and (E) gives `(1,1)`.
- At time `r+2`, `b = 1` is odd:
  - if `y` is even, then `x` is odd and (D) gives `(0,0)`;
  - if `y` is odd, then (U) gives `(2,2)`. ∎

**Corollary 1.2 (a single flip on an undecided orbit).** Let `n` satisfy `T^j(n) > n` for `1 <= j <= L`. Consider the orbit that follows `T` up to `v = T^k(n)` (odd, `k <= L−1`), then goes to `(v−1)/2`, and then follows `T`. This is the single-flip member, as long as the orbit does not return to the pair of `v`; that holds for all but finitely many `n` per class (Lemma 5).

Assume:
- the Collatz word of `T^(k+1)(n)` begins `1^r00`;
- `T^(k+j)(n) >= 3n+2` for `1 <= j <= min(r+2, L−k)`.

Then this orbit does not fall below `n` within `L` steps.

*Proof.*
- For `j <= r+2`: `x_j = (y_j − b_j)/3` with `b_j ∈ {1,2}` and `y_j = T^(k+j)(n) >= 3n+2`, so `x_j >= n`.
- For `j >= r+3`: `x_j = y_j = T^(k+j)(n) > n`. ∎

**Probabilities.**
- At a stopping time `k` (the class mod `2^(k+1)` fixed), the letters of `T^(k+1)(n)` are uniform (Terras).
- So `P(1^r00 for some r) = Σ_r 2^(−r−2) = 1/2`.
- Under `Q`, the analogous sum is `Σ_r c^r (1−c)² = 1 − c = 0.369`.
- The runner checks the exact count `2^14 − 1` over the odd classes mod `2^16` (§1 of the output).

**Data.** The single flip is placed at the highest odd point among the first `L` points (the peak lane's first candidate):

| `L` | rescued | re-merged | other failure |
|---|---|---|---|
| 16 | 0.334 | 0.395 | 0.271 |
| 32 | 0.250 | 0.548 | 0.202 |
| 48 | 0.199 | 0.636 | 0.165 |

The peak candidate's next letter is forced to be `0` (the peak is even). It therefore merges when the letter after that is `0`, which for undecided orbits happens more often than `1/2`.

This explains the peak lane's observation that the accepted flip moves down the candidate list (mean rank `1.3 → 4.2`). The fix is not a better single flip but many flips.

## 2. Prefix-determined flip rules: bijection and tilt

A **flip rule** decides, at each time `t`, whether to flip the current point. The decision is a function of `t` and the *modified* parity prefix `(u''_0, …, u''_t)`, and a flip is allowed only if `u''_t = 1`. The modified orbit is:
- `x_0 = n`;
- `x_(t+1) = (x_t − 1)/2` if flipped, `T(x_t)` otherwise;
- `u''_t = x_t mod 2`;
- `U_t` = number of unflipped odd steps before `t`;
- `J_t` = number of flips before `t`;
- slope level `S'_t = U_t − ct`, which may be used by the rule.

**Lemma 2 (bijection).** For every `L`, the map `n mod 2^L ↦ (u''_0, …, u''_(L−1))` is a well-defined bijection `Z/2^L → {0,1}^L`. On the class of a word, `2^t x_t = 3^(U_t) n + B_t` with an integer `B_t` that depends only on the word.

*Proof.* We argue by induction on `t`.
- **Inductive hypothesis.** On the class `K` of `n mod 2^t` with a given prefix, the decisions before `t` are fixed, and `2^t x_t = 3^(U_t) n + B_t`.
- **Splitting the class.** For `n = r + 2^t s`, `x_t = x_t(r) + 3^(U_t) s`. The factor `3^(U_t)` is odd, so the parity of `x_t` alternates with `s`. Hence `K` splits into its two subclasses mod `2^(t+1)`, carrying `u''_t = 0` and `u''_t = 1`.
- **The next step.** On each subclass the decision at `t` is fixed, and `2^(t+1) x_(t+1)` equals:
  - `3^(U_t) n + B_t` (even step),
  - `3^(U_t+1) n + 3B_t + 2^t` (unflipped odd step), or
  - `3^(U_t) n + B_t − 2^t` (flipped odd step).
- **Conclusion.** Distinct classes separate at their first differing bit, and both sides have `2^L` elements. ∎

The decision may use `S'_t`, since `S'_t` is a function of the prefix. It may not use the actual size of `x_t` or any future letter.

**Lemma 3 (tilt).** For every word, `e(ω) − cL = S'_L + J_L`, where `e(ω) = U_L + J_L` is the number of ones. Hence, for every event `A` determined by the modified word,

    |{n mod 2^L : ω(n) ∈ A}| / 2^L = 2^(−(1−H)L) E_Q[ λ^(S'_L + J_L) ; A ],

and each flip multiplies the tilted weight by `λ = 0.585`. (Apply the peak note's identity to the modified word.) ∎

This is the mechanism that replaces "freshness". The flipped continuation need not be independent of the source's undecidedness. A surviving modified word with `J` flips simply carries `J` extra factors `λ`.

## 3. The reflected barrier (Theorem A)

`Π^(m)` (`m >= 1`) flips at time `t` iff `u''_t = 1` and `S'_t >= m−1`. **Survival** `F_L = {S'_t > 0 for 1 <= t <= L}`; `N_m(L) := |F_L|` (words). Put `R_t = e_t − ct = S'_t + J_t` (the walk of the letters).

**Lemma 4.**
- (a) `S'_t < m − c` for all `t`.
- (b) At a *zone* state (`S'_t >= m−1`) both letters give `S'_(t+1) = S'_t − c`.
- So the exact count is a DP over `U`: zone steps are counted twice, every other step branches. (`N_m(L)` is computed exactly for `L <= 180`, and in floating point beyond.)

*Proof.*
- (a) Induction: from `S' < m−1` an up-step gives `S' + 1 − c < m − c`; from the zone the step is down.
- (b) An odd point in the zone is flipped; an even point is halved. ∎

**Theorem A.** For all integers `m >= 1` and `L >= 1`, with `n = n_m = ceil(2(m+1)²/σ²)`,

    N_m(L)/2^L  <=  2^(−(1−H)L) ρ_b^floor(L/n),        ρ_b = 1 − (1 − λ²)/14 = 0.953013.

*Proof.*
1. **Tilt.** By Lemma 3 and `S'_L > 0` on `F_L`, `N_m(L)/2^L <= 2^(−(1−H)L) E_Q[λ^(J_L) 1_(F_L)]`.
2. **Blocks.** Put `a = m+1` and consider the blocks `[kn, (k+1)n]` for `0 <= k < floor(L/n)`. On `F_L`, `S'` lies in `(0, m−c)` at the block ends (at `k = 0` the start is `0`), so `Δ_k S' ∈ (−(m−c), m−c)`. Since `Δ_k R = Δ_k S' + Δ_k J` with `Δ_k J >= 0`:
   - `Δ_k R > −a`;
   - if `Δ_k R >= a`, then `Δ_k J > a − (m−c) = 1 + c`, so `Δ_k J >= 2`.
3. **Pointwise bound.** Hence `λ^(J_L) 1_(F_L) <= Π_k φ(Δ_k R)`, where `φ(x) = 1{|x| < a} + λ² 1{x >= a}`.
4. **Expectation.** Under `Q` the block increments of `R` are iid, with the law of `S_n`. So
   `E φ(S_n) = 1 − Q(S_n <= −a) − (1−λ²) Q(S_n >= a) <= 1 − (1−λ²) Q(|S_n| > a) <= 1 − (1−λ²)/14`
   by Lemma C(4) of the peak note (`a >= 1`, `nσ² >= 2a²`; PROVED and audited there). ∎

No property of the transducer (merge probabilities, the parity of `b`, the growth of `a`) is used. A merge simply means that the reflected orbit returns to the barrier, where it is flipped again and pays another `λ`.

**Checks.**
- The DP equals the integer enumeration (Lemma 2 check).
- Theorem A holds on `1 <= L <= 180` for every `m` (exact; max ratio `0.518`), and at `L = 256, 512, 1024` for `m <= 40` (float; ratios `2.4e−3, 8.3e−4, 3.1e−4`).
- The bound is far from sharp; for example, at `L = 180`:

| `m` | 2 | 3 | 4 | 5 | 6 | 8 | (`ρ_L`) |
|---|---|---|---|---|---|---|---|
| `N_m(180)/2^180` | `2.0e−23` | `1.7e−11` | `2.0e−8` | `3.6e−7` | `1.5e−6` | `4.7e−6` | `7.4e−6` |

## 4. The private pairing price (Theorem B)

**Definition.**
- For `n >= 3`, a *private certificate* is a finite set `F` of pair indices such that, in the pairing member with exactly the pairs of `F` flipped, the orbit of `n` falls below `n` within `L` steps.
- Its cost is `c(F) = Σ_(i∈F) n/i`, and `c*(n) = min_F c(F)`. We have `c*(n) = 0` if `n` falls under `T`, and `c*(n) < 2` always, since odd `n` can flip its own pair.
- The **private price** is `π_L = lim sup_(X→∞) X^(−1) Σ_(3<=n<=X) c*(n)`.

Why this normalisation: a consistent construction whose flips are each owned by one source has flip density `X^(−1) Σ_(owners n <= X) Σ_(i ∈ F_n) n/i` (a flip at relative index `i/n = r` is counted below `X` iff `n <= X/r`). So `π_L` is the price with no sharing of flips (riders) and no interference.

**Lemma 5 (realization and descent).** Fix a flip rule as in §2 and a word `ω`. For all but finitely many `n` in the class of `ω`:
- (i) the pairing member whose flipped pairs are the pairs of the flipped points of `n`'s modified orbit reproduces that orbit for `L` steps;
- (ii) `x_t < n ⟺ S'_t < 0` for `1 <= t <= L`;
- (iii) a flip at time `t` has `n/(pair index) <= 2·3^(−S'_t) (1 + o(1))` as `n → ∞`.

*Proof.*
- **Affine form.** `x_t = 3^(S'_t) n + β_t` with `β_t = B_t/2^t` fixed (Lemma 2, `3^(U_t)/2^t = 3^(S'_t)`). The slopes `3^(S'_t)` are pairwise distinct, because `S'_s = S'_t` with `s ≠ t` would make `c` rational.
- (i) A deviation needs `x_t − x_s ∈ {0, 1}` for a flip time `s ≠ t` (the pair of an odd `x_s` is `{x_s, x_s+1}`). This is a linear equation in `n` with nonzero slope, so it has at most two solutions per `(s,t)`.
- (ii) Let `δ = min_(1<=t<=L) |S'_t| > 0`. If `S'_t <= −δ` then `x_t <= 3^(−δ)n + |β_t| < n` for large `n`, and symmetrically in the other case.
- (iii) The pair index is `(x_t+1)/2 >= (3^(S'_t)n − |β_t|)/2`. ∎

**Theorem B.** For every `L >= 1`:

1. `π_L >= 2ρ^peak_L`.
2. For every `m >= 2`: `π_L <= 2·3^(1−m) L ρ_L + 2 N_m(L)/2^L`. For every `m_top >= 3`:
   `π_L <= 2L [3^(1−m_top) ρ_L + Σ_(m=2)^(m_top−1) 3^(1−m) N_(m+1)(L)/2^L] + 2 N_2(L)/2^L`.
   In the single-scale bound, `L ρ_L` may be replaced by `E_m`, the expected number of flips that `Π^(m)` makes before its slope walk dies (over all classes, by the same DP), because the flips of an undecided source before its descent are among them. The table below uses `min(L ρ_L, E_m)`.
3. `lim sup_(L→∞) L^(−1/3) ln(2^((1−H)L) π_L) <= −γ`, where `γ = (3/2)(ln(1/ρ_b) σ²)^(1/3) (ln 3)^(2/3) = 0.3574`.
4. Hence `exp(−(C_3+o(1)) L^(1/3)) <= π_L/ρ_L <= exp(−(γ−o(1)) L^(1/3))`, with `C_3` the peak note's elementary constant, and `π_L = ρ_L exp(−Θ(L^(1/3)))`.

*Proof.*
1. **Lower bound.**
   - For `n` in a bad class, a certificate must flip the pair of some point `T^k(n)` with `k <= L−1`: otherwise `n` follows `T` for `L` steps.
   - The first such pair has index `<= (T^k(n)+1)/2`. So `c*(n) >= 2n/(max_(j<L) T^j(n) + 1) >= 2n/(w*(u)(n + L/3) + 1)` by the peak note §2.1.
   - Averaging over `n <= X` gives `2·2^(−L) Σ_(Bad_L) 1/w* = 2ρ^peak_L`. The exceptional set is finite.
2. **Upper bounds.**
   - For `n` in a bad class (odd, and off a finite set), use the flips of `Π^(m)` made before the first descent if its slope walk dips below `0` by time `L`. By Lemma 5 this is a certificate of cost `<= 2·3^(1−m) J (1+o(1)) <= 2·3^(1−m) L (1+o(1))`, since flips happen only at `S' >= m−1`.
   - Otherwise use the own pair, at cost `< 2`.
   - By Lemma 2 the classes on which `Π^(m)` survives number `N_m(L)`. This gives the single-scale bound.
   - **Multi-scale.** For bad `n` let `m(n)` be the largest `m ∈ [2, m_top]` for which `Π^(m)` brings `n` down. Then either `m(n) = m_top`, or `Π^(m(n)+1)` survives on its class (at most `N_(m(n)+1)(L)` classes), or no `m` works (at most `N_2(L)` classes, cost `< 2`).
3. **Exponent.**
   - Insert Theorem A with `m_top = L`.
   - The top term is `2L·3^(1−L)ρ_L`. The bottom term is `2 N_2(L)/2^L <= 2·2^(−(1−H)L) ρ_b^floor(L/78)`. Both are `e^(−Ω(L))` relative to `2^(−(1−H)L)`.
   - In the middle, the scales `m <= log L` contribute `exp(−Ω(L/log² L))`.
   - For `m > log L`, `n_(m+1) = (2/σ²)(m+2)²(1+O(log^(−2)L))`. Writing `m+2 = θL^(1/3)`, the term is `exp(−L^(1/3)[θ ln 3 + ln(1/ρ_b)σ²/(2θ²)](1+o(1)))`.
   - Minimising over `θ` gives `γ`, and there are at most `L` terms.
4. **Ratios.** Use `ρ_L >= 2^(−(1−H)L)/(2.72L(L+1))` (`L >= 10`, peak note Consequence (b)) and, for the lower side, part 1 together with the peak note's Theorem 3. ∎

**Numbers.** Float DP; the lower bound is `2ρ^peak_L`, and the upper bounds are those of part 2 with exact `N_m(L)`:

| L | `ρ_L` | `2ρ^peak_L` | single-scale bound | /`ρ_L` | multi-scale bound | /`ρ_L` |
|---|---|---|---|---|---|---|
| 32 | 9.63e−3 | 6.02e−4 | 1.45e−2 | 1.51 | 1.65e−1 | 17.1 |
| 64 | 1.49e−3 | 2.60e−5 | 1.86e−3 | 1.25 | 1.77e−2 | 11.9 |
| 128 | 6.98e−5 | 2.25e−7 | 7.38e−5 | 1.06 | 3.62e−4 | 5.18 |
| 256 | 3.28e−7 | 1.12e−10 | 3.06e−7 | 0.936 | 4.13e−7 | 1.26 |
| 512 | 1.62e−11 | 2.85e−16 | 4.96e−12 | 0.307 | 2.31e−12 | 0.143 |
| 768 | 1.33e−15 | 2.81e−21 | 1.75e−16 | 0.132 | 3.57e−17 | 0.0269 |
| 1024 | 1.15e−19 | 4.31e−26 | 7.17e−21 | 0.0622 | 7.48e−22 | 0.0065 |

**How good the bounds are.**
- The rigorous bounds carry the crude factor `L` (flips per source `<= L`), so they beat `ρ_L` only from `L ≈ 256` on.
- The elementary exponent is reached slowly: `ln(2^((1−H)L) B_L)/L^(1/3) = −0.160, −0.303, −0.343` at `L = 10^6, 10^8, 10^10`.
- The actual private cost of the multi-scale barrier (flips before descent only) was sampled over uniformly random bad classes (`1500, 1500, 700` samples):

| `L` | private cost `/2ρ^peak` | `/ρ_L` | fallbacks | mean flips |
|---|---|---|---|---|
| 32 | 6.48 | 0.405 | 2/1500 | 1.43 |
| 64 | 8.49 | 0.148 | 0 | 1.58 |
| 128 | 9.06 | 0.029 | 0 | 1.74 |

**Reading.**
- On the tested range the private price tracks `ρ^peak`, not `ρ_L`.
- The constant `6.5 → 9.1` grows slowly. It comes from integer scales, from flipping at the first barrier point rather than at the peak, and from about 1.6 flips per source.

## 5. The Robin inequality and the conditional private HYP-9140 (Theorem C)

Let `A_M(L) = #{u : S_t > 0 (1<=t<=L), S_t < M (0<=t<=L−1)}`. These are the bad words with `w* < 3^M`, so `3^(−M) 2^(−L) A_M(L) <= ρ^peak_L`.

As sets of letter sequences, both `A_M(L)` and the survivors of `Π^(m)` are contained in `Bad_L`: the unreflected walk is `R = S' + J >= S'`. The Robin inequality asks that the reflected (partially absorbing) top boundary cost at most a constant compared with a hard wall one unit higher.

**Conjecture R.** `sup_(L>=1, m>=1) N_m(L)/A_(m+1)(L) < ∞`.
- FINITE-EXACT: the supremum over `L <= 180` is `1.003255`, at `(50,7)`.
- VERIFIED: `<= 1.0005` at `L = 256, 512, 1024`, `m <= 60`.
- For `K = 2` the maximum is `1.0000889`, so `K_0 = 1` fails for `K = 1, 2`.

**Heuristic.**
- In the tilted measure a zone visit has weight `c·λ + (1−c) = 2(1−c) = 0.738`, a Robin boundary with absorption `0.262` per visit.
- In the scaling limit such a boundary is Dirichlet, shifted by `O(1)`. So `N_m(L) ≈ A_(m+O(1))(L)`.
- A proof needs a sharp comparison of the two lattice problems. The elementary block constants lose a constant factor in the exponent (§8).

**Theorem C.** Suppose that for some `L`, some `K >= 0` and some `K_0 >= 1`, `N_m(L) <= K_0 A_(m+K)(L)` for `2 <= m <= L`. Then

    π_L  <=  ( K_0 3^K (27L + 18) + 4L 2^(−L) ) ρ^peak_L.

*Proof.* Use Theorem B(2) with `m_top = L`, and write `M(u) = max_(t<L) S_t(u)`, so `3^(−M(u)) = 1/w*(u)`.
1. **Middle scales.** A word `u` is counted in `A_(m+1+K)(L)` iff `M(u) < m+1+K`. Hence
   `Σ_(m=2)^(L−1) 3^(1−m) A_(m+1+K)(L) <= Σ_(u∈Bad_L) Σ_(m > M(u)−1−K) 3^(1−m) <= Σ_u (3/2)·3^(2+K−M(u)) = 13.5·3^K·2^L ρ^peak_L`.
2. **Bottom scale.** `N_2(L) <= K_0 A_(2+K)(L) <= K_0 Σ_u 3^(2+K−M(u)) = 9K_0 3^K 2^L ρ^peak_L`.
3. **Top scale.** Since `w* <= (3/2)^(L−1)`, `ρ_L <= (3/2)^(L−1) ρ^peak_L`. So `2L·3^(1−L)ρ_L <= 4L·2^(−L) ρ^peak_L`. ∎

So Conjecture R implies HYP-9140 for the private price, with polynomial factor `O(L)`, and it would give the sharp constant of `π_L`: `ln π_L = −(1−H)L ln 2 − κ_3 L^(1/3)(1+o(1))` modulo Mogul'skii (THM-4480 Theorem 4).

## 6. Consistency I: the domination lemma

Harm to a Collatz-good `n` is the first consistency issue: an earlier flip deflects its orbit and the deflected orbit comes down later. The following lemma gives a sufficient condition that is local to each flip point.

Say that `v` satisfies **1FD(v, L)** (one-flip domination) if `T^(j−1)((v−1)/2) <= T^j(v)` for `1 <= j <= L`.

**Lemma 6 (domination).** Let `E` be a set of flipped pairs, and let `F_E` be the member. Assume:
- every flipped odd point `v` is an up-image, `v = T(p)` with `p` odd and `pair(p) ∉ E`;
- every flipped odd point `v` satisfies 1FD(v, L).

Then:
- (a) **Isolation.** An `F_E`-orbit that starts at `s` and visits the even member of a flipped pair at a time `t >= 1` before falling below `s` must have started at such an even member.
- (b) **Domination.** For every `n` that is not such an even member, `F_E^t(n) <= T^t(n)` for `t <= min(L, τ(n))`, where `τ(n)` is the first descent time under `F_E`.
- In particular every `n` that falls under `T` within `L` also falls under `F_E` within `L`.

*Proof.*
- **(a) Isolation.**
  1. *Preimages of `y` under `F_E`.* They are:
     - `2y`, if pair `y` is unflipped;
     - `2y+1`, if pair `y+1` is flipped;
     - `(2y−1)/3`, if it is odd with unflipped pair;
     - `2y/3`, if it is even with flipped pair.
  2. *An even member of a flipped pair.* Such a `u = v+1` is `3(p+1)/2 ≡ 0 (mod 3)`, and its preimages reduce to `2u`:
     - `(2u−1)/3` is not integral;
     - `2u+1 ≡ 1 (mod 3)` is not an up-image, so it is never a flipped point;
     - `2u/3 = p+1` has `pair(p) ∉ E`.
  3. *The doubles `2^j u`.* Each is again `≡ 0 (mod 3)`. Its preimages are `2^(j+1) u`, or `2^(j+1) u/3` when that number is itself an even member of a flipped pair. The other candidates are excluded: `2^(j+1) u + 1 ≡ 1 (mod 3)` is never flipped, and `(2^(j+1) u − 1)/3` is not integral.
  4. *The chain.* Let `t >= 1` be the first time the orbit is at an even member `u`. By minimality no earlier point (except possibly `s` itself) is an even member. So, going backwards from `t`, the orbit is at `2u, 4u, …, 2^t u` at times `t−1, t−2, …, 0`, unless `s` is itself an even member. In the first case `s = 2^t u > u`, so the orbit is below `s` at time `t`.
- **(b) Domination.** We prove by induction on `k`: for every flipped odd point `w` and every `i <= k` such that the orbit of `w` visits no even member of a flipped pair at times `1..i`, `F_E^i(w) <= T^i(w)`.
  1. Put `x_1 = (w−1)/2`. The orbit of `w` follows `T` from `x_1` until it meets a flipped odd point `w'` at some time `s`.
  2. For `i < s`: `F_E^i(w) = T^(i−1)(x_1) <= T^i(w)` by 1FD.
  3. For `i >= s`: `F_E^i(w) = F_E^(i−s)(w') <= T^(i−s)(w')` by induction, and `T^(i−s)(w') = T^(i−1)(x_1) <= T^i(w)`.
  4. For a source `n`: its orbit follows `T` up to the first flipped odd point `w` (at time `t_0`). Then apply the claim to `w`. By (a), no even member is met before descent. ∎

1FD is not rare to fail. It fails for `0.200, 0.292, 0.361` of random `v ≡ 11 (mod 12)` (`L = 16, 32, 64`). The first violations have a heavy tail (the transducer's `a_t` is a lazy symmetric walk under Haar).

A construction that demanded 1FD at every barrier flip (MODE 3 in the scratch runs) rejected many barrier points. Its fallback flips violate 1FD, so it was still not valid (§8). Lemma 6 therefore explains the mechanism of harm, but it is not the tool that removes it.

## 7. Consistency II: the consistent barrier construction (FINITE-EXACT) and interference

**Construction.** `procgen_pairpeak_20260926_greedy.c`, MODE 0; THM-4475's framework.
- **Processing.** Process `n = 3, 4, …`. The default path uses frozen bits and reads free bits as Collatz. If it descends within `L`, freeze its pairs; otherwise `n` is rescued.
- **Barrier rescue.** Try the scales `W = 12, 11, …, 2`. At scale `W`, follow frozen bits, and flip a FREE odd point `x` when:
  - `x` is an up-image, so its partner `x+1 ≡ 0 (mod 3)` is isolated;
  - `x ≡ 3 (mod 4)`, so the partner goes `x+1 → 3(x+1)/2 → 3(x+1)/4` and descends in 2 steps;
  - `x >= n·3^(W−1)`;
  - `pair(T(x))` is protected at 0.
  Accept the rescue if `n` descends within `L`, and freeze the certificate.
- **Fallback.** THM-4475's A/F/B, then P, then a depth-first search.

This construction realizes the barrier of §3 at the actual heights, with the partner rules of THM-4475. It is designed to be analysable, not optimal: the peak lane's greedy is about 2× cheaper at `L = 32`.

**Validity (FINITE-EXACT).**
- For `L = 8, 10, …, 40, 48, 56, 64` and every `3 <= n <= 10^6`, there is no stuck `n`, and a final re-run shows that every `n` falls below itself within `L`.
- The same holds for `L = 16, 32` at `N = 10^7`, where the density moves from `0.013534 → 0.013578` and `0.002050 → 0.002045`.
- Every rescued `n` is Collatz-bad, i.e. no harm was observed.
- That it never gets stuck for all `n` is not proved. THM-4475's lemmas G5–G9 would have to be re-proved in the presence of barrier flips at level `>= 1`, which can be hit late within the first 9 steps.

**Data** (`N = 10^6`, density of flipped pairs `<= N/2`):

| L | density | `2ρ_L` | `2ρ^peak` | /`2ρ^peak` | /`2ρ_L` | peak-lane greedy | rescues | barrier | fallback | degraded |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 0.052052 | 0.148438 | 0.026206 | 1.99 | 0.351 | 0.031274 | 46430 | 14585 | 31845 | 0 |
| 12 | 0.032246 | 0.110352 | 0.013363 | 2.41 | 0.292 | 0.017886 | 32257 | 14812 | 17445 | 0 |
| 16 | 0.013534 | 0.064514 | 0.005322 | 2.54 | 0.210 | 0.006902 | 16918 | 11393 | 5525 | 231 |
| 20 | 0.009512 | 0.052124 | 0.003414 | 2.79 | 0.183 | 0.004658 | 13055 | 9626 | 3429 | 320 |
| 24 | 0.005008 | 0.034163 | 0.001667 | 3.00 | 0.147 | 0.002284 | 8027 | 6587 | 1440 | 442 |
| 28 | 0.003288 | 0.026260 | 0.001032 | 3.19 | 0.125 | 0.001542 | 5906 | 5047 | 859 | 410 |
| 32 | 0.002050 | 0.019254 | 0.000602 | 3.40 | 0.107 | 0.000908 | 4118 | 3626 | 492 | 313 |
| 40 | 0.000928 | 0.011647 | 0.000252 | 3.68 | 0.080 | − | 2294 | 2085 | 209 | 208 |
| 48 | 0.000432 | 0.007075 | 0.000109 | 3.98 | 0.061 | − | 1321 | 1214 | 107 | 135 |
| 56 | 0.000272 | 0.004536 | 0.000052 | 5.26 | 0.060 | − | 809 | 743 | 66 | 93 |
| 64 | 0.000138 | 0.002986 | 0.000026 | 5.30 | 0.046 | − | 527 | 495 | 32 | 57 |

(The full table, including `L = 10, 14, …, 38` and the rider counts, is in the output. At `L >= 56` the density rests on fewer than 150 flipped pairs.)

**What the table shows.**
- **Sharing.** Only 35–52% of the undecided `n <= N` need their own rescue (`L = 16..64`); the rest ride earlier certificates. For example, at `L = 32` there are 4118 rescues against `ρ_L N = 9627`. This makes the consistent price cheaper than the private sum over the rescued set up to `L = 48`: the "private cost" column of the output is `0.002474` against a density of `0.002050` at `L = 32`.
- **Interference.** A rescue is *degraded* when its actual best scale is below its private best scale (pure Collatz plus its own barrier rule).
  - Degraded rescues: 1.4% (`L = 16`), 5.5% (24), 7.6% (32), 9.1% (40), 10.2% (48), 10.8% (64).
  - At least 95% of them (100% at `L >= 48`) are **blocked**: the desired flip point lies on an earlier frozen certificate.
- **The expensive fallbacks** (A/F/B/P/dfs) per undecided source:
  - MODE 0: `0.171, 0.084, 0.051, 0.036, 0.030, 0.021` (`L = 16, 24, 32, 40, 48, 64`).
  - At `L = 64` all 32 are interference-caused: the private process never needs a low rescue there.
- **Diagnostic MODE 1** (do not freeze default paths of Collatz-good `n`; not a valid member):
  - The fallbacks drop to `0.174, 0.065, 0.028, 0.013, 0.008, 0.0027`, and the degraded rescues to 9–86.
  - But 55, 16, 3, 6, 5, 1 earlier `n <= 10^6` are broken by later flips (harm).
  - So freezing is what blocks, and not freezing harms. Lemma 6 is the natural repair, but 1FD fails for 20–36% of points.

**Reading.**
- If the MODE 0 interference fraction (expensive rescues per undecided source, `0.021` at `L = 64` and falling slowly) decays only polynomially, this construction's price is `ρ_L/poly(L)`, not `ρ^peak`.
- The MODE 1 fraction falls faster (`0.0027` at `L = 64`, below `ρ^peak/ρ_L = 0.0087`).
- The data cannot decide the asymptotic question for `δ_L`.

## 8. What is proved about HYP-9140, and what is missing

**Proved.**
- The pairing *move* is not an obstruction to peak discounting.
- Theorems A–B: flips placed by any prefix-determined rule pay a factor `λ` each in the tilted measure. Coupling of the flipped and unflipped continuations (Lemma 1) is thereby sidestepped, not controlled.
- So an undecided source can be sent down privately at cost `ρ_L exp(−Θ(L^(1/3)))` on average, and at `poly(L)·ρ^peak` given Conjecture R.

**Missing.** A single member of the pairing family that realizes the private certificates of most undecided sources. The two observed interference mechanisms are:
1. **Blocking by merging.** An undecided orbit merges, before reaching its barrier point, into an earlier frozen certificate. It arrives late, since otherwise it would ride that certificate down. From then on it cannot flip.
2. **Harm.** Without freezing, a later flip on an earlier Collatz-good orbit can make that orbit come down too late (§6; 1–55 per `10^6` in MODE 1).

**Why the natural proof routes fail.**
- **Union bounds over late merges.** Counting the `n` whose orbit meets a given set `Z` of certificate points within `L` steps needs equidistribution of `Z` modulo powers of 3. For a general set, the backward trees can be exponentially large (THM-4475's `W_L`).
- **Bijection plus tilt for the actual construction.** The frozen environment makes the flip decisions depend on the exact value of `x_t`, not on its class. Lemma 2 then fails, and with it the `λ`-per-flip accounting.
- **Horizon slack.** Absorbing late merges by certifying within `L + Δ` costs `2^((1−H)Δ)`. Hence `Δ = O(L^(1/3))` at most. Late merges can be late by `Θ(L)`.
- **Domination (Lemma 6).** 1FD fails for a constant fraction of points. Its induction also needs every flip, fallbacks included, to be 1FD.

**Precise open statements.**
- **(O1)** HYP-9140: `δ_L <= poly(L) ρ^peak_L`.
- **(O2)** The weak form: `δ_L <= ρ_L e^(−cL^(1/3))` for some `c > 0`. A proof would follow from any consistent construction in `P_L` in which the undecided sources that pay more than `3^(−W)`, with `W ≍ L^(1/3)`, have density `<= ρ_L e^(−cL^(1/3))`.
- **(O3)** Conjecture R. It gives HYP-9140 for `π_L` and the sharp constant `κ_3` for `π_L`.

## 9. Failures and caveats

- **The freshness heuristic.** The HYP-9140 text ("fresh continuation, falls with probability tending to 1") is wrong for undecided sources: Lemma 1, and the rescue rate `0.33 → 0.20` of the peak candidate. The peak lane's greedy works because it tries many candidates. The provable replacement is the reflected barrier.
- **First attempt at a proof.** It tried to show that for most undecided sources *some* single high flip succeeds, via a probabilistic analysis of the transducer's `(a_t, b_t)` under `Q`. This needs the parity of `b_t` to be equidistributed, which I could not prove. The barrier plus tilt argument avoids it.
- **Frozen versus unfrozen.**
  - MODE 0 (freeze everything) is valid, but interference (blocking) grows to about 10% of rescues.
  - MODE 1 (freeze only rescued certificates) has little blocking, but breaks 1–55 earlier `n` per `10^6`.
  - MODE 3 (MODE 1 plus 1FD at every barrier flip; scratch only, not in the deliverable) was still invalid at `L = 16, 24, 32` (66, 5, 1 broken `n`), because the THM-4475 fallback flips are not 1FD. It also had 3–10× more fallbacks.
- **Constants.** The elementary exponent `γ = 0.357` is far from `κ_3 = 2.108`. The multi-scale bound carries a factor `L` (flips per source `<= L`), so it beats `ρ_L` only from `L ≈ 256`.
- **Theorem C is formal at small L.** With `K_0 = 1.0033` it gives `π_L <= (81.3L + 54.2)ρ^peak_L` for `L <= 180`, which is weaker than the trivial `2ρ_L` there. Its value is the conditional asymptotic statement.
- **Statistics at large L.** The consistent-construction densities at `L = 56, 64` rest on 136 and 69 flipped pairs.
- **The private cost column in the output** is the private cost of the construction's own rule, summed over the consistent construction's rescued set. It is not `π_L`.
- **What is and is not proved.** The float DP checks (`L = 256…1024`) are not integer-exact, but all DP sums are sums of positive terms (relative error `~1e−13`). Theorem A itself is a proof and does not rest on them.

## 10. Reproduction

```
python3 04-computation/experiments/procgen_pairpeak_20260926_run.py > 05-knowledge/results/procgen_pairpeak_20260926.out
```

- **Environment.** Python 3.10.0, numpy 2.2.6, Apple clang 17.0.0. The runner compiles the C program into a temporary directory.
- **Dependency.** It imports `procgen_peak_20260926_lib.py` read-only: exact floors, `ρ_L`, `ρ^peak_L`, the Mogul'skii constant. Its sha256 is unchanged from the peak note.
- **Final run.**
  - Wall time 28.2 s on the shared 8 GB machine.
  - Peak RSS 306 MiB by `/usr/bin/time -l`, which includes the C program at `N = 10^7`; the runner itself reports 216 MB.
  - 222 output lines, 100 checks, ending in `ALL CHECKS PASSED`.
- **Determinism.** A second run gave identical output apart from the `[time]` lines.

| file | sha256 (raw bytes) |
|---|---|
| `04-computation/experiments/procgen_pairpeak_20260926_lib.py` | `aa915a333d913f2ada76a90ee1dc6602a17179f38b24b7f23267053f39b4ec56` |
| `04-computation/experiments/procgen_pairpeak_20260926_greedy.c` | `e142a6bd23fa35f3a3b6b64a070cb87e5a68bb4478c510190439eb2a344b18e6` |
| `04-computation/experiments/procgen_pairpeak_20260926_run.py` | `e4c52557d7ffb12b5ae8890e5511721c98f59dbf5780813e1b3e407297b1a108` |
| `05-knowledge/results/procgen_pairpeak_20260926.out` | `e742b603b1ba928ab86802108bbe89ea2ab227da1f7c6e7d8e270160f6b0f54d` |
| same output without `[time]` lines (`grep -v '^\[time\]'`) | `9fbc95e525a31e11b52c9da930e5b596c5fd0f50ad0ec54589913611b53de1cc` |
| `04-computation/experiments/procgen_peak_20260926_lib.py` (imported, unchanged) | `ecc5f8a27c13ca40fad95ae8f6486a88a982d4cfa386d5b3808ac307ec778070` |

**Single runs.** Build with `cc -O2 -o greedy procgen_pairpeak_20260926_greedy.c -lm`, then run `./greedy L N WMAX MODE`:
- MODE 0: the valid construction;
- MODE 1: the no-freeze diagnostic;
- MODE 2: THM-4475's `G_L`.

**Runner sections.**
- 0 constants;
- 1 coupling;
- 2 bijection and realization;
- 3 Theorem A;
- 4 the private price;
- 5 Robin;
- 6 the consistent construction;
- 7 1FD.

Temporary files lived only under `scratch/procgen_pairpeak/` (not part of the deliverable).

## 11. Hypothesis bookkeeping (no files created)

- **HYP-9140.**
  - Remains OPEN.
  - Its "Why it is plausible" paragraph should be corrected per §1: one flip re-merges with probability 1/2, and the working mechanism is a reflected barrier with `λ` per flip.
  - It now splits into Conjecture R (the private, pairing-move part; §5) and consistency/interference (§§6–8).
- **New candidates**, for the coordinator to number:
  - **(R)** Conjecture R, the Robin inequality `sup N_m(L)/A_(m+1)(L) < ∞`. FINITE-EXACT `<= 1.003255` for `L <= 180`. It implies HYP-9140 for the private price.
  - **(I)** Interference: there is a consistent member of `P_L` whose expensive rescues have density `<= ρ^peak_L poly(L)`. Equivalently for the purpose of HYP-9140: blocking by late merges can be avoided at polynomial cost.
- **Promotion candidate** after audit: Theorems A and B (the private price `π_L = ρ_L exp(−Θ(L^(1/3)))`), with Lemmas 1, 2, 5 and 6.
