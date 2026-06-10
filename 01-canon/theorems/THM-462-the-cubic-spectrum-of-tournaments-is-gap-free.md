# THM-462 — The cubic spectrum of tournaments is gap-free: every c3 ∈ [0, max] is realized

**Status:** PROVED (self-contained elementary proof below: classical identity + Landau + a
four-square near-regular perturbation covering the top window + induction) + VERIFIED
(exhaustively, exact integer DP over Landau sequences, **all n ≤ 60**: zero gaps; the proof's
constructive step machine-checked for **every (n,t), n ≤ 80** — 21,379 exact constructions).
The identity and the max formula are **classical** (Kendall & Babington Smith 1940); the
gap-free statement at general n is **plausibly folklore** — our searches found small-n
verifications but no explicit general-n proof in print (see Attribution). The proof here is
self-contained either way.
**Source:** kind-pasteur-2026-06-10-S1 (Thread C, HYP-2368 — which this REFUTES as posed).

## The theorem

> **THM-462.** Let `A(n)` be the set of values of `c3(T)` (number of directed 3-cycles =
> cyclic triples) over all tournaments `T` on `n` vertices. Then for every `n ≥ 1`
>
> `A(n) = {0, 1, 2, …, M(n)}`,  where  `M(n) = (n³−n)/24` (n odd), `(n³−4n)/24` (n even).
>
> **The degree-3 channel has NO forbidden levels.** Contrast: the repo's degree-`n` invariant
> `H(T)` has the forbidden values 7 and 21 (THM-029/THM-079). Impossibility begins higher up
> the invariant hierarchy than the cubic layer.

Throughout: `s_1 ≤ … ≤ s_n` is the (sorted) score sequence, `C(a,b)` binomial,
`f(s) = Σ_i C(s_i,2)`, and `c3 = C(n,3) − f(s)`.

## Lemmas

**L1 (cyclic-triple identity — classical).** For every tournament,
`c3(T) = C(n,3) − Σ_v C(s_v,2)`. In particular c3 depends only on the score sequence.
*Proof.* In any 3-subset the three within-triple out-degrees sum to 3, so they are
`(2,1,0)` (transitive) or `(1,1,1)` (cyclic). A triple is transitive iff exactly one of its
vertices beats the other two, and `Σ_v C(s_v,2)` counts the pairs (v, {2 vertices v beats}),
i.e. each transitive triple exactly once. ∎
[Kendall & Babington Smith 1940; reproduced in Moon, *Topics on Tournaments* (1968) and
David, *The Method of Paired Comparisons*. Sanity-verified below.]

**L2 (Landau 1953).** A nondecreasing integer sequence `s_1 ≤ … ≤ s_n` is the score sequence
of some tournament iff `Σ_{i≤k} s_i ≥ C(k,2)` for all k, with equality at `k = n`.
Hence by L1: `A(n) = { C(n,3) − f(s) : s Landau }`.

**L3 (upper bound — classical).** `c3(T) ≤ M(n)`, with equality for (near-)regular scores.
*Proof.* Minimize `f` over integer sequences with sum `C(n,2)`: if `s_j ≥ s_i + 2`, the trade
`s_i → s_i+1, s_j → s_j−1` changes `f` by `s_i − s_j + 1 < 0`; so the minimum has all scores
within 1 of each other: all `h = (n−1)/2` (n odd), or `h−1, h` each `n/2` times (`h = n/2`,
n even). These are Landau (prefix sums `k·(h or h−1) ≥ C(k,2)` directly). The minimum values
are `m(n) = n·C(h,2)` (odd), `m(n) = (n/2)(C(h,2)+C(h−1,2))` (even), and exact algebra gives
`C(n,3) − m(n) = M(n)` in both parities. ∎

**L4 (deviation identity).** (odd `n = 2h+1`) If `s_i = h + e_i` with `Σ e_i = 0`, then
`f(s) = m(n) + (Σ e_i²)/2`  [since `C(h+e,2) = C(h,2) + (e(2h−1) + e²)/2`].
(even `n = 2h`) If the `h` lower-half scores stay `h−1` and the upper half is `h + d_i` with
`Σ d_i = 0`, then `f(s) = m(n) + (Σ d_i²)/2`.

**L5 (window length).** `W(n) := M(n) − M(n−1) = h(h+1)/2` for `n = 2h+1`, `= h(h−1)/2`
for `n = 2h`. (Direct algebra; asserted exactly for n = 8..80 in the script.)

**L6 (dominant-vertex embedding).** `A(n−1) ⊆ A(n)`: add a vertex beating all others —
every triple containing it is transitive, so c3 is unchanged.

**L7 (window construction).** For `h ≥ 8` (i.e. `n ≥ 16`) and every integer `t ∈ [0, W(n)]`,
there is a Landau sequence with `f = m(n) + t`, hence a tournament with `c3 = M(n) − t`.
So `[M(n−1), M(n)] ⊆ A(n)`.

*Proof.* By Lagrange, write `t = a₁² + a₂² + a₃² + a₄²` with `0 ≤ a_j ≤ √t`. Note
`a_max ≤ √t ≤ √(h(h+1)/2) ≤ h−2` for `h ≥ 8` (⟺ `(h−1)(h−8) ≥ 0`; even case
`√(h(h−1)/2) ≤ h−2` already for `h ≥ 6`), and `A := Σ a_j ≤ 2√t` (Cauchy–Schwarz, 4 terms).

*Odd `n = 2h+1`:* take 4 scores `h − a_j`, 4 scores `h + a_j`, the other `n−8` equal `h`.
By L4, `f = m(n) + t`. Landau check (`P_k` = sum of k smallest):
- `k ≤ 3`: every score `≥ h − a_max ≥ 2`, so `P_k ≥ 2k ≥ C(k,2)`.
- `4 ≤ k ≤ n−4`: `P_k ≥ kh − A` (total negative deviation is A), and
  `kh − A ≥ C(k,2)` ⟺ `A ≤ k(2h+1−k)/2`, whose min over `[4, 2h−3]` is `4h−6` (both
  endpoints); `A ≤ 2√t ≤ √(2h(h+1)) ≤ 4h−6` ⟺ `14h² − 50h + 36 ≥ 0` ⟺ `h ≥ 3`. ✓
- `k = n−j, 1 ≤ j ≤ 3`: the top j scores sum `≤ j(h + a_max) ≤ j(2h−2) ≤
  j(2h+1) − j(j+1)/2 = C(n,2) − C(n−j,2)` (⟺ `j ≤ 5`), i.e. `P_{n−j} ≥ C(n−j,2)`. ✓
- `k = n`: equality. ✓

*Even `n = 2h`, `h ≥ 8`:* keep the lower half at `h−1`; inside the upper half (h ≥ 8 slots)
take 4 scores `h − a_j`, 4 scores `h + a_j`, the rest `h`. By L4, `f = m(n) + t`. Landau:
- `k ≤ 3`: every score `≥ min(h−1, h−a_max) ≥ 2`, so `P_k ≥ 2k ≥ C(k,2)`.
- `4 ≤ k ≤ n−4`: for any k-set S, `Σ_S s ≥ Σ_S base − A`, so `P_k ≥ B_k − A` where `B_k` is
  the near-regular prefix sum; the base slack `B_k − C(k,2)` is piecewise concave in k with
  minimum `4h−10` on `[4, 2h−4]` (attained at both ends); `A ≤ 2√t ≤ √(2h(h−1)) ≤ 4h−10`
  ⟺ `14h² − 78h + 100 ≥ 0` ⟺ `h ≥ 4`. ✓
- `k = n−j, 1 ≤ j ≤ 3`: top j sum `≤ j(h + a_max) ≤ j(2h−2) ≤ 2jh − j(j+1)/2`
  (⟺ `j ≤ 3`). ✓   `k = n`: equality. ✓
All scores lie in `[2, n−1]` (`h + a_max ≤ 2h−2 = n−2`). ∎

## Proof of THM-462

`A(n) ⊆ [0, M(n)]` by L3. For the reverse: induction on n. **Base** `n ≤ 16`: exact
computation (bitset DP over Landau sequences, gap-free for all `n ≤ 60`; see Verification).
**Step** `n ≥ 17`: by induction `[0, M(n−1)] = A(n−1) ⊆ A(n)` (L6), and
`[M(n−1), M(n)] ⊆ A(n)` (L5 + L7, `h ≥ 8`). Union: `[0, M(n)] ⊆ A(n)`. ∎

The realizing tournaments are fully explicit: a transitive tower (L6) of height `n − n₀` on
top of a four-square-perturbed near-regular core on `n₀` vertices (L7).

## Verification (all exact integer arithmetic, no floats)

- `04-computation/c3_spectrum_kpc1.py` (+ `.out` in 05-knowledge/results/):
  **L1 sanity** — identity vs direct triple counting on ALL labeled tournaments n = 4, 5
  (64 + 1024) and 500 random tournaments each at n = 6, 7, 8 (seed 20260610): 0 mismatches.
  **Brute spectrum** — full backtracking enumeration of Landau sequences n = 3..13
  (counts 2, 4, 9, …, 48107 = A000571(3..13) exactly): gap-free, max = formula at every n.
- `04-computation/c3_spectrum_dp_kpc1.py` (+ `.out`): exact bitset DP over score
  *values with multiplicities* (state = (#scores, partial sum) → bitset of achievable f;
  in-block Landau prefixes reduce to block boundaries because `w + jv − C(c+j,2)` is concave
  in j). Cross-checked = brute enumeration for all n ≤ 13. **Result: n = 3..60 all gap-free**,
  `|A(n)| = M(n) + 1`, max = M(n), min f = m(n), at every n. (38 s total; A000571(60) ≈ 10^31
  sequences would be unenumerable — the DP is what makes n = 60 exact.)
- `04-computation/c3_window_proof_check_kpc1.py` (+ `.out`): the L7 recipe rebuilt and
  checked **from the Landau definition** for every n in 8..80 and every t in [0, W(n)]
  (21,379 constructions): valid for **all n ≥ 16** (and odd n = 11, 13, 15), `f = m(n)+t`
  exact, `a_max ≤ h−2` always; W(n) algebra asserted for n = 8..80. The failures at
  n ∈ {8,10,12,14} and t-tails of n = 9 are exactly the `h < 8` cases excluded by L7 —
  covered by the computational base.

## Attribution / honesty / literature (searches 2026-06-10)

- **Classical:** L1 and the max formula M(n) are due to **Kendall & Babington Smith,
  "On the method of paired comparisons", Biometrika 31 (1940) 324–345** (their "circular
  triads" d); reproduced in Moon (1968) and David (1988). The brief's "Goodman" attribution
  apparently conflates this with **Goodman's 1959 formula** for monochromatic triangles in
  2-edge-colorings of K_n — the *undirected* analogue. L2 is Landau (1953). L3 is classical
  (Kendall–Babington Smith). **New here:** the explicit four-square window construction L7
  and the assembled gap-free proof + the n ≤ 60 exact spectrum.
- **Gap-freeness in the literature:** K&BS (1940) tabulated the distribution of d for small n
  (all values occur there); Alway (Biometrika 49, 1962) computed distributions further;
  McShane & Harris, *Counting tournaments with a specified number of circular triads*,
  J. Integer Sequences 27 (2024) (OEIS A357242/A357257/A357248/A357266 = exactly 2/3/4/5
  circular triads) count tournaments per level. We could NOT locate an explicit general-n
  statement or proof that every value in [0, M(n)] occurs (searched: "circular triads
  spectrum / every value / all integer values attainable", "cyclic triples tournament
  possible values", "number of 3-cycles tournament spectrum", Kendall/Moon/Goodman variants).
  Very plausibly folklore among paired-comparison statisticians; we claim the *proof record*,
  not necessarily priority.
- **OEIS:** max-c3 values 1, 2, 5, 8, 14, 20, 30, 40, 55, … = **A006918(n−2)** (confirmed by
  value-search; A006918 carries the comment "Maximal number of inconsistent triples in a
  tournament on n+2 nodes [Kac]"). The spectrum-size sequence `M(n)+1` (n ≥ 3):
  2, 3, 6, 9, 15, 21, 31, 41, 56, 71, 92, 113, 141, … has **no OEIS match** (candidate
  submission; trivially A006918 + 1, so low priority). #Landau sequences = A000571 (matched
  exactly n = 3..13).

## Punchline for the repo

The cubic invariant c3 — the first nontrivial layer of the cycle hierarchy — is
**spectrally complete**: every level [0, M(n)] is realized, by an explicit transitive tower
over a four-square-perturbed near-regular core. Forbidden values are NOT a generic feature
of tournament invariants; they first appear higher up (H ∉ {7, 21}, THM-029/079). For the
FAST-channel/hallucination-detector calibration: every 3-cycle contradiction count is a
realizable calibration point — the c3 dial is continuous (integer-dense), unlike the H dial.

**Artifacts:** `04-computation/c3_spectrum_kpc1.py`, `c3_spectrum_dp_kpc1.py`,
`c3_window_proof_check_kpc1.py` (+ `.out` files in `05-knowledge/results/`).
Resolves **HYP-2368** (REFUTED as posed: no forbidden zones below max). Relates: THM-029,
THM-079 (forbidden H), THM-005 (3-cycle bijection), Landau/Moon canon.
