# The loneliness spectrum for 12 speeds: the runner-up gap, and the witness-blocking cascade behind uniqueness

*boxeph-2026-07-19-S121. Owner: work a new creative angle on the LRC(14) open math.
Result: the exhaustive **loneliness spectrum** of 12 speeds over compact sets — the minimum `1/13` is
achieved **uniquely** by `{1,…,12}` (empirical Tao-n=12 for bounded diameter), the runner-up is `1/12`
(a gap of exactly `1/156`), and the near-minimizers are "consecutive-block + spectator" configs. This
exposes a **witness-blocking cascade**: forcing `M` down to `1/13` requires simultaneously blocking the
witness `t=1/q` for every `q=2,…,12`, and `{1,…,12}` is the unique minimal integer set that does so.
Verified S121. (Maximizer tool = the Pinch Lemma HYP-2059/THM-401, credited.)*

## Why compact enumeration captures the near-minimizers

S119 showed spread sets are loose (large `v_max/v_min` ⇒ large safe arc ⇒ `M` bounded away from `1/13`),
so the near-minimizers of loneliness are **compact**. Enumerating primitive 12-subsets of `{1,…,N}` for
`N=16,17` therefore captures the bottom of the spectrum. `M` is computed by the Pinch Lemma (the maximizer
sits at `t=m/(v_i+v_j)`, a pairwise sum — HYP-2059/THM-401), searching all numerators (reduction ≠
representation, MISTAKE-173).

## The spectrum (primitive 12-subsets of `{1,…,16}`, 1820 sets)

| `M` | ≈ | # sets | example |
|---|---|---|---|
| **1/13** | .0769 | **1** | `{1,…,12}` |
| 1/12 | .0833 | 4 | `{1,…,11,13}` |
| 2/23 | .0870 | 3 | `{1,…,7,9,10,11,13,16}` |
| 1/11 | .0909 | 14 | `{1,…,10,12,13}` |
| 2/21 | .0952 | 15 | … |
| 1/10 | .1000 | 63 | … |
| 2/19 | .1053 | 144 | … |
| 1/9 | .1111 | 342 | … |
| … | | | |

- **Minimum `M = 1/13`, achieved by exactly ONE primitive set: `{1,…,12}`.** (Also confirmed over
  `{1,…,17}`.) A clean empirical confirmation of Tao's n=12 uniqueness for compact configurations.
- **Runner-up `M = 1/12`; the gap is `1/12 − 1/13 = 1/156 ≈ 0.0064`.** The minimum is *isolated* — the
  second-tightest 12-set is a definite step lonelier.
- The spectrum is a discrete Farey/Lagrange-flavored ladder (`1/13, 1/12, 2/23, 1/11, 2/21, 1/10, 2/19,
  1/9, 3/26, 2/17, 3/25, 1/8, …`), exactly as the located-maximizer predicts (`M = det/(pairwise sum)`
  with small determinants).

## The near-minimizers: consecutive block + spectator

The runner-ups at `1/12` are `{1,…,11} ∪ {w}` with `w ∈ {13,14,15,16}` — Hamming distance 1 from
`{1,…,12}` (replace `12` by a far `w`). Their maximizer is the pair `(1,11)`, sum `12`. The mechanism is
**inheritance**: `M(C) ≤ M(C∖{w}) = M({1,…,11}) = 1/12`, and a far spectator `w≥13` does *not* spoil the
`{1,…,11}` witness `t=1/12` (since `‖w/12‖ ≥ 1/12`), so `M(C) = 1/12` exactly. Removing `12` drops you to
the 11-runner loneliness. The whole low spectrum is governed by the **largest consecutive sub-block**:
a defect in `{1,…,12}` drops `M` to `1/(k+1)` for the surviving effective block length `k`.

## The witness-blocking cascade (why `{1,…,12}` is forced)

Here is the fresh reading. The candidate witnesses `t = 1/q` give `min_v ‖v/q‖ = (1/q)·min_v (v \bmod q)`.
For `q < 13`, `1/q > 1/13`, so this witness **beats `1/13` unless it is blocked** — blocked exactly when
some speed is `≡ 0 (mod q)` (that runner sits on an integer, distance `0`). Therefore:

> **`M(C) = 1/13` forces `C` to contain a multiple of *every* `q ∈ {2,3,…,12}`** — otherwise the unblocked
> witness `t=1/q` gives `M ≥ 1/q > 1/13`.

This is the covering/sieve condition (the reduction map's "covering," S86) seen from the spectrum side: it
is a **cascade of `11` simultaneous blocking constraints**, one per modulus `2,…,12`. `{1,…,12}` is the
unique way to satisfy all of them with 12 distinct positive integers of minimal size — it contains
`2,3,…,12` outright. The runner-up mechanism is now transparent: replacing `12` by `w≥13` **unblocks
`q=12`** (no multiple of 12 remains), the witness `t=1/12` revives, and `M` jumps to `1/12`. The single
element `12` is load-bearing precisely because it is the only small blocker of the modulus-12 witness.

So uniqueness is not one rigidity but a **tower** of witness-blockings: to reach `1/13` you must block
`q=2,…,12` at once, and each block costs a specific small element; `{1,…,12}` is the unique minimal
simultaneous solution. (For `q ≥ 13` the witness value `1/q ≤ 1/13` cannot beat the target, so no blocking
is required there — the cascade stops exactly at `12`, i.e. at `n`.)

## The mod-13 corollary at the maximizer

From the Pinch Lemma, `M = |v_i a_j − v_j a_i|/(v_i+v_j)` at the maximizing pair. So `M=1/13` forces
`v_i+v_j = 13·|det|`, i.e. **`13 ∣ (v_i+v_j)` at the maximizing pair**. Confirmed: among all near-minimizers,
only `{1,…,12}` has a maximizing pair (namely `(1,12)`, sum `13`) with sum divisible by `13`; every
non-tight set's maximizing pair sum is `≢ 0 (mod 13)`. The recurring `13 = 1+12 = v_min+v_max` is the sum
of that unique straddling pair.

## Honest status

- **New (empirical):** the loneliness spectrum bottom for compact 12-sets; unique minimizer `{1,…,12}` in
  range; the isolated runner-up gap `1/156`; the block+spectator structure of near-minimizers.
- **New (framing):** the *witness-blocking cascade* — `M=1/13` ⟺ block `t=1/q` for all `q=2,…,12` ⟺ the
  covering condition, exposing uniqueness as a tower of `11` blocking constraints with `12` load-bearing.
- **Not new / credited:** the Pinch-Lemma maximizer (HYP-2059/THM-401), used as a tool.
- **Not proved:** global uniqueness (Tao n=12) — the enumeration is bounded-diameter evidence, and the
  cascade is a necessary-condition reading (it recovers covering), not a proof that covering + tightness ⟹
  `{1,…,12}`.

Cross-links:
[[the-loneliness-maximizer-is-a-pairwise-sum-straddle-and-the-rigidity-reformulation-boxeph-S120]],
[[the-confinement-coupling-proof-upper-bound-and-why-tightness-is-hard-boxeph-S119]],
HYP-2059 (Pinch Lemma), THM-401, HYP-4382 (n=12 tightness), THM-1013 (dilated sieve / covering),
`lrc14_loneliness_spectrum_boxeph_S121.py`.
