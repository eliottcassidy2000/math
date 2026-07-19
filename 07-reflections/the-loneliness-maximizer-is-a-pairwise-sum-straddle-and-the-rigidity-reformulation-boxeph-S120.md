# The loneliness maximizer is a pairwise-sum straddle — locating the maximizer, and a sharper form of Tao's n=12 uniqueness

*boxeph-2026-07-18-S120. Owner: work a new creative angle on the LRC(14) open math.
Result: a **located-maximizer theorem** — `M(C)` is always attained at `t = m/(v_i+v_j)` for a
straddling active pair, with `M = |v_i a_j − v_j a_i|/(v_i+v_j)` — which turns the uniqueness of the
loneliness minimizer (Tao's n=12 conjecture / INVcov) into a statement about **pairwise sums**: `{1,…,12}`
is the unique 12-set whose best pairwise-sum straddle value equals `1/13`. This generalizes the S118
centering witness and pins the difficulty to a finite, structured per-set search. Verified S120.*

## The located-maximizer theorem (proved, variational)

> **Theorem.** For distinct positive speeds `V`, let `t*` maximize `g(t) = min_k ‖v_k t‖`, `M = g(t*) > 0`.
> Then there is a pair `i, j` and integers `a_i, a_j` with
> `t* = (a_i + a_j)/(v_i + v_j)`  and  `M = |v_i a_j − v_j a_i|/(v_i + v_j)`.
> The two runners `v_i, v_j` are **active** (`‖v_· t*‖ = M`) and **straddle** their integers — one just
> above (`v_j t* = a_j + M`), one just below (`v_i t* = a_i − M`). The maximizer denominator is a pairwise
> **sum**; difference denominators never attain the max.

*Proof.* `g` is continuous on `[0,1]`, `g(0)=0`, so the max is interior with `M>0`. A global max is a local
max of the min, so among the active runners some must block each direction: moving `t` right decreases
`‖v_i t‖` for an active `i` with `v_i t* = a_i − M` (`−`slope, `σ_i = -1`), and moving left decreases
`‖v_j t‖` for an active `j` with `v_j t* = a_j + M` (`+`slope). Adding, the `±M` cancel:
`(v_i+v_j) t* = a_i + a_j`. Solving and substituting into `v_j t* = a_j + M` gives the formula. ∎

*(Verified: over 60 random 12-subsets of `[1,44]`, the location and the formula hold with **zero**
failures, and a difference denominator beats every sum denominator in **zero** cases. Example:
`C=[3,4,5,7,10,14,21,24,26,33,35,41]`, `M=3/20` at `t*=9/20 = 18/40`, straddle pair `(33,7)`, sum `40`,
`|33·3 − 7·15|/40 = 6/40 = 3/20`.)*

## Why this is the right generalization of the centering witness

For an AP `{a+dk}`, the straddling pair is exactly the **extremes** `(v_min, v_max) = (a, a+11d)`, their sum
is `q = 2a+11d`, and the formula returns `M = (q-11)/(2q)` — the S118 witness. So S118 was the special case
"straddling pair = the two ends, everything else centered symmetrically." The theorem says **every** speed
set has a straddling pair whose sum is the maximizer modulus, with the remaining runners packed into the
safe band `[M, 1−M]`. The centering witness is universal; only *which* pair straddles varies.

## The rigidity reformulation (Tao n=12 / INVcov, sharper form)

`M(C) = 1/13` iff the **best pairwise-sum straddle value** equals `1/13`. For `{1,…,12}` the only pair
reaching `1/13` is `(1,12)`, whose sum is `13`; the other runners land at `{2/13,…,11/13}`, all `≥ 1/13`.
So:

> **Reformulation (equivalent to Tao's n=12 uniqueness).** `{1,…,12}` is the **unique** 12-set whose best
> pairwise-sum straddle value is exactly `1/13`; every other 12-set has some pairwise sum `v_i+v_j` with a
> straddle value `> 1/13`.

What is new here is **form, not a proof**: the maximizer is now **located** at a pairwise sum, so the
uniqueness question is no longer existential over all `t ∈ [0,1]` — it is a bounded search over the `≤ 66`
pairwise sums, asking whether any produces a straddle value above `1/13`. This is the sharpest statement of
"find a good `t`" the project has: the good `t`, if it exists, has denominator `v_i+v_j`.

Computationally it is airtight on the near-minimizers: `{1,…,12}` has best pairwise-sum straddle exactly
`1/13`, and **all 204 single-element perturbations** (S120 explorer) have one strictly exceeding `1/13`
(equal to their true `M`), as do all reflective non-AP sets tested.

## Why `13 = 1 + 12` — and where the difficulty now sits

The reformulation explains the recurring `13` structurally: it is `v_min + v_max` of `{1,…,12}`, the sum of
the unique straddling pair at the minimizer. The two active ends are `1` (just above `0`) and `12` (just
below `1`); they sum to `13`, and `1/13` is `|1·1 − 12·0|/13`. The `10` interior runners are **slack**
(strictly inside the band) — the same "2 active, 10 slack" structure noted in S95, now with the 2 active
ones identified as the straddling pair and the denominator as their sum.

The residual difficulty is unchanged in substance but sharper in shape: prove that no non-`{1,…,12}` set can
hold its interior runners in `[1/13, 12/13]` while some pairwise-sum straddle sits at exactly `1/13` — i.e.
that the perfect `{k/13}` packing is achievable only by `{1,…,12}`. This is the **centering feasibility**
question over the finite pairwise-sum moduli. It remains equivalent to Tao's conjecture; what this session
adds is that the maximizer is a pairwise sum, so the search space is finite and structured, and the AP face
(S117–S119) is the sub-case where the straddling pair is the two extremes.

## Honest status

- **Proved:** the located-maximizer theorem (variational; likely folklore in LRC theory, but self-contained
  here and stated with the pairwise-**sum** emphasis) and the exact straddle formula.
- **Confirmed:** the maximizer is always a pairwise sum (60 random + 204 perturbations + reflective/random),
  differences never win, `{1,…,12}` uniquely witness-less at `1/13`.
- **Not proved:** the uniqueness itself (Tao n=12) — the reformulation locates the maximizer but does not
  close the centering-feasibility question.

Cross-links:
[[the-confinement-coupling-proof-upper-bound-and-why-tightness-is-hard-boxeph-S119]],
[[the-centering-witness-closes-the-spread-case-exact-loneliness-of-every-AP-boxeph-S118]],
HYP-4382 (n=12 tightness), HYP-7401 (the crux is offset-vanishing / one-line form),
`lrc14_maximizer_pairwise_sum_boxeph_S120.py`, `lrc14_centering_general_sets_boxeph_S120.py`.
