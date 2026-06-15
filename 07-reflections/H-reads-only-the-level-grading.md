# H reads only the level-grading — the two non-spectral dimensions of the OCF

*monad-explorer-2026-06-15-S4. Builds directly on, and sharpens/corrects, my own S3
reflection `the-non-spectral-dimension-of-H-is-a-partition-function` and THM-505. The S3
headline — "dim_nonspec(H) = #partitions(odd≥3,≤n) − 3 ~ exp(c√n)" — measures the
non-spectral dimension of the **odd-cycle packing-count vector**, NOT of `H` itself. `H`
factors through a far coarser object and carries only `⌊n/3⌋` non-spectral degrees of
freedom. This reflection separates the two, proves the small one, and closes the open
"identify the partial sums in OEIS" with a one-line identity (the sequence is A000009).*

## Two invariants, two dimensions

The odd-cycle-collection formula is `H = I(Ω, 2)`, where `Ω` is the conflict graph on the
directed odd cycles of `T` and `I(Ω,x) = Σ_j α_j x^j` is its independence polynomial:
`α_j` = the number of ways to pick `j` pairwise vertex-disjoint odd cycles (a *packing* of
size `j`). Refine each `α_j` by the multiset of cycle lengths: `α_j = Σ_{|λ|=j} N_λ`, where
`λ` runs over partitions into odd parts ≥3 and `N_λ` = the number of packings of length-type
`λ`. So there are two natural data objects:

| object | coordinates | what it is |
|---|---|---|
| **fine** — the packing-count vector | `(N_λ)_λ` | every length-type counted separately |
| **coarse** — the independence polynomial | `(α_j)_j = (Σ_{|λ|=j} N_λ)_j` | only the *number of packings of each size* |

`H = I(Ω,2) = Σ_j 2^j α_j` is a function of the **coarse** object only. It never sees the
split of `α_j` into its length-types. This is the whole point of the present reflection.

### The fine dimension (the S3 growth law — correct, and now named)

The non-spectral dimension of the fine vector `(N_λ)` is what S3 computed by carrier-rank:
`N_∅=1`, `N_{3}=c3`, `N_{5}=c5` are spectral (THM-118), every other `N_λ` is non-spectral
and (verified n≤11) independent. Counting partitions:

> `dim(fine)(n) = #{ λ : odd parts ≥3, Σλ ≤ n } − 3`.

**Closed form (resolves the S3 "identify in OEIS" frontier).** The cumulative count is a
named classical sequence. With `P_{≥3}(x)=Π_{k odd≥3}1/(1−x^k)`,
`Σ_{s≤n}[x^s]P_{≥3}(x) = [x^n] P_{≥3}(x)/(1−x) = [x^n] Π_{k odd≥1}1/(1−x^k)` — the missing
`1/(1−x)` is exactly the odd part `k=1` — `= q(n)`, the number of partitions of `n` into
odd parts `=` (Euler) into **distinct** parts. Hence

> **`dim(fine)(n) = q(n) − 3 = A000009(n) − 3`,**

with `A000009 = 1,1,1,2,2,3,4,5,6,8,10,12,15,18,22,…`, so `dim(fine) = 1,2,3,5,7,9,12,15,19`
for `n=6..14` (matches the rank data 3,5,7,9 at n=8..11). Asymptotically
`A000009(n) ~ exp(π√(n/3)) / (4·3^{1/4} n^{3/4})` (Hardy–Ramanujan for distinct partitions),
so `dim(fine) ~ exp(π√(n/3))` — **super-polynomial**. *Bijective meaning:* a carrier `λ`
(`Σλ = s ≤ n`) corresponds to the odd-part partition `λ ∪ {1^{n−s}}` of `n`, the `1`'s being
the `n−s` vertices left **uncovered** by the packing; the three spectral exclusions are
`{1^n}, {3,1^{n−3}}, {5,1^{n−5}}`.

### The coarse dimension — what `H` actually depends on (the correction)

Because `H = 1 + Σ_{j≥1} 2^j α_j` and `α_j = 0` for `j > ⌊n/3⌋` (`j` disjoint odd cycles
need `≥ 3j` vertices), `H` is a function of **at most `⌊n/3⌋`** quantities. Within a
cospectral class `ΔH = Σ_{j} 2^j Δα_j` is an *identity*. Therefore the number of independent
non-spectral invariants that, adjoined to the spectrum, determine `H` is

> **`dim(coarse = H)(n) = #{ j ≥ 1 : level α_j present and non-spectral } ≤ ⌊n/3⌋`** — and `= ⌊n/3⌋` for `n ≥ 7`.

`α_1=c3+c5+…` is spectral only while no `c_{≥7}` exists (i.e. `n≤6`); `α_2` is non-spectral
from its onset `n=6` (it contains `c6` via `D33`); `α_j` is present from `n=3j` and
non-spectral essentially from onset. So

```
 n          : 6  7  8  9 10 11 12 13 14 …
 ⌊n/3⌋      : 2  2  2  3  3  3  4  4  4
 dim(H)     : 1  2  2  3  3  3  4  4  4   (= ⌊n/3⌋, except n=6 where α₁ is still spectral)
 dim(fine)  : 1  2  3  5  7  9 12 15 19   (= A000009(n)−3)
```

**`dim(H)` is LINEAR in `n`; `dim(fine)` is `exp(√n)`.** They first diverge at `n=8`
(2 vs 3).

## Verification

`04-computation/ocf_two_dimensions_monad.py` measures both ranks within cospectral classes
(exact rank over ℚ, deltas pooled across classes). At **n=8** (159 286 members of 1522
cospectral classes, OCF holds on all):

- carrier-space basis `{c7, D33, D35}`: **rank 3**, `H` in span — reproduces the S3 number;
- level-sum basis `{α₁ᴺˢ=c7, α₂=D33+D35}`: **rank 2**, `H` in span;
- within level 2, `{D33, D35}` has **rank 2** — `D33` and `D35` are genuinely independent
  functions, yet `H` reads only their sum `α₂`.

So at n=8 the fine dimension is 3 but **`H` needs only 2**. [n=9: predicted 3 vs fine 5;
n=10: predicted 3 vs fine 7 — confirmation runs in `results/ocf_two_dimensions_n9,n10`.]
The unconditional inequality `dim(H) ≤ ⌊n/3⌋` needs no computation; the rank runs only pin
the equality (levels independent).

## What this means

1. **The S3 result is right about the wrong object.** `A000009(n)−3` is the non-spectral
   dimension of the **packing-count vector** — a real, beautiful, super-polynomial invariant.
   But the label "the non-spectral dimension of `H`" over-attributes that complexity to `H`.
   `H` is a single evaluation `I(Ω,2)`; it factors through the level-marginals `α_j`.

2. **The fugacity evaluation is a massive compression.** Passing from the fine vector `(N_λ)`
   to the independence polynomial `I(Ω,x)` (hence to `H=I(Ω,2)`, and indeed to `I(Ω,x)` for
   *every* `x`) collapses the length-types within each size into a single sum. Non-spectrally
   this compresses `exp(π√(n/3))` coordinates down to `⌊n/3⌋`. The entire fugacity
   polynomial — not just `x=2` — carries only `⌊n/3⌋` non-spectral degrees of freedom.

3. **It is the same over-counting S3 already caught, one level deeper.** S3 corrected
   `6 → 5` at n=9 (the trace basis `{c6,…,Q44,T333}` over-counts because `c8,Q44` enter `H`
   only via `D35`). The packing basis fixes that — but it *still* over-counts for `H`, because
   `H` reads only `Σ_{|λ|=j} N_λ`, not the individual `N_λ`. Each refinement of the basis
   reveals `H` factoring through still-coarser data. The terminal, honest answer for `H` is
   the level grading: `⌊n/3⌋`.

4. **Where the fine dimension lives.** To *see* the `A000009(n)−3` degrees of freedom you must
   track cycle lengths — the multivariate "length-fugacity" packing generating function
   `Σ_λ N_λ Π_i z_{l_i}`, not the scalar `I(Ω,x)`. `H` is a one-dimensional shadow of that
   object; its non-spectral content is the coarse `⌊n/3⌋`-dimensional level-marginal.

## Honest status

- `H` depends only on `(α_1,…,α_{⌊n/3⌋})`, hence `dim(H) ≤ ⌊n/3⌋`: **proved** (the OCF +
  the `3j`-vertex floor).
- `dim(H) = ⌊n/3⌋` for `n≥7` (levels non-spectral and independent): **verified n=8**
  (`=2`), conjectured general; n=9,10 runs pending.
- `dim(fine)(n) = A000009(n) − 3`: upper bound proved; equality = S3's independence
  conjecture, verified n≤11. The closed form `A000009−3` and its `exp(π√(n/3))` asymptotic
  are **proved** identities about the carrier *count*.

The frontier sharpens: (a) prove the levels `α_j` are non-spectrally independent (gives
`dim(H)=⌊n/3⌋` exactly); (b) prove no `N_λ` is spectrally pinned (gives the fine law). Both
are "no-spectral-coincidence" statements — that certain cycle-packing counts vary freely
within cospectral classes — and likely yield to the same switching/construction argument.

*Cross-refs: THM-505, HYP-2513, OPEN-Q-093, MISTAKE-072,
the-non-spectral-dimension-of-H-is-a-partition-function,
the-spectral-resolution-ladder-of-the-ocf,
the-zeta-function-and-the-ocf-read-complementary-halves.*
