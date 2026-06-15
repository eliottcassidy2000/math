# The non-spectral dimension of H is a partition function (and the n=9 "6" was a basis artifact)

> **SCOPE CORRECTION (monad-explorer-S4, see `H-reads-only-the-level-grading`).** The
> partition law below correctly counts the non-spectral dimension of the **packing-count
> vector** `(N_λ)` — and that cumulative count is the named sequence **A000009** (partitions
> into distinct=odd parts), so `dim(packing) = A000009(n)−3 ~ exp(π√(n/3))`. But it is NOT
> the dimension of `H` itself: `H = I(Ω,2) = 1+Σ_j 2^j α_j` factors through the **level-sums**
> `α_j = Σ_{|λ|=j}N_λ`, of which there are only `⌊n/3⌋`. So `dim_func(H)(n) ≤ ⌊n/3⌋` (LINEAR,
> PROVED), `< A000009(n)−3` for n≥8 (verified: dim(packing)=3 but dim(H)=2 at n=8). Read
> "the non-spectral dimension of `H`" below as "of the packing vector." The fugacity-2
> evaluation compresses `exp(√n) → n/3`.


*monad-explorer-2026-06-15-S3. Builds directly on THM-505 and my own two prior reflections
`the-overlaps-stop-being-shadows-the-correlation-tower` and
`the-zeta-function-and-the-ocf-read-complementary-halves`. It sharpens the correlation-tower
picture into a closed form, and corrects one number in it: the intrinsic non-spectral
dimension of `H` at n=9 is **5, not 6** — the 6 was an over-count caused by working in a
redundant (trace) basis.*

## The slogan, made exact

The project's organizing slogan is **"the spectrum is mean-field, the OCF is correlation."**
The eigenvalue spectrum of `A` (equivalently the trace vector `tr A^3,…,tr A^n`, equivalently
the Bowen–Lanford zeta `1/det(I−uA)`) is a *single-closed-walk* invariant. The
odd-cycle-collection formula `H = I(Ω,2) = Σ_j 2^j α_j` is a *correlation* invariant: `α_j` is
the number of ways to pick `j` pairwise vertex-disjoint odd cycles. The non-spectral content
of `H` is precisely the part of this correlation data that the mean field cannot see.

The question I have been chasing across three sessions: **how big is that part?** Define
`dim_nonspec(H)(n)` = the number of independent non-spectral invariants that must be adjoined
to the full spectrum to determine `H`. Last session I reported `0,1,2,3,6` for `n=5..9` and
called the jump `3→6` a "dimension break above `n−5`." This session I found the clean law —
and that the right number at n=9 is 5.

## The law

Expand the OCF by the **length-multiset** of each packing. A packing of disjoint odd cycles
has a multiset of lengths `λ = {l_1,…,l_j}`, every `l_i` odd and `≥3`. Let `N_λ(T)` = the
number of packings with that exact multiset. Then, grouping the sum,

```
H(T) = Σ_λ  2^{|λ|} · N_λ(T) ,     over all partitions λ into odd parts ≥3 with Σλ ≤ n.
```

The map `λ ↦ N_λ` is the natural basis. Three of these counts are **spectral or trivial**:

- `N_∅ = 1`               (the empty packing — the constant term),
- `N_{3} = c3 = tr A³/3`   (a single triangle — spectral, THM-118),
- `N_{5} = c5 = tr A⁵/5`   (a single pentagon — spectral, THM-118).

**Every other `N_λ` is non-spectral, and — verified exhaustively-by-sampling at n=8,9,10 —
they are all linearly independent within cospectral classes.** Hence

> **dim_nonspec(H)(n) = #{ partitions of s into odd parts ≥ 3, with 0 ≤ s ≤ n } − 3.**

Generating function: with `P(x) = Π_{k odd, k≥3} 1/(1−x^k)`,
`dim_nonspec(H)(n) = ( Σ_{s≤n} [x^s] P(x) ) − 3`. The values:

```
 n :    6  7  8  9  10  11  12  13  14
dim :    1  2  3  5   7   9  12  15  19
              ^  ^   ^   ^   <- VERIFIED by rank (3,5,7,9 at n=8,9,10,11)
```

The increment `dim(n) − dim(n−1)` is exactly `p_{odd≥3}(n)` = the number of partitions of `n`
itself into odd parts ≥3 (for `n ≥ 6`). New carriers switch on at `n` = the sum of a new
length-multiset:

```
 {3,3}@6  {7}@7  {3,5}@8   {9},{3,3,3}@9   {3,7},{5,5}@10   {11},{3,3,5}@11   …
```

## Why `n−5` looked right, and where it really breaks

For `n ≤ 8` there is *at most one* odd-≥3 partition of each integer, so the cumulative count
grows by exactly 1 per step and `dim = n−5` coincidentally. The first integer with **two**
such partitions is `9 = {9} = {3,3,3}`, and that is exactly where the growth accelerates:
`dim(9)−dim(8) = p_{odd≥3}(9) = 2`. The "break above `n−5`" is not a pathology — it is the
moment the restricted partition function first exceeds 1. Asymptotically `dim_nonspec(H)`
grows like `exp(c√n)` (Hardy–Ramanujan for a restricted partition function), not linearly:
the non-spectral complexity of `H` is sub-exponential but super-polynomial.

## The correction: the n=9 dimension is 5, not 6

Last session's chain worked in the **trace basis** `{c6,c7,c8,c9,Q44,T333}` and counted 6.
But `c8` and `Q44` enter `H` **only through their sum**: the closed form has `+4c8 + 4Q44`,
and both terms come from `4·D35` via `D35 = c3c5 − W8 + c8 + Q44` (W8 spectral). So `H` sees
only the single quantity `c8 + Q44 = D35` (mod spectrum). Adding `c8` then `Q44` as two
separate steps in the chain double-counts one degree of freedom.

This session's rank computation makes it basis-independent: the rank of the within-class
delta-matrix of the trace carriers `{c6,c7,c8,c9,Q44,T333}` is 6, but the rank of
`{c6,c7,(c8+Q44),c9,T333}` is 5 and **already contains `H`** (rank does not increase when the
`H` column is appended). So the intrinsic dimension is 5, realized by the OCF carriers
`{c7, c9, D33, D35, T333}`, each independently needed (drop-one drops the rank). The trace
basis is over-complete; the OCF packing basis is the honest one. (This does not contradict
the n=9 *computation* — the witnesses and `ΔH = 4ΔQ44 + 8ΔT333` are all correct — it
re-coordinatizes it. The headline number changes 6→5.)

## What this says, beyond the count

1. **The carriers are indexed by partitions, and the weight is `2^{#parts}`.** This is the
   independence-polynomial fugacity `x = 2` acting level-by-level (the `2^level` rule of
   THM-505), now seen to run over *all* odd-≥3 partitions, not just the simple cycles. `H` is
   the fugacity-2 partition function of the odd-cycle gas; its non-spectral coordinates are
   the gas's configuration counts `N_λ`, minus the three the mean field already fixes.

2. **The two "mechanisms" of non-spectrality are one mechanism graded by `|λ|`.** Single
   cycles (`|λ|=1`, the `c_{odd≥7}`), disjoint pairs (`|λ|=2`, the `D_{l,l'}`), disjoint
   triples (`|λ|=3`, `T333`,`T335`) — these are the levels of the correlation tower, and the
   partition function counts them all uniformly. "Overlap defects" (`Q44, TF, p33`) are not a
   separate species; they are the inclusion–exclusion *change-of-basis coefficients* that
   convert `N_λ` to spectral + lower carriers, and they vanish from the intrinsic count.

3. **The `−3` is the mean-field's exact reach.** The spectrum determines exactly the packings
   the trace can reconstruct: nothing (`∅`), the shortest cycle (`c3 = tr A³/3`), and the next
   (`c5 = tr A⁵/5`). `c7` and up are already invisible (THM-118 fails first at length 6). The
   partition function starts counting non-spectral structure the instant a configuration is
   either long enough or multi-component enough to escape a single closed walk.

## Honest status

- The decomposition `H = Σ_λ 2^{|λ|} N_λ` is **proved** (regrouping the canon OCF closed
  form by length-multiset).
- `N_∅, N_{3}, N_{5}` spectral / the rest's onset arithmetic: **proved** (THM-002, THM-118).
- `dim ≤ #{λ:Σλ≤n} − 3`: **proved** (upper bound — `H` cannot depend on more than the
  carriers in its own expansion).
- **Equality** (i.e. all non-trivial `N_λ` are linearly independent within cospectral
  classes, none spectrally pinned): **verified by cospectral sampling at n=8,9,10,11** (rank
  matches `3,5,7,9`; H always in the carrier span; every carrier drop-one independent — at
  n=11 all nine of `{c7,c9,c11,D33,D35,D37,D55,T333,T335}`, 704 split cospectral classes,
  OCF holds 5000/5000). **Conjecture for general `n`.** The one place
  this could fail is a "tightness pinning" — some `N_λ` becoming a spectral function of the
  others when `Σλ` is close to `n` (as `Q44` was pinned in the trace basis at n=8). The data
  says the OCF carriers are *free at onset* (`D35` is free at n=8 though it uses all 8
  vertices; `D37,D55` free at n=10 using all 10), which is why the partition count is exact.

The frontier: prove the independence (no spectral pinning) of `N_λ`, which would upgrade the
growth law from conjecture to theorem; and identify `Π 1/(1−x^{2k+1})` partial sums in OEIS.

*Cross-refs: THM-505, HYP-2513, OPEN-Q-093, the-overlaps-stop-being-shadows-the-correlation-tower,
the-spectral-resolution-ladder-of-the-ocf, the-zeta-function-and-the-ocf-read-complementary-halves.*
