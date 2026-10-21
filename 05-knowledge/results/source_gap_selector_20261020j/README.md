# A positive complement register and same-source debt repair

20 October 2026. **PROVED:** a positive companion to discarded debt for
every actual positive odd history; a half-descent payment theorem; a
three-letter native repair after any descent; an exact rational cap selector;
and the optimal affine state barrier in a precisely stated relaxation.
**FINITE-EXACT:** the inherited175 supplied first exits and four named
source histories. Existence of a qualifying descent or repair run for every
source, first-exit payment without restriction, and Collatz convergence remain
**OPEN**. No literature-priority claim is made.

Files: [checker](check.py), [four supplied histories](inputs.json),
[receipt](receipt.json). Production does not search for a first exit or ROOT.

## 1. Inheritance, the missing coordinate, and the new target

The closest mechanism is the [unit-cap and source-gap controller](../grounded_block_children_20261020i/README.md):
a terminal word with multiplier μ and discarded debt D preserves [0,c)
exactly when μ+D/c≤1. Combined with an ordinary cone prefix and an actual
source gap at least c+1, this pays the complete word. Its cap3 extension
completed the declared5000-source finite bank; it did not prove entry into
such a controller for an arbitrary source.

The [earlier first-exit payment packet](../first_exit_complement_games_20261020/README.md)
already has a sharp all-length coarse-cylinder bound, a pre-step last-letter
threshold, and finite-depth and logarithmic-depth regions. We reuse their
scope rather than present the same threshold as a new universal argument.
Its source55 later descent to53 is the canonical hostile: actual descent
alone does not pay discarded debt.

The incoming [THM-4735, valuation-dependent cones and the debt-two obstruction](../../../01-canon/theorems/THM-4735-valuation-dependent-debt-cones-optimal-terminal-splits-and-the-debt-two-obstruction-to-universal-cone-coverage.md)
records its pre-audit status separately from exact finite data. Its source4575
and437247 controls show why a uniform threshold2 is insufficient. We do not
promote the large census or a pre-audit label; the elementary identities used
here are proved below and reproduced independently by exact arithmetic.

The corrected near miss is the quantifier switch

    for every head, some inflated new source completes it
                 ⇒ every original source reaches a good terminal.

That implication is unavailable. This packet does not inflate sources.
The least-used sidecar is the **positive complement x+1−D**, rather than
another upper estimate on D alone.

Live board: **original source / actual odd state / ordered valuation /
positive complement / normalized reset sum / strict ROOT stop**.
The anchor is fixed-source payment; the niche is an invariant positive cone;
the wildcard is a sharp limitation on all affine state-only debt barriers.

Connection contract: an authenticated history is mapped to (x,λ,K), retaining
its original source and word. Payment, its exact gap, and the monotone
normalized register can be reconstructed. This map alone does not retain
the next native valuation or establish that a future stopping event exists.
The cheapest hostile is source55's unpaid later descent; the positive control
is the next actual reset at the same source.

## 2. The companion coordinate that never changes sign

Let n>1 be positive and odd, and let w=(a₁,…,aₘ) be an actual finite word
of positive odd Collatz valuations. Write x₀=n and

    xᵢ=(3xᵢ₋₁+1)/2^aᵢ,     aᵢ=v₂(3xᵢ₋₁+1).

Stop at the first1; ROOT padding is not a permitted word. Set λ₀=1,D₀=0,
and for μᵢ=3/2^aᵢ define

    λᵢ=μᵢλᵢ₋₁,
    Dᵢ=μᵢDᵢ₋₁+δ(aᵢ),
    δ(1)=δ(2)=0,    δ(a)=1−4/2^a for a≥3.             (1)

Thus D is exactly the inherited nonnegative discarded debt, not a signed
ROOT charge or an annulus-energy increment. Define

    Kᵢ=xᵢ+1−Dᵢ,             Eᵢ=xᵢ−1+Dᵢ=2xᵢ−Kᵢ.     (2)

**Positive complement theorem.** Every such history satisfies

    Kᵢ>0,                   0≤Dᵢ<xᵢ+1.               (3)

More precisely,

    Kᵢ=3Kᵢ₋₁/2                         if aᵢ=1,
    Kᵢ=(3Kᵢ₋₁+2)/2^aᵢ                 if aᵢ≥2,      (4)

and

    Kₘ/λₘ = n+1 + (2/3) Σ_(i:aᵢ≥2) 1/λᵢ₋₁.        (5)

Proof: substitute (1) and the literal state equation into xᵢ+1−Dᵢ.
The initial value K₀=n+1 is positive; both maps in(4) preserve strict
positivity. Divide by λᵢ and telescope to obtain(5).

This proves a source-specific lower bound that survives arbitrary native
refinement: K/λ never decreases. A letter1 leaves it unchanged; every
letter at least2 adds exactly 2/(3λ_previous). The increment is independent
of the size of that non-one valuation. Nevertheless the **unnormalized**
K also carries the changing multiplier λ; a lower bound on K/λ is not
automatically a time-uniform lower bound on K.

The positive envelope has the complementary recurrence

    Eᵢ=3Eᵢ₋₁/2+1                  if aᵢ=1,
    Eᵢ=3Eᵢ₋₁/2^aᵢ                if aᵢ≥2.           (6)

Equations(5) and(6) split the inverse-slope sum into non-one and one
positions. They retain the order and actual clocks, not just a reset count.
The raw carrier and an independent inverse-word reader authenticate the
original integer source in every saved receipt.

## 3. Payment from half descent, with no first-exit premise

The complete word pays the original source exactly when

    Dₘ<n−h,     equivalently Eₘ<n−1,
    equivalently Kₘ>2h−n+1,             h=xₘ.         (7)

**Any-history half-descent theorem.** If an actual nonempty history from
odd n>1 reaches h with 2h+1≤n, then it pays:

    Dₘ<h+1≤n−h.                                     (8)

This allows arbitrarily many earlier resets and earlier descents. It does
not require a cone cut, a uniform debt cap, or first-exit minimality.
In particular every supplied strict ROOT word pays: h=1 and D<2≤n−1.
This last statement is conditional on possessing the actual ROOT word;
it is not a way to produce one for an arbitrary source.

The boundary cannot be sharpened to D<h. The supplied strict ROOT word of
source123027327 has125 odd letters, valuation cost225, and

    D=210645860432040711441183989615078067661202342765484255223476419889
      /210624583337114373395836055367340864637790190801098222508621955072
      >1=h,

while D<h+1 remains strict. No floating comparison is used.

Among the inherited175 named first-exit words,124 satisfy(8) directly.
This is a new proof route for those supplied records, not175 new sources
or a claim that the remaining51 are unpaid. The older complete controller
certificates for all175 remain intact.

## 4. A bounded repair after any actual descent

Suppose a supplied head takes n>1 to h<n. If h=1, it already pays by(8).
Otherwise retain the **actual next letters**, with no freedom to tune them.

**Three-letter repair theorem.** If the next three valuations are all
at least2, the concatenated word pays the same original source n. If ROOT
occurs earlier within that run, stop there and it already pays.

Proof: at the descending checkpoint E<2h by(2)–(3). A run of r valuations
at least2 multiplies E by at most(3/4)^r, by(6). Since h≤n−2,

    E_after3 < (27/32)h ≤ (27/32)(n−2) < n−1.         (9)

No new source, ROOT padding, or repeated budget is introduced. More
generally the explicit sufficient run condition is

    2h(3/4)^r≤n−1.                                  (10)

The required native run is a remaining guard, not a control action we may
impose on the fixed orbit. The26 immediate repair examples in the receipt
are taken only from the inherited175 heads and their supplied child suffixes.
The count includes shorter runs that actually reach ROOT.

Three is sharp **given only the relaxed information D<h+1**. For the
synthetic incoming state n=131,h=129,D=129, the literal native tail is

    129 →97 →73 →55,       valuations(2,2,2).

Here E=257, so two twos leave E=2313/16≥130, while three leave
E=6939/64<130. This incoming debt is not asserted to arise on a history
from131. The example establishes the limit of the history-forgetting
invariant, not the existence of a native first-exit obstruction requiring
exactly three letters.

The old hostile is repaired without changing its source: source55's later
39-letter endpoint53 is unpaid, but its next actual valuation5 takes it
to5 and pays. The saved strict ROOT word then reaches1 with valuation4.
Its earlier first descent was at step3, which is why the unpaid39-letter
prefix does not refute first-exit payment.

## 5. Exact selectors and the obligation they do not discharge

For a retained non-ROOT head, the least integer a≥2 satisfying

    2^a(n−1)>3E_head                                (11)

is the exact threshold for the next non-one letter to pay. It is computed
before taking that step. This is the inherited last-letter threshold,
now expressed in the companion coordinate. An actual valuation below
this threshold is not replaced by a larger one. Existence of a future
head where the actual valuation meets(11) is not proved.

There is also an exact cap-architecture selector. For a nonempty terminal
word v of multiplier μ and debt D_v, a finite rational cap c≥1 exists
only if μ<1, and its least value is

    c_min(v)=max(1,D_v/(1−μ)).                        (12)

A prescribed actual gap g≥2 admits such a cap with c+1≤g exactly when

    μ<1 and c_min(v)≤g−1.                            (13)

The proof is simply solving μ+D_v/c≤1 for c; nonempty positive valuation
words cannot have μ=1 because3^m≠2^A. At each proposed split the ordinary
cone-prefix condition must still be checked. Minimizing(12) over eligible
splits is necessary and sufficient **for this architecture**, not for
positive payment in general. `cap_selector` checks this iff on supplied
histories; it does not call a first-exit finder or certify future arrival.

On the four supplied controls and all175 prior words, the checker retains
each complete state, source, first-exit boundary and strict ROOT suffix.
The source437247 first exit has92 letters and debt exceeding3; (3) remains
valid at every prefix and(12) supplies a source-gap certificate without
forcing the debt below2. No larger census is rerun.

The exact universal remaining obligation can be stated in two distinct
ways, which must not be conflated:

1. **Conditional payment:** if a first exit exists, prove (7) there.
   The invariant proves it on the half-descent region but leaves the band
   (n+1)/2≤h<n.
2. **Fixed-source termination:** prove that every n>1 reaches a paid
   descent, a half-descent, or another authenticated terminating controller.
   Such an all-source paid-descent theorem would imply convergence by
   well-founded descent. Neither (3) nor a finite test bank proves it.

In particular `for every head there exists a new inverse source` from the
previous factory does not provide the second quantifier.

## 6. The affine barrier is optimal in its stated relaxation

Consider the local region

    x>1 positive odd,       0≤D<αx+β,

where α,β are nonnegative rationals and α+β>0. Permit every incoming
rational debt in this region, even one not generated by a prior source.
Require preservation under every actual native step at x, including a
step terminating at1. Then:

    the region is invariant  iff  α≥β≥1.             (14)

Here invariant means the output obeys0≤D'<αx'+β; the terminal1 need not
have an outgoing step. Thus α=β=1 is the pointwise smallest member of
this entire affine family, and gives exactly(3).

Proof: the slack at the upper boundary is

    (α−β)/2                              at a=1,
    (α+β)/4                              at a=2,
    β−1+(α−3β+4)/2^a                    at a≥3.      (15)

The native valuation1 at x=3 forces α≥β. Arbitrarily large even a occur
at x=(2^a−1)/3, so the last expression forces β≥1. Conversely these
two inequalities make every slack nonnegative: for a≥3 its minimum is
positive already at a=3 when its numerator is negative, and its limit
β−1 is nonnegative otherwise. Multiplication by μ>0 preserves the strict
incoming inequality even when the boundary slack is zero.

If a condition fails, take an actual edge with negative boundary slack
and choose rational D just below αx+β; the output violates the barrier.
The checker supplies exact rational examples. The excludedα=β=0 has
empty debt domain and would be vacuously invariant.

This does not rule out a stronger history-aware nonlinear invariant.
It explains why merely optimizing a single affine function of the current
state cannot improve the half-descent threshold derived from(3). Retaining
the normalized sum(5), source residue, or a history restriction can still
sharpen the bound; those extra coordinates are precisely what(14) forgets.

## 7. Reproduction and limits

The finite universe is the175 inherited supplied heads, their113 retained
child ROOT words, and four further supplied histories at55,63,437247 and
123027327. The latter were discovered once with a1000-odd-step cap and
stored in inputs.json; production only replays them. It checks every prefix
of these four histories, and independently reverses every complete word.

Independent symbolic controls cover all words of lengths0 through5 over
valuations1 through4; those rational affine words are never mislabeled
as native integer sources. Rational cap comparisons and inverse-slope
sums are computed by separate formulas. The affine-barrier grid has
α,β∈{0,1/4,…,2}, with its zero-domain case excluded. Type, changed-source,
changed-register, missing-half-gap, wrong-tail and ROOT-padding hostiles
remain active under optimization.

```text
python -B -X utf8 05-knowledge/results/source_gap_selector_20261020j/check.py
python -B -O -X utf8 05-knowledge/results/source_gap_selector_20261020j/check.py --verify
```

The intrinsic operations are ordered affine updates and exact stopping
guards. There is no pairwise binary relation here that benefits from being
forced into a tournament. The signed ROOT-current/annulus ledgers retain
different weights and clocks; this discarded-debt result does not settle
their uniform energy inequality.

Normal generation and optimized immutable replay both pass 15,631 checks.
The LF-normalized receipt SHA256 is
`b6f67b7af556cc7479d7029b7084fcd6f03df95425bd491d5c6d6088deef9b02`.
Independent peer review by `tiling_fold` passed the full proof, production
interfaces, and optimized saved-receipt comparison. In particular it checked
the original-source half-descent and three-letter repair, the synthetic-state
scope of the affine-barrier iff, and the inherited architecture required by
the minimal-cap selector. No correction was requested.
