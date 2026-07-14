---
id: THM-766
title: The first-window tooth ladder for tight sets — every tight n-speed set has penultimate span a_(n-1)/a_1 at least n^2/(n+2); below span n its first safe window must occupy one explicit danger tooth, giving exact projective cones and a finite k-ladder in the settled-core range
status: PROVED (elementary first-window and one-tooth geometry). The bound 1 <= k <= n-1 additionally uses M(A minus {a_n}) >= 1/n through THM-759; it is unconditional for n=12 by settled LRC(<=13)
source: codex-2026-07-14-S3 (uniform audit of the n=12 sporadic tight branch)
depends_on:
  - THM-759   # only for the finite upper range k <= n-1; Part D independently refines its cap
related:
  - THM-619   # the general component-tooth band system in the LRC14 confinement setting
  - THM-668   # core maxima occur at rational pair-sum events
  - THM-757   # the LRC13 tight-locus rigidity target
  - THM-763   # uniform finite height; complementary global compression
  - THM-765   # hereditary primitivity and full component-tooth deck
  - HYP-6820 # the two-part uniformity audit
external: LRC(<=13) SETTLED (for M(P) >= 1/12 in the n=12 corollary)
---

# THM-766 — The first-window tooth ladder for tight sets

## Statement

Let `n >= 2` and let

`A = {a_1 < a_2 < ... < a_n}`

be distinct positive integer speeds with

`M(A) = 1/(n+1)`.

Write

`m = a_1`, `b = a_(n-1)`, `w = a_n`, `P = A minus {w}`, and
`delta = 1/(n+1)`.

Then the following hold.

### A. Reciprocal lower bound for the top speed

> **`w >= n m`.**

### B. The first-window tooth cone

If `b < n m`, define the nonempty open interval

> **`J = (1/((n+1)m), n/((n+1)b))`.**

Every speed in `P` is strictly `delta`-safe throughout `J`. Tightness therefore
forces `J` into one closed danger tooth of `w`. Consequently there is an integer
`k >= 1` such that

> **`((n+1)k - 1)m <= w`,**
>
> **`n w <= ((n+1)k + 1)b`.**

Equivalently, `(b/m,w/m)` belongs to the `k`-th projective tooth cone. In
particular,

> **`b/m >= n((n+1)k - 1)/((n+1)k + 1) >= n^2/(n+2)`,**

and comparison of interval widths gives the independent cap

> **`w <= 2mb/(nm-b)`.**

If instead `b >= nm`, the first inequality `b/m >= n^2/(n+2)` is automatic.
Thus every tight `n`-speed set satisfies the unconditional penultimate-span
bound

> **`a_(n-1)/a_1 >= n^2/(n+2)`.**

### C. Finite tooth range when the core floor is available

Assume additionally that `M(P) >= 1/n`. Then THM-759 gives `w <= nb`. In the
first-window case `b < nm`, the tooth label satisfies

> **`1 <= k <= n-1`.**

For `n=12`, the elementary penultimate-span bound specializes to

> **`a_11/a_1 >= 72/7`.**

Settled LRC for the eleven-speed core, through THM-759, additionally makes the
tooth range finite: if `a_11 < 12a_1`, there is a label
`k in {1,...,11}` for which

> **`(13k-1)a_1 <= a_12`,**
>
> **`12a_12 <= (13k+1)a_11`.**

The corresponding span thresholds form an increasing ladder

> **`a_11/a_1 >= 12(13k-1)/(13k+1)`.**

For example, the first three thresholds are `72/7`, `100/9`, and `57/5`.
Hence `a_11/a_1 < 100/9` forces `k=1`, in which case

> **`12a_1 <= a_12 <= 7a_11/6`.**

### D. Exact center-alignment refinement

Assume `mu = M(P) > delta`, let `t_0` maximize the core, and put
`r = (mu-delta)/b`. Then the open arc `|t-t_0| < r` is contained in the
strict safe set of `P`. Tightness forces this whole arc into one `w`-tooth,
so

> **`||w t_0|| <= delta - wr`.**

In particular `wr <= delta`, which recovers THM-759's cap

> **`w <= b/((n+1)mu-1)`.**

Unlike the scalar cap, the displayed center inequality retains the arithmetic
alignment. By THM-668, `t_0` can be chosen at a rational pair-sum event, so it
is an exact residue-band condition.

## Proof

### Step 1: `w >= nm`

Take `t = 1/(m+w)`. For every `a in A`,

`m/(m+w) <= at <= w/(m+w) = 1-m/(m+w)`.

This interval lies in `[0,1]`, and its two endpoints have the same circular
distance `m/(m+w)` from the nearest integer. Therefore

`min_(a in A) ||at|| >= m/(m+w)`.

Since `M(A)=1/(n+1)`, we must have

`m/(m+w) <= 1/(n+1)`,

which rearranges to `w >= nm`.

### Step 2: the first window is strictly safe for the core

Suppose `b<nm`. Then

`1/((n+1)m) < n/((n+1)b)`,

so `J` is nonempty. For `p in P` and `t in J`, using `m <= p <= b`,

`pt > p/((n+1)m) >= 1/(n+1) = delta`,

while

`pt < np/((n+1)b) <= n/(n+1) = 1-delta`.

In particular `0<pt<1`, so no reduction modulo one is hidden here and

`||pt|| > delta`

for every `p in P`. Thus `J` lies in

`G_P = {t : min_(p in P) ||pt|| > delta}`.

### Step 3: endpoint-safe one-tooth containment

Let

`D_w = {t : ||wt|| <= delta}`.

If any `t in J` lay outside `D_w`, then all speeds in `A` would be strictly
farther than `delta` at `t`, contradicting `M(A)=delta`. Hence `J subset D_w`.

On the real line,

`D_w = union_(j in Z) [(j-delta)/w,(j+delta)/w]`.

Because `delta<1/2`, these closed teeth are pairwise disjoint and separated by
open gaps. The connected interval `J` must therefore lie in one tooth, say the
one labelled `k`. The zero tooth cannot contain `J`: its right endpoint is
`delta/w`, whereas the left endpoint of `J` is `delta/m > delta/w` because
`w>m`. Thus `k>=1`.

Although `J` is open, containment in the closed `k`-th tooth gives weak
endpoint inequalities by taking one-sided limits:

`(k-delta)/w <= 1/((n+1)m)`,

`n/((n+1)b) <= (k+delta)/w`.

Substituting `delta=1/(n+1)` and clearing the positive denominators gives

`((n+1)k-1)m <= w`,

`nw <= ((n+1)k+1)b`,

as claimed. This limit argument is why the theorem uses weak inequalities even
though the core-safe window is open.

### Step 4: the span ladder and the width cap

Combining the two cone inequalities yields

`b/m >= n((n+1)k-1)/((n+1)k+1)`.

The right-hand side is strictly increasing for `k>=1`; its value at `k=1` is

`n*n/(n+2) = n^2/(n+2)`.

If `b>=nm`, then `b/m>=n>n^2/(n+2)`, so the same universal penultimate-span
bound holds in both cases.

Also,

`length(J) = (nm-b)/((n+1)mb)`.

One `w`-tooth has length `2/((n+1)w)`. Since `J` fits in one tooth,

`(nm-b)/((n+1)mb) <= 2/((n+1)w)`.

All factors are positive in the first-window case, and rearrangement gives

`w <= 2mb/(nm-b)`.

### Step 5: the finite range for `k`

Under `M(P)>=1/n`, THM-759 gives `w<=nb`. If `b<nm`, then `w<n^2m`.
If `k>=n`, the first cone inequality would instead give

`w >= ((n+1)n-1)m = (n^2+n-1)m > n^2m`,

a contradiction. Hence `k<=n-1`. The `n=12` formulas follow by direct
substitution. The second threshold is

`12(2*13-1)/(2*13+1)=12*25/27=100/9`,

so a ratio below `100/9` leaves only `k=1`; its upper cone inequality is
`12a_12<=14a_11`, or `a_12<=7a_11/6`.

### Step 6: center alignment

For `|h|<r=(mu-delta)/b` and every `p in P`, the triangle inequality for
distance to the integers gives

`||p(t_0+h)|| >= ||pt_0||-p|h| >= mu-b|h| > delta`.

Thus the open radius-`r` arc around `t_0` lies in `G_P` and hence in one
closed `w`-tooth. Choose a real lift and let `j/w` be that tooth's center.
Taking limits at both ends of the symmetric open arc gives

`|t_0-j/w|+r <= delta/w`.

Multiplication by `w` yields

`|wt_0-j| <= delta-wr`.

The right side is at most `delta<1/2`, so `j` is a nearest integer and the
left side is `||wt_0||`. This proves the center-alignment inequality. Its
right side must be nonnegative, giving `wr<=delta` and therefore

`w <= delta*b/(mu-delta) = b/((n+1)mu-1)`.

This completes the proof. ∎

## Exact component-containment interpretation

For any core `P` and proposed top speed `w`, define `G_P` and `D_w` as above.
Then

`G_P subset D_w  implies  M(P union {w}) <= 1/(n+1)`.

Conversely, tightness implies this containment. Whenever the matching LRC
floor `M(P union {w})>=1/(n+1)` is known, containment is therefore equivalent
to tightness. Componentwise, every connected component `[a,b]` (with endpoint
closure understood) must satisfy, for some integer tooth label `j`,

`wa >= j-delta` and `wb <= j+delta`.

This is the theorem-facing object. A quotient that remembers only total safe
measure, runner scores, or pairwise overlaps loses simultaneous containment.

## Measure-only guardrail, independently checked

The natural sufficient strategy `meas(G_P)>meas(D_w)=2/13` cannot hold
uniformly even on nonextremal eleven-speed cores. For

`P={1,2,3,4,5,6,7,8,9,10,12}`,

one has `M(P)=1/11>1/12` without a numerical maximization: adding `12` to
`{1,...,10}` cannot raise its standard AP value `1/11`, while `t=1/11`
gives every speed in `P` distance at least `1/11`. An exact breakpoint
partition at level `1/13` gives ten open safe components:

```
(7/78, 6/65)       (7/39, 12/65)      (7/26, 25/91)
(14/39, 19/52)     (53/117, 6/13)     (7/13, 64/117)
(33/52, 25/39)     (66/91, 19/26)     (53/65, 32/39)
(59/65, 71/78).
```

Their lengths sum to

`461/8190 = 0.056288... < 2/13`,

with exact deficit `2/13-461/8190=799/8190`. The same value was independently
recovered as one minus the union length of all eleven danger combs
(`1-7729/8190`). Thus scalar measure is too coarse: the uniform sporadic-branch
proof must retain component addresses and their simultaneous tooth bands.

## Tournament Analysis and assumption challenge

The natural exact vertices here are **safe components / tooth obligations**,
not runners. A useful pairwise observable is whether two components admit a
common compatible tooth label after the projective `(m,b,w)` gauge; cyclic
endpoint order supplies the tie Hamiltonian path. Score histograms, cycles,
SCCs, edge flips, and Hamiltonian-path counts can diagnose compatibility.

This tournament is still a lossy shadow: pairwise compatibility does not imply
simultaneous compatibility of every component. The theorem's exact predicate
is the component-tooth incidence hypergraph with endpoint widths, center
residues, divisor pins, and the rational core-maximum sidecar retained. The
challenged assumption is that runners should be the vertices; proof
obligations preserve the tightness predicate more faithfully.

## Scope

THM-766 is a uniform reduction, not a proof that the `n=12` sporadic branch is
empty. It removes every proposed tight twelve-speed set with
`a_11/a_1<72/7`, stratifies the remaining sub-`12` projective range into eleven
explicit cones, and strengthens THM-759's scalar cap to an exact alignment
band. The residual is to prove that every primitive nonextremal eleven-speed
core has at least one safe component whose band is incompatible with every
top speed in its allowed cone.
