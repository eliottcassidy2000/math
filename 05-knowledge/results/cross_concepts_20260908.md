# Past concepts that change the live questions

**Status: PROVED scoped results + FINITE-EXACT controls + independent audits.
LRC(14), the uniform coin optimum, the named Mahler/Rule 30 questions, and
general smooth-fluid regularity/singularity frontiers remain OPEN.**

This pass recovered old knot continuation, finite-state realizations,
checksum orientations, support/multiplicity distinctions, and integral graph
repair. Four audited theorem files now turn those connections into exact
statements. The most productive question was: **after an operation has
forgotten information, what must be retained to make its next step lawful?**
The answers differ by problem: a continuation state, an integer carry,
a repair capacity, and spatial matrix minors.

This continues the [Euler bridge](euler_bridge_20260908.md), including its
[pinned external-source audit](../reference/OPENAI-EULER-AUDIT-2026-09-08.md).
It does not upgrade the external repository's full PDE verification status.
The inherited base is `47bf7e622`; the four namespaces were explicitly
reserved in `8b25fa628` before promotion. Current canon, corrections and
reproducible computation were used over historical prose.

## What was recovered and what changed

| Recovered concept | Exact new object | Change to the live question |
|---|---|---|
| Gordian continuation and common-context localization, THM-2176/2191 | [THM-4460](../../01-canon/theorems/THM-4460-adverse-budget-continuation-state-and-shear-composition-collapse.md): the adverse budget's complete future state | Ask which observation supports the actual next operation; favorable present cost alone cannot be discarded |
| Cyclic coin checksum and moving donors, THM-2225/3340/3344 | [THM-4461](../../01-canon/theorems/THM-4461-checksum-orientation-transcendence-and-nonlinear-dyadic-carry.md): transcendental orientation germs with a three-component nonlinear scale recurrence | Replace a stationary linear row state by an explicit scale action, while retaining its causal decoder |
| Integral graph repair quotient, THM-3990 | [THM-4462](../../01-canon/theorems/THM-4462-cyclic-repair-capacity-frontier-and-reflection-torsion-collapse.md): a sharp cyclic cancellation capacity and complete variance curve | Once symmetry kills the obstruction group, measure the size of the allowed repair |
| Support versus multiplicity, THM-2352; unrestricted realization versus physical transport, THM-3163 | [THM-4459](../../01-canon/theorems/THM-4459-pure-shear-mixtures-have-exact-nuclear-moment-cost.md): exact shear-mixture cost and sharp spatial-minor inequality | Solve the arbitrary-law relaxation, then enforce the spatial gradient identities it loses |

The linked new files give the full paths of the older theorem IDs, including
the ambiguous legacy THM-853 slug used in the continuation lane. None of the
generic old mechanisms is advertised as a newly discovered theorem.

## 1. LRC: the current scalar budget is not the continuation state

For the fixed centered observer `r=R-mean R`, THM-4458 defined

```
B_r(L)=min_b mean (-r_s(L_s-b))_+.
```

THM-4460 identifies this exactly as optional transport from positive-role
capacity to negative-role capacity, with profit `L_negative-L_positive`.
It proves

```
B_r(L+K)=B_r(L'+K) for every real K
iff L-L' is constant on the coordinates where r is nonzero.
```

Thus the complete future state keeps the visible vector modulo constants.
The favorable cone `B_r(L)=0` is much larger. For example, `0` and `r`
both cost zero, but the same future addition `-r` separates them.

The unrestricted-addition theorem was then tested against genuinely legal
operations in the auxiliary indicator-density model. Activating any of
12 labelled phase densities is a commuting idempotent operation. Among
all 4,096 masks there are only 19 current budget values, but their classes
under continuations of length at most d have exact counts

```
d:       0     1     2     3     4     5     6     7
classes:19  2434  3148  3688  3973  4073  4094  4096.
```

Seven is necessary and sufficient. Two independent finite classifiers and
an analytic shortest separating pair agree. The exact adverse budget
certifies all 4,095 nonconstant masks in this family; the inherited norm
gate certifies 34. This is an improvement in a fully declared auxiliary
universe, not a count of solved canonical LRC configurations.

The **next problem** is to determine the observation closed under the
actual canonical owner/word action. The finite example shows why labels
can be needed even when the present scalar has a large permutation symmetry.
Unrestricted real additions are not a substitute for that owner action.

## 2. LRC: repairability becomes a sharp capacity problem

On a cycle, an integral zero-mean profile `a` has an integral Laplacian
repair exactly when `sum s*a_s=0 mod n`. Reflection symmetry forces that
class to be two-torsion, so it vanishes for every odd n. This rules out
the proposed torsion obstruction for the symmetric danger-overlap profile.

THM-4462 retains the size of the repair. For `div f_s=f_s-f_(s-1)` and
uniform edge capacity `|f_s|<=k`, exact cancellation occurs iff

```
k >= (max_s A_s-min_s A_s)/2,   A_s=sum_(j<=s) a_j.
```

For the canonical auxiliary 13-phase role, after clearing denominator
1183, `a=(144,53,-25,...,-25,53)`. Its sharp cancellation threshold is
`k=125`, or **125/1183 in the original density units**. Its minimum
remaining variance, over every flux within the capacity, is

```
[(144-2k)^2+2*53^2+10*(-25+k/5)^2]/(13*1183^2),  0<=k<=91/2;
(250-2k)^2/(30*1183^2),                         91/2<=k<=125;
0,                                             k>=125.
```

Explicit extremizers and a direct variational certificate prove the entire
curve. An integral repairing potential exists, but its unavoidable range
is 447. The result generalizes to the specified overlap family at all odd
`p>=5`; it does not require p prime.

The **next problem** is an a priori capacity bound for flux constructed
from the actual physical owner action on the same cyclic ordering. Every
centered leak has some flux, so mere flux existence adds no obstruction.
The new sharp threshold makes the required missing estimate concrete.

## 3. Fair coins: nonrationality does not exclude a nonlinear scale state

The old checksum extractor has an orientation transfer that total fairness
hides. Let `p=P(0)` and let `F(p)` be its heads probability with initial
zero. Let `G(u)` be initial-one heads probability in `u=P(1)`. THM-4461
recovers the literal constant-half words and obtains

```
L(t)=sum_(r>=0) t^(2^r),
F(p)=(3p-2p^2-L(p-p^2))/2,
G(u)=(2u^2-u+L(u-u^2))/2.
```

Both germs are transcendental. An elementary logarithmic-growth argument
for L proves this, without a finite-prefix extrapolation. Yet the exact
finite shell sums have an autonomous three-component polynomial recurrence
at the dyadic scale, explicitly supplied in the theorem.

The lost coordinate is particularly visible modulo two:

```
L(p-p^2)=p in F_2[[p]],
C(p)=(L(p-p^2)-p)/2 in Z[[p]],
C(p)=-p^3+p^4-2p^5+3p^6-... .
```

Parity reduction has erased a transcendental series; the divided carry
recovers it. The component count is fixed, but the degree, scale index and
integer state range are unbounded. This is compatible with the existing
obstruction to fixed finite-dimensional linear row realizations.

The **next problem** remains a cross-shell donor decoder with improved
deadlines. This exact nonlinear state is a positive control that such a
proposal can be tested against. It changes no deadline by itself. The
Rule 30 connection is an exact warning about parity information loss;
no map to its physical temporal dynamics is supplied. The Mahler `3/2`
reset comparison stops at its distinct native arithmetic state.

## 4. Fluids: arbitrary shear mixtures pass, gradient mixtures must pay

THM-4460 first shows that every nonnegative subadditive cost below pure-shear
distance collapses to zero. Two pure shears already sum to a matrix with
nonzero self-stretching. Thus the catalytic metric construction cannot be
applied to that distance as though it were a metric on matrix differences.

THM-4459 solves the stronger arbitrary-mixture problem in every dimension.
For any real trace-free M, over laws on pure shears S with `E S=M`,

```
min E||S||_F^2=||M||_*^2,
min E||S-M||_F^2=||M||_*^2-||M||_F^2.
```

Here `||M||_*` is the sum of singular values. A minimizing law uses exactly
rank M nonzero atoms of equal norm. The elementary proof combines polar
factorization with the classical orthogonal zero-diagonal basis lemma.

The initial cheap hostile led to a second, more physical pull. In an
affine-periodic cell `A=M+grad v`, `v in W^(1,2)`, the mean of every
second minor equals the corresponding minor of M. This is a classical
null-Lagrangian identity. It yields the **sharp** three-dimensional bound

```
||cof M||_F^2 <= mean||A||_F^2 * mean d_F(A,pure shears)^2.
```

For `M=E12+E23+E31`, the arbitrary equal mixture of `3E12,3E23,3E31`
has second moment 9 and zero shear distance at every atom, while M has
stretching rate one. Any affine-periodic gradient field with that same
mean and second moment at most 9 must instead have mean squared shear
distance at least **1/3**. This distinguishes an arbitrary moment law from
a spatial gradient using a precise invariant and the same energy budget.

The theorem credits the classical hollow-basis and null-Lagrangian facts,
supplies their needed proofs, and makes no global priority claim for the
derived bounds. These results concern squared velocity gradients, not fluid
kinetic energy. A globally periodic velocity has mean gradient zero; the
nonzero-mean cell is affine plus periodic. Pressure evolution, localization,
boundary conditions and temporal ordering still have to be established for
any proposed Euler construction.

## Cross-check of the concept board

| Live concept | What the other lanes now require it to retain |
|---|---|
| Present observable versus future action | LRC's 19-to-4096 experiment makes loss of labels explicit; coin carries likewise cannot be dropped just because a quotient resets |
| Allowed repair and its cost | Symmetry kills cycle torsion but leaves capacity; shear support is algebraically universal but leaves moment cost and spatial minors |
| Scale state versus fixed row state | The coin recurrence survives at dyadic scales; no implication transfers to native Mahler addresses or a changed LRC observer |
| Realization versus aggregate feasibility | THM-3163's old automatic realization boundary reappears as a solved shear moment problem and a failed gradient realization |
| Equality and failure boundaries | Full activation is the unique cancellation mask; flux capacity125 is exact; shear mean variance has an attaining law; periodic minors exclude that law in the named cell |

The Anchor was LRC continuation/cancellation, the Niche was coin scale
transport, and the Wildcard was shear aggregation. The wildcard produced
a spatial obstruction after its arbitrary-law relaxation proved too broad.
No artificial tournament was introduced: the genuine carriers were a
labelled activation monoid, a cyclic flux graph, a dyadic recurrence and
a matrix moment law.

The exact connection maps and discarded coordinates are stated separately
in each theorem. There is no asserted coin-to-Mahler, coin-to-Rule30, or
auxiliary-density-to-canonical-LRC intertwiner. Stopping at those specific
missing maps preserves useful progress without converting analogy into proof.

The final incoming-work pass read `1fa733c37` and `ed6517c0d`: the
[fixed-DG quartic class](planar_jc48_sep08_quartic_closure.md) is now
excluded in its stated global ring, and the
[all-m boundary calculation](planar_jc48_sep08_all_m_constant_d.md)
uses moving inverse indices `2m-3` and `6m-7`. Its fixed-index hostile
reinforces the scale-state distinction above: the observer must change
with the native scale. The corrected full numerator and volume factor
were retained in that source. These are incoming proved suppliers with
their own audits, not dependencies of the four theorems here; no map from
their rational Jacobian equation to the coin recurrence is asserted.

## Evidence and handoff

Each new theorem includes its reproduction command, explicit universe,
hostile controls and independent proof audit. The scripts use exact
standard-library arithmetic and always-active exceptions; normal and
optimized Python replays agree with the stored transcripts. No Lean
formalization is claimed.

- [Shear moment controls](cross_concepts_shear_mixture_20260908.out): 399 rational polar certificates, rank-deficient cases, and composition hostiles.
- [Gradient minor controls](cross_concepts_gradient_minors_20260908.out): independent exact Fourier convolution, matched mean/energy comparison, and sharpness controls.
- [Continuation controls](cross_concepts_catalysis_20260908.out): 7,020 profile pairs, 18,198 scanner queries, all 4,096 activation masks by two classifiers, and 6,561 matrix decompositions.
- [Coin controls](cross_concepts_donor_20260908.out): all 548 legal words through shell endpoint16, polynomial recursion through degree512, and non-dyadic/parity hostiles.
- [Cyclic repair controls](cross_concepts_tournament_20260908.out): 1,089 torsion vectors, 2,430 capacity instances, 140 exact optimum certificates, and 3,369 brute flux vectors.

The next high-value work has a precise gate in each lane: canonical labelled
owner transitions, a physical flux bound, a better causal donor decoder, or
an affine-cell/pressure-compatible shear construction paying both the spatial
minor budget here and the trajectory budget of THM-4457. Repeating an
unrestricted aggregate feasibility computation would not discharge them.
