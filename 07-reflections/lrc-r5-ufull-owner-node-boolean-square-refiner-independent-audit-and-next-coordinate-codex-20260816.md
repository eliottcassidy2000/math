# A second-path audit of the Boolean square and its cheapest ancestry crossing

**Status: FINITE-EXACT SECOND-PATH HOSTILE AUDIT ACCEPTS THE SCOPED
BOOLEAN-SQUARE REFINER; THE INVERSE-OWNER CROSSING IS AN OPEN, EXACTLY TYPED
NEXT TEST.  LRC(14) remains open.**  This audit does not import the submitted
Boolean-square candidate.  It starts from the separately audited common-base
parent, independently reconstructs the source-state partition, and checks the
result by exact marginal, rank, guard, erasure, and fixed-class gates.  Its
tables and all seven comparison digests agree with the fully independent audit
at `ca7e8935c`.

## The corrected hierarchy is load-bearing

The root spine is

```text
R={0,6,12}.
```

Its six nonempty proper subsets form a useful *abstract cut completion*.  Each
orients the two cross-cut pairs and leaves the within-side pair missing, hence
has census

```text
(missing, one-way, two-way)=(1,2,0).
```

There are no both-way edges.  The completion is not the physical source
carrier.  Reconstructing the THM-2471/THM-2594 source profiles gives exactly

```text
{0} -> {0,6} -> {6} -> {6,12} -> {12},
toggle word (6,0,12,6),
measure vector (1,12,2,12,1)/28.
```

The sixth abstract cut `{0,12}` is globally absent.  `OWNER` then removes only
the physically realized centre `{6}`.  It does not remove `{0,12}`, because
there is no such source state to remove.  This is the exact hierarchy:

```text
six abstract cuts -> five physical Gray-path states -> four owner states.
```

The four surviving states are

```text
(left,single)  {0}
(left,double)  {0,6}
(right,single) {12}
(right,double) {6,12}.
```

Each has measure `1/28`; complement toggles both intrinsic bits, hence acts by
XOR with `11`.  The four objects are states of partial three-root tournaments,
not vertices of a size-four tournament.

## Exact audit result

The new implementation imports only the hash-pinned independent parent audit,
not the candidate.  It rebuilds the source partition from the parent's exact
THM-2594 profiles and splits the independently desheeted THM-3514 endpoint
integrand before integration.  It then checks:

- the five-state path, measures, toggle word, abstract cut census, and absence
  of `{0,12}`;
- all four owner-state measures and complement involution;
- pointwise same-root zero;
- five literal delete-and-restore guard controls;
- exact summation back to every row of the audited rank-one parent;
- coupled and source-weight-erased ranks `4`, and centered rank `3`;
- both intrinsic bit marginals of rank `2`;
- full Walsh by relation support `(52,1,3,12,36)` and centered support
  `(36,0,0,0,36)`; and
- four nonzero state values and four nonzero Walsh values at fixed relation
  class `(1,0,6)`.

Normal and optimized runs are byte-identical.  The semantic digest is
`af0d543232869e82ee8d0191478ba7a833954cb19b387dedb6fb6f44a6fa272c`.

The source-erasure control remains an essential boundary: after replacing the
numerical source weights by ones, the state partition still has rank `4` and
all `36` mixed modes.  Thus the typed support partition is load-bearing, but
the computation does not prove that the numerical source weights cause the
nonseparability.

## Cheapest lawful ancestry coordinate

The least expensive canonical coordinate destroyed by the present marginal is
the last inverse-owner branch from THM-2461 and THM-2471, not the missing
abstract cut.  On the current leg of the THM-2471 stalk,

```text
d=13^5,
X_(u,a)=(w_u+a)/13^5,
T^4 X_(u,a)=(w_u+(a mod 13))/13.
```

Therefore

```text
r_owner := a mod 13
```

is the canonical last inverse branch above `w_u`.  It is already a genuine
Boolean sheet label before the `a`-sum.  It is not the deep-probe label
`r_deep`: THM-2461 explicitly distinguishes their base variables, maps, and
observables.  Equal cardinality thirteen is not an intertwiner.

This is cheaper than grouped exact `C(a;X,m)`.  The latter needs the complete
address orbit and its grouping/reselection semantics; `r_owner` retains only
one canonical radix digit and marginalizes all older ancestry.  A binary
partition of the thirteen branches would be cheaper numerically but would not
be intrinsic without a further theorem, so it is not the preferred first
test.

### Exact proposed crossing

Reopen the finite Boolean stalk in THM-2471 equation (31), retain
`r_owner=a mod 13` on the current `X_(u,a)` leg, and marginalize the older
`a`, `b`, and opposite-leg `e'` data only after forming the same common-base
endpoint product.  The output is a typed

```text
four square states x thirteen owner branches x thirteen relation residues
```

array.  Its mandatory first gate is

```text
sum_(r_owner) table(square,r_owner,relation)
  = audited Boolean-square table(square,relation)
```

row by row.  The cheapest decisive positive is not formal Fourier fullness;
it is a nonzero three-way ANOVA component, or equivalently a strict conditional
rank increase after all square, owner-branch, and relation marginals have been
removed.  Literal guards, source anchor, pointwise same-root zero, Lucas field,
and fixed `(1,0,6)` controls must be rerun on the retained-sheet object.

## Connection contract

| field | exact answer |
|---|---|
| source | audited four-state source-support square and THM-2471 Boolean stalk |
| target | `V_4 x F_13(r_owner) x F_13(relation)` exact array |
| map | retain `a mod 13` on `X_(u,a)` before ancestry marginalization |
| preserved | common base, source state, current root, literal endpoint guard, last inverse branch |
| destroyed | older `a` digits, packet `b`, opposite sheet `e'`, grouped exact address, chronology |
| required marginal | summing `r_owner` recovers the audited `4 x 13` square exactly |
| cheapest hostile | conditional rank/three-way ANOVA after all lower-order marginals |
| forbidden identification | `r_owner=r_deep`, source time = arrival time, or owner branch = exact `C(a;X,m)` |

This identifies a lawful next coordinate, not a completed transport.  No exact
`C(a;X,m)`, `U_clock` chronology, uniform-row statement, row exclusion, or
LRC(14) consequence is claimed.

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_independent_audit_20260816.py
```
