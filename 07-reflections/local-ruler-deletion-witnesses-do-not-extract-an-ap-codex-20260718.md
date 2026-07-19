# Local ruler deletion witnesses do not extract an AP

*Codex 2026-07-18.  Exact obstruction to a tempting local surrogate; this is
not a counterexample to crown collapse, the `n=12` equality conjecture, or
LRC(14).*

Boxeph-S119's confinement--coupling proof for an actual arithmetic
progression is global: a hypothetical better time confines **all twelve**
phases to one arc, and the identity `a(dt)=d(at)` couples that confinement to
the AP step.  It is tempting to replace that global condition at the crown
interface by the finite pointwise data already visible at the twelve shallow
points:

```text
every j/13 is a two-owner local deletion witness;
all active pairs use one pair-sum ruler Q;
the deleted runner Q kills every one of those witnesses.
```

That replacement is false even after adding primitivity, Cover14, and the
compact top-ratio hypothesis.

## The general local-ruler mechanism

Let a twelve-speed packet be labelled by its nonzero residues,

```text
W={w_r:1<=r<=12},             w_r=r (mod 13),
```

and suppose `Q` is divisible by `13` and

```text
w_r+w_(13-r)=Q                 for 1<=r<=6.              (1)
```

At `t=j/13`, multiplication by `j` permutes the twelve nonzero residues.
Therefore every runner is at distance at least `1/13`, and exactly the two
runners whose multiplied residues are `+1` and `-1` bind.  If they are
`w_+` and `w_-`, then (1) gives

```text
w_++w_-=Q.                                                   (2)
```

The point is genuinely a local maximum.  For sufficiently small positive
`epsilon`, the exact one-sided germs are

```text
phi_W(j/13-epsilon)=1/13-w_+ epsilon,
phi_W(j/13+epsilon)=1/13-w_- epsilon.                       (3)
```

All other runners start at distance at least `2/13`, so they do not alter
this first-order calculation.  Finally, adjoining the runner `Q` kills every
one of the twelve points because

```text
Qj/13 is an integer.                                         (4)
```

Thus the complete shallow local-ruler/deleted-runner pattern follows from
the complementary-lift law (1).  It does not require an arithmetic
progression.

## A primitive compact Cover14 counterexample to the surrogate

Take

```text
W={14,28,29,56,70,84,98,112,126,153,154,168},
V=W union {182}.                                             (5)
```

This is obtained from `14{1,...,12}` by the complementary lift-neutral
replacement

```text
{42,140} -> {29,153}={42-13,140+13}.                         (6)
```

Consequently the residues remain `1,...,12`, and every complementary pair
still sums to `182`.  The replacement destroys both the AP shape and the
common factor.  Directly:

```text
gcd(V)=1,
max(V)/second(V)=182/168=13/12<13,
V contains a multiple of every q=2,...,14.                  (7)
```

One convenient carrier table is

```text
q:       2  3  4  5  6  7  8   9  10  11  12  13  14
carrier:14 84 28 70 84 14 56 126  70 154  84 182  14.
```

For `j=1,...,12`, the active pairs are

```text
{14,168}, {98,84}, {126,56}, {153,29}, {112,70}, {154,28},
{28,154}, {70,112}, {29,153}, {56,126}, {84,98}, {168,14},
```

and every displayed pair sums to `182`.  The extra runner is zero at all
twelve points.

Nevertheless the global values are

```text
M(W)=2/13,                 M(V)=2/17.                       (8)
```

In particular `t=3/119` gives every runner of `V` clearance at least `2/17`;
the binders are `84` and `154`.  The twelve killed shallow points simply miss
this other Farey chart.  Exact enumeration of every self-cusp and every
sum/difference pair intersection proves the upper halves of (8), so these
are exact maxima rather than sampled lower bounds.

## What this removes from the proof menu

The row (5) refutes the implication

```text
primitive + Cover14 + compact
+ full shallow residue packet
+ common active pair-sum ruler
+ deleted runner kills all twelve shallow local maxima
    => AP deletion, tight deletion, or strict 1/13 cover.    (9)
```

It does **not** refute THM-1149's crown-collapse target.  That target starts
with the global premise `M(V)<1/13`; here `M(V)=2/17>1/13`.  Nor does it
touch the open `n=12` statement `M(W)=1/13 => W=d{1,...,12}`; here
`M(W)=2/13`.  THM-1171 remains exactly right: once a tight deletion is
already known to be an AP, its exact AP loneliness forces the homogeneous
ray.  The counterexample shows why neither "AP" nor "tight" can be extracted
from the shallow local spine.

The viable pointwise target must retain the **whole** deletion-safe set

```text
E_Q={t:min_(w in W)||wt||>=1/13}
```

and prove `E_Q` lies in the deleted runner's danger comb.  Twelve canonical
points, even with owner slopes, a common ruler, and the correct `13*14`
carrier, do not determine that containment.  This is precisely the global
confinement missing from the local surrogate.

## Tournament and carrier audit

Taking the twelve shallow points as vertices and comparing their active-pair
sums makes all `66` comparisons ties: every sum is `182`.  Resolving ties by
the order `1,2,...,12` produces the transitive score histogram
`0,1,...,11`, zero directed triangles, twelve singleton SCCs, and one
Hamiltonian path.  It preserves the common ruler and loses the global escape
at `3/119`.  Replacing the vertices by the six complementary owner pairs has
the same failure: all pair observables tie.

The challenged assumption is that canonical witnesses, or their tournament,
are the underlying object.  The smallest faithful carrier at this interface
is instead

```text
the full edge/curtain event atlas
  + every component of E_Q
  + owner and tooth addresses
  + the deleted-runner incidence on those components.       (10)
```

That is the global deletion-containment predicate used by THM-1149.  The
S119 confinement--coupling argument succeeds on APs because AP algebra
compresses (10) without losing it; the complementary-lift row proves that
the compression is not available before AP extraction.

## Exact replay

```text
04-computation/lrc13_local_ruler_deletion_surrogate_referee_codex_20260718.py
SHA-256 2d0e4404fd2c1706b79a95db0bf1b498e95d12e46f0599b54a751c03f46b06fe

05-knowledge/results/lrc13_local_ruler_deletion_surrogate_referee_codex_20260718.out
SHA-256 8fac327be940c20960175c84b47e58a56422f4f67d9e84bc6510a5dce531f96c
```

Normal and optimized runs are byte-identical.  The global `n=12` rigidity
problem, Cover14 crown collapse, and LRC(14) remain open.
