# The selector cut is a parity carrier, not a two-state quotient

**Status: FINITE-EXACT PARTIAL SCOUT / NON-CANONICAL SYNTHESIS.**  The exact
companion is
[`gmc_selector_paired_half_transfer_partial_scout_20260803.py`](../04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py),
with frozen
[`output`](../05-knowledge/results/gmc_selector_paired_half_transfer_partial_scout_20260803.out).
It rebuilds the three static witness-decorated graphs of
[THM-3288](../01-canon/theorems/THM-3288-maximizing-witness-lifted-walk-rational-series.md)
from the pinned relation artifact of
[THM-3287](../01-canon/theorems/THM-3287-weighted-backbone-dominance-witness-section-and-selector-cut.md).
It introduces no theorem ID and changes no claim about chronology, factorial
moments, tournaments, or any open problem.

## Inheritance and question

THM-3287 gives an exact selector-origin sign cut on the twelve response rows:

```text
small/positive = {2,11,16,18,22},
full/negative  = {3,7,10,13,17,19,21}.
```

Every core row edge crosses this cut.  THM-3288 consequently notes that its
active `(row,static witness)` graphs are bipartite and obtains recurrence
polynomials involving only even powers of `z`, except for the full core's
extra factor `z`.  The new questions were:

1. does the selector bit itself define a lawful two-cell quotient;
2. what is the coarsest equitable refinement retaining the all-ones walk
   observer;
3. can the exact orders `10/14/15` be reconstructed from the two half-step
   transfers; and
4. what vector, rather than only what scalar, causes the full-core residual
   `-1392`?

The answers are respectively **no**, `11/14/16` cells, **yes**, and one
explicit selector-small zero mode.

## Exact block structure and the first hostile

Order the active vertices with selector-small rows first and selector-full
rows second.  For all three families the adjacency matrix is literally

```text
          [ 0   X ]
A       = [       ].                                      (1)
          [ X^T 0 ]
```

The exact dimensions are:

| family | active vertices | directed arrows | selector split | `X` |
|---|---:|---:|---:|---:|
| backbone | 18 | 32 | `13+5` | `13 x 5` |
| selected tree | 23 | 40 | `14+9` | `14 x 9` |
| full core | 26 | 68 | `14+12` | `14 x 12` |

Thus the selector cut is a genuine parity carrier.  It is not an equitable
two-cell lumping.  Already the degree histograms within its two sides are:

```text
backbone      small ((1,10),(2,3))
              full  ((1,1),(2,1),(4,2),(5,1))
selected tree small ((1,10),(2,3),(4,1))
              full  ((1,5),(2,1),(4,2),(5,1))
full core     small ((1,5),(2,5),(3,1),(4,1),(5,1),(7,1))
              full  ((1,5),(2,3),(5,1),(6,3)).          (2)
```

For example, in the full core `(11,Q+2)` has degree one while `(22,Q+4)`
has degree seven, although both have the same selector bit.  This is the
cheapest exact hostile to a binary state compression.

## Equitable versus observer quotients

Starting with the two selector cells and repeatedly splitting vertices by
their complete neighbour-count vectors gives the coarsest equitable
refinement of that supplied partition.  The exact result is:

| family | equitable cells | all-ones Krylov order | observer-unseen equitable modes |
|---|---:|---:|---:|
| backbone | 11 | 10 | 1 |
| selected tree | 14 | 14 | 0 |
| full core | 16 | 15 | 1 |

The companion prints every cell and every quotient-matrix row.  Their typed
cell/quotient digests are

```text
backbone      de4dc6a0bc1e204036c33ac99e5590eadc719141eea143c503df95b33d11abdb
selected tree 137cd64446a6985deaac8e8e070bb94ad29b7ca4474633b131f2813783c2861d
full core     cc2f113e46d1137e28eaf1753514a68d10bdc1b70849ec9b476c72b70a7ba375.
```

If `P` is the cell-indicator matrix and `B` the printed quotient, the scout
checks exactly

```text
A P = P B,                                             (3)
a_n = cell_sizes^T B^n 1                              (4)
```

for all tested powers, together with the vector Krylov relation that proves
the continuation for every power.  The backbone quotient has characteristic
polynomial `z q_B(z)` but its all-ones observer is orthogonal to that zero
mode.  The full-core quotient has `z^2 q_C(z)`, while the observer retains
exactly one of the two quotient zero modes.  This explains why an equitable
quotient can be slightly larger than the minimal scalar observer without
contradicting THM-3288's Hankel minimality.

## The selector-paired half transfer

Squaring `(1)` gives

```text
A^2 = diag(XX^T,X^TX).                                (5)
```

Let `M=A^2`.  Then

```text
a_(2m)   = 1^T M^m 1,
a_(2m+1) = 1^T M^m A1.                                (6)
```

Exact vector Krylov closure, not recurrence fitting, gives:

| family | even-half order | odd-half order |
|---|---:|---:|
| backbone | 5 | 5 |
| selected tree | 7 | 7 |
| full core | 8 | 7 |

For the first two families the even and odd halves have the same minimal
polynomial `p(t)`, and the full observer polynomial is literally

```text
m(z)=p(z^2).                                           (7)
```

For the full core, the odd half has

```text
p_C(t)=(t-4)
 (t^6-30t^5+284t^4-1149t^3+2105t^2-1676t+452),       (8)
```

while the even half has minimal polynomial `t p_C(t)`.  Hence the interlaced
sequence has minimal observer polynomial

```text
m_C(z)=z p_C(z^2),                                     (9)
```

of degree fifteen; THM-3288's degree-fourteen tail factor is
`q_C(z)=p_C(z^2)`.  Equations `(5)--(9)` explain the `z^2` denominator
structure and the orders `10/14/15` directly from the selector-paired
transfer, while retaining separate start and output vectors.

For a fixed family, the two smaller companion matrices may be powered in
logarithmically many exact ring-arithmetic steps.  This is an arithmetic-
operation observation only; no coefficient bit-complexity bound is claimed.

## The full-core head is an exact zero projection

Put

```text
h=p_C(A^2)1.                                           (10)
```

The exact vector is nonzero, lies entirely on the selector-small side, and
satisfies

```text
A h=0,
1^T h=-1392,
h^T h=2516736,
p_C(0)=-1808.                                         (11)
```

Its twelve nonzero coordinates are

```text
(2,Q+4)   -384
(11,Q+2)   192     (11,Q+4)  -384     (11,Q+5)   192
(18,Q+1)  1008     (18,Q+2)  -336
(18,Q+4)  -336     (18,Q+5)  -336
(22,Q+1)  -464     (22,Q+2)  -464
(22,Q+4)   384     (22,Q+5)  -464.                    (12)
```

Moreover

```text
pi_0(1)=h/p_C(0)
```

is exactly the orthogonal projection of the all-ones vector onto the kernel
of `A` inside the full graph, with spectral mass

```text
1^T pi_0(1)=87/113.                                   (13)
```

The degree-seven even recurrence therefore fails at its first possible index
`m=7` by

```text
1^T p_C(A^2)1=-1392,                                  (14)
```

and succeeds at every later index because `A^2h=0`.  In the original
interlaced sequence this is exactly the `n=14` residual and the coefficient
of `x^14` in THM-3288's generating numerator.  The trailing zero in the
order-fifteen recurrence is thus the scalar shadow of the concrete vector
`h`, not a formatting artifact.

## Operation boundary and hostile paths

The active node sets are strictly nested.  Passing from the backbone to the
selected tree creates five vertices,

```text
(3,Q-8),(13,Q-1),(16,Q+1),(17,Q-1),(17,Q-8),
```

and passing to the full core creates three more,

```text
(10,Q-8),(13,Q-8),(21,Q-8).                            (15)
```

Consequently an edge-addition formula must retain the changing active-start
vector; padding by isolated vertices and silently keeping the old all-ones
observer changes `a_0`.

The canonical path `2-17-16` remains a direct structural hostile.  Its first
edge requires middle fibre `Q-1`, its second requires `Q-8`, and the fibres
are disjoint.  Therefore the decorated graph has no length-two walk
projecting to that row path.  Matrix powers count exactly the surviving
static witness walks; they do not repair the missing section.

## Connection contract and scope

```text
source:
  the three active THM-3287 maximizing-witness relation graphs;

target:
  selector-paired half transfers, equitable refinements, and all-ones
  observer/Krylov companions;

map:
  row selector -> bipartition -> A^2 half transfer, together with stable
  neighbour-count refinement and the all-ones start/output vectors;

preserved:
  every static witness-decorated walk count and its exact rational series;

destroyed:
  individual walk identity, row-path simplicity, response magnitude,
  physical one-pole chronology, and continuation after reset;

needed sidecars:
  the active-node indicator, both parity boundary vectors, and for the full
  core the head vector h;

hostiles:
  selector-only degree heterogeneity, active-node creation, and 2-17-16.
```

This is not a tournament: the arrows are symmetric static compatibility
arrows, and no tie-breaking orientation is introduced.  It is not response
composition, a positive current, a factorial-moment argument, or an
`LRC(14)` mechanism.

## Reproduction

Run

```text
python3 04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py
python3 -O 04-computation/gmc_selector_paired_half_transfer_partial_scout_20260803.py
```

and compare LF-normalized bytes with the frozen output.  Current hashes are

```text
script da47f9db80558343bd6442196c65150417241acbed65385b876ca09a66f0490f
output ef276b6a99a3ec1ac59d5cc352a3aa9a50fa64440a6615695d684a68a6683208.
```

The source pins the promoted THM-3287 and THM-3288 primary artifacts, replays
the relation source against its stored transcript, rebuilds all three graphs,
prints the complete equitable cells and quotient matrices, and has no Python
assert node or floating literal.
