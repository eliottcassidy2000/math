# LRC(14): a concave cubic OR majorant closes the fixed-pool rank-four census

**Status: session reflection.  Mathematical truth is the frozen exact packet;
LRC(14) remains open.**

## Inheritance and portfolio

- **Closest proved mechanism:** THM-4333 turns low-rank pair-safe cells into a
  weighted failure hypergraph and optimizes its covered mass by decreasing
  marginals.  THM-4336 proves the desired rank-four floor on `(50,70)` and
  `(509,640)` by exact quartic-cell coverage.
- **Hostile:** `(50,70)` defeats the raw rank-four degree screen, so another
  pass through first moments cannot be uniform.
- **Corrected near miss:** extrapolating the two THM-4336 controls was barred;
  the other `181,192` pairs needed their own exact certificates.
- **Least-used sidecar:** the intersection count `t` of one rank-at-most-four
  cell with a partial body.  Keeping `t`, rather than only its hit/no-hit bit,
  exposes a low-degree concave envelope.
- **Anchor:** all residual pairs at rank four.
- **Niche:** concave polynomial majorants of the Boolean OR cost.
- **Wildcard:** treat truncation order as an interpolation parameter and ask
  how much incidence order can be discarded while retaining submodularity.

The live concept board was: literal wall cells; failure hyperedges; cardinality
eight coverage; low-degree incidence tensors; concave discrete envelopes; and
the one-appender discrepancy transfer.

## The move

On a rank-at-most-four failure cell, exact coverage is `1_(t>0)` for
`t in {0,1,2,3,4}`.  The cubic

```text
g(t)/14 = t-(9/14)binom(t,2)+(3/14)binom(t,3)
```

majorizes that OR cost.  Its values are `(0,1,19/14,9/7,1)` and its
increments are `(1,5/14,-1/14,-2/7)`.  Two features matter simultaneously:

1. it returns exactly to one at `t=4`, preventing the four-hit cells from
   creating uncontrolled slack; and
2. its increments decrease, so the induced body objective remains
   submodular even though the last two marginals are negative.

Thus the union coverage is bounded using degrees, codegrees, and triple
degrees only.  A quartic tensor is unnecessary.  At every partial body, the
sum of the largest required current marginals upper-bounds every completion;
an exact branch-and-bound is therefore cheap enough across the full pair
universe.

This is the real reason the computation scales.  The raw degree bound passes
only `29,595` pairs.  The overlap-aware cubic passes the other `151,599`, and
an exact-all replay passes all `181,194`.  The previous hostile `(50,70)` is
also the unique weakest normalized cubic row, so the audit lands on the known
hard geometry rather than an accidental new corner.

## What changed

The stated rank-four frontier is gone:

```text
for every Q in the THM-4231 remainder and every K in binom(P,8),
L_4(K;Q)/D>7/81.
```

The finite statement yields a cubic-only uniform one-appender cutoff `13,737`.
Inheriting THM-4336's sharper exact controls lowers the uniform cutoff to
`12,274`, uniquely at `(50,140)`.  This replaces THM-4333's rank-three
`3,370,132,808` cutoff on the same frozen pair universe.

The computation is not the LRC(14) entry theorem.  It proves safety for a
large structured fixed-pool family after a sufficiently large appender.  It
does not map an arbitrary thirteen-speed row into that family.

## Connection contract

```text
source:       rank-at-most-four weighted failure hypergraph
target:       cubic degree/codegree/triple-incidence objective
map:          replace OR(t) by g(t)/14
preserved:    a uniform lower bound on every cardinality-eight survivor mass;
              decreasing-marginal branch-and-bound
destroyed:    exact union weight and the identity of omitted rank>=5 cells
sidecar:      optimizer body plus direct L4 evaluation; literal-wall controls
test:         flat-enumerate all C(30,8) bodies on stratified hostile pairs
```

Seven literal-wall flat enumerations, including both previous controls, the
first and last universe rows, two interior indices, and the former cutoff
hostile `(721,746)`, reproduce the exact optimizer rewards and least masks.

## Reusable method and next barrier

The reusable move is a **concave OR-envelope compiler**:

1. freeze the maximum hyperedge rank `m`;
2. find a low-degree polynomial in binomial coordinates that majorizes
   `1_(t>0)` on `t=0,...,m`;
3. require decreasing discrete increments so cardinality branch-and-bound has
   sound current-marginal bounds;
4. optimize using only the corresponding low-order incidence tensors;
5. keep direct cell evaluation as a sidecar because the envelope loses exact
   coverage.

For rank four, the cubic envelope is enough.  At higher ranks, the cheapest
hostile probe is a small linear program for the lowest-degree concave envelope,
followed by a comparison between envelope slack and tensor cost.  The method
is contraindicated when the target needs the exact minimizing body or when
envelope slack is of the same size as the desired margin.

For LRC(14), the honest barrier has moved back to entry: either connect
arbitrary `2+12` rows to the fixed pool without losing owner/address data, or
use the now-complete rank-four chart as a terminal for a new reduction.  More
fixed-pool optimization is no longer the leading need.
