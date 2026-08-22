# U-spine mixed-order rank break, cover hypergraphs, and the harmonic boundary

**Research reflection / provenance, not a truth source.**  Exact claims are
routed to THM-3416, THM-3455, THM-3461, and the independently audited
THM-3464.

## The mechanism behind the q=123 surprise

The first seven parabolic labels and provisional exact ranks are

```text
q:    11  27  51  83  123 171 227
rho:   6   4   7   9    8   4   9.
```

The tempting square-root pattern at `51` and `83` fails at `123`.  The useful
fact is not merely that the answer is eight.  The witness

```text
(1,40,42,81,82,83,117,122)
```

has quotient orders

```text
123,123,41,41,3,123,41,123.
```

No proper divisor supplies a rank-seven cover, but `q=41` does supply an exact
rank-eight cover:

```text
(3,5,11,19,28,33,37,39).
```

Its threefold pullback covers at `q=123` with active gcd three.  The displayed
mixed-order packet is a second realization of the same rank: it uses one
order-three backbone, three order-41 blocks, and four genuinely order-123
blocks, while retaining active gcd one.  Thus the numerical grade forgets
which primitive quotient layer attains it: both `Q=41` and `Q=123` attain
eight.  The divisor-minimum layers remain disjoint; coexistence is equality of
their minima, not a failure of independent layer pricing.

This suggests a concrete search heuristic for later composite U-spine labels:
factor `q`, enumerate quotient-order profiles before owner residues, and price
mixed profiles by their possible sheet-density and residue intersections.
For the next unresolved composite label

```text
q_8=291=3*97,
```

the first profile to test is the direct analogue

```text
order 3 backbone + several order 97 blocks + primitive order 291 repairs.
```

That is a structural hypothesis and not yet a rank claim.

## Why a tournament is the wrong compression

At fixed `q`, the exact object has:

- sheet vertices `Z/qZ`;
- owner-labelled danger hyperedges;
- a two-valued common-centre layer;
- a cover predicate on the full union; and
- mode widths, centres, and endpoint scales as sidecars.

There is no intrinsic orientation between two owners.  One can manufacture an
arc from asymmetric private-sheet counts, but ties, missing arcs, and reversed
choices change under harmless relabelling and still do not determine whether
eight edges cover.  At `q=123` the occupancy profile is

```text
1^106 2^6 4^11.
```

The eleven fourfold sheets are an explicit witness that pairwise comparison
data is not the full predicate.  The correct finite carrier is the labelled
incidence hypergraph.  A tournament can be a diagnostic shadow only after its
orientation gauge and lost higher intersections are declared.

## Ternary words help only with the right labels attached

For one chosen family, each sheet has the ternary occupancy symbol

```text
0 = uncovered,  1 = private,  2 = multiply covered.
```

This word immediately certifies coverage and exposes redundancy.  It does not
reconstruct which owners meet at a multiply covered sheet, and it collapses
all multiplicities two and above.  A faithful compiler therefore stores

```text
binary owner-by-sheet incidence matrix
  + ternary occupancy projection
  + owner/mode/centre sidecars.
```

This is the same loss boundary isolated by the Rule-30 ternary-intersection
work: low-order overlap codes are useful observables, not substitutes for the
labelled source object or its current.

## Natural subsets and harmonic sampling

Every property of the U-spine indices defines a subset of the natural numbers.
For example,

```text
S_7={t:rho_ZMC(q_t)<=7}.
```

THM-3455 proves `S_7` periodic and therefore gives both its natural density
and the coefficient of

```text
sum_(t in S_7, t<=N) 1/t.
```

The new exact values at `t=5,7` refine the complementary high-rank set but do
not change `S_7`.  A finite rank prefix cannot determine a harmonic
coefficient for exact rank eight or nine.  To obtain one, the next theorem
would need a finite automaton, a genuine period, or an analytic distribution
law for the higher-rank mask banks.

The sequence of **labels** has a different harmonic question:

```text
sum_t 1/q_t
```

converges because `q_t` is quadratic.  It must not be confused with harmonic
sampling by the index `t`, where positive-density subsets diverge
logarithmically.

## Connection to the current transplant

All fixed-half-twist packets lie on one THM-3410 projective ray, so their
pairwise wedges vanish.  The new relation-residue transplant shows that this
does not preclude endpoint activity after owner values are lifted while mask
residues and the physical clock are preserved.  The map is therefore

```text
mask hypergraph at q
  -> CRT-decorated owner tuple
  -> delayed endpoint word
  -> grouped target coefficient.
```

THM-3464 advances only the first object.  The live obstruction on the last
arrow is cancellation in a grouped coefficient, not cover rank.  Mixing those
two predicates would erase precisely the progress obtained on both fronts.
