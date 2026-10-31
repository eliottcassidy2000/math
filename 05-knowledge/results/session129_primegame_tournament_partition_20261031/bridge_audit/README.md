# Independent audit: the user's PRIMEGAME partition really contains a weighted triangle

**Status: PROVED elementary statements, with FINITE-EXACT independent checks.**
Scope: the fourteen stated PRIMEGAME fractions, their four stated groups, and
the declared control/data split. This report proves no new tournament census,
prime-generating-program correctness theorem, or LRC consequence. No publishing
was performed by this auditor.

## 1. Main result and strongest useful repair

Write the user's group products as `a,b,c,d`, reserving capital letters below
for individual instructions in Conway's usual program order. Exact arithmetic
gives

```
a = (11/13)(13/11)                         = 1,
b = 55(1/7)(1/17)(17/91)                  = (5/49)(11/13),
c = (19/51)(77/19)(77/29)(29/33)          = (49/9)(11/17),
d = (78/85)(23/38)(95/23)(15/2)           = (45/2)(13/17).
```

Declare the six control primes `11,13,17,19,23,29`, with data primes `2,3,5,7`.
Taking control valuations, the three nonunit products have precisely the arcs

```
17 --d--> 13 --b--> 11,
17 --------c------> 11.
```

This is the **transitive tournament on three vertices**, carrying gains
`45/2,5/49,49/9`. Its orientations come from negative-to-positive control
valuation, not from a cosmetic tie breaker. The unit `a` is a zero-boundary
chain, not an extra tournament arc. Thus there is a precise positive result
inside the proposed partition, even though the originally hoped-for node
counts and interpretation need correction.

Moreover, `b,c,d` are multiplicatively independent, with an **integral decoder**.
For any integers `x,y,z`, put `q=b^x c^y d^z` and

```
T(q) = nu_2(q) + nu_5(q) + nu_23(q).
```

Then

```
x = T(q),
y = nu_11(q) - T(q),
z = nu_13(q) + T(q).
```

Consequently the user's own three block products can faithfully encode any
three-coordinate integral vector, including a tournament face-incidence vector
after its three coordinates have been explicitly associated with `b,c,d`.
Prime factorization recovers the vector exactly. This requires no new bank of
pairwise coprime integers. If another proved result classifies ten tournament
objects by ten distinct primitive vectors on three named coordinates, this
codec transports that classification to ten distinct rational numbers.

**Limit:** the association of the three tournament-node types with the three
blocks remains a declared labeling sidecar. This report does not derive that
association from PRIMEGAME dynamics, and the codec does not explain why there
are ten vectors. It is a valid transport of an independently proved classifier,
not a new classifier or a universal arithmetic/tournament equivalence.

## 2. Exact original graph

For a fraction `f`, write its control part as `p_target/p_source`; use boundary
vertex `0` when no control prime occurs on that side. Every stated fraction has
at most one control prime on each side, with exponent one. Its remaining
positive rational factor is its gain. The graph is therefore uniquely defined
**relative to this declared control/data split**.

| Program letter | Fraction | Control arc | Data gain |
|---|---:|---|---:|
| A | 17/91 | 13 → 17 | 1/7 |
| B | 78/85 | 17 → 13 | 6/5 |
| C | 19/51 | 17 → 19 | 1/3 |
| D | 23/38 | 19 → 23 | 1/2 |
| E | 29/33 | 11 → 29 | 1/3 |
| F | 77/29 | 29 → 11 | 7 |
| G | 95/23 | 23 → 19 | 5 |
| H | 77/19 | 19 → 11 | 7 |
| I | 1/17 | 17 → 0 | 1 |
| J | 11/13 | 13 → 11 | 1 |
| K | 13/11 | 11 → 13 | 1 |
| L | 15/2 | 0 → 0 | 15/2 |
| M | 1/7 | 0 → 0 | 1/7 |
| N | 55 | 0 → 11 | 5 |

Its underlying simple graph, after deleting loops, is `K_(2,3)` on bipartition

```
{11,17} | {0,13,19},
```

with a pendant `29` at `11` and a pendant `23` at `19`. Its shape is a theta
graph with two attached leaves, not a tournament. The directed multigraph also
has opposite pairs on `11—13`, `13—17`, `11—29`, `19—23`, and two loops at `0`.

There are seven vertices and fourteen edges. The graph is connected, so its
integral incidence matrix has rank six and its cycle lattice has rank eight.
The complete ten-prime valuation matrix has rank ten; its multiplicative
relation lattice therefore has rank four. Equivalently, data gains detect four
independent directions in the eight-dimensional control-cycle space.

The unweighted simple graph has just the identity and the interchange `0↔13`
as automorphisms. The directed multigraph, even before gains are included, has
only the identity. These are finite graph checks, not evidence of an arithmetic
identification with independently defined tournament classes.

## 3. An integral chain map, not an ordinary graph quotient

Let `C1=Z^14`, `C0=Z^7`, and let `∂` be the directed incidence matrix above.
The user's blocks in program letters are

```
user A = JK,   user B = NMIA,   user C = CHFE,   user D = BDGL.
```

Let `P:Z^4→C1` send each block basis vector to its instruction-count vector.
Let `J:Z^3→C0` include the vertices `11,13,17`. In this vertex order, define

```
        a  b  c  d
d = [   0  1  1  0
        0 -1  0  1
        0  0 -1 -1 ].
```

Then the exact integral chain identity is

```
∂ P = J d.
```

**Proof.** The control valuation of a product is the sum of the control
valuations of its factors. The four explicitly factored products in §1 give
the displayed columns. The boundary row is minus the sum of all six control
rows, so the equality holds in all seven rows. ∎

Both `P` and `J` are injective; for `P` this follows because the four instruction
sets are disjoint and nonempty. Hence the map preserves every formal cycle of
the macro chain complex, rather than merely preserving its scalar product.
Its integral cycle lattice is exactly

```
ker(d) = Z a ⊕ Z(b-c+d).
```

Indeed the row equations force the `b,d` coefficients equal and the `c`
coefficient their negative; the `a` coefficient is free. Of the original eight
independent cycle directions, this map retains two. The `a` cycle is the actual
two-edge cycle `J+K`. The second cycle is a signed circulation, not a directed
cycle in the transitive tournament.

The zero-boundary generator `a` can be represented by a loop at any macro
vertex for homological bookkeeping, but the chain complex itself does not
select such a base vertex. It must not be presented as a canonical fourth arc
of a tournament.

This is not an ordinary graph contraction. In block B, the trail
`13→17→0→11` passes through `17`, itself a retained macrovertex. In block D,
the arc `17→13`, the two-cycle `19→23→19`, and the loop at `0` form **three
disconnected components**. The block cannot be the edge set of a single walk,
even though its total incidence equals that of the arc `17→13`.

## 4. Gain holonomy is precisely the information missed by the unweighted triangle

The route `17→13→11` has gain `db=225/98`, while the direct arc has gain
`c=49/9`. Therefore

```
h = db/c = 2025/4802 = 3^4 5^2 / (2·7^4) ≠ 1.
```

Control balance does not imply multiplicative balance. If vertex-dependent
positive rational gauges transform an edge gain by
`w(u,v)↦w(u,v)g(u)/g(v)`, this ratio is unchanged, because all vertex factors
cancel around the signed cycle. Two arc gains can be normalized to one, but
all three can be so normalized **iff** `h=1`. Here the triangle is genuinely
weighted, despite having no directed cycle.

The macro valuation matrix has rank three: its control part has rank two,
and `T` evaluates to one on `b` and zero on `a,c,d`, so it detects the remaining
control-balanced direction. Thus the only multiplicative relation among the
four block products is the unit generator `a`.

For the codec of §1, the valuation rows `T,nu_11,nu_13` on the columns `b,c,d`
are

```
[ 1 0 0
  1 1 0
 -1 0 1 ],
```

which has determinant one. Its integral inverse gives the decoder. This proves
injectivity over all integer triples, not merely the finitely checked examples.
It also proves that these three valuation columns generate a primitive
rank-three sublattice of the full ten-prime lattice.

## 5. Minimal no-go for replacing a block by its product

The reduced block-D fraction is `585/34`. It is integral at the **least possible
positive input**, `N=34`, and formally sends `34→585`. Yet no permutation of
its four fractions can execute at `34`, even if the usual program priority rule
is discarded.

**Proof.** Initially there is no factor `19` or `23`. The fractions `78/85` and
`15/2` neither create nor remove these primes. The other two instructions,
`23/38` and `95/23`, each require one of the two primes in their denominator
before they can produce the other. Hence the first of these two instructions
can never execute. This empty-place obstruction is invariant under every
choice of the first two kinds of instruction. ∎

The checker separately tests all 24 instruction permutations. Under the
actual full priority list, the very first step from `34` is instead
`I=1/17`, giving `34→2`.

Positive controls prevent overclaiming that aggregation never corresponds to
a real finite execution:

```
user A at11:  KJ,   11 → 13 → 11;
user B at637: AIMN, 637 → 119 → 7 → 1 → 55;
user C at153: CHEF, 153 → 57 → 231 → 203 → 539.
```

Every arrow here is the first applicable instruction in the full conventional
fourteen-instruction ordering. These establish particular legal macro words,
not universal applicability domains for their reduced product fractions.

The minimal failed implication is therefore:

```
nonnegative final valuations / integral block product
  DOES NOT IMPLY
an executable ordering of the block instructions.
```

Restoration requires the ordered word, intermediate nonnegative valuations,
and, for priority semantics, proofs that every earlier competing fraction is
inapplicable at every step. Flow balance alone omits this information.

## 6. Independent audit of the full nonnegative multiplicative kernel

Use lower-case instruction counts `a,...,n` in program order, only in this
section. The valuation functional `T=nu_2+nu_5+nu_23` changes by zero on the
first thirteen instructions and by one on `N=55`. Therefore any nonnegative
multiplicative relation must have `n=0`. The total number of control-prime
tokens changes only on `I` and `N`, so balance gives `i=n=0`.

The remaining prime-valuation equations give

```
f=e, g=d, h=c, a=b+c, k=c+j,
m=e-b, l=c+e-b, 2b=c+d+e.
```

Solving without division or congruence conditions gives the complete
parametrization by four arbitrary nonnegative integers `c,d,m,j`:

```
a=2c+d+m, b=c+d+m, e=f=c+d+2m,
g=d, h=c, k=c+j, l=c+m, i=n=0.
```

The omitted coordinates are the parameters themselves. Conversely these
formulas satisfy every prime-valuation equation and are nonnegative. Hence the
nonnegative relation monoid is **free on exactly four generators**:

```
JK,
ABDEFG,
A^2 B C E F H K L,
A B E^2 F^2 L M.
```

This agrees with the separate partition worker's result, independently
derived here. It describes commutative multiplicative relations; it does not
assert that any chosen ordering of a generator is executable or that the
priority program follows such an ordering from every suitable state.

## 7. Connection contract and boundaries

| Item | Exact content |
|---|---|
| Source | The user's partition of the fourteen fractions; declared control primes |
| Target | A rank-two incidence complex of a weighted transitive triangle, plus one zero-boundary generator |
| Map | Each block maps to the sum of its instruction valuation columns |
| Preserved | Integral boundary, formal cycles, exact rational block products, gain holonomy |
| Lost | Six source-cycle directions, intermediate control occupancy, ordering, priorities, legal execution domains |
| Restoring sidecar | Full instruction word and prefix/priority certificates; labeling for any external tournament-class codec |
| Decisive positive test | Three factorized products, integral chain identity, determinant-one decoder |
| Decisive hostile | Block D at34 is integral but unschedulable under all24 orderings |

There is no claim that the full PRIMEGAME control graph is a tournament, that
its graph contraction canonically produces the user's order-four merged-node
atlas, or that the group `A` equals one order-three merged tournament node.
The graph-level, gain-level, arithmetic-codec, and execution-level statements
are distinct. The strongest bridge established here is exact and useful
without identifying them.

## 8. Reproduction and evidence

```
python3 check_bridge.py > results.json
python3 -O check_bridge.py > results.optimized.json
cmp results.json results.optimized.json
# Both commands below must fail (exit1), showing optimization cannot drop gates:
python3 check_bridge.py --test-failing-gate
python3 -O check_bridge.py --test-failing-gate
```

The standard-library checker verifies exact rational products and valuations,
both matrix ranks, the chain square, both macro-cycle generators, every
positive-control priority trajectory, all24 hostile permutations, all graph
automorphisms, the four multiplicative generators,625 independent parameter
samples, and729 signed-exponent codec roundtrips. Generality of the two integer
parametrizations is supplied by the proofs above, not inferred from sampling.

`results.json` is the matching output; `SHA256SUMS` covers this report, source,
and both identical normal/optimized outputs. The checker uses always-active
`require` gates rather than optimizable bare assertions. No borrowed tournament
census appears in this checker.
