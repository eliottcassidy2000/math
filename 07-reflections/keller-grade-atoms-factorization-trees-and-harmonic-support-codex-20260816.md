# The Keller degree functor preserves products and forgets atoms

**Status: proved synthesis from THM-3438 and THM-3541.**  This note replaces
the proposed “degree four is Keller's seven” analogy by the exact
factorization theory of the degree grading.  It classifies no individual map
up to equivalence and gives no `JC(2)` or LRC consequence.

## 1. Inheritance pass

The closest proved mechanism is THM-3438:

```text
KDeg(r)={1} union {n>=3},
AtomDeg(r)={n>=3},
DecompDeg(r)={uv:u,v>=3}.                              (1)
```

The canonical hostile to grade-only reasoning is already at degree nine:
the weighted `S_9` map is a composition atom, while `F_3 o F_3` is a
composite of the same degree.  The corrected near miss is therefore not a
missing construction.  It is the assumption that irreducibility of an
integer in the grade monoid is equivalent to composition-irreducibility of a
map.

The least-used sidecar is the complete factorization word of the integer
grade.  THM-3541 computes it exactly.

## 2. Two atom notions and one lossy functor

Let `K_r` be the composition monoid of Keller maps in dimension `r>=3`, up
to whatever equivalence is chosen only after composition is defined.  Degree
is a multiplicative map

```text
deg:K_r -> M={1} union {n>=3}.                         (2)
```

It has two asymmetric properties.

1. If `deg(F)` is an atom of `M`, then `F` is forced atomic.
2. If `deg(F)` is reducible in `M`, `F` may still be atomic.

Thus degree reflects forced atomicity but does not preserve or detect map
atoms.  The restoring object is the intermediate-field/block lattice.  An
`S_n` cover is primitive even when `n` factors; a composition carries an
ancestry block system even when its total degree is the same.

The atoms of `M` are

```text
odd primes, 4, 8, and twice an odd prime.              (3)
```

The exceptional composites `4,8,2p` are exactly the crevasses cut by the
missing grade two.  Degree four is not an isolated sporadic value; every
grade from three through eight is forced atomic.

## 3. Factorization trees have an exact depth interval

Write `N=2^a m` with `m` odd and `b=Omega(m)`.  A leaf decomposition of the
numerical grade into atoms `(3)` has every possible leaf count in

```text
[b+ceil(max(a-b,0)/3), b+floor(a/2)].                  (4)
```

The endpoints have a transparent tree meaning.

- Pair as many powers of two as possible with odd prime leaves (`2p`), then
  package the remaining twos mostly in blocks of three (`8`): this minimizes
  leaves.
- Split powers of two mostly into pairs (`4`), using one `2p` or one `8` to
  absorb odd parity: this maximizes leaves.

Every intermediate leaf count occurs.  The global ratio is at most `3/2`,
and `36` attains it:

```text
36=6*6=4*3*3.                                          (5)
```

The first nonunique tree already occurs at

```text
24=3*8=4*6,                                            (6)
```

but both presentations have two leaves.  Equation `(5)` is the first time
the leaf count itself changes.

These are numerical factorization trees.  They do not prove that one
polynomial map has two decomposition trees.  Conversely, THM-3438 lets one
realize every expression depth from one to the upper endpoint of `(4)` by
choosing atomic weighted maps at grouped leaf degrees and composing them.
This varies the map inside the grade.

## 4. Why the ternary tree is one fibre, not the whole family

The fixed sporadic cubic tower has a literal ternary ancestry tree because
every factor degree is three.  Its depth-`k` grades are `3^k`, a single ray
inside `(2)`.  THM-3541 shows what changes outside that ray:

- a factor `4`, `8`, or `2p` is one numerical leaf, not a forced cluster of
  binary leaves;
- mixed grades admit different arities and different leaf counts; and
- an atomic `S_N` map collapses the entire numerical factor tree to one
  primitive cover.

So the family is indeed degree-graded and multiplicative, but it is not one
ternary tree.  Ternary ancestry is retained data of the fixed cubic
submonoid.  The global grade monoid is a forest of multiplicative
factorizations together with primitive `S_N` points in every grade.

## 5. The exact strong-atom-floor survivor

The H-spectrum strategy uses a lower floor on allowable strong atoms to
exclude a candidate integer from an additive or tiling semigroup.  Its
Keller analogue is only the necessary statement

```text
F=G o H nontrivially
=> deg(F)=deg(G)deg(H), with both factors >=3.           (7)
```

Equation `(7)` explains why grades `3` through `8` are forced atomic and why
`9` is the first mixed grade.  It does not exclude an `S_9` atom.  The
analogy survives as a **floor on block sizes**, not as a classification by
the total integer.

Accordingly, degree four is the first ordinary composite integer whose only
factorization `2*2` is illegal in the Keller grade monoid.  It is one member
of the infinite crevasse family `(3)`, not a counterpart of a uniquely
forbidden H-spectrum value.

## 6. Harmonic support is thick on both sides

The grade-atom support is the disjoint Boolean union

```text
A=P_odd disjoint-union 2P_odd disjoint-union {4,8}.     (8)
```

For `Re(s)>1`, its Dirichlet series is exactly

```text
sum_(n in A)n^(-s)
 =(1+2^(-s)) P_odd(s)+4^(-s)+8^(-s),                  (9)
```

where `P_odd` is the odd-prime zeta sum.  At the harmonic boundary `s=1`,
Euler's prime-reciprocal divergence makes `(8)` harmonic-thick.  The
reducible grades are also harmonic-thick because they contain `9N`.

This supplies a sharp example for the user's subset principle: a set can be
arithmetically sparse yet carry divergent reciprocal mass, and both sides of
a factorization dichotomy can diverge.  Harmonic convergence therefore
cannot replace the factorization word or block lattice.

## 7. Tournament and XOR boundary

The alternatives in `(6)` do not orient one another.  They are two
factorization witnesses with no intrinsic pairwise comparison.  Likewise,
the four numerical packet types in other lanes do not become `V4` unless a
Boolean operation is supplied.

One may build a tournament on proof obligations—shorter word, smaller
largest leaf, simpler monodromy—but that gauge is external.  It does not
turn numerical factorizations into edges of the map-composition monoid.  The
lawful relation is a divisibility/refinement poset; its sidecar is the actual
map and intermediate-field lattice.

## Connection contract

| field | exact answer |
|---|---|
| source | composition monoid of Keller maps |
| target | numerical monoid `M={1} union {n>=3}` |
| map | generic field degree |
| preserved | units and multiplication; map decomposition implies grade factorization |
| destroyed | map atom status at reducible grades, formulas, block lattice, Jelonek data |
| restoring sidecar | monodromy group and intermediate-field/block lattice |
| first mixed grade | `9`: atomic `S9` and composite `3*3` both exist |
| first numerical nonuniqueness | `24=3*8=4*6` |
| first length nonuniqueness | `36=6*6=4*3*3` |

The exact companion and theorem carry the full proof.  Nothing here
classifies tame conjugacy classes, fixed-map decomposition uniqueness,
`JC(2)`, a physical current, or LRC(14).
