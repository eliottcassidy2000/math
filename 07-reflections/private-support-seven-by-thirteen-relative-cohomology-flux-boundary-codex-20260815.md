# The private-support graph is literally 7 by 13, but its static currents are exact

**Research reflection / finite-exact structural sidecar, not a truth source.**
The support graph is supplied by proved THM-3473.  The exact linear
algebra is frozen by
`04-computation/lrc_private_support_7x13_incidence_h1_probe_20260815.py`.
No LRC(14), bispectrum, or Jacobian conclusion is claimed.

## The unexpected exact carrier

THM-3473's proved Boolean atlas has coactivity packets

```text
A={1,4,5,6},       B={2,3,5,8},       C={5,7}.
```

Their undirected two-section is two `K4`s and one edge, all glued at owner
five.  It has

```text
|V|=8,      |E|=6+6+1=13.
```

Deleting the hub row from an oriented incidence matrix gives a canonical
shape

```text
B: Q^13 -> Q^7.                                      (1)
```

This is a literal `7 by 13` carrier, not numerology imported after the fact.
The exact probe finds

```text
rank(B)=7,
det(B B^T)=256=16*16*1,
nonzero maximal minors=256,
dim ker(B)=13-7=6.                                   (2)
```

Every nonzero maximal minor is `+/-1`.  By Cauchy--Binet, the determinant in
`(2)` is the sum of their squares.  Graphically, the `256` terms are exactly
the spanning trees: `16` choices in the first `K4`, `16` in the second, and
the bridge is forced.

Thus the static support geometry already has an unweighted spectral
nonvanishing theorem.  It does **not** prove the LRC `7 tensor 13` bispectrum
nonvanishing because no LRC phase or current has yet been placed on the
thirteen edges.

## Seven relative owner classes are also canonical

Let `C_q` be the sheet cycle and let

```text
P=P_1 disjoint-union ... disjoint-union P_8
```

be the nonempty private-sheet set, coloured by owner.  The eight colour
indicators span an eight-dimensional subspace of `H^0(P)`.  The diagonal
constant is the restriction of `H^0(C_q)`.  The connecting map for the pair
`(C_q,P)` therefore injects the reduced colour space

```text
Q^8/<(1,...,1)> ~= Q^7                              (3)
```

into `H^1(C_q,P)`.  These are honest relative classes: an owner mask equals
one on its own private colour and zero on the other seven, so its boundary
primitive is not allowed to vanish on `P`.

The same reduced owner-potential space is the source of the graph coboundary

```text
B^T: Q^7 -> Q^13.                                    (4)
```

Equations `(3)--(4)` give a concrete candidate for the first half of a D5
map: relative private-owner classes become thirteen labelled coactivity-edge
currents.

## The exact H1 hostile

The candidate stops one step short of flux.  For the coactivity graph `G`,

```text
H^1(G;Q)=Q^13 / im(B^T),       dim H^1=6.             (5)
```

Every current produced by `(4)` is a coboundary, so it pairs to zero with all
six cycle vectors in `ker(B)`.  Full row rank of the `7 by 13` incidence
matrix proves injectivity of owner potentials; it proves no nonzero absolute
`H^1(G)` class.

This is the precise distinction between three statements that had been easy
to blur:

```text
relative owner class:        nonzero, seven-dimensional;
edge response/coboundary:    nonzero, thirteen coordinates;
absolute graph flux class:   zero without an added holonomy sidecar.
```

The static Boolean atlas therefore supplies a lawful source and an exact
no-go boundary for the D5 program.

This no-go concerns **absolute graph flux only**.  It does not imply that an
edge cochain used as the diagonal weight of a Laplacian has zero determinant.
MISTAKE-409 records that the first version of this reflection blurred those
two predicates.

## Weighted spectral closure reduces to two tetrahedra

Put a weight `w_e` on every coactivity edge and form

```text
L_w=B diag(w_e) B^T.                                  (6)
```

The matrix-tree theorem gives

```text
det(L_w)=w_(5,7) * Tau_A(w) * Tau_B(w),               (7)
```

where each `Tau` is the sixteen-term spanning-tree polynomial of a `K4`.
At unit weights `(7)` is `1*16*16=256`.

This factorization makes the next spectral test cheap and exact.  A proposed
relation-current or theta-slaved phase assignment must provide:

1. a nonzero bridge weight on `(5,7)`;
2. a nonzero weighted tree sum on the first tetrahedron;
3. a nonzero weighted tree sum on the second tetrahedron.

Nonzero individual edge weights are insufficient: complex cyclotomic terms
can still cancel inside either sixteen-term `Tau`.  Conversely, positivity
or a separated character orbit on each tetrahedron would prove weighted
nonvanishing immediately.

Proved THM-3482 now supplies a sharper exact control.  Give each owner
its private-sheet count and orient edges by owner order.  The resulting edge
gradient has zero absolute `H^1` class, but

```text
Tau_A<0,       Tau_B<0,       w_(5,7)=-4k<0
```

for every `k>=1`; hence the weighted determinant is strictly negative.  Four
or six edge weights are actually zero, depending on `k mod 3`.  Therefore
neither nonzero holonomy nor all-edge support is necessary for this
graph-level spectral predicate.  What remains missing for LRC is a lawful
map from a physical relation current to these thirteen labelled weights.

## Connection to incoming current work

The THM-3472 repair taught that a weaker morphism can preserve the target
predicate after an isomorphism fails.  Here the opposite warning applies:
incidence rank is strong enough for a static spectral Gram determinant but
too weak for `H^1` flux.  The missing coordinate is phase/holonomy.
That coordinate is missing **for absolute flux**; THM-3482 shows it is not a
formal prerequisite for weighted determinant nonvanishing.

THM-2334 provides relation-residue currents as character-twisted finite
transforms.  The incoming relation-current transplant and the nonzero
theta-slaved `Q(zeta_91)` contraction are therefore the right candidate
sources for the weights in `(7)`.  A lawful transplant must state which
thirteen coefficients map to the graph edges and why the six cycle pairings
or the two tree polynomials carry the desired physical current.

## Immediate frontier

- Obtain a canonical edge-weight map from the relation-residue current, not
  from graph aesthetics.
- Evaluate the two sixteen-term `K4` tree polynomials exactly in the relevant
  cyclotomic field.
- Use the bridge edge `(5,7)` as the sharp hostile: it is the unique repair
  packet, so any transplant killing it makes the full determinant zero.
- Compare the six graph-cycle holonomies with the seven relative owner
  classes.  Their `7` versus `6` mismatch is information, not an error: one
  owner-potential direction feeds the forced bridge/tree response but no
  independent cycle.

The result is both encouraging and restrictive.  The repo now has a genuine
`7 by 13` combinatorial carrier, an exact absolute-`H^1` no-go for static
owner gradients, and a provisional nonzero weighted determinant on one
canonical gradient family.  Physical LRC closure still requires the missing
relation-current-to-edge map, not non-coboundarity by itself.
