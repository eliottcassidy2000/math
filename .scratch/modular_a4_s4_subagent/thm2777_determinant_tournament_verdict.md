# THM-2777 determinant-tournament verdict

## Verdict

The proposed object is **not cosmetic while the `A2` chamber is retained**.
The marked tournament is rigid and chiral, so global reversal is a genuine
observable.  A single unweighted tournament is nevertheless chamber-gauge
dependent.  Its marked preserved content is:

1. the quartic sign character, realized as global reversal under odd
   `W(D3)` elements;
2. the absolute determinant weights, which recover the `A2` triple and the
   three opposition pairs.

There is a sharper quotient loss.  After chamber reorientation is forgotten,
the **weighted switching class is self-converse**: switching the `a12`
vertex and relabelling `a23<->a13`, `s13<->s23` identifies it with its global
reverse while preserving every determinant weight.  Therefore quartic sign
does *not* survive the chamber-free quotient.  What remains there is the
weighted `A2`/opposition incidence and its projective oriented-matroid class.

The construction supplies no new affine Kummer, Jelonek, or LRC coordinate.
In particular, the orientation character cannot distinguish the `000` and
`110` even boundary words.  The fully labelled six-line action does
distinguish the identity from the `110` face square, but that is the already
faithful six-edge action of THM-2753; the tournament does not say whether the
nontrivial element occurs as divisor inertia.

Promotion is mathematically reasonable if the theorem is scoped this way.
It should not be promoted as a canonical tournament on six unmarked roots.

## 1. Exact marked construction

Take

```text
h=(1,1,1),
a12=e1-e2,  a23=e2-e3,  a13=e1-e3,
s12=e1+e2,  s13=e1+e3,  s23=e2+e3.
```

The `a` roots are the positive `A2=h^perp` roots for one chamber; the `s`
roots are oriented by positive dot product with `h`.  For distinct roots
`alpha,beta`, orient

```text
alpha -> beta  iff  det(alpha,beta,h)>0,
w(alpha,beta)=|det(alpha,beta,h)|.
```

The complete determinant table, with rows interpreted in the displayed
vertex order, is

| pair | det | pair | det | pair | det |
|---|---:|---|---:|---|---:|
| `a12,a23` | 3 | `a12,a13` | 3 | `a23,a13` | -3 |
| `a12,s12` | 2 | `a23,s23` | 2 | `a13,s13` | -2 |
| `a12,s13` | -1 | `a12,s23` | -1 | `a23,s12` | -1 |
| `a23,s13` | -1 | `a13,s12` | 1 | `a13,s23` | 1 |
| `s12,s13` | -1 | `s12,s23` | 1 | `s13,s23` | -1 |

Thus there are no ties and

```text
weight spectrum:  1^9, 2^3, 3^3.
```

The weight-three graph is the complete graph on the three `A2` roots.  Each
of its three edges is a pair of positive `A2` roots spanning the subsystem.
The weight-two graph is exactly

```text
{a12,s12}, {a23,s23}, {a13,s13},
```

the three same-support orthogonal/opposition pairs.  All other pairs have
weight one.

The resulting tournament has score multiset

```text
(1,2,2,3,3,4),
```

is strongly connected, contains six cyclic triples, and has 23 labelled
Hamiltonian paths.  It has trivial automorphism group and is not isomorphic
to its converse.  The determinant table itself is the unambiguous
isomorphism-class certificate.

## 2. Gauge and covariance

Changing the `A2` chamber changes signs on some of the three `a` vertices.
Negating one representative reverses every incident arc, i.e. performs
tournament vertex switching.  The six chambers give six distinct marked
tournaments in one switching class.  Moreover the switch at `a12`, followed
by

```text
a23 <-> a13,       s13 <-> s23,
```

takes the tournament to its converse and preserves every weight.  Therefore:

```text
marked h + marked chamber + ambient orientation -> one tournament;
forget chamber                                -> switching class.
```

For every `g in W(D3)`,

```text
det(g alpha,g beta,g h)=det(g)det(alpha,beta,h).
```

Hence even elements transport the weighted tournament orientation
faithfully, while odd elements globally reverse every arc; all absolute
weights are preserved.  With the chamber retained, this realizes the
quartic sign character on a genuine pairwise observable.  After quotienting
by chamber switching, self-converseness erases that character.

## 3. Binary/ternary spectrum and sharp loss

For THM-2775's generators, the induced permutations on the six root lines
are

```text
X_S: (a12 s13)(a13 s12),          type 2^2 1^2,
Y_S: (a12 a13 a23)(s12 s13 s23), type 3^2,
D=(X_SY_S)^2:
     (a23 s23)(a13 s13),          type 2^2 1^2.
```

Thus `X_S` and the nonzero `V4` face square collide on unlabelled six-line
cycle type exactly as in THM-2753.  The determinant orientation separates
their characters:

```text
det(X_S)=-1  -> global reversal,
det(Y_S)=+1  -> preservation,
det(D)=+1    -> preservation.
```

This restores quartic sign, but not affine boundary placement.  Both the
identity (`000`) and `D` (`110`) have determinant `+1`.  If the labelled
line permutation is retained they are different, but the static tournament
does not provide the missing statement that `D` is or is not inertia at a
particular divisor.

## 4. Required promotion loss ledger

```text
source:
  oriented Euclidean D3 root system, marked even body diagonal h,
  marked A2 chamber.

target:
  six-vertex tournament with edge weights 1,2,3.

map:
  sign and magnitude of det(alpha,beta,h).

preserved:
  marked oriented-matroid chirotope; quartic sign as global reversal
  while the chamber is retained;
  A2 triple; same-support opposition matching; binary/ternary line spectra.

forgotten:
  after chamber quotient, quartic sign/global reversal;
  absolute quartic roots and common translation; Kummer valuation row;
  divisor at which a V4 element occurs; units and Cl[2]; polynomial
  realization; endpoint/owner/factor allocation.

sidecars needed:
  for JC, THM-2685 normalized boundary rows plus unit/class-group data;
  for LRC, a lawful action on THM-2763's address carrier preserving W,
  L_13, carrier allocation, and endpoint current.

cheapest hostile:
  THM-2769 realizes 110 while preserving the same finite W(D3) frame.
```

So THM-2777 can canonize a clean finite representation theorem, but it does
not advance `JC(2)`, `DC(2)`, or `LRC(14)` without a new realization map.
