# The relation-role transplant lands in `B^1`, not `H^1`

**Status: research synthesis around a FINITE-EXACT unnumbered sidecar; not a
truth source and not an LRC(14) result.**  Reproduce the exact claims with
`04-computation/lrc_relation_role_chart_weighted_closure_probe_20260815.py`.

## 1. Inheritance pass

The closest proved mechanism is THM-3482: an owner-potential gradient can
have zero absolute graph cohomology and nevertheless have a nonzero signed
weighted tree determinant.  The canonical hostile is the flat potential,
whose edge gradient and determinant both vanish.  The corrected near miss is
MISTAKE-409: nonzero `H^1` flux is not necessary for this spectral predicate.

The least-used sidecar is the exact relation now proved in audited THM-3479,

```text
a=(-27,-27,-27,20110798,-41,-27,-27,-27,38)          (1)
```

in order `(c1,c2,c3,H,q1,q2,q3,q4,q5)`.  It contains a role packet that was
not visible when one tried to biject a native thirteen-character orbit with
the thirteen graph edges.

## 2. The six-slot packet matches the two `K4` wings

Exactly six coefficients in (1) equal `-27`:

```text
{c1,c2,c3} union {q2,q3,q4}.                          (2)
```

The proved THM-3473 private-support carrier has exactly six outer vertices,
split into the two triples of its two `K4` blocks.  It also has a unique
degree-seven hub `u5` and a unique degree-one repair leaf `u7`.

This gives the explicit role contract

```text
delete q1 as the relation gauge;
H  -> u5 (hub);
q5 -> u7 (repair leaf);
{c1,c2,c3} and {q2,q3,q4} -> the two outer triples.    (3)
```

There are

```text
2*(3!)*(3!)=72                                       (4)
```

labelled maps respecting (3).  Exact enumeration gives `|Aut(G)|=72`, and
the automorphism orbit of one role map is all 72 maps.  Thus (3) determines
one chart **up to target graph automorphism**.  This is much smaller than the
`8!/72=560` arbitrary graph-labelling classes, but it uses the extra semantic
choices in (3); it is not forced by the native `C13` action.  THM-3479's
no-`C13`-equivariant-bijection obstruction remains intact.

## 3. The explicit D5 linear map

Let

```text
Lambda_a=ker(a:Q^9 -> Q).
```

Because `a_q1=-41`, deleting `q1` is an isomorphism

```text
pi_q1: Lambda_a -> Q^8.                               (5)
```

Choose any representative of the unique role chart (3), use it to identify
the eight coordinates with vertex potentials, quotient by constant
potentials, and apply the oriented graph coboundary.  The resulting map is

```text
D_role:
Lambda_a / Q*s  --~-->  Q^8/Q*1  --d-->  B^1(G;Q)
                                      subset C^1(G;Q), (6)
```

where `s` is the unique relation vector whose eight retained coordinates are
all one.  Its deleted coordinate is

```text
s_q1=20110674/41.                                     (7)
```

The map in (6) has rank seven and is an isomorphism onto the seven-dimensional
coboundary space.  Since the graph has eight vertices and thirteen edges,

```text
dim H^1(G;Q)=13-8+1=6,                                (8)
```

but every value of `D_role` has zero pairing with an exact six-cycle basis.
Consequently

```text
[D_role(U)]=0 in H^1(G;Q) for every U in Lambda_a.     (9)
```

This is the explicit map the earlier D5 discussion was reaching for, with a
crucial correction: its natural target is `B^1`, not a nonzero class in
`H^1`.

There is also a sharp integral boundary.  Projection (5) identifies the
integer relation lattice with an index-41 congruence sublattice of `Z^8`, not
all of `Z^8`; (7) is the cheapest witness.  Hence (6) is a rational
isomorphism, while an integral formulation must retain the mod-41 sidecar.

## 4. Weighted closure survives every role chart

Feed the two exact THM-3479 tuples into (6) as relation-coordinate
potentials, orient every edge by owner order, and form the signed weighted
reduced Laplacian.  For the canonical displayed chart, the determinant
factors through the forced bridge:

```text
U_full:
  T_A=-408230882586462649,
  T_B=140608,
  bridge=-312,
  det=17908964716879810126984704;

U_clock:
  T_A=-408094818615760385,
  T_B=-282469436237523288,
  bridge=-286,
  det=-32968453616912573701504956325888921680.         (10)
```

Here `T_A,T_B` are the two signed `K4` tree sums.  More robustly,

```text
                         U_full       U_clock
nonzero role charts      72/72        72/72
distinct determinants    30           36
v_13(det)                 4            4   (every chart). (11)
```

The exact valuation has a transparent source.  The hub and the three
coordinates on one wing are congruent modulo 13, forcing `13^3` in that
`K4` tree sum; the hub and repair leaf are congruent modulo 13 but not 169,
forcing one more 13 in the bridge.  The other `K4` tree sum is a 13-adic
unit.  Thus the relation's 13-adic decoration is visible as a universal
`13^4` spectral factor without causing cancellation.

Equations (9) and (11) coexist.  The determinant is a nonlinear polynomial
on exact cochains, not a function of the absolute cohomology class.  This is
the same mechanism as THM-3482, now reached from the THM-3479 relation
coordinates rather than private-sheet counts.

## 5. What this does not transplant

The potentials in (10) are the nine integer coordinates of `U_full` and
`U_clock`.  They are not THM-2334 endpoint amplitudes, Fourier phases,
ancestry currents, or THM-2512 physical word currents.  Therefore (10) does
not prove that an endpoint bank lands on the weighted carrier.

The exact connection contract is:

| field | exact content |
|---|---|
| source | `Lambda_a/Q*s`, formed from relation-coordinate potentials |
| target | `B^1(G;Q)` on the THM-3473 thirteen-edge carrier |
| map | delete `q1`, apply role chart (3), take vertex differences |
| preserved | rational linear rank, role packet, mod-13 congruence, exactness |
| destroyed | endpoint phase, clock, masks, ancestry, grouped address, all-unit projector |
| needed sidecar | a proof that physical endpoint responses factor through the eight role potentials |
| cheapest hostile | flat potential: zero edges and zero determinant |

The next decisive test is not another arbitrary edge bijection.  It is to
construct eight role-labelled endpoint potentials on one THM-3479 tuple and
verify that their edge differences reproduce a lawful endpoint/current
observable.  Only then may (11) be read as spectral closure for that physical
row.

## 6. The Jacobian `H^1` class is a different type

For the fixed sporadic Jacobian map, the separate discriminant deduction is

```text
[-L] in H^1_et(Spec Q(a,b,c),mu_2),                   (12)
```

the common Kummer/sign class of its three coordinate cubics.  Equation (12)
is étale/Galois cohomology of a function field.  Equation (9) is ordinary
cohomology of a finite graph.  The common notation `H^1` does not provide a
map between them.

Thus a literal LRC-current-to-JC-flux identification would be ill-typed
without a realization functor or sheaf carrying graph edges to boundary
valuations/Kummer cocycles.  The lawful present diagram is a cospan:

```text
Lambda_a/Q*s  --D_role--> B^1(G;Q) --> H^1(G;Q) [zero]

JC cubic discriminants -----------------------> H^1_et(K,mu_2) [[-L]]. (13)
```

No cross-arrow in (13) has been constructed.  What has been gained is a
precise target for one: it must transport the role/valuation sidecars, not
merely identify two six- or seven-dimensional vector spaces.

## 7. Tournament and recurrence interpretation

The useful finite carrier is two four-vertex complete packets joined by one
forced edge.  Its `1+6+6` edge-orbit split and `1+1+6` vertex-role split are
intrinsic; a tournament orientation is only a gauge for the signed tree
polynomial.  The 72 charts forming one automorphism orbit show exactly which
label permutations are harmless.

THM-3484 then explains what happens when the corresponding degree-seven tree
determinant varies periodically across three residue states: its minimal
recurrence has order 22, with the common degree-seven leading term removed
from both nontrivial cubic Fourier colours.  This produces a clean hierarchy:

```text
role packet -> exact gradient -> nonzero tree determinant
            -> ternary determinant word -> order-22 recurrence.            (14)
```

Only the first three arrows are instantiated by the present relation-role
sidecar, and the first arrow is semantic rather than physical.  Ancestry,
bispectrum nonvanishing, scalar-row exclusion, and LRC(14) remain open.
