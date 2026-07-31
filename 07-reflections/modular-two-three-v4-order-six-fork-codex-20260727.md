---
title: "Two and three on one object: the modular order-six fork and its missing sidecars"
date: 2026-07-27
status: RESEARCH SYNTHESIS; current proved corrections route through THM-2632, THM-2633, THM-2641, and THM-2646; older frontier paragraphs below preserve provenance
source: codex-2026-07-27-modular-v4-session
---

# The answer is a fork, not a coincidence

> **Current correction (2026-07-28).**  THM-2632 proves that the Farey
> six-state object is the `S3` permutahedral hexagon with three intrinsic
> `V4` theta cuts, whereas the simultaneous integral `2/3` lift has cyclic
> abelianization `C6`; their uncoloured hexagons are not canonically the same.
> THM-2641 strengthens this to a no-go: `T^1274` preserves the live mod-2,
> mod-49, and mod-91 residue skeleton while traversing the missing `C3`
> fibre, and the three conjugates of `T` keep one abelianized class while
> exhausting the theta directions.  THM-2633 independently excludes
> `C4,V4,D4` quartic Keller monodromy, so later sentences treating `D4` as
> live are historical and superseded.  The surviving quartic types are
> `A4,S4`; the cubic resolvent remains action data, not an affine Keller map.

## The common object is the braid helix

THM-2646 supplies the literal object that the earlier synthesis was missing:

```text
B3=<x,y | x^2=y^3>,       z=x^2=y^3,
B3/<z>=C2*C3=PSL2(Z).                                      (0)
```

Thus the binary and ternary generators are roots of one common full twist
*before* quotienting.  Exponent sum makes this exact:

```text
B3 = PSL2(Z) x_(C6) Z.                                    (0a)
```

The modular hexagon records only height modulo six; the braid lives on its
integer helix.  This is consequential, not decorative: `(sigma1 sigma2)^-1`
and `(sigma1 sigma2)^2` occupy the same exact modular point but close to the
unknot and trefoil, and their central fibre contains torus knots of unbounded
genus.  The standard affine `V4:S3` shadow kills the full twist, so Feuerbach
origins, theta cuts, and uncoloured graceful repairs cannot recover it.  A
monic cubic discriminant winding does recover the missing integer height,
but only together with the full modular braid class; finite `S3` inertia is
still too coarse.  This is the sharp form of the recurring rule: a quotient
kernel is a gauge only for consumers that literally annihilate it.

The recurring `2,3,4,6` pattern has a rigorous core, but it is not one
universal tree.  There are two natural order-six quotients of

```text
PSL_2(Z)=C2*C3:

  abelianization                  mod-two quotient
  C2*C3 -> C2 x C3 = C6          C2*C3 -> GL2(F2)=S3.        (1)
```

The first forces the two torsion generators to commute.  The second keeps
their noncommuting action on the three nonzero directions of `F2^2`.
Cardinality six hides this distinction.

The smallest separator is the commutator.  It dies in `C6` and is a
nontrivial three-cycle in `S3`.  Consequently any carrier built only from
commuting diagonal insertions can see at most the first shadow.  This is
exactly the information boundary of THM-2602: a vertex calculus is
commutative, while the live LRC obligation is an ordered transition kernel.

## The three exact incarnations

### 1. Farey and Pythagorean trees

THM-2596 separates the group action from the rooted grammar.

- The Bass--Serre tree is `(2,3)`-biregular.  Contracting its degree-two
  vertices gives a trivalent tree; rooting makes the later branching binary.
- The positive Farey children are parabolic words in the torsion
  generators, not the torsion generators themselves.
- The three Berggren/Pythagorean branches form a PGL reduction partition of
  the same rational interval.  With

```text
lambda(x)=x/(1+x),  rho(x)=1/(2-x),  j(x)=1-x,
```

  their exact prefix code is

```text
{rho, lambda rho j, lambda lambda}.                         (2)
```

Thus the ternary grammar retains one binary branch and refines the other,
using one orientation reversal.  It is not the `C3` torsion action.

The faithful quantitative state for the THM-2056 flank is not a tree label
but

```text
(G,ell)=(U^T U,U^T w),
F_w(Uz)=z^T G z-91 ell^T z.                                (3)
```

The cross Gram entry is load-bearing: two acute unimodular flanks with the
same two endpoint defects `(-90,-89)` have mediant defects `-177` and `+1`.
This is a reusable warning: a quotient that preserves endpoints need not
preserve insertion.

### 2. The six-vertex bicycle

THM-2597 identifies the first nonzero tile bicycle as a `C6`.  Its three
partial-cube Theta classes have two edges each.  The order-two half-turn
swaps the two edges in every class; the order-three rotation cycles the
classes.  These commuting rotations are exactly the `C6` side of (1).

The same cycle is `K3,3` minus one perfect matching.  It retains only the
two derangement matchings.  Restoring the missing matching restores all six
`S3` matchings and the permanent-versus-determinant sign obstruction of
THM-2290.  Completion changes the matching family; it is not a relabeling.

This also types the Krenn--Gu--Soltesz connection.  The alternating
matchings of the even cycle give the two monochromatic inherited vertex
colorings.  The selected pair tensor and hafnian contraction are the
faithful data; the uncolored support or a tournament is not.

The `C6` support is a partial cube and is not graceful, while deleting one
edge gives a graceful `P6`.  Partial-cube and graceful predicates therefore
do not recover each other.

### 3. Quartic roots and their cubic resolvent

THM-2595 and THM-2598 expose the other side of (1).  Four quartic sections
form an affine `V4` torsor.  Choosing an origin identifies

```text
S4=AGL2(F2)=V4 semidirect GL2(F2),
S4/V4=GL2(F2)=S3.                                          (4)
```

The cubic resolvent remembers the three nonzero directions / `2+2`
pairings and forgets the origin.  Its discriminant equals the quartic
polynomial discriminant identically, and its cube-minus-square cusp is
universal invariant theory.  Neither fact transfers Keller anatomy.

The field-theoretic reason is decisive.  For generic `S4` closure `L/K`,

```text
quartic root field = L^S3,       cubic resolvent field = L^D8. (5)
```

The subgroups are incomparable.  The resolvent is a Galois-closure
correspondence, not an intermediate cover of the quartic source.  A raw
double root can also be only an order-index collision: in the exact hostile
family the normalization chart is

```text
R(-2+tY)=t^2[tY^3-4Y^2+4sY-1],                             (6)
```

whose residual quadratic is separable away from `s^2=1`.

THM-2606 makes the partial-cube content exact.  The three nonzero parities
of an affine `V4` torsor are simultaneously:

```text
three omitted translation directions;
three 2+2 quartic-resolvent channels;
three translation-invariant C4 partial cubes.              (7)
```

Full affine `S4` symmetry leaves only the empty graph and `K4`, so it cannot
choose a connected partial cube.  Choosing an origin reduces to `S3` and
selects the unique invariant connected partial cube `K1,3`.  On the six
unordered pairs, the intrinsic object is the incident/disjoint association
scheme, not a tournament.  Its four antipode-compatible `C6` structures are
equivalent to the four possible origins.

### 4. A nonabelian seven-edge quotient at the LRC prime

The modular fork has a second exact finite shadow which is much closer to
the live ordered-transition problem than either order-six quotient.  Work in
`PSL2(F13)` with

```text
U^a=[[1,a],[0,1]],             g_t=[[0,1],[-1,t]],
t in {3,5,6}.                                                (7a)
```

Each `g_t` has projective order seven.  Put

```text
A_k=g_t^k U^a g_t^(-k).
```

The forward and reverse ordered conjugation products telescope before any
entrywise calculation:

```text
A_0 A_1 ... A_6=(U^a g_t)^7 g_t^(-7),
A_6 A_5 ... A_0=g_t^7(g_t^(-1)U^a)^7.                       (7b)
```

The two moving traces are respectively `t-a` and `t+a`.  The noncentral
projective-order-seven trace set in `F13` is

```text
{plus_or_minus 3, plus_or_minus 5, plus_or_minus 6}.         (7c)
```

Consequently the forward product is central exactly when `t-a` lies in
(7c), and the reverse product is central exactly when `t+a` lies there.
For `t=3,5,6`, the two orientations together hit every
`a in F13^*`.  THM-2603 proves the resulting exact finite atlas: the
coverage is a trace-intersection theorem, not a
coincidence of twelve matrix products.

This is also the cleanest demonstration yet of why order cannot be
marginalized.  Reversing the same seven conjugates changes `t-a` to `t+a`.
The construction is an exact noncommutative design target for the seven-edge
LRC clock, but it is not yet a Boolean carrier, semantic owner path, positive
kernel, or relation current.  A lawful use must realize the conjugating
`g_t` motion and the parabolic `U^a` insertion on the same physical ancestry.

The `84=14*6` frame count also has a stronger affine law.  Removing the six
frames over infinity leaves `78` states which form a regular torsor for

```text
B=C13 semidirect C6,
(c,kappa)(d,lambda)=(c+kappa d,kappa lambda).              (7d)
```

Thus proved THM-2611's thirteen chronological choices and the six tangent
frames cannot simply be checked independently.  Their cyclic obstruction is
a frame product together with a frame-weighted root sum.  The missing
physical object is a lawful principal-B transition carrier; the abstract
regular action does not identify the chronological `C13` with the root
translation subgroup.

# Geometry clues: what survives

## Feuerbach

The incenter and three excenters have four projective sign classes

```text
{+/-1}^3 / simultaneous sign = V4.                          (8)
```

Internal versus external tangency to the nine-point circle singles out the
incenter class.  This supplies an origin and hence the abstract `K1,3` in
THM-2606.  Sign multiplication is not a motion between circles, and no
metric or Keller dynamics transfers.

## Four-body coordinates

For centred scalar positions `sum x_i=0`, the three Hadamard/Jacobi pair
coordinates

```text
X=x1+x2-x3-x4,  Y=x1-x2+x3-x4,  Z=x1-x2-x3+x4             (9)
```

have squares equal to four times the three squared-pair resolvent roots.
The double transpositions act by even sign flips.  This is an exact
kinematic identity, but it imports no force law, collision regularization,
or N-body integrability.

## Napkin ring

The radius-forgetting phenomenon is a useful no-go model for controlled
forgetting.  In dimension `d`, bore radius `a`, height `h`, and axial
coordinate `z`, the cross-section is

```text
kappa_(d-1)[(a^2+h^2/4-z^2)^((d-1)/2)-a^(d-1)].           (10)
```

Dependence on `a` disappears exactly when `(d-1)/2=1`, hence only in
dimension three.  Then integration gives `pi h^3/6`.  This is a quadratic
codimension-two accident, not a generic license to discard a parameter.
The repo analogue is: eliminate a nuisance coordinate sectionwise and
prove the cancellation before integrating; never infer it from the final
aggregate alone.

## Super parity and gracefulness

A nontrivial `Z2` grading on pointed `V4` is one of the three parity
functionals in (7).  `S3` permutes them, so no nonzero grading is invariant.
This captures the super-parity analogy but supplies neither a Lie
supergroup nor a `Z3` grading.

`K4`, `C4`, and `K1,3` are all graceful.  Gracefulness therefore cannot
choose the missing affine origin or resolvent channel.

# Frontier consequences

## LRC(14)

The six-state resemblance is now a stopping rule.  THM-2587's
`C2 x {low,middle,high}` wall is ordered/reflected, not cyclic.  THM-2600's
six middle digit is a digit value, not an order-six group action.  A
commutative gate algebra collapses free words to the `C6` shadow; it cannot
manufacture the noncommuting adjacent-root transition demanded by
THM-2602.  The next useful object is more sharply one positive composable
`C13 semidirect C6` kernel

```text
K_ell((q,kappa),(q',kappa'))
```

with the adjacent root, clock, owner, and order retained before
marginalization.

## Degree-four Keller / G1

The cubic resolvent gives a normalized branch sub-divisor and an exact
index tax, not a degree-three Keller map.  THM-2612 now separates the three
live monodromy shadows exactly:

```text
D4 -> C2,        A4 -> C3,        S4 -> S3=S4/V4.          (11)
```

Only `A4,S4` are quotients of the modular `C2*C3`; `D4` has no element of
order three and cannot occur as such a quotient.  Its replacement detector
is geometric: the unique punctured-source deck involution cannot extend
polynomially, so some source coordinate has a genuine Jelonek-boundary pole.
For a depressed quartic its deck denominator has norm `S'(w)^2`, while a
tame double-transposition/four-cycle index prime divides both `q` and
`p^2-4r`.  This still does not exclude `D4`, but it turns the old cardinality
analogy into two exact computational gates.

# Reusable meta-pattern

When `2,3,4,6` reappear, ask which of the following was chosen or forgotten:

```text
free word -> commutative exponent sums -> C6;
free word -> mod-two linear action     -> S3;
four affine sections -> forget origin -> three channels;
choose channel -> C4;
choose origin  -> K1,3 or one of four edge-state C6 cycles.
root torsor + frame torsor -> C13 semidirect C6, not C13 x C6.
```

The preserved cardinality is never enough.  The required sidecar is one of
word order, Gram cross term, affine origin, selected channel, matching
family, or physical transition correspondence.
