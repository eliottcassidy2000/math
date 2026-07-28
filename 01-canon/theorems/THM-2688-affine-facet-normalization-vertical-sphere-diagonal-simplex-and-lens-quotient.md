---
id: THM-2688
title: "Affine facet normalization, vertical sphere, diagonal simplex, and lens quotient"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every prime affine missing-label graph, the disconnected facet
  normalization has boundary-sphere vertical image but filled-simplex
  diagonal quotient.  Its mapping cones and mapping torus are explicit.
  At p=13 the free label quotient is L^11(13;1,2,3,4,5,6), and the
  THM-2657 physical lift is its nonzero coefficient Bockstein class 7.
  These are abstract/coarse invariants only: no physical cyclic gluing,
  endpoint transition, thirteenfold overlap, row exclusion, or LRC(14)
  conclusion follows.
source: root-2026-07-28-affine-facet-holotopy
depends_on:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
related:
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
script: 04-computation/affine_facet_normalization_quotient_holotopy_thm2688.py
output: 05-knowledge/results/affine_facet_normalization_quotient_holotopy_thm2688.out
script_sha256: 16685599b678e1e77e0a9a4c984cad70281e7edec02cb83ce9a715de7874aa7b
output_sha256: fecc676276c90dc21111cb2ef1524f64d5830ddd03a90dadb2f034f913f55fab
independent_script: 04-computation/affine_facet_normalization_quotient_holotopy_thm2688_referee.py
independent_output: 05-knowledge/results/affine_facet_normalization_quotient_holotopy_thm2688_referee.out
independent_script_sha256: 63ed1445d06abed8b736cfe1a9291310e248945f9de40691c19f89a0f2a7e3ba
independent_output_sha256: 4add8f23084dd238b4fbbe99c87176892e978fd54a9b7657459d25920457f76a
hash_basis: LF-normalized bytes
---

# THM-2688 -- affine facet normalization and quotient holotopy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The theorem separates four objects that were previously easy to conflate:
the disconnected carry normalization, its vertical label image, its diagonal
orbit quotient, and the physical odometer lift.  The first three have exact
homotopy types; the last is obstructed by the same nonzero class that appears
as torsion on the lens quotient.

## The general prime theorem

Let `p>=3` be prime, let `C=D=F_p`, and let

```text
m(c)=a+b c,                    b!=0.
```

Define a simplicial complex `K_m` with vertex set

```text
V(K_m)={(c,delta):delta!=m(c)}
```

and with a set of vertices declared a simplex exactly when all its vertices
have the same first coordinate `c`.  Let

```text
pi:K_m -> Delta(D),            pi(c,delta)=delta.
```

For any `s!=0`, put `r=bs` and

```text
g(c,delta)=(c+s,delta+r).
```

Then:

1. `K_m` is the disjoint union of `p` filled `(p-2)`-simplices.
2. The **simplicial image** of `pi` is `boundary Delta^(p-1)`.  This is not an
   orbit quotient.
3. `g` generates a free `C_p`-action on `K_m`.  The nonzero offset
   `t=delta-m(c)` is invariant, and it identifies

   ```text
   K_m/<g> = Delta(F_p^*) = Delta^(p-2).
   ```

4. With the standard one-apex topological mapping cone,

   ```text
   Cof(pi)  ~= S^(p-2) wedge (wedge^(p-1) S^1),
   Cof(q)   ~= wedge^(p-1) S^1,
   Tor(g)    = Delta^(p-2) x S^1,
   ```

   where `q:K_m->K_m/<g>` is the orbit map and `Tor(g)` is the mapping torus.
   Consequently

   ```text
   H~_1(Cof(pi);Z)=Z^(p-1),     H~_(p-2)(Cof(pi);Z)=Z,
   H~_1(Cof(q);Z)=Z^(p-1),
   ```

   with the two displayed groups in `Cof(pi)` added when `p=3` (so its
   `H~_1` has rank `3`).

5. Let `rho(delta)=delta+r` act on `boundary Delta^(p-1)`.  Then `pi` is
   equivariant.  Because `p` is prime and `r!=0`, `rho` acts freely on the
   geometric boundary.  After radial projection from the barycenter, this is
   the unit sphere in the reduced real regular representation of `C_p`.
   Hence, for `p=13`, keep the physical generator fixed and orient each real
   Fourier two-plane by its preferred complex coordinate to obtain

   ```text
   boundary Delta^12 / <rho> = L^11(13;1,2,3,4,5,6).
   ```

   Its integral homology is

   ```text
   H_0=H_11=Z,
   H_1=H_3=H_5=H_7=H_9=Z/13,
   H_i=0 otherwise.
   ```

   The equivariant map descends to a surjection

   ```text
   Delta^11 = K_m/<g> -> L^11(13;1,2,3,4,5,6).
   ```

   Its source is contractible, so its mapping cone is homotopy equivalent to
   the lens space.  Thus the odd-degree `13`-torsion is a precise invariant of
   the vertical facet gluing that is invisible in the diagonal orbit simplex.

## Proof

For fixed `c`, the allowed vertices are exactly the `p-1` labels in
`D\{m(c)}`, and different first-coordinate fibres have no common vertex.
This proves the first assertion.

A label set `S subset D` is in the simplicial image of `pi` precisely when
there is a `c` with `S subset D\{m(c)}`, equivalently `m(c) notin S`.  Since
`m` is bijective, this holds precisely for the proper subsets of `D`.  Those
are the faces of `boundary Delta^(p-1)`.

Compatibility is exact:

```text
m(c+s)=m(c)+bs=m(c)+r.
```

Thus `g` preserves the deleted graph.  Every nonidentity power changes the
first coordinate because `s!=0` and `p` is prime, so the action is free.  Also

```text
(delta+r)-m(c+s)=delta-m(c)=t.
```

Every row contains exactly one vertex for every `t in F_p^*`; hence the
vertex orbits are indexed by these `p-1` offsets and every subset of them is
a quotient simplex.  This proves `K_m/<g>=Delta^(p-2)`.

The complex `K_m` is homotopy equivalent to `p` discrete points, while the
vertical image is `S^(p-2)`.  Because every component of `K_m` is
contractible and the target is path connected, `pi` is homotopic to a
constant map (contract each component to one common target basepoint).
Therefore

```text
Cof(pi) ~= S^(p-2) wedge Sigma(K_m).
```

The unreduced suspension of `p` contractible components is a wedge of `p-1`
circles.  The same argument with contractible target `Delta^(p-2)` proves the
formula for `Cof(q)`.  In `(c,t)` coordinates the generator is
`(c,t)->(c+s,t)`; the mapping torus of one cyclic permutation of the `p`
components is therefore `Delta^(p-2) x S^1`.

For the equivariant quotient, a nonidentity power of `rho` is a single
`p`-cycle.  Its only fixed barycentric point in the full simplex has all
coordinates `1/p`, hence lies in the interior and not on the boundary.
Radial projection about this barycenter identifies the boundary
equivariantly with the sphere in

```text
{(x_delta) in R^p: sum_delta x_delta=0}.
```

Complex Fourier coordinates split that representation into the six conjugate
pairs of nontrivial characters.  Keep the physical generator fixed; within
each real two-plane choose whichever conjugate character has exponent in
`{1,2,3,4,5,6}`, then order the planes by that exponent.  The displayed
weights are therefore `1,2,3,4,5,6` without changing the generator.  This is
the stated lens space.  Its standard one-cell-in-each-degree cellular complex has
boundary multiplication by `13` in positive even degrees and zero in odd
degrees, giving the displayed homology.

## Exact THM-2672 specialization

In THM-2672,

```text
p=13,       m(c)=12-2c,       s=7,       r=-1,
g(c,delta)=(c+7,delta-1),      t=delta-m(c).
```

Thus the theorem gives exactly

```text
K = 13 disjoint Delta^11,
image(pi)=boundary Delta^12,
K/<g>=Delta^11,
Cof(pi) ~= S^11 wedge 12 S^1,
Cof(q) ~= 12 S^1,
Tor(g)=Delta^11 x S^1.
```

The abstract generator has `g^13=1`.  Its canonical physical lift is circle
translation by `7/13^6`, and THM-2657 gives

```text
(g_physical)^13: x -> x+91/13^6=x+7/13^5 != x.
```

The speed-one guard has no nontrivial rotation stabilizer.  Therefore the
abstract orbit simplex, mapping torus, lens quotient, and mapping cones do not
supply a physical `C_13` transition or a gluing of THM-2672's carry facets.
Present factors need not covary under the kernel translation.  No common
component map, endpoint current, mixed-configuration thirteen-fold overlap,
row exclusion, or LRC decrement follows.

## The lens torsion is the odometer Bockstein

Let `L=L^11(13;1,2,3,4,5,6)` be the label quotient above, and let

```text
alpha in H^1(L; C_13)
```

classify the principal `C_13` cover `S^11 -> L`.  Use the exact cyclic
extension from THM-2657,

```text
0 -> C_(13^5) -> C_(13^6) -> C_13 -> 0.                 (1)
```

Freeze the additive gauges from THM-2657: all coefficient actions are
trivial, the kernel embedding is `i(a)=13a`, the quotient is
`pi(k)=2k mod 13`, and the deck/root generator

```text
g:(c,delta)->(c+7,delta-1)
```

is identified with `1 in C_13`.  Use the cellular orientation
`boundary(e^2)=+13e^1` for the connecting map.

The connecting homomorphism sends `alpha` to

```text
beta_(1)(alpha)=7 in H^2(L;C_(13^5)) ~= C_13,           (2)
```

where the kernel coordinate is the canonical translation by `1/13^5` and
the physical slope-seven lift is translation by `7/13^6`.

Indeed, `pi_1(L)=C_13` and the universal cover `S^11` has vanishing first and
second cohomology.  The Cartan--Leray comparison therefore gives in degree two

```text
H^2(L;C_(13^5)) = H^2(C_13;C_(13^5))
                 = C_(13^5)/13 C_(13^5) ~= C_13.        (3)
```

The class `alpha(g)=1` lifts set-theoretically by `s(1)=7`.  Thirteen turns
translate by

```text
13*(7/13^6)=7/13^5.                                     (4)
```

Since `91=i(7)`, its extension cocycle represents `+7` modulo `13`, exactly
as in THM-2657.  Naturality of the connecting homomorphism identifies this
group extension class with `(2)`.  Reversing the deck generator or cellular
orientation rescales it by a unit; only nonvanishing is gauge invariant.  In
particular `alpha` has no lift to
`H^1(L;C_(13^6))`; equivalently, the principal label cover has no compatible
physical `C_(13^6)` lift whose quotient is the abstract label action.

This makes the two post-THM-2672 obstructions one object:

- the odd integral torsion `H_1(L)=Z/13` is the coarse quotient footprint;
- its coefficient Bockstein `(2)` is the physical odometer clutch preventing
  cyclic facet gluing.

The conclusion is only an obstruction to the specified cyclic lift.  It does
not exclude a noncyclic enlarged carrier, a groupoid-valued transition, an
affine handoff using changing digits, or any positive endpoint packet.

### Sharp boundaries

1. For a split extension `C_(13^5) x C_13`, the connecting class is zero and
   a lift exists.
2. Dimension one is exceptional: `L^1(13)=S^1` is not a `K(C_13,1)` through
   degree two, `H^2(S^1;-)` vanishes, and this Bockstein formulation does not
   detect the cyclic extension.  The `11`-dimensional THM-2688 quotient is
   safely above that boundary.
3. Replacing the generator by a unit rescales `7` but cannot make the class
   zero.

## Dowker/Čech/hocolim verdict

The relation

```text
R(c,delta) <=> delta!=m(c)
```

is `K_(p,p)` minus a perfect matching.  Its two Dowker complexes are both
`boundary Delta^(p-1)`: a proper set on either side has a common neighbour,
and the full set does not.  Likewise, the images of the `p` row facets form a
finite simplicial good cover of the boundary whose nerve is the same boundary
(equivalently, thicken the facets to open stars).  This is a useful typing
explanation, but it adds no invariant beyond the vertical image.

There is one categorical trap.  The hocolim over the **discrete** carry index
is the disjoint normalization `K_m`; only the full intersection/Čech diagram
recovers the boundary.  Saying merely “the hocolim is the sphere” silently
inserts the missing overlap morphisms.  The mapping cone of `pi` is the
cheapest exact ledger for that insertion.  The equivariant lens quotient is
the only genuinely additional invariant here, and it remains coarse and
nonphysical.

## Positive and hostile controls

The exact companion

```text
python 04-computation/affine_facet_normalization_quotient_holotopy_thm2688.py
python -O 04-computation/affine_facet_normalization_quotient_holotopy_thm2688.py
```

checks:

- all face counts and all `13` diagonal vertex orbits at `p=13`;
- the vertical and orbit quotients separately;
- chain-level mapping-cone Betti numbers for `p=3,5,7` over
  `F_2,F_3,F_5,F_7,F_11`;
- the free label action and lens weights;
- a nonbijective constant `m` hostile, where the vertical image collapses to
  one filled simplex rather than a boundary; and
- a non-generator step on `Z/6`, where the diagonal quotient has two filled
  simplex components rather than one.

These hostiles show that bijectivity of `m`, compatibility of the two shifts,
and transitivity of the component step are all load-bearing.

An independent direct referee reconstructs the relation and diagonal orbits
without importing the primary mapping-cone code.  It verifies the prime
cases `3,5,7,13`, the reduced-regular character split
`{+/-1,...,+/-6}=F_13^*`, the lens cellular homology, and all `169` values of
the section defect for `s(a)=7a`.  The latter has `91` zero entries and `78`
wrap entries equal to `7`, so the coefficient Bockstein is nonzero.  Its
additional incompatible-label-shift and split-extension hostiles show,
respectively, that `r=bs` and nonsplitting are essential.  Normal and
optimized executions byte-match

```text
05-knowledge/results/affine_facet_normalization_quotient_holotopy_thm2688_referee.out
```

with the frontmatter hashes above.

QED.
