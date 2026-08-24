---
id: THM-3951
title: "Affine-plane boundary incidence is a forest and excludes the equianharmonic survivor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On a normal surface containing A2, the incidence multigraph of prime
  boundary curves and the exceptional paths between their branches embeds
  in the tree boundary of a common resolution. In particular two distinct
  boundary primes cannot meet at two distinct smooth points. Applied to the
  minimal-debt A1 internal-split row of THM-3950, the graph ramification
  prime and its irreducible equianharmonic residual meet over at least two
  finite color fibres whenever gcd(c,RUV)=1. Both primes would have to be
  deleted by a same-function-field Keller chart, contradicting the forest
  lemma. The explicit R=1,S=Y survivor has two transverse smooth
  intersections at Y=omega and Y=-omega^2. Thus its reduced one-place
  component, normal quadratic surface, and two independent Cardano classes
  still cannot produce a planar Keller map.
source: jc-degree6-one-place / post-THM-3950 boundary-incidence audit, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-zero-debt-lift and jc-cohn3709,
  2026-08-24). Both audits reconstructed the boundary-tree multigraph and
  connected-fibre paths, the monic-cubic domain proof, residual
  irreducibility, the three-color finite-point count, generic ramification
  after normalization, and the Zariski-Main source/target bridge. Tangencies
  are correctly allowed and gcd(c,RU0V0)=1 is load-bearing. Normal and
  optimized runs byte-match the frozen 61-gate output, all hashes agree,
  documentation checks pass, and no repair was required.
depends_on:
  - THM-3950-a1-internal-split-denominator-debt-and-equianharmonic-shadow
related:
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
script: 04-computation/jc2_affine_plane_boundary_incidence_forest_thm3951.py
output: 05-knowledge/results/jc2_affine_plane_boundary_incidence_forest_thm3951.out
script_sha256: 3e021a8b9591bc37e634c021cb34881aef67cf2a29b410decf0a1e35d5e900c8
output_sha256: 0bd0fb1f4eccb6946358bc2d829b34705905ff65c9ad13ed058ce6a1d06bc44f
semantic_sha256: fda10e083ea8bb512edaca2275da87745bafe289b24b5ff6ee5eddab546b3a76
hash_basis: raw LF bytes
---

# THM-3951 -- the equianharmonic survivor fails the boundary-incidence forest

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero. Fix

```text
omega^2+omega+1=0,       delta=omega-omega^2,       delta^2=-3.       (1)
```

This theorem separates three surfaces that must not be conflated:

```text
A2_(P,t)       the polynomial target;
X              the finite normal cubic surface over that target;
U ~= A2        a hypothetical Keller source open inside X.            (2)
```

The discriminant curves are finite curves downstairs in `A2_(P,t)`. Their
ramification primes upstairs in `X` would be boundary curves of `U`. The
obstruction concerns the incidence of those **upstairs** primes; it does not
move a finite target branch curve to infinity.

## 1. Boundary incidence on an affine-plane completion is a forest

Let `X` be a normal integral surface containing a dense open

```text
U ~= A2.                                                            (3)
```

Choose a normal projective completion `Xbar` and resolve the birational map,
equal to the identity over `U`, between `Xbar` and the standard completion
`(P2,L_infinity)`. Thus there is a smooth projective surface `W` and proper
birational morphisms

```text
                         W
                       /   \
                    Xbar   P2                                      (4)
```

which are isomorphisms over `U`, and whose reduced boundary
`D=W minus U` is SNC. The map `W -> P2` factors into point blowups supported
over `L_infinity`. Starting from the one-vertex graph of `L_infinity`, a
blowup at a smooth boundary point adds a leaf and a blowup at a crossing
subdivides an edge. Therefore the dual **multigraph** of `D` is a tree.

This gives the following usable form of the boundary-incidence forest lemma.

> Let `D1,D2 subset X minus U` be distinct prime curves. They cannot meet at
> two distinct points `x1,x2` at which `X`, `D1`, and `D2` have smooth
> distinct germs. Tangency is allowed in this statement.

Indeed, take a sufficiently high relative log resolution in `(4)` and let
`D1tilde,D2tilde` be the strict transforms. If the map to `Xbar` is already
an isomorphism at `xi`, the corresponding intersection gives an edge between
the two vertices. Otherwise its proper birational fibre is connected because
`Xbar` is normal. The two branch attachments are consequently joined by a
path of exceptional boundary components. Fibres over `x1` and `x2` are
disjoint, so the two paths are internally disjoint. Their union is a cycle.
If both paths are direct edges, the two parallel edges are a multigraph
two-cycle. Both possibilities contradict the tree property.

The same proof gives THM-3920's unibranchness lemma by taking two branches of
one prime curve. The present two-prime consequence is genuinely different:
each curve may be smooth and unibranch, yet their repeated incidence is still
forbidden.

## 2. The universal minimal-debt cubic

Let `R,S,c in k[t]` satisfy

```text
gcd(R,S)=1,       R/S nonconstant,       c!=0,
U0=S+omega^2 R,  V0=S-omega R,          gcd(c,R U0 V0)=1.             (5)
```

Define the minimal denominator-debt packet

```text
g=c U0 V0,                  r=c R U0 V0,
A=c(S-R)V0^2,               B=c(S+R)U0^2,
E=(B-A)AB,
C=(1-omega^2)B-(1-omega)A=3r+delta c S^3.             (6)
```

The internal-split identities are

```text
AB=(cS U0 V0)^2-r^2,
A[(cS U0 V0)^2-omega r^2]=(cS U0 V0)^3-r^3,
B[(cS U0 V0)^2-omega^2 r^2]=(cS U0 V0)^3+r^3.          (7)
```

Put

```text
F=T^3-3PT-(E+CP),
X0=Spec k[P,t,T]/(F).                                  (8)
```

This cubic is a domain. If the monic cubic `F` were reducible over the UFD
`k[P,t]`, it would have a root in `k[P,t]`: a root in the fraction field is
integral, and the UFD is normal. Comparing `P`-degrees first makes the root
independent of `P`; comparing the linear and constant rows then gives

```text
root=-C/3,                    C^3+27E=0.                (9)
```

After division by `B^3`, `(9)` is a nonzero algebraic equation over `k` for
`A/B`. But

```text
A/B=((S-R)(S-omega R)^2)/((S+R)(S+omega^2 R)^2)         (10)
```

is the composition of two nonconstant rational maps, hence is nonconstant.
This contradiction proves the claim.

Let `X` be the normalization of `X0`. It is finite over `A2_(P,t)`, because
`X0` is finite there and finitely generated algebras over `k` are excellent.

## 3. The ramification divisor has two repeatedly incident primes

The relative derivative of `(8)` is

```text
F_T=3(T^2-P).                                             (11)
```

On this scheme put `h=-T`, so `P=h^2`. The remaining equation factors as

```text
E+C h^2-2h^3=(h-r)Q,
Q=-2h^2+(C-2r)h+r(C-2r).                                (12)
```

The first factor gives the graph ramification prime

```text
EGamma: P=r(t)^2,                 T=-r(t),               (13)
```

whose normalization is `A1_t`. The residual quadratic is irreducible. Its
discriminant has square class

```text
(R-S)(R+S)(R+delta S)(3R+delta S),                       (14)
```

the pullback under the nonconstant map `R/S` of a smooth genus-one double
cover of `P1`. If `(14)` were a square in `k(t)`, it would give a nonconstant
map `P1 ->` that genus-one curve, impossible by Riemann--Hurwitz. Thus `Q`
defines one irreducible residual ramification prime `ER`.

The exact incidence identity is

```text
Q(r)=2r(C-3r)=2delta c r S^3.                            (15)
```

Consider the three colors

```text
R=0,                 U0=0,                 V0=0,
R/S=0,               R/S=-omega,           R/S=omega^2. (16)
```

At any finite point in one of these fibres, `(5)` gives

```text
r=0,       S!=0,       c!=0,
F_P=3h-C=-delta cS^3!=0,
Q_h=delta cS^3!=0.                                      (17)
```

Consequently `X0`, `EGamma`, and `ER` have smooth distinct germs there; the
normalization `X -> X0` is an isomorphism near the point. The nonconstant
map `R/S:P1_t -> P1` hits all three colors in `(16)`. Only one color fibre
can contain the single point `t=infinity`, so at least two of these clean
intersection points lie in the affine parameter line. They are distinct
because the three color fibres are disjoint. Hence `EGamma` and `ER` meet at
at least two distinct smooth points of the affine normal surface `X`.

This argument allows ramification of the color map: the two germs may be
tangent. The forest lemma needs two distinct incidences, not transversality.
The hypothesis `gcd(c,R U0 V0)=1` is the stated boundary of this assertion;
additional common-debt collisions are not classified here.

## 4. Why both curves would be Keller boundary

Suppose, for contradiction, that the same function field admitted source
coordinates `x,z` such that

```text
k(x,z)=Frac(X),          P,t in k[x,z],
Jac_(x,z)(P,t) in k*.                                      (18)
```

Then

```text
f:A2_(x,z) -> A2_(P,t)                                    (19)
```

is etale and quasi-finite. Zariski Main factors `(19)` through the finite
normalization as

```text
A2_(x,z) = U  --open immersion-->  X  --finite--> A2_(P,t). (20)
```

At the generic point of `EGamma`, equation `(17)` with `t` generic reads
`F_P=-delta cS^3`, which is nonzero. At the generic point of `ER`, the
quadratic `Q` cannot divide the linear polynomial `3h-C`, so `F_P` is again
nonzero generically. Thus `X0` is normal at both generic points, `(11)`
shows that the finite map is ramified there, and both primes lift uniquely
to ramification primes of `X`.

The restriction of the finite map to the open `U` in `(20)` is the etale
Keller map. Therefore

```text
EGamma union ER subset X minus U.                         (21)
```

They are boundary **upstairs** even though their images are ordinary finite
branch curves **downstairs**. Their two smooth intersections now contradict
Section 1. Hence no same-function-field planar Keller chart `(18)` exists.

## 5. The smallest survivor is already excluded

For THM-3950's explicit row take

```text
R=1,                 S=Y,                 c=1,
s=(Y+omega^2)(Y-omega)=Y^2-delta Y-1,
A=(Y-1)(Y-omega)^2,
B=(Y+1)(Y+omega^2)^2.                                  (22)
```

Here `r=s` and

```text
C=3s+delta Y^3,
E=-s^2(s+delta Y^3),
Q=-2h^2+(s+delta Y^3)h+s(s+delta Y^3).                  (23)
```

The residual discriminant is

```text
disc_h(Q)=(s+delta Y^3)(9s+delta Y^3)
 =delta(Y-1)(Y+1)(delta Y+1)(Y-delta)^3,                (24)
```

so `Q` is irreducible. Moreover

```text
Q(s)=2delta sY^3.                                       (25)
```

The two simple roots of `s` give the distinct points

```text
(P,T,Y)=(0,0,omega),             (0,0,-omega^2).         (26)
```

At each point

```text
F_P=-delta Y^3!=0,       Q_h=delta Y^3!=0,              (27)
```

and `(25)` has a simple zero. Thus these are transverse intersections of
two smooth ramification curves on a smooth surface. The normalization cannot
separate or remove them. Sections 1 and 4 exclude the explicit survivor from
every same-function-field planar Keller chart.

The conclusion is deliberately narrow. It does not refute the normal
quadratic surface or the two independent Cardano classes proved in THM-3950;
it shows that those successful algebraic invoices are incompatible with the
topology of an affine-plane source boundary. It also does not classify rows
where `c` shares a factor with `R U0 V0`, distributions across several cube
factors, non-`A1` primary branches, or arbitrary cubic fields with different
branch data. `JC(2)` remains open.

Run

```bash
python3 04-computation/jc2_affine_plane_boundary_incidence_forest_thm3951.py
python3 -O 04-computation/jc2_affine_plane_boundary_incidence_forest_thm3951.py
```

for the assertion-only exact companion.
