---
id: THM-2598
title: "Quartic V4 resolvent torsor, universal cubic cusp, and the exact transfer boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED + CORRECTION.
  The integral cubic resolvent of a general quartic has
  exactly the quartic polynomial discriminant; its depressed cubic and
  cube-minus-square cusp identity are universal invariant theory, not
  Keller constraints.  On the nondegenerate depressed open, reconstruction
  from the three squared pair sums is an exact four-point V4 torsor.  The
  homogeneous leading-coefficient boundary has a forced 1+2 shadow and an
  explicit residual-quadratic valuation gate.  The complete transitive
  branch list C4,V4,D4,A4,S4 has matching images C2,1,C2,C3,S3; only
  S4 has a full S3 quotient.  The D4 branch has a generic deck C2 whose
  polynomial extension is not proved, correcting THM-1375/2465's deletion
  of D4 by that route.  THM-2633 later excludes D4 by the independent affine
  point-stabilizer abelianization gate, without extending the deck map.  The
  S4/A4/D4 field lattice,
  tame inertia table, and local order-index tax prove that the cubic
  resolvent is a Galois-closure correspondence, not an intermediate cover
  of the quartic source.  The Zariski-main model further pins the D4
  extension gap to an antipodal present/omitted unramified divisor pair in
  the finite normalization.  Consequently no grade-three Keller anatomy
  descends without an additional affine/Keller realization and V4-origin
  sidecar.  This theorem alone excludes no D4/A4/S4 branch; after THM-2633
  only A4/S4 remain possible.  THM-3438 later realizes the S4 branch and
  settles global G1 positively; A4 realization, planar degree four, JC(2),
  and DC(2) remain open.
source: codex-2026-07-27-quartic-resolvent-transfer
depends_on:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-1365-galois-reduction-jc-fixed-point-bridge
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
related:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-1375-reduced-jacobian-lattice-reflection-driven
  - THM-1310-conic-pair-fibers-and-design-equations
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2566-two-chart-saturated-cusp-atlas-and-parasitic-plane-ledger
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
  - MISTAKE-297
script: 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
output: 05-knowledge/results/jacobian_quartic_v4_resolvent_thm2598.out
script_sha256: 62683c2ee2a1b4f123b155cc22ea9ae52da707c43c7af0b8de16e7dc2fd1f888
output_sha256: 780d663948ea30bee8f65fe4887ea7107b07a6380bb804b75e60576e9fa1eea5
hash_basis: working-tree bytes (LF)
---

# THM-2598 -- what the quartic resolvent remembers, and what it cannot carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED + CORRECTION.**

> **GLOBAL G1 UPDATE (THM-3438 / MISTAKE-396).**  The weighted quartic is an
> explicit `S_4` 2-jet Keller map.  This realizes one of the two surviving
> monodromy rows without altering any resolvent, `V_4`, or boundary statement
> proved here.  The `A_4` and planar lanes remain open.
All algebraic identities and finite permutation claims below
are reproduced by the exact companion.  The Galois/valuation assertions are
proved here.  This theorem strengthens the classical layer of THM-2455; it
does not promote the open G1 verdict of THM-2465.

Throughout, the coefficient field has characteristic different from `2,3`
when division or tame inertia is used.  The polynomial identities themselves
hold integrally over `Z` after denominators in (7) are cleared.

## 1. The integral resolvent and the universal cusp identity

For

```text
f(T) = A T^4 + B T^3 + C T^2 + D T + E                         (1)
```

define the monic integral resolvent

```text
R_f(W) = W^3 - C W^2 + (B D - 4 A E) W
         + 4 A C E - B^2 E - A D^2.                            (2)
```

If `alpha_1,...,alpha_4` are the roots of `f`, the three roots of (2)
are

```text
beta_1 = A(alpha_1 alpha_2 + alpha_3 alpha_4),
beta_2 = A(alpha_1 alpha_3 + alpha_2 alpha_4),
beta_3 = A(alpha_1 alpha_4 + alpha_2 alpha_3).                   (3)
```

For example,

```text
beta_1-beta_2 = A(alpha_1-alpha_4)(alpha_2-alpha_3),             (4)
```

and the other two differences give the other four root differences.
Squaring their product proves the exact identity

```text
Disc_W(R_f) = Disc_T(f).                                        (5)
```

Set

```text
I = C^2 - 3 B D + 12 A E,
J = 2 C^3 - 9 B C D - 72 A C E + 27 B^2 E + 27 A D^2.          (6)
```

The shift `W=Z+C/3` changes (2) into

```text
Z^3 - (I/3) Z - J/27,                                          (7)
```

and hence

```text
27 Disc_T(f) = 4 I^3 - J^2.                                    (8)
```

Thus the depressed cubic law, equality of discriminants, and the
cube-minus-square cusp are automatic for every quartic.  They are not
consequences of the Keller equations.  This is the first transfer boundary:
THM-2473's depressed law and its particular physical coefficients cannot be
inferred merely from (7)--(8).

There is also no leading-factor transfer hidden here: (2) is monic even when
the quartic leading coefficient `A` varies.  Dividing to an unscaled root
resolvent introduces rational denominators, while clearing them changes the
visible leading factor by gauge and square-index terms.

## 2. Depressed reconstruction is exactly a `V_4` torsor

Let

```text
f(T)=T^4+pT^2+qT+r                                             (9)
```

be separable, with `q != 0`, and label its roots `x_1,...,x_4` so that
their sum is zero.  Put

```text
s_1=x_1+x_2,  s_2=x_1+x_3,  s_3=x_1+x_4,  u_i=s_i^2.           (10)
```

The `u_i` are the three roots of

```text
S_f(U)=U^3+2pU^2+(p^2-4r)U-q^2,                               (11)
```

and `s_1 s_2 s_3=-q`.  Conversely, order the three roots of (11),
choose square roots satisfying

```text
s_i^2=u_i,             s_1 s_2 s_3=-q,                         (12)
```

and set

```text
x_1=( s_1+s_2+s_3)/2,   x_2=( s_1-s_2-s_3)/2,
x_3=(-s_1+s_2-s_3)/2,   x_4=(-s_1-s_2+s_3)/2.                  (13)
```

Equations (13) recover the four roots of (9).  There are exactly four
sign triples in (12).  Flipping an even number of signs acts on (13) by

```text
1, (12)(34), (13)(24), (14)(23),                               (14)
```

so the reconstruction fibre is a free `V_4` torsor.  At `q=0` one of the
`u_i` is zero and this free torsor degenerates; that boundary is not included
in the claim.

The smallest marked-coordinate control is

```text
f_+(T)=T^4-7T^2+6T,       f_-(T)=T^4-7T^2-6T,
S(U)=(U-1)(U-4)(U-9).                                        (15)
```

The two quartics have the identical normalized resolvent but opposite product
orientation and are exchanged by `T -> -T`; thus (15) is a marked-gauge
hostile, not a pair of nonisomorphic covers.  For `f_+`, the four choices

```text
(1,2,-3), (1,-2,3), (-1,2,3), (-1,-2,-3)                      (16)
```

recover `(0,1,2,-3)` and its three double-transposition relabellings.  This
is a literal example of four different `V_4` sections over one ordered
resolvent triple.

At field-algebra level the loss is stronger.  The distinct biquadratic fields

```text
Q(sqrt(2),sqrt(3)): T^4-10T^2+1,    S=U(U-8)(U-12),
Q(sqrt(2),sqrt(5)): T^4-14T^2+9,    S=U(U-8)(U-20)              (17)
```

both have completely split cubic resolvent algebra `Q^3`.  Thus the quotient
algebra can be identical while the `V_4` torsor/field is not.  This is a
boundary control (`q=0`), separate from the free-torsor open above.

### 2.1 Universal label torsor versus specialized monodromy

The torsor in (12)--(14) is the universal ordered-label quotient

```text
S4 -> Sym{12|34,13|24,14|23}=S3,       kernel V4.              (17a)
```

It must not be confused with either the kernel or the deck group of a
specialized quartic field.  For transitive monodromy `G<=S4`, the invisible
matching kernel is `G intersect V4`; for quartic point stabilizer `H`, the
generic deck group is `N_G(H)/H`.  The exact companion enumerates all `30`
subgroups of `S4`; its `9` transitive subgroups form exactly these five
conjugacy types:

| quartic `G` | matching image | matching algebra | `|G intersect V4|` | `|H|` | `|N_G(H)/H|` |
|---|---|---|---:|---:|---:|
| `C4` | `C2` | `1+2` | 2 | 1 | 4 |
| `V4` | `1` | `1+1+1` | 4 | 1 | 4 |
| `D4` | `C2` | `1+2` | 4 | 2 | 2 |
| `A4` | `C3` | irreducible cyclic cubic | 4 | 3 | 1 |
| `S4` | `S3` | irreducible full-`S3` cubic | 4 | 6 | 1 |

Thus a cubic resolvent polynomial is not always a cubic **field**.  The
`D4` row has one rational matching plus a quadratic conjugate pair; the
`A4` row has a cyclic cubic; only `S4` has a full-`S3` cubic.  Likewise a
nontrivial `N_G(H)/H` is initially a function-field deck group over the
finite-etale locus, not automatically a polynomial automorphism group of
the affine source.

The `D4` row has an additional exact square model.  Its unique fixed
matching `c` is the rational resolvent root.  Declare the other two nonzero
`V4` translations to be edges on the four sheets; this is a `C4` whose two
diagonals are the `c`-pairs, and

```text
D4=Aut(C4)=Stab_S4(c).                                        (17b)
```

The central half-turn of this square is translation by `c`; it represents
the nontrivial element of `N_D4(H)/H`.  Thus the fixed resolvent channel and
the generic deck involution share one antipodal label direction, although
the matching quadratic field and the root-field deck quotient remain the
distinct algebraic constructions typed in (26a).  In particular this
half-turn lies in the matching kernel `V4^0` and therefore acts trivially on
the auxiliary matching quadratic field; the common label is not a field
inclusion.

## 3. The homogeneous leading drop is a universal `1+2` shadow

For the nonmonic depressed quartic

```text
f_A(T)=A T^4+pT^2+qT+r                                        (18)
```

translate (2) to the squared-pairing coordinate.  The integral polynomial is

```text
S_A(U)=U^3+2pU^2+(p^2-4Ar)U-Aq^2,                             (19)
Disc_U(S_A)=Disc_T(f_A).                                       (20)
```

On the leading-coefficient divisor,

```text
S_0(U)=U(U+p)^2.                                               (21)
```

For `p != 0`, this is one simple root plus one double root; for `p=0`, all
three roots coincide.  This `1+2` shadow is therefore universal coefficient
geometry, not a proved one-plus-two fibre law for a Keller map.

The exact discriminant expansion is

```text
Disc(f_A) = A [ 4p^3(4pr-q^2)
  + A(-128p^2r^2+144pq^2r-27q^4) + 256A^2r^3 ].               (22)
```

Consequently, at a reduced divisor `A=0`, the generic valuation is exactly
one whenever `p(4pr-q^2)` is a unit.  Higher valuation can occur only on

```text
p=0                   or                   q^2=4pr             (23)
```

in the residue field (or when `A` itself is nonreduced).  Formula (22) is a
useful residual-quadratic gate for HYP-9027, but no current Keller theorem
excludes that projection collision.

## 4. The cubic is a Galois-closure correspondence, not a source quotient

Let `L/K` be the Galois closure of the generic quartic field.  In the `S_4`
case,

```text
quartic root field       K_4 = L^{S_3},
cubic resolvent field    K_3 = L^{D_8},
resolvent splitting field K_6 = L^{V_4},  Gal(K_6/K)=S_3.      (24)
```

The subgroups `S_3` and `D_8` are incomparable and generate `S_4`, so

```text
K_3 is not a subfield of K_4,        K_3 intersect K_4 = K.     (25)
```

In the `A_4` case,

```text
K_4=L^{C_3},             K_3=L^{V_4},                          (26)
```

and the same incomparability/intersection statement holds.  Hence a dominant
rational map from the quartic source normalization to the cubic resolvent
normalization would induce a field inclusion forbidden by (25) or its
`A_4` analogue after (26).

The `D4` row has an additional distinction.  If `H=C2` is the quartic point
stabilizer and `V4^0=G intersect V4` is the normal matching kernel, then

```text
quartic root field K_4=L^H,
matching quadratic M_2=L^{V4^0},
root-field quadratic intermediate E_2=L^{N_D4(H)}.             (26a)
```

Here `H` is not contained in `V4^0`, the two order-four subgroups
`V4^0` and `N_D4(H)` are distinct, and `H,V4^0` generate `D4`.
Consequently `M_2` is not a subfield of `K_4` and is **not** the quadratic
intermediate `E_2` inside `K_4`.  The matching quadratic and the generic
deck quotient `N_D4(H)/H=C2` are also different constructions.

The distinction occurs in the smallest exact field-theory hostile:

```text
f(T)=T^4-2,                    S(W)=W(W^2+8).                (26b)
```

The quartic has transitive Galois group `D4`; its root field
`Q(2^(1/4))` contains the proper intermediate `Q(sqrt(2))`, whereas the
nontrivial matching factor gives `Q(sqrt(-2))`.  This is not a Keller
example.  It proves that the root-field intermediate and matching quadratic
cannot be identified from their shared degree or `C2` label.

The three functions in (3) combine four conjugate sheets; they are not
functions on one quartic source sheet.  `V_4` acts only after passage to `L`.

This is the exact lost-origin sidecar.  Even an abstract degree-three cover
with the right discriminant is not automatically a polynomial source, an
affine space, an etale/Keller map, or a chosen physical coordinate.

## 5. Tame inertia and the resolvent index tax

At a tame divisorial valuation, the permutation discriminant exponent is
`degree - number of inertia orbits`.  The quotient action on the three
pairings gives the complete table

| quartic inertia | quartic exponent `d_4` | pairing action | cubic exponent `d_3` |
|---|---:|---|---:|
| identity | 0 | identity | 0 |
| `(12)` | 1 | transposition | 1 |
| `(12)(34)` | 2 | identity | 0 |
| `(123)` | 2 | 3-cycle | 2 |
| `(1234)` | 3 | transposition | 1 |

Thus double transpositions are completely invisible after the `V_4`
quotient, and a 4-cycle retains only its order-two shadow.

Let `i_4,i_3` be the valuation lengths of the quartic and resolvent polynomial
orders inside their maximal orders.  Since order discriminant equals field
discriminant times index squared, (5) forces

```text
i_3-i_4 = (d_4-d_3)/2.                                        (27)
```

The tax is zero for a transposition or 3-cycle and exactly one for a double
transposition or 4-cycle.  Equality of the raw polynomial discriminants can
therefore record an index collision where the normalized cubic cover is
unramified.

An exact `S_4` hostile is

```text
f_{s,t}(T)=(T^2-1)^2-t(T+s),
R_{s,t}(W)=W^3+2W^2+(4st-4)W+(8st-t^2-8),                     (28)
Disc=-t^2(256s^3t-256s^2-288st+27t^2+256).
```

The specialization `(s,t)=(2,1)` is irreducible of type `4` modulo `2`
and type `1+3` modulo `3`; hence its Galois group, and therefore the generic
group in (28), is `S_4`.  At `t=0` and `s^2 != 1`, the two roots near `+1`
swap and the two near `-1` swap: inertia is `(12)(34)`.  Yet

```text
R_{s,0}=(W+2)^2(W-2),                                         (29)
```

The normalization is visible without appealing only to the inertia table.
Put `W=-2+tY`.  Exact division gives

```text
R_{s,t}(-2+tY)=t^2[tY^3-4Y^2+4sY-1].                         (29a)
```

At `t=0` the residual quadratic `4Y^2-4sY+1` has discriminant
`16(s^2-1)`, so it is separable when `s^2!=1`.  The normalized cubic
extension is therefore unramified there.  The raw double root in (29) is
precisely the one-unit nonmaximal-order tax in (27).

Two minimal controls are

```text
T^4-t             -> W(W^2+4t)                 (4-cycle tax),
(T^3-t)(T-1)      -> W^3-3tW-t-t^2             (visible 3-cycle). (30)
```

Their discriminants are respectively `-256t^3` and
`-27t^2(t-1)^2`.

## 6. The exact Keller/Jelonek survivor and stopping boundary

There is first a correction to the live monodromy list.  The cited Campbell
Galois criterion excludes the regular groups `C4,V4`.  THM-1375 also deleted
`D4` by treating `N_G(H)/H=C2` as a polynomial deck group, but THM-1365's
fixed-point argument explicitly assumes that generic deck transformations
extend polynomially across the Jelonek divisor.  No such extension is proved
for `D4`.  The honest currently proved non-Galois list is therefore

```text
D4, A4, S4,                                                   (30a)
```

not only `A4,S4` (MISTAKE-297).  This restores a quadratic intermediate to
the quartic root field and invalidates any unconditional “no intermediate
field” statement outside the `A4/S4` primitive branches.

> **SUPERSEDING UPDATE (THM-2633).**  The correction above remains the right
> repair of the old Smith argument, but it is no longer the final monodromy
> list.  THM-2633 proves that the point stabilizer of any affine-space Keller
> map surjects onto monodromy abelianization.  `D4` fails that gate, so the
> live degree-four list is now `A4,S4`.  The `D4` normalization and split-owner
> analysis below remains a correct conditional atlas on an empty Keller lane.

The missing extension arrow has an exact Zariski-main formulation.  Write a
Keller map as

```text
X=A^3 --j--> Xbar --Fbar--> Y=A^3,                            (30b)
```

where `j` is the canonical open immersion and `Fbar` is the finite
normalization of `Y` in `K(X)`.  Every generic deck automorphism `sigma` of
`K(X)/K(Y)` preserves the integral closure of `C[Y]`, so it extends uniquely
to a `Y`-automorphism `sigmabar` of `Xbar`.  It extends to a polynomial deck
automorphism of `X` **if and only if**

```text
sigmabar(j(X))=j(X),   equivalently
sigmabar(D_infty)=D_infty,    D_infty=Xbar\j(X).              (30c)
```

Let `R` be the ramification support of `Fbar`.  Since `F` is etale,
`R subseteq D_infty`; since `sigmabar` is over `Y`, the different (hence
`R`) is invariant.  Therefore `D_infty=R` would force (30c).  In the live
`D4` row, the resulting polynomial involution is impossible by the
polynomial-deck fixed-point argument of THM-1365 together with the
prime-order Smith fixed-point theorem.  Consequently any `D4` Keller
survivor must have an **unramified missing boundary divisor** in
`D_infty\R`, and its involution must fail to preserve the ownership of the
whole missing set.  Indeed, if every divisorial component of `D_infty` lay
in `R`, then `sigmabar` would preserve the boundary in codimension one.  On
the normal affine source, the induced rational self-map and its inverse would
then be regular off codimension two, so their coordinate functions extend;
they would give a polynomial deck automorphism after all.  Thus some
unramified boundary divisor `D_0` is exchanged with an unramified divisor
inside `j(X)`, and the two lie over the same target divisor because
`sigmabar` is over `Y`.  This included/missing sheet pair is the precise
remaining sidecar, not an abstract normalizer calculation.

Set-theoretically `B_4=Fbar(R)` and `A_F=Fbar(D_infty)`.  The obstruction in
(30c) lives upstairs: even equality `B_4=A_F` would not by itself exclude an
unramified missing sheet lying over the same target divisor as a ramified
one.

Equivalently, a surviving `D4` branch must exhibit a **split owner pair**.
For some irreducible unramified boundary divisor `D_0` and its target image
`C`, the involution gives another divisor `D_1=sigmabar(D_0)` with

```text
Fbar(D_0)=Fbar(D_1)=C,       D_0 subset D_infty,
D_1 subset j(X),             sigmabar: D_0 <--> D_1.           (30d)
```

At the generic point of `C`, these are two isomorphic unramified local
branches with opposite present/omitted ownership bits.  By (17b), they are
an antipodal pair in the canonical `C4`, not an arbitrary pair of quartic
sheets.  Thus the first decisive `D4` computation is not another
discriminant: it is a normalization and boundary-owner calculation testing
whether any antipodal deck `C2` pair is split as in (30d).  If no such pair
exists, the `D4` branch is excluded.

For a hypothetical degree-four Keller map, let `B_4` be the branch divisor
of the finite quartic normalization and `B_M` that of the nontrivial
normalized matching factor: degree two for `D4`, degree three for `A4/S4`.
The inertia table gives

```text
B_M subseteq B_4 subseteq A_F.                                 (31)
```

The first inclusion is strict in the double-transposition control (28): that
quartic branch component is invisible in the normalized cubic cover.  The
second inclusion is a logically different boundary loss: an open source may
omit a divisor above an unramified normalization point (the degree-one model
is an open immersion with a deleted target point).  No theorem says this
abstract strictness occurs in a Keller witness.  Thus (31) refines, but does
not replace, THM-2465's tower constraints.

In an `S_4` lane, `B_M` has `S_3` monodromy.  Under the same
Zariski--Lefschetz/Deligne--Fulton hypotheses used in THM-1330, a generic
plane section of `B_M` cannot be nodal-only: its complement group must
surject onto nonabelian `S_3`.  This worse-than-nodal selection rule for the
non-`V_4` branch sub-divisor is the sharp positive transfer.  In an `A_4`
lane the quotient is cyclic `C_3`, and in a `D4` lane it is `C2`, so this
nonabelian conclusion disappears in both branches.

No stronger grade-three transfer follows.  The normalized resolvent map is
finite (hence proper), ramified, and its affine normal source is not generally
an affine space; it is not a Keller map.  THM-2473/2546 describe the fixed sporadic map, not a
classification of all degree-three Keller maps.  In particular none of the
following is preserved by (2): a Jelonek leading factor, a physical
trace-zero coordinate, the `-4 square times L` factorization, the cusp
cylinder, a survivor/empty-fibre law, or sheet ownership at infinity.

The two-parameter control

```text
f=T^4+uT+v,               R=W^3-4vW-u^2,
Disc=256v^3-27u^4                                             (32)
```

makes the distinction literal: both monic root covers are finite over the
coefficient base, so both have empty nonproper-value sets, although their
common discriminant is cuspidal.  A discriminant cusp is not a Jelonek cusp.

The next valid G1 test is therefore to compute the maximal resolvent order,
its normalized branch divisor and index tax, and then its generic
plane-section singularities.  Testing only depression or the raw
discriminant identity cannot constrain G1.

## 7. The exact `2,3,4` interpretation (scoped analogy)

The even sign-flip group `V_4` acts simply transitively on the four sections
in (12).  Thus the sections form an affine `V_4` torsor with no distinguished
zero.  Choosing one quartic section supplies a reconstruction origin.  The
three nonidentity translations are the three perfect matchings of four
labels, and

```text
Aut(V_4)=GL_2(F_2)=S_3.                                       (33)
```

After choosing an origin, the full relabelling action is

```text
S_4 = AGL_2(F_2) = V_4 semidirect GL_2(F_2).                  (33a)
```

This is the exact algebraic reason that a four-sheet binary sign torsor has a
three-object permutation shadow.  Moreover `sign:S_4 -> C_2` kills `V_4`
and therefore factors through `S_4/V_4=S_3`; this is the representation-
theoretic reason the discriminant square class survives the resolvent.
THM-2596/2597 meet the same abstract group through
`PSL_2(Z) -> PSL_2(F_2)=S_3` (and its order-six bicycle frame), but there is
no canonical physical identification between that modular quotient and this
resolvent monodromy.

There is also an exact partial-cube refinement, developed systematically in
THM-2606.  A nonzero covector on the translation space of the affine
`V_4` torsor has a two-element nonzero level set; using those two directions
as Cayley generators makes the four sections a square `C4`.  Its omitted
direction is one of the three perfect matchings.  Conversely every connected
translation-invariant partial cube on four torsor points is one of these
three squares.  Since `GL_2(F_2)` is transitive on the three omitted
directions, full `AGL_2(F_2)=S4` monodromy selects none of them.  Thus a
binary square is another way to encode a chosen resolvent channel, not a
monodromy-invariant structure and not a new Keller constraint.

The same three objects are the `2+2` clusterings
`12|34,13|24,14|23` used by a four-body Jacobi tree.  The quotient remembers
the clustering and forgets particle/sheet ownership.  This is an exact
combinatorial identification only: no N-body dynamics, Hamiltonian structure,
or integrability statement transfers through it.

Likewise one may label the incircle/excircles by four sign choices, but
Feuerbach's internal-versus-external tangency distinguishes the incircle;
the geometric configuration is not a `V_4` translation symmetry.  These are
interpretive correspondences, not dependencies or Keller consequences.

Krenn--Gu--Soltesz, [*Questions on the Structure of Perfect Matchings inspired
by Quantum Physics*](https://arxiv.org/abs/1902.06023), use the decomposition
of `K4` into three disjoint perfect matchings as a four-vertex/three-colour
carrier.  That is exactly the carrier indexing (3).  Their amplitudes multiply
edge weights and permit cancellation between matchings, whereas (3) adds two
root products; no amplitude theorem transfers from the common carrier alone.

## 8. Relation to the existing frontier

- THM-2455 proves the monic depressed discriminant/cusp scaffold.  This
  theorem adds the general integral formula, explicit reconstruction torsor,
  homogeneous leading drop, residual valuation gate, and normalization loss.
- THM-2465's maximality/no-intermediate-field statement applies to both live
  branches `S4/A4`.  This theorem correctly restored `D4` after MISTAKE-297;
  THM-2633 later excludes it by a different gate.  The exact `D4` field
  lattice and distinction between matching and root-field quadratics remain
  valid conditional algebra.
- HYP-9027 remains open.  Formula (22) identifies the first residual
  collision its odd-valuation strategy must control; it does not control it.
- No G1 witness is constructed and no live `A4,S4` branch is excluded here.

## 9. Reproduction

```bash
python 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
python -O 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
```

Both modes must byte-match the stored transcript and end in
`FAILED CHECKS: NONE`.
