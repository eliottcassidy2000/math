---
id: THM-2598
title: "Quartic V4 resolvent torsor, universal cubic cusp, and the exact transfer boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The integral
  cubic resolvent of a general quartic has
  exactly the quartic polynomial discriminant; its depressed cubic and
  cube-minus-square cusp identity are universal invariant theory, not
  Keller constraints.  On the nondegenerate depressed open, reconstruction
  from the three squared pair sums is an exact four-point V4 torsor.  The
  homogeneous leading-coefficient boundary has a forced 1+2 shadow and an
  explicit residual-quadratic valuation gate.  The S4/A4 field lattice,
  tame inertia table, and local order-index tax prove that the cubic
  resolvent is a Galois-closure correspondence, not an intermediate cover
  of the quartic source.  Consequently no grade-three Keller anatomy
  descends without an additional affine/Keller realization and V4-origin
  sidecar.  No degree-four Keller map, G1, JC(2), or DC(2) is excluded.
source: codex-2026-07-27-quartic-resolvent-transfer
depends_on:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2546-integral-coordinate-dichotomy-and-parity-lens-scope
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
output: 05-knowledge/results/jacobian_quartic_v4_resolvent_thm2598.out
script_sha256: c198b1b969f3284b259771bfdd612ceb019add3d7fa9c7ea5bfe1375492499b0
output_sha256: 1c1a899eca2062f8fc29066e410679544e9f306823f45e27bdf518381c841a3a
hash_basis: working-tree bytes (LF)
---

# THM-2598 -- what the quartic resolvent remembers, and what it cannot carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  All algebraic
identities and finite permutation claims below
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
`16(s^2-1)`, so it is separable when `s^2 != 1`.  The normalized cubic
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

For a hypothetical degree-four Keller map, let `B_4` be the branch divisor
of the finite quartic normalization and `B_3` that of the normalized cubic
resolvent correspondence.  The inertia table gives

```text
B_3 subseteq B_4 subseteq A_F,                                 (31)
```

The first inclusion is strict in the double-transposition control (28): that
quartic branch component is invisible in the normalized cubic cover.  The
second inclusion is a logically different boundary loss: an open source may
omit a divisor above an unramified normalization point (the degree-one model
is an open immersion with a deleted target point).  No theorem says this
abstract strictness occurs in a Keller witness.  Thus (31) refines, but does
not replace, THM-2465's tower constraints.

In an `S_4` lane, `B_3` has `S_3` monodromy.  Under the same
Zariski--Lefschetz/Deligne--Fulton hypotheses used in THM-1330, a generic
plane section of `B_3` cannot be nodal-only: its complement group must
surject onto nonabelian `S_3`.  This worse-than-nodal selection rule for the
non-`V_4` branch sub-divisor is the sharp positive transfer.  In an `A_4`
lane the quotient is cyclic `C_3`, so this nonabelian conclusion disappears.

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

After one reconstruction section is chosen, the even sign-flip group `V_4`
acts simply transitively on the four sections in (12).  Thus the sections
form an **affine** `V_4` torsor: they have no distinguished zero until a
quartic sheet/origin is supplied.  Its three nonidentity translations are
the three perfect matchings of four labels, and

```text
Aut(V_4)=GL_2(F_2)=S_3.                                       (33)
```

After choosing an origin this is the standard affine decomposition

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

The same three objects are the `2+2` clusterings
`12|34,13|24,14|23` used by a four-body Jacobi tree.  The quotient remembers
the clustering and forgets particle/sheet ownership.  This is an exact
combinatorial identification only: no N-body dynamics, Hamiltonian structure,
or integrability statement transfers through it.

Likewise one may label the incircle/excircles by four sign choices, but
Feuerbach's internal-versus-external tangency distinguishes the incircle;
the geometric configuration is not a `V_4` translation symmetry.  These are
interpretive correspondences, not dependencies or Keller consequences.

## 8. Relation to the existing frontier

- THM-2455 proves the monic depressed discriminant/cusp scaffold.  This
  theorem adds the general integral formula, explicit reconstruction torsor,
  homogeneous leading drop, residual valuation gate, and normalization loss.
- THM-2465 already records maximality/no-intermediate-field and the broad
  `S_4/A_4` tower constraints.  The new content here is the exact field
  lattice comparison, inertia-to-index ledger, and the proof that the cubic
  is not a cover of the quartic source.
- HYP-9027 remains open.  Formula (22) identifies the first residual
  collision its odd-valuation strategy must control; it does not control it.
- No degree-four Keller witness is constructed or excluded.

## 9. Reproduction

```bash
python 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
python -O 04-computation/jacobian_quartic_v4_resolvent_thm2598.py
```

Both modes must byte-match the stored transcript and end in
`FAILED CHECKS: NONE`.
