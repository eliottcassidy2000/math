---
id: THM-2864
title: "Quartic edge/orientation sextic resolvents and the D8 radicand product"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A depressed
  quartic has two explicit sextic lifts of its cubic
  matching resolvent.  The edge lift has generic stabilizer the nonnormal
  edge V4 and square discriminant 64 q^2 Delta^2.  The oriented-four-cycle
  lift has generic stabilizer C4 and discriminant
  2^66 q^12 J_or^4 Delta^3, so it retains the quartic discriminant square
  class away from its separate collision wall J_or=0.  Over one matching,
  its radicand F(u), the edge radicand u, and the quartic discriminant obey
  C(u)^2 u F(u)=16 q^2 Delta.  This realizes the three D8 characters of
  THM-2862 by explicit polynomial formulas.  It is not a graph-quartic
  Keller realization or a Jacobian-conjecture exclusion.
source: root/quartic-sextic-resolvent-triangle-2026-07-28
audit: >
  thm2864-formula-audit-2026-07-28 (independent matching/edge/orientation
  derivation, elementary-symmetric and discriminant factorization, D8
  radicand, stabilizer, specialization, and sharp-boundary audit: ACCEPT)
depends_on:
  - THM-2862-modular-level-three-four-congruence-ladder-and-inequivalent-six-lifts
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-2862-modular-level-three-four-congruence-ladder-and-inequivalent-six-lifts
script: 04-computation/quartic_sextic_resolvent_triangle_thm2864.py
output: 05-knowledge/results/quartic_sextic_resolvent_triangle_thm2864.out
script_sha256: fb499ea79bd806e2e00ec94d31d2f5232dd3dfe7c17afb4b60100f1353be471c
output_sha256: 2ba7da5a16497e6d4c92766127e0c84a8258219bae8aa457c66b11e87c00aeb4
hash_basis: LF-normalized bytes
---

# THM-2864 -- two sextic lifts of the quartic cubic resolvent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2598 constructs the cubic matching resolvent of a quartic and proves
that it has the quartic discriminant.  THM-2753 shows that the six-edge
action erases ambient parity.  THM-2862 then finds another six-point lift,
the oriented-four-cycle action, which retains parity and identifies the
three abstract quadratic characters over one matching.

The present theorem makes those two six-point lifts polynomial.  More
importantly, it finds the exact radicand identity relating them.  The
identity is universal.  The field and stabilizer interpretations require
the separability hypotheses stated below.

Work first over a characteristic-zero field `K`.  Let

```text
f(X)=X^4+pX^2+qX+r=prod_(i=1)^4 (X-x_i),       sum_i x_i=0,    (1)
```

and write

```text
Delta=256r^3-128p^2r^2+144pq^2r-27q^4
      +16p^4r-4p^3q^2.                                      (2)
```

All displayed coefficient identities lie in `Z[p,q,r]`.  Thus they remain
polynomial identities in every characteristic, but statements involving
sign, square classes, separability, or the `D8` subgroup lattice are made
only in characteristic zero (and have the usual characteristic-not-two
versions).

## 1. The matching cubic and the six edges

Fix a matching

```text
M={{i,k},{j,l}}.
```

Put

```text
s=x_i+x_k=-(x_j+x_l),                    u=s^2.              (3)
```

Writing

```text
f(X)=(X^2-sX+a)(X^2+sX+b)                                (4)
```

gives

```text
p=a+b-u,                    q=s(a-b),                    r=ab. (5)
```

Eliminating `a,b` from `(5)` gives the matching cubic

```text
S(U)=U^3+2pU^2+(p^2-4r)U-q^2.                              (6)
```

Its three roots are the squared sums attached to the three matchings.
This is the depressed-coordinate resolvent of THM-2598, and

```text
Disc_U(S)=Delta.                                             (7)
```

Remembering which of the two opposite edges in a matching was selected
amounts to choosing the sign of `s`.  Therefore the edge sextic is simply

```text
E(Y)=S(Y^2)
    =Y^6+2pY^4+(p^2-4r)Y^2-q^2.                             (8)
```

If the roots of `S` are `u_1,u_2,u_3`, the roots of `E` are
`+-sqrt(u_i)`.  Multiplying the three within-pair differences and the
four cross-pair differences for each pair of indices gives

```text
Disc_Y(E)
 =4^3 (u_1u_2u_3) prod_(i<j)(u_i-u_j)^4
 =64 q^2 Delta^2.                                           (9)
```

Consequently `q Delta !=0` is exactly the separability condition for this
particular edge generator.  On the universal labelled-root field, or on
an `S4` specialization satisfying it, a root `s=x_i+x_k` has stabilizer

```text
Stab(s)=< (i k),(j l) > ~= V4_edge,                         (10)
```

the nonnormal Klein four that preserves the selected edge and its
opposite setwise.  This is the edge stabilizer of THM-2862, not the normal
translation `V4` which is the kernel of the three-matching quotient.

## 2. The oriented-four-cycle sextic

An oriented four-cycle is a cyclic ordering

```text
gamma=(i j k l)
```

modulo cyclic rotation.  There are six such orderings; inversion pairs
them over the three matchings.  Define

```text
d_1=x_i-x_k,                  d_2=x_j-x_l,
Omega_gamma=d_1 d_2 (d_1^2-d_2^2).                         (11)
```

Rotation sends `(d_1,d_2)` to `(d_2,-d_1)` and fixes `Omega_gamma`.
Reversal fixes the matching and negates `Omega_gamma`.  Thus `(11)` has
the correct `C4`-versus-reflection transformation law without choosing
root square signs by hand.

Equations `(4)--(5)` give

```text
d_1^2+d_2^2 = -4p-2u,
d_1^2d_2^2  = 16r-4pu-3u^2,
(d_1^2-d_2^2)^2
              =16(u^2+2pu+p^2-4r).                         (12)
```

Hence the squared orientation value over a matching is

```text
F(U)=16(16r-4pU-3U^2)(U^2+2pU+p^2-4r),                    (13)
Omega_gamma^2=F(u).                                        (14)
```

Taking the norm from the cubic algebra gives the orientation sextic

```text
O(Y)=Res_U(S(U),Y^2-F(U))
    =prod_(S(u_i)=0)(Y^2-F(u_i)).                           (15)
```

There is a short coefficient proof which also exposes the collision wall.
Exact reduction in `Z[p,q,r][U]/(S)` gives

```text
F(U)=256rU^2+(512pr-48q^2)U
     +(256p^2r-64pq^2-1024r^2).                            (15a)
```

The trace, second elementary symmetric function, and norm of `(15a)` are

```text
Tr(F)=32(8p^2r-3pq^2-32r^2),
e_2(F)=-256q^2(p^2+12r)(32pr-9q^2),
Norm(F)=4096q^4Delta.                                      (15b)
```

Consequently `(15)` expands to

```text
O(Y)=Y^6
 -32(8p^2r-3pq^2-32r^2)Y^4
 -256q^2(p^2+12r)(32pr-9q^2)Y^2
 -4096q^4 Delta.                                           (16)
```

Define the orientation-collision polynomial

```text
J_or=1024r^3+768p^2r^2-288pq^2r+27q^4.                    (17)
```

This `J_or` is not the binary-quartic invariant named `J` in THM-2598.
For `{i,j,k}={1,2,3}`, subtracting the three values of `(15a)` and using
the Vieta relations for `S` gives

```text
F(u_i)-F(u_j)
 =-16(u_i-u_j)(16r u_k+3q^2),
prod_k(16r u_k+3q^2)=q^2 J_or.                             (17a)
```

It follows first that

```text
Disc_Z prod_i(Z-F(u_i))=2^24q^4J_or^2Delta,                (17b)
```

and then, after adjoining both signs of every square root, that the exact
discriminant of `(16)` is

```text
Disc_Y(O)=2^66 q^12 J_or^4 Delta^3
         =Delta (2^33 q^6 J_or^2 Delta)^2.                 (18)
```

Thus `q Delta J_or !=0` makes all six orientation values distinct.  Their
generic stabilizer is the cyclic group generated by `gamma`:

```text
Stab(Omega_gamma)=<gamma> ~= C4.                            (19)
```

One direct generic witness is

```text
f=X(X-1)(X-2)(X+3)=X^4-7X^2+6X.                            (20)
```

Its six edge values are `+-1,+-2,+-3`, while its six orientation values
are

```text
+-24, +-96, +-120.                                         (21)
```

So no extra labelled permutation fixes `(11)` identically.  Equations
`(18)` and `(21)` together prove the generic orbit and stabilizer claim.
For a specialized quartic with smaller monodromy, the sextic can of course
factor even when it remains separable.

## 3. The exact `D8` radicand product

For the same matching define its cross-difference product

```text
C(U)=(x_i-x_j)(x_i-x_l)(x_k-x_j)(x_k-x_l)
    =3U^2+4pU+p^2-4r.                                      (22)
```

In fact `C(U)=S'(U)`, so

```text
Norm_(K[U]/(S)/K)(C(u))=-Delta.                            (22a)
```

This makes the role of quartic separability in the radicand identity
visible before any root calculation.

There are two short proofs of the central identity.  First, `(5)` gives

```text
s(d_1^2-d_2^2)=-4q.                                        (23)
```

Second, the quartic Vandermonde factors into the two within-pair
differences and the four cross differences:

```text
Delta=d_1^2d_2^2 C(u)^2.                                   (24)
```

Multiplying `(23)^2` by `(24)` and using `(11)--(14)` yields

```text
C(u)^2 u F(u)=16q^2 Delta.                                 (25)
```

Equivalently, direct polynomial division gives the integral quotient
identity

```text
C(U)^2 U F(U)-16q^2Delta = 0             in Z[p,q,r][U]/(S). (26)
```

Assume `q Delta !=0`, and let `A=K[U]/(S)`, interpreted factor by factor
if the cubic is reducible.  Then `u,C(u),F(u)` are units in `A`, and `(25)`
gives the square-class law

```text
[u] [F(u)] = [Delta]                    in A^times/(A^times)^2. (27)
```

In the universal `S4` splitting field, fix the matching `M`.  Its
stabilizer is `D8`, and the three quadratic lifts of the cubic matching
field are obtained by adjoining

```text
sqrt(u)=s                    edge choice, stabilizer V4_edge;
sqrt(F(u))=Omega_gamma       orientation, stabilizer C4;
sqrt(Delta)                  quartic sign, stabilizer V4_normal. (28)
```

Equation `(27)` is therefore the coefficient-level realization of
THM-2862's character identity

```text
chi_disc=chi_edge chi_or.                                  (29)
```

It is important that `(27)` is a product law, not an independence claim.
On special monodromy strata one or more of the three square classes can
collapse.  Likewise `J_or=0` can make the chosen orientation generator
nonprimitive even though the abstract subgroup field still exists in the
Galois closure.

## 4. The two sextics retain different parity information

Equation `(9)` says that the edge sextic discriminant is a square whenever
it is nonzero:

```text
Disc(E)=(8qDelta)^2.                                        (30)
```

This is the polynomial shadow of THM-2753's ambient six-edge parity
erasure.  In contrast, `(18)` gives

```text
[Disc(O)]=[Delta]                         in K*/K*2          (31)
```

away from the displayed zero divisors.  The orientation sextic therefore
retains the quartic sign character, exactly as the projective/four-cycle
action in THM-2862 does.

Neither statement follows from the number `6` or from the common
three-matching quotient.  The square class changes only because one sextic
uses the nonnormal `V4_edge` lift and the other uses the `C4` lift.

## 5. Sharp collision boundaries

There are three distinct degeneration mechanisms.

1. `Delta=0`: the quartic and matching cubic are inseparable, so both
   six-object constructions collide.
2. `q=0`: one matching radicand `u` is zero and the free edge-sign torsor
   degenerates.  The edge sextic has a repeated root; the orientation
   sextic also degenerates, though not necessarily on the same matching.
3. `J_or=0` with `q Delta !=0`: the edge sextic remains separable but
   the orientation generator collides.

The third boundary is genuine.  At

```text
(p,q,r)=(1,4,-3),                       Delta=-22000,
J_or=0,                                                     (32)
```

one has

```text
E(Y)=(Y-1)(Y+1)(Y^4+3Y^2+16),                              (33)
O(Y)=(Y^2-1280)^2(Y^2+14080).                              (34)
```

Thus `E` is squarefree, whereas `O` has

```text
gcd(O,O')=Y^2-1280.                                        (35)
```

This prevents replacing the orientation hypothesis by quartic
separability alone.  It also explains the fourth power of `J_or` in
`(18)`: `J_or` is a primitive-generator collision wall, not a second copy
of the quartic discriminant.

## 6. What transfers to the Keller frontier

The theorem supplies a concrete replacement for the vague instruction
"use the quartic resolvent."  A monic depressed graph quartic now has three
explicit objects over one matching cubic:

```text
matching cubic S,
edge sextic E with square discriminant,
orientation sextic O with quartic discriminant class.       (36)
```

The useful possible transfer is local and conditional.  If a graph-quartic
Keller chart realizes `(1)` without losing its leading-coefficient,
normalization, or boundary-owner data, then `(18)` separates the odd
quartic discriminant divisor from the even `q` and `J_or` collision walls.
The orientation maximal order may therefore retain a parity/different
sidecar which the edge maximal order cannot see.

None of those affine hypotheses is automatic:

- a general graph quartic is not already monic and depressed;
- depressing it can introduce leading-coefficient denominators;
- a raw polynomial discriminant can differ from the maximal-order
  discriminant by an index square;
- `q=0` and `J_or=0` are coordinate/generator boundaries until their
  graph-chart meaning is proved;
- the cubic, edge-sextic, and orientation-sextic fields live in the
  Galois closure and are not automatically functions on one source sheet;
- no Jelonek owner, polynomial deck map, or Keller source is constructed.

The next exact Keller test is therefore to homogenize `(16)--(18)` for

```text
A X^4+B X^3+C X^2+D X+E
```

and compute the leading-divisor valuations, order indices, and normalized
different of the orientation sextic.  Merely observing the square-class
identity `(31)` does not exclude either live `A4` or `S4` monodromy.

There is likewise no LRC consequence.  The three quadratic characters in
`(28)` are algebraic sheet labels, not physical owner currents, endpoint
phases, or carry coordinates.

## 7. Exact companion and status ledger

Run

```text
python 04-computation/quartic_sextic_resolvent_triangle_thm2864.py
python -O 04-computation/quartic_sextic_resolvent_triangle_thm2864.py
```

Both modes must LF-normalized byte-match

```text
05-knowledge/results/quartic_sextic_resolvent_triangle_thm2864.out.
```

The standard-library companion implements sparse multivariate polynomials
over `Z`.  It checks the matching elimination, both sextics, the resultant
norm, both discriminants, the quotient identity `(26)`, all six root-level
edge/orientation values, and the generic and boundary examples.  It also
enumerates all `24` root relabellings and verifies the exact
`V4_edge/D8/C4` stabilizer triangle.  It uses explicit exceptions and no
truth-bearing Python `assert`, floating-point decision, optional CAS, or
scratch dependency.

```text
PROVED IN THE CANDIDATE: matching cubic and edge sextic formulas;
                         orientation invariant and sextic formula;
                         exact edge and orientation discriminants;
                         separate J_or collision wall;
                         D8 radicand and square-class product;
                         generic V4_edge/C4 stabilizers;
                         sharp q, Delta, and J_or boundaries.

NOT PROVED:              nonmonic/homogeneous graph-quartic realization;
                         maximal-order divisor or different comparison;
                         a Jelonek or polynomial-source lift;
                         exclusion of A4 or S4 Keller monodromy;
                         a physical LRC carrier;
                         JC(2), DC(2), G1, or LRC(14).             (37)
```

An independent hostile audit rederived `(6)--(9)` from the quadratic
factorization, reconstructed `(12)--(18)` from the elementary symmetric
functions of `F(u_i)`, and proved `(25)` both from the Vandermonde and by
quotient-ring reduction.  It separately checked the `V4_edge`, `D8`, and
`C4` stabilizers, the smaller-monodromy factorization boundary, and all three
collision loci.  It also reproduced the exact wall example `(32)--(35)`.
Independent normal and optimized runs LF-normalized byte-match the stored
transcript, both declared hashes match, and the documentation gate passes.

**QED.**
