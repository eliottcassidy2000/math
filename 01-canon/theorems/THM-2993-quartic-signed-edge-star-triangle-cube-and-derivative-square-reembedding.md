---
id: THM-2993
title: "Quartic signed-edge star/triangle cube and derivative-square re-embedding"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a separable
  depressed quartic, the eight simultaneous choices
  of its three signed matching roots split into four coherent root-stars and
  four complementary triangles.  Their block-discriminant products are the
  roots of two explicit reciprocal quartics A_star and A_triangle.  The star
  quartic is the derivative-square resultant and, off q*K*Delta=0, re-embeds
  the same unpointed quartic etale algebra.  This is not a quartic Keller map,
  affine-origin section, or integral owner transfer.
source: codex-quartic-c3-parity-atlas-2026-07-30
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary
  - THM-2987-binary-clock-ledger-gl2-orbit-rank-and-first-v4-corner
script: 04-computation/quartic_signed_edge_star_triangle_reembedding_thm2993.py
output: 05-knowledge/results/quartic_signed_edge_star_triangle_reembedding_thm2993.out
script_sha256: 9ad6834d99c9d0ff083a021322e2a44ac8ccad5c1948d016835ff739be4773b0
output_sha256: 60d38b0a78ab519d7e8f3cb3b78f2aa0bb6985369ad1c1a0c713ff48d02f45b5
hash_basis: LF-normalized bytes
---

# THM-2993 -- the three signed edges assemble into a star/triangle cube

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and main statement

[THM-2992](THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary.md)
decodes the intrinsic fixed two-sheet block on one signed matching edge, but
does not assemble the three matchings.  [THM-2968](THM-2968-quartic-edge-and-oriented-cycle-s4-complements.md)
identifies the physical signed-pair image as the even-flip complement
`H0 < F2^3 semidirect S3`.  [THM-2595](THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go.md)
proves that the binary and ternary free factors are separately gauge-trivial
while their common affine origin can still be obstructed.

The simultaneous object is exact.  Let `k` have characteristic zero and let

```text
f(T)=T^4+pT^2+qT+r                                      (1)
```

be separable.  Put

```text
S(U)=U^3+aU^2+bU+c,
a=2p,                 b=p^2-4r,             c=-q^2,     (2)
D=Disc_U(S)=Disc_T(f).
```

For labelled roots `x0,x1,x2,x3` with sum zero, order the matchings as

```text
01|23,                  02|13,                  03|12.  (3)
```

The three coherent pair sums based at `x0` are

```text
s1=x0+x1,              s2=x0+x2,              s3=x0+x3,
si^2 a root of S,       s1*s2*s3=-q.                        (4)
```

When `q!=0`, every `si` is nonzero.  Choose independently `si` or `-si`.
The resulting cube `F2^3` has two intrinsic halves:

```text
even parity:  000,011,101,110 = four root-stars;
odd parity:   111,100,010,001 = four complementary triangles. (5)
```

For a signed edge sum `s`, THM-2992's selected quadratic block has

```text
d(s)=-2p-s^2-2q/s.                                     (6)
```

For the coherent star based at a root `x`, the three selected block
discriminants are `(x-y)^2`, `y!=x`, and therefore

```text
D_star(x)=product_(y!=x)(x-y)^2=f'(x)^2.                (7)
```

The complementary triangle uses the other three edges, hence

```text
D_triangle(x)=D/f'(x)^2.                                (8)
```

Thus all three signed matchings produce two canonical **unpointed** quartic
covariants.  They do not choose one of their four roots.

## 2. The two quartic covariants and every wall

Define

```text
E=4a^3c-3a^2b^2+18abc+4b^3-135c^2,
K=b^3-a^3c=Res_U(S,aU+b).                               (9)
```

Then the star and triangle products are the roots of

```text
A_star(Z)=Z^4+4(ab+7c)Z^3-2E Z^2
                 +4(ab-9c)D Z+D^2,                     (10)

A_triangle(Z)=Z^4+4(ab-9c)Z^3-2E Z^2
                 +4(ab+7c)D Z+D^2.                     (11)
```

They obey

```text
A_star(Z)=Res_T(f(T),Z-f'(T)^2),                        (12)
A_triangle(Z)=Z^4 A_star(D/Z)/D^2,                      (13)
A_star-A_triangle=64c Z(Z^2-D).                         (14)
```

The coefficient of `Z` in `(10)` also has the compact quartic form

```text
D(8p^3-32pr+36q^2).                                    (15)
```

Their individual discriminants are identical:

```text
Disc_Z(A_star)=Disc_Z(A_triangle)=2^24 c^2 K^2 D^3.     (16)
```

Consequently each factor is a separable four-sheet cover exactly on the open
`cKD!=0`.  This is not yet an eight-sheet etale statement.  Put

```text
H=a^6c^2-2a^4bc^2-2a^3b^3c-26a^3c^3+29a^2b^2c^2
  -2ab^4c-18abc^3+b^6-26b^3c^2+189c^4.                (17)
```

Then

```text
Res_Z(A_star,A_triangle)=2^32 c^4 D^4 H,                (18)
Res_Z(A_star,Z^2-D)=2^8 D^2 H.                         (19)
```

Off `cD=0`, `H=0` says exactly that one star product collides with its own
complementary triangle.  Off-diagonal star/triangle collisions force `c=0`.
Thus the full eight-value decoder is etale only on

```text
c*K*D*H != 0.                                           (20)
```

Equations `(16)--(20)` keep four logically different failures visible:

```text
c=0: signed-product/parity collapse;
K=0: collision inside one four-sheet factor;
D=0: original quartic branch collision;
H=0: a star meets its own reciprocal triangle.          (21)
```

## 3. Why the cube has no canonical origin

Use bit `0` for the edge in `(3)` incident with vertex `0`.  A sheet
permutation sends a bit vector by an affine signed-pair map

```text
epsilon |-> delta+sigma(epsilon),
delta in F2^3,             sigma in S3.                 (22)
```

Exact enumeration gives

```text
image(S4)=H0={(delta,sigma): |delta| is even},           (23)
```

of order `24`.  It preserves each half in `(5)` and acts regularly through
its normal `V4` on that half.  Translation by the central word `111` exchanges
the two halves but is not a physical sheet permutation.  The four subgroups
splitting `H0 -> S3` are precisely the four point stabilizers; equivalently,
they are the four choices of one root-star as affine origin.

Therefore `(10)` is the same unpointed `S4` cover in a new coordinate, not a
canonical affine-origin section.  The cubic coefficients distinguish the
star **half** from the triangle half when `q!=0`, because `(14)` contains
`c=-q^2`; they do not point one star inside that half.  Replacing

```text
(T,q) by (-T,-q)                                        (24)
```

fixes `S,D,A_star,A_triangle` and reverses the missing quartic coordinate.

## 4. Exact derivative-square re-embedding

In the quartic algebra `B=k[T]/(f)`, let `z=f'(T)^2`.  Reduction modulo `f`
gives

```text
z=-8qT^3+4(p^2-4r)T^2+4pqT+q^2.                        (25)
```

The basis transition is

```text
det(1,z,z^2,z^3 relative to 1,T,T^2,T^3)
       =2^12 q^2 K D.                                  (26)
```

Thus on `qKD!=0`, the map `Z |-> z` is an isomorphism

```text
k[Z]/(A_star)  ~=  k[T]/(f).                            (27)
```

This is stronger than equality of splitting fields.  An exact inverse is
also polynomial.  With `P_i=P_i(a,b,c)` given below,

```text
32 q K D T=P0+P1 z+P2 z^2+P3 z^3,                      (28)
```

where

```text
P0=(4a^3c-a^2b^2-18abc+4b^3+27c^2)
   (8a^3bc-a^2b^3-12a^2c^2+2ab^2c-4b^4-9bc^2),

P1=40a^5c^2-14a^4b^2c+3a^3b^4+64a^3bc^2-87a^2b^3c
   -378a^2c^3+4ab^5+117ab^2c^2+100b^4c-297bc^3,

P2=-(4a^3bc-3a^2b^3+56a^2c^2+2ab^2c-4b^4+57bc^2),
P3=-(2a^2c-ab^2+3bc).                                  (29)
```

Only the sign of `q=sqrt(-c)` is absent from the cubic coefficient field;
choosing it chooses between `T` and `-T` in `(28)`.  More importantly for
integral geometry, `(26)` contains a full factor `D`.  Indeed

```text
Disc(A_star)=(2^12 q^2 K D)^2 Disc(f).                  (30)
```

Along a simple transposition divisor with `qK` a unit, the original quartic
order has discriminant valuation `1`, the derivative-square order valuation
`3`, and their index has valuation `1`.  The generic isomorphism `(27)` is
therefore non-unimodular exactly at the branch boundary.  It cannot transport
a Keller/Jelonek owner sheet through that divisor.

## 5. The `C2*C3` compatibility and the degree-four cusp word

Fix one affine involution lift `s` and one affine order-three lift `t` above
the standard matching quotient.  The involution fixes a two-point affine
line; the ternary element fixes one star.  Across all `48` such lift pairs,

```text
fixed star in fixed line
  <=> <s,t> is a point-stabilizer S3
  <=> order(st)=2,                                      (31)

fixed star outside fixed line
  <=> <s,t>=S4
  <=> order(st)=4.                                      (32)
```

There are `24` pairs of each kind.  With the convention that `st` means first
`t`, then `s`, exact controls are

```text
s=(01),       t_split=(013),       order(st_split)=2,
s=(01),       t_full =(123),       order(st_full) =4,    (33)
```

and the two ternary elements induce the same quotient `3`-cycle.  Hence the
binary and ternary faces are each translation-gauge trivial, while their
product/cusp word detects the compatibility class.  The order `4=|V4|` in
the transverse class is the precise finite shadow of the user's proposed
`C2*C3`/forbidden-degree mechanism.

This is a physical realization of THM-2595's unique `H^1(C2*C3,V4)` bit, not
a new invariant beyond it.  In particular, the quotient actions and the two
separate local gauges still do not determine the common origin.

## 6. Grade-three branch-character corollary

Suppose the matching cubic `(2)`--or another cubic conditionally identified
with that matching cubic--has the THM-2473 discriminant form

```text
D=-4W^2 L.                                               (34)
```

Then `(16)` becomes

```text
Disc(A_star)=-2^30 (c K W^3)^2 L^3.                     (35)
```

On `c*K*W!=0`, the local reduced branch divisor is still `L`.  Globally the
raw discriminant support is `V(W) union V(L)`, while the even multiplicity of
`W` means that the odd-multiplicity discriminant character remains `[-L]`,
the same square class as for `D`.  The exponent of `L` has risen from one to
three because the derivative-square basis index contributes `D^2` in `(30)`.
Thus the grade-three discriminant anatomy transfers as an odd branch
character while its raw support and integral owner/Jelonek sheet do not.

There is one useful **conditional coefficient comparison**, not a physical
quartic construction.  THM-2473's core cubic is

```text
L0 X^3+B X-2c0,                 B=4-3b0c0.              (36)
```

On `L0!=0`, if one additionally identifies its monic normalization with the
matching cubic `(2)`, then

```text
a=0,       b=B/L0,       c=-2c0/L0,       K=B^3/L0^3.   (37)
```

So the derivative-square `K` wall is exactly the old `B=0` coefficient wall;
the exceptional empty curve `L0=B=0` is where the normalization itself
fails.  No such physical resolvent identification is proved here, and
[THM-2681](THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion.md)
warns that a normalized `S3` cubic need not carry the required quartic `V4`
kernel globally.

Writing THM-2473's other coefficient as

```text
Q=27a0c0^2-9b0c0+8,
```

its monic discriminant is `D=-4Q^2/L0^3`.  Under the same conditional
coefficient identification, `(26)` and `(16)` pull back exactly to

```text
det(1,z,z^2,z^3)=-2^15 c0 B^3 Q^2/L0^7,
Disc(A_star)=-2^32 c0^2 B^6 Q^6/L0^17.                 (38)
```

At the actual target `(a0,b0,c0)=(0,4/3,1)`, one has

```text
L0=16/27,       Q=-4,       B=0,
monic cubic X^3-27/8,
A_star=(8Z-243)^3(8Z-27)/4096.                          (39)
```

The cubic is separable while the derivative-square quartic has a triple
collision.  Thus the canonical generic re-embedding introduces a genuine
divisorial `B=0` degeneration.  This is compatible with, rather than a
counterexample to, THM-2681: on its full normalization the corresponding
coefficient has the form `B=rs`, so the `r=0` or `s=0` boundary is precisely
where a connected quasi-etale `V4` torsor cannot be extended.

For this depressed coefficient lane `a=0`, the cross wall has the sharp real
form

```text
H=(b^3-13c^2)^2+20c^4.                                 (40)
```

Hence `H>0` for real `c!=0`; its nontrivial roots are genuinely complex
(`b^3/c^2=13+-2i sqrt(5)`).

## 7. Sharp controls and first failed implications

### 7.1 Split positive control

For roots `(-6,1,2,3)`, one has

```text
(p,q,r)=(-25,60,-36),       D=1016064,
star products     =64,196,324,254016,
triangle products =4,3136,5184,15876.                  (41)
```

The star based at `-6` has signed sums `(-5,-4,-3)`, block discriminants
`(49,64,81)`, and product `254016=f'(-6)^2`.

### 7.2 The `K` wall is load-bearing

For

```text
(p,q,r)=(1,1,3/4),             D=125,             K=0,
```

```text
A_star=(Z-25)^2(Z^2+6Z+25).                              (42)
```

Thus `qD!=0` does not make the star cover separable.

### 7.3 The cross wall is independent

For

```text
f=T^4-7T^2+7T,       D=2401,       K=-16807,       H=0,
```

both four-sheet factors are separable but

```text
gcd(A_star,A_triangle)=Z-49.                             (43)
```

Thus `qKD!=0` does not make the full eight-value decoder etale.

### 7.4 The `q=0` and sign hostiles

For `f=(T^2-1)(T^2-4)`, `q=0` while `KD!=0`; nevertheless

```text
A_star=A_triangle,                 det(1,z,z^2,z^3)=0.   (44)
```

One pair sum is zero, so `(6)` also loses its signed denominator.  Separately,
`f_+(T)=T^4-7T^2+6T` and `f_-(T)=f_+(-T)` have the same matching cubic and
the same two covariants but opposite quartic coordinates.  Positivity or a
real coefficient presentation does not restore the missing origin.

The first failed stronger implications are therefore

```text
three matching roots + signed parity half
    does not imply a pointed V4 origin;

generic quartic-algebra isomorphism
    does not imply an integral branch-owner isomorphism;

same binary and ternary quotient gauges
    does not imply a common affine gauge.                (45)
```

## 8. Connection contract and exact evidence

The exact bridge is

```text
source:     the three signed matching roots of one depressed quartic;
map:        choose one edge per matching, multiply its block discriminants,
            then pass to z=f'(T)^2;
target:     A_star and its discriminant-reciprocal A_triangle;
preserved:  the unpointed four-sheet S4 cover, matching quotient,
            discriminant character, and C2/C3 compatibility bit;
destroyed:  q-sign/T-orientation, a pointed V4 origin, integral branch index,
            present/omitted owner, and global Keller realization;
sidecars:   one affine-origin section and one integral owner/Jelonek section;
cheapest hostile tests:
            K=0 at (1,1,3/4), H=0 at (-7,7,0), and q=0 at (-5,0,4).
```

The standalone companion verifies `(2)--(44)` over exact polynomial rings.
It enumerates all `24` physical sheet permutations, the complete signed cube,
all four quotient splittings, and all `48` binary/ternary lift pairs.  It
checks the direct resultant, reciprocal factor, both discriminants, the
cross-resultant and own-complement factor, the derivative-square basis
determinant, the explicit cubic inverse `(28)--(29)`, the grade-three formula,
and all stated hostiles.  No floating-point decision or optimization-sensitive
truth gate is used.

Reproduce with

```text
python 04-computation/quartic_signed_edge_star_triangle_reembedding_thm2993.py
python -O 04-computation/quartic_signed_edge_star_triangle_reembedding_thm2993.py
```

Normal and optimized transcripts are LF-identical to

```text
05-knowledge/results/quartic_signed_edge_star_triangle_reembedding_thm2993.out.
```

Frozen LF-normalized hashes (candidate evidence) are

```text
script  9ad6834d99c9d0ff083a021322e2a44ac8ccad5c1948d016835ff739be4773b0
output  60d38b0a78ab519d7e8f3cb3b78f2aa0bb6985369ad1c1a0c713ff48d02f45b5
```

The independent hostile audit reproduced every exact identity and transcript,
and caught the raw-support versus odd-character distinction repaired in
Section 6.
