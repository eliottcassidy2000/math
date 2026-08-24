---
id: THM-3924
title: "Decic cubic index-five ramification-class obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The normal
  rational cubic completion B of THM-3921 has class group Z. Its two primes
  above A=0 have residue degrees one and two and generate Cl(B) modulo their
  diagonal principal relation. The unique ramification prime E has pole
  valuations -6 and -1 in its localized equation, so [E] is five times a
  primitive generator. THM-3922 requires every affine-plane boundary prime
  to be primitive; because E must be deleted and Cl(B) has rank one, no
  additional divisor can repair the failure. This independently excludes
  the genus-collapsed decic cubic atlas. The same five is the A^5
  power-basis index exponent; JC(2) remains open.
source: jc_degree6_one_place / post-THM-3922 primitive-boundary lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (incoming_thm3923_3924_audit/root and
  jc_zero_debt_lift, 2026-08-23).  They reconstructed the radial UFD chart and
  full unit group, inclusion B subset k[u,z], simple-plus-double A=0 prime
  packet and residue degrees, diagonal Nagata quotient Cl(B)=Z, irreducible
  ramification prime, complete divisor relation with valuations (-6,-1), and
  THM-3922 primitive-basis contradiction.  They checked that no extra divisor
  repairs rank one and that the conclusion is confined to the actual
  THM-3921 completion, not THM-3915.  The companion verifies the basis
  determinant, residual irreducibility, valuation gap, class divisibility,
  and power-index exponent in 20 gates.  LF-normalized normal and optimized
  streams match the frozen LF output; raw script/output and semantic hashes
  and documentation checks pass.  No repair was required.
depends_on:
  - THM-3921-quintic-genus-collapse-decic-degeneration-packet
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3917-quintic-parameter-rational-collapsed-cubic
  - THM-3915-rational-decic-cube-resolvent-index-debt-euler-tariff
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
script: 04-computation/jc2_decic_cubic_index_five_ramification_class_thm3924.py
output: 05-knowledge/results/jc2_decic_cubic_index_five_ramification_class_thm3924.out
script_sha256: 3c0af6ff50dffb74336732367f98451aa108fb1d501a5c8daa4f268bf48e15cd
output_sha256: cf93d07d262c30ca5d19a203e95a33a7f75e1534af53eef947193f129bfa785b
semantic_sha256: a0c6d59e67172193dfe81763d33589dd46ed1eca185e0e60ba32f01c6202e18b
hash_basis: raw LF bytes
---

# THM-3924 -- the A5 index debt is a five-divisible boundary class

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let `b` be any root
of the irreducible quintic

```text
K(b)=2304b^5+10176b^4+4064b^3+996b^2+84b+5.               (1)
```

Retain the THM-3921 polynomials

```text
p=(u^2-1)(u^2+b)^2,
h=[p^(3/2)]_+,
F=z^3-3pz+2h,                                               (2)

A=F/4,                         C=uF/4,
Z=A^3z.                                                     (3)
```

Let `R=k[A,C]` and let `B` be the explicit normal finite-free cubic order
of THM-3921. Let `E` be its ramification prime above the irreducible decic
discriminant. Then

```text
Cl(B)=Z<g>,                         [E]=5g.                  (4)
```

In particular `[E]` is nonzero but nonprimitive. There is no open immersion
`A2 -> Spec B` on which the finite cubic map to `Spec R` is etale. Thus the
THM-3921 field cannot give a planar Keller counterexample for either of two
independent reasons:

```text
topology:  six normalization addresses meet one boundary point (THM-3917),
classes:   the required ramification boundary has class 5g (this theorem).
                                                                    (5)
```

The integer five is not an unrelated afterthought. It is simultaneously

```text
the power-order index exponent A^5,
the difference 6-1 of the two pole orders of z^2-p,
the divisibility of the ramification class [E].             (6)
```

## 1. Localizing at A recovers a factorial radial chart

Put `X=Spec B`. Since the two extra integral basis elements in THM-3921 have
only powers of `A` in their denominators, inverting `A` removes the
overorder distinction:

```text
B_A=R_A[Z].                                                  (7)
```

In this localization `u=C/A` and `z=Z/A^3`. The homogenized cubic relation
is exactly

```text
f(A,C,A^3z)=A^9(F(u,z)-4A).                                 (8)
```

Eliminating `A=F/4` therefore gives

```text
B_A isomorphic to k[u,z,F^(-1)].                            (9)
```

THM-3917 proves that `F` is irreducible on the quintic locus `(1)`. Hence
the right side of `(9)` is a localization of a polynomial UFD and

```text
Cl(B_A)=0,                    B_A^*=k^* A^Z.                (10)
```

The units of `B` itself are scalar. Indeed every element of `B` is integral
over `R subset k[u,z]`, hence integral over the normal ring `k[u,z]`; because
both rings have fraction field `k(u,z)`,

```text
B subset k[u,z].                                             (11)
```

A unit of `B`, together with its inverse, is therefore a unit of the
polynomial ring, so

```text
B^*=k^*.                                                     (12)
```

## 2. The two primes over A=0 and the Nagata quotient

At the generic point of the target line `A=0`, the power polynomial reduces
to

```text
f(0,C,Z)=(Z-C^3)^2(Z+2C^3).                                (13)
```

The discriminant of the normal order is a nonzero constant times the
irreducible decic `Delta`, and exactly

```text
Delta(0,C)=C^8(512C+96b^3+48b^2+18b+3)/128 !=0.           (13a)
```

Thus the normal fibre over the generic point of `A=0` is etale. The simple
root in `(13)` gives one prime `P_-` of residue degree one. THM-3921 resolves
the doubled power root by

```text
c=b-1/2,                       d=(4b+1)/8,
L=C^3+cA^2C,                   m=L-dA^4/C,
Z=m+A^5x,                                                  (14)

384C^4x^2-512C-(96b^3+48b^2+18b+3)=0 modulo A.            (15)
```

The right side of the corresponding square equation has a simple zero as a
rational function of `C`; its constant is nonzero modulo `K`. Hence `(15)`
is irreducible over `k(C)`. It gives one prime `P_2` of residue degree two.
Both primes are unramified, so

```text
div_B(A)=P_-+P_2.                                          (16)
```

Nagata localization for `(9)` now reads

```text
B_A^*/B^* -> Z[P_-] direct_sum Z[P_2] -> Cl(B) -> 0.       (17)
```

By `(10),(12)`, the first term is `Z`, generated by `A`, and `(16)` sends
its generator to `(1,1)`. Consequently

```text
Cl(B)=(Z[P_-] direct_sum Z[P_2])/<P_-+P_2>
     isomorphic to Z,

g=[P_-],                         [P_2]=-g.                  (18)
```

This calculation uses the actual cubic completion. The three-torsion on
the quadratic resolvent in THM-3921 neither proves nor enters `(18)`.

## 3. The ramification equation has pole vector (-6,-1)

On the UFD chart `(9)`, the Jacobian calculation of THM-3920 shows that the
ramification divisor is

```text
H_ram=z^2-p(u)=0.                                           (19)
```

The polynomial in `(19)` is irreducible because `p` is not a square in
`k(u)`. It remains prime after `F` is inverted: `F` and `H_ram` have no
common factor because their resultant is the nonzero residual `-4r`.
Therefore `(19)` cuts out the restriction of one prime `E` of `B`, with
multiplicity one.

Homogenize the square term by

```text
Pi=A^6p(C/A),
H_ram=A^(-6)(Z^2-Pi).                                      (20)
```

At `P_-`, equation `(13)` gives `Z=-2C^3 modulo A`, and hence

```text
(Z^2-Pi) modulo A=3C^6 !=0,
v_(P_-)(H_ram)=-6.                                         (21)
```

At `P_2`, substitute `(14)`. Exact expansion gives

```text
(m+A^5x)^2-Pi
 =A^5(2C^3x+A(2b+1)/8+O(A^2)).                            (22)
```

The residue of `x` is nonzero by `(15)`, so the leading coefficient in
`(22)` is nonzero. Therefore

```text
v_(P_2)(H_ram)=5-6=-1.                                     (23)
```

There are no other horizontal zeros of `H_ram` on `(9)` and no other primes
above `A=0`. Taking the divisor of the rational function `H_ram` on `X`
therefore gives the complete relation

```text
div_X(H_ram)=E-6P_--P_2.                                   (24)
```

Using `(18)` in `(24)` proves the central identity

```text
[E]=6[P_-]+[P_2]=5g.                                       (25)
```

## 4. Why no extra boundary divisor can pay the debt

Suppose there were an affine-plane Keller atlas

```text
j:A2 -> X                                                   (26)
```

for the finite cubic map. Its image must omit `E`, because the finite map is
ramified at the generic point of `E` while its restriction to the source
plane is etale. THM-3922 says that the prime divisorial components of
`X minus j(A2)` must form a `Z`-basis of `Cl(X)`.

But `(18)` has rank one. Hence a proper affine-plane open can omit exactly
one prime divisor, and that divisor's class must be primitive. The already
mandatory divisor is `E`, whose class is `5g` by `(25)`. It is not
primitive. Deleting an additional prime cannot help: two boundary primes
would have to be a two-element basis of a rank-one group, contradicting
THM-3922 before any relation is chosen. Thus `(26)` is impossible.

This argument is independent of the six-address topology. It would still
exclude the order if a deformation separated every normalization address
while retaining the class packet `(18),(25)`.

## 5. Why the same five appears twice

For the THM-3921 integral basis `(1,e_1,e_2)`, the determinant relative to
the power basis `(1,Z,Z^2)` is

```text
det((1,e_1,e_2)/(1,Z,Z^2))=-(4b+1)/(8A^5).                 (27)
```

Thus the power order has index ideal `(A^5)`. Independently, the two
infinity blocks see the ramification equation with pole orders six and one;
their difference is five and becomes the coefficient in `(25)`. Equations
`(25),(27)` prove the exact identity of integers claimed in `(6)`:

```text
power-index debt 5 = infinity valuation imbalance 5
                   = boundary-class divisibility 5.        (28)
```

The theorem does not claim that every index-`A^m` cubic order has an
`m`-divisible ramification class. Here the equality comes from the explicit
simple-plus-resolved-double infinity packet `(13)--(15)`. Nor does it
transfer the conclusion to the differently scaled THM-3915 order without a
separate coordinate audit. What is proved is the complete class obstruction
for the quintic degeneration of THM-3921. No planar Jacobian counterexample
is obtained, and `JC(2)` remains **OPEN**. **QED.**
