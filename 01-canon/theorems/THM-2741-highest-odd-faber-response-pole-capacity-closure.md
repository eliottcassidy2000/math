---
id: THM-2741
title: "Highest odd-Faber response pole-capacity closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the full
  chosen-sheet split degree-22 Faber response
  family, every geometrically integral member with at least one nonzero odd
  Faber coefficient is physically empty.  If a_j is the highest odd seed and
  r=22-j, the G2 infinity face has 3*gcd(r,6) coarse normalization branches.
  The homogenized third response has a pole on every branch, with order at
  least eight.  Pullback would give the source rational primitive at least
  three poles, whereas THM-2723 permits at most one.  Reducible/nonreduced
  response members, the all-even zero-flux boundary, JC(2), and DC(2) remain
  open.
source: root/highest-odd-faber-response-pole-capacity-2026-07-28
audit: thm2705-2709-audit-2026-07-28 and thm2694-full-lift-fibre-scout-2026-07-28 (independent recurrence/multinomial, weighted-branch, quotient-parity, response-pole, pullback, replay, and docs audits)
depends_on:
  - THM-2719-full-split-odd-faber-generic-normalization-genus-four-hundred-nineteen
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
related:
  - THM-2725-split-even-faber-nonzero-first-flux-unified-closure
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
script: 04-computation/jc2_degree22_highest_odd_response_poles_thm2741.py
output: 05-knowledge/results/jc2_degree22_highest_odd_response_poles_thm2741.out
script_sha256: 6f0a940c71b8b0e890e7c9cb51fc4a2b55a9b2eabbea6a4ecad946269672c4dc
output_sha256: 28ba0ec2df8bf9300e84d37b44241672db73a9f686ec271d92c50c7faca52beb
hash_basis: LF-normalized bytes
---

# THM-2741 -- every highest odd seed creates forbidden response poles

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2726 uses the physical coefficient `q` to remove the `a_21!=0`
geometrically integral locus.  The same infinity germ contains a stronger
object: the third response is transverse for **every** possible highest odd
seed.  This includes `a_3` and `a_1`, where `q` itself no longer has an
infinity pole.

## 1. Statement

Use THM-2719's full chosen-sheet split response family

```text
F_23=Phi_22+sum_(j in J) a_j h^(22-j)Phi_j-lambda h^23,
G_24=Psi_22+sum_(j in J) a_j h^(22-j)Psi_j-W h^24,     (1)

J={1,2,3,5,6,7,9,10,11,13,14,15,17,19,21}.
```

Let `C_a` be a geometrically integral member.  If at least one coefficient
with odd index is nonzero, then `C_a` supports no physical polynomial Keller
trajectory arising from the split polynomial exact-square prefix and the
constant-coefficient Faber gauge in `(1)`.

No condition is imposed on `lambda`, `W`, or the four even coefficients.
No generic smoothness or genus hypothesis is used.

## 2. The highest odd seed and its weighted initial system

Let `j` be the largest odd index with `a_j!=0`, and put

```text
r=22-j in {1,3,5,...,21},             g=gcd(r,6).       (2)
```

Work on the `d=1` index cover at the unique infinity point.  Give the local
variables the integral weights

```text
v(h)=6/g,                    v(q)=v(s)=r/g.             (3)
```

The top degree-22 faces both have weight `6r/g`:

```text
Phi_(22,6)=-(231/128) q s(q^2-3s^2)(3q^2-s^2),

Psi_(22,6)=-(231/256)(q^2-s^2)
                         (q^4-14q^2s^2+s^4).           (4)
```

For every odd reduced degree in `(1)`, exact Faber recurrence gives

```text
Phi_j(1,0,0)!=0,            R_j(1,0,0)!=0,
ord_(q,s) Psi_j(1,q,s)>=1.                               (5)
```

Thus the chosen `a_j h^r Phi_j` has the same weight as `(4)`.  Every lower
odd seed has gap strictly larger than `r`.  The four even columns have the
following common local orders for `Phi`, `Psi`, and `R`:

```text
degree j        14   10    6    2
gap 22-j         8   12   16   20
local order       4    3    2    1.                    (6)
```

For an even column with gap `e` and local order `o`, its suppressed common-
`1/g` weight is `6e+or`.  In every row of `(6)` and for every odd
`r<=21`,

```text
6e+or>6r.                                             (7)
```

The `lambda h^23` and `W h^24` terms are also strictly higher.  The exact
companion checks all three observables and every one of the eleven choices
of `j`; its smallest strict margin is positive in every row.  Consequently
the complete weighted initial system is exactly

```text
Psi_(22,6)(q,s)=0,
Phi_(22,6)(q,s)+a_j Phi_j(1,0,0)h^r=0.                (8)
```

## 3. Three or nine coarse normalization branches

The second polynomial in `(4)` has six distinct lines

```text
q=alpha s,              alpha^2 in
{1,7+4sqrt(3),7-4sqrt(3)}.                            (9)
```

The first polynomial in `(4)` is nonzero on every line in `(9)`.  On each
line, `(8)` therefore reduces at its initial face to

```text
h^r=c_alpha s^6,                    c_alpha!=0.         (10)
```

Newton--Puiseux, or weighted Hensel lifting after the substitution in `(3)`,
gives exactly `g=gcd(r,6)` normalization branches over each line.  Hence the
index cover has `6g` branches.  Its central involution is

```text
(h,q,s)|->(-h,-q,s).                                  (11)
```

It exchanges the nonzero tangent lines `alpha` and `-alpha`, so it acts
freely on these branches.  The coarse response curve therefore has

```text
3g in {3,9} distinct normalization points at infinity. (12)
```

On each lifted branch a normalization parameter `t` may be chosen with

```text
ord_t(s)=ord_t(q)=r/g,             ord_t(h)=6/g.       (13)
```

## 4. The third response is transverse on every branch

Homogenize the third Faber observable to weight `25`:

```text
R_25=R_22+sum_(k in J) a_k h^(22-k)R_k.               (14)
```

The affine response is the involution-invariant rational function

```text
R_aff=R_25/h^25.                                      (15)
```

Indeed THM-2129's exact parity is
`R_k(d,-q,s)=(-1)^(k+1)R_k(d,q,s)`.  Multiplication by
`h^(22-k)` makes every summand of `R_25` anti-invariant under `(11)`, as is
`h^25`; their quotient therefore descends to the coarse curve.

Restricting the weight-`6r/g` face of `(14)` by the first equation of `(8)`
gives

```text
R_(22,6)-[R_j(1,0,0)/Phi_j(1,0,0)]Phi_(22,6)
 =c_j q s(q^2-3s^2)(3q^2-s^2),                       (16)
```

where the exact nonzero constants are

```text
j       21       19       17         15       13
c_j  1771/1024 441/256 4389/2560 1309/768 3465/2048

j       11       9       7          5        3       1
c_j   429/256 847/512 2079/1280 1617/1024 385/256 693/512. (17)
```

The transversality is an exact quadrature identity, not an accidental gcd.
For `z=q+i s`, the two faces are proportional to

```text
Psi_(22,6)=-(231/256) Re(z^6),
q s(q^2-3s^2)(3q^2-s^2)=(1/2) Im(z^6).               (17a)
```

Thus the response zeros interlace the six `G2` tangent lines.  The central
`C2` quotient pairs each six-line family into three directions; it does not
erase their transverse phase.

The polynomial in `(16)` is coprime to `Psi_(22,6)`.  Hence it is nonzero
on every line in `(9)`, and on every normalization branch

```text
ord_t(R_25)=6r/g,
pole_order(R_aff)=(150-6r)/g>0.                       (18)
```

The pole orders, in descending odd degree, are

```text
144,44,120,108,32,84,72,20,48,36,8.                  (19)
```

Thus even the weakest row `j=1` has nine coarse response poles of order
eight.  The response, unlike `q`, remains decisive in the two final odd
rows.

## 5. Source rational-primitives have at most one pole

A physical trajectory gives a rational map from the source line to `C_a`.
It is nonconstant because THM-2723 gives the exact third-flux identity

```text
U R_source'=kappa,                 kappa!=0.            (20)
```

Since `C_a` is geometrically integral, the map lifts to its normalization
and extends to a nonconstant morphism from `P1`.  It is therefore finite and
surjective.  Every one of the at least three distinct pole points in `(12)`
has a nonempty, pairwise disjoint source fibre, so `R_source` would have at
least three distinct poles.

But the rational-primitive classification used in THM-2723 is sharper for
the response itself.  Either

```text
U in C* and R_source is affine-linear,
```

or, after translating `x`,

```text
U=u_0 x^m,            R_source=s_0+s_1 x^(1-m), m>=2. (21)
```

The first case has only the pole at infinity; the second has only the finite
pole at `x=0`.  Hence `R_source` has at most one pole on `P1`, contradicting
the pullback of `(18)`.  This proves the stated exclusion.

## 6. Boundary and inheritance

This theorem strictly extends THM-2726 from `a_21!=0` to every nonzero odd
Faber direction, and it includes arbitrary values of the four even seeds and
both flux constants.  Combined with THM-2725, it makes the all-even
nonzero-first-flux chart and every geometrically integral odd deformation
empty.

It does **not** exclude a reducible or nonreduced odd response member.  Its
infinity branches can lie on different components, while a physical
trajectory sees only one image component.  It also does not remove the
all-even zero-first-flux boundary, justify this degree-22 normal form for an
arbitrary Keller pair, or prove `JC(2)` or `DC(2)`.  The sharp next target is
therefore componentwise: prove that any component carrying the third-flux
trajectory inherits at least two response branches, or classify the
factorization locus on which the `3g` branches split into one-pole pieces.

There is already a sharp arithmetic constraint on the latter escape.  If a
physical component contains exactly one response pole of order
`P=(150-6r)/g`, then the constant-`U` case is impossible: an affine-linear
source response has pole order one, whereas `P>=8`.  Hence `(21)` holds, the
unique target pole has one totally ramified source preimage, and

```text
degree(P1_source -> normalization(component))=e,
m-1=eP.                                                (22)
```

In particular `P` divides `m-1`.  This divisibility law does not prove that
the physical component contains one of the forced branches; it records the
strongest surviving one-branch boundary rather than hiding it inside the
reducible caveat.

## 7. Exact companion

Run

```text
python 04-computation/jc2_degree22_highest_odd_response_poles_thm2741.py
python -O 04-computation/jc2_degree22_highest_odd_response_poles_thm2741.py
```

and compare with

```text
05-knowledge/results/jc2_degree22_highest_odd_response_poles_thm2741.out.
```

The companion reconstructs all Faber rows from the recurrence, checks every
weighted inequality rather than assuming the table, proves each response
coefficient in `(17)` nonzero, verifies every bivariate gcd, and records the
branch and pole invoices for all eleven odd seeds.  It uses explicit
exceptions rather than optimization-sensitive `assert` statements.

## 8. Independent hostile audit

The audit independently reconstructed the load-bearing local data rather
than replaying the companion.  At the square point `P=(w^2+1)^2` it obtained,
for every odd `j`,

```text
c_1=binom(j/2,(j+1)/2)!=0,       c_3/c_1=-1/(j+3),
R_j/Phi_j=(j+1)/(2(j+3))!=0,
```

while `Psi_j` has no constant term by parity.  This gives
`c_j=231(j+2)/(128(j+3))` and reproduces all eleven constants in `(17)`.
For the four even columns it reduced every strict face inequality to
`r<24`, so no omitted even seed or flux term reaches the initial system.

After eliminating the simple transverse coordinate in `Psi_(22,6)`, each
of its six lines has equation `h^r=c s^6+higher`.  The audit therefore
recovered exactly `gcd(r,6)` branches per line, checked that the central
involution freely pairs opposite lines, and verified that `R_25/h^25`
descends.  Its pole calculation gives `(150-6r)/g` on every coarse branch,
including the minimum order eight.  Finally it separately checked that the
source classification has exactly one pole in either case of `(21)` and
that a finite surjective normalization map pulls back distinct pole fibres.
Normal, optimized, and stored transcripts agree byte for byte and both
declared hashes match.  No additional hypothesis was found.
