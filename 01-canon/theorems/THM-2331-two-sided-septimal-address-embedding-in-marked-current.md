---
id: THM-2331
title: "Two-sided septimal address embedding in a marked current"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. Let
  w be a primitive nine-coordinate scalar word with at least two
  coordinates nonzero modulo seven. Given an exact relation r with
  r.w=0, any integer frequency triangle X+m w_3=Y, and any requirement
  that both endpoint rectangle harmonics avoid the septimal zero set,
  there are at least 3*5^7=234,375 residue choices u modulo seven such
  that u.w=X and v=u+m e_3-r has v.w=Y, with every coordinate of u and
  v nonzero modulo seven. An exact Bezout lift preserves these conditions.
  Consequently, under THM-2325 and THM-2327, every exact all-91-unit
  relation address in every prescribed nonzero target-vector fibre occurs
  as a nonzero Abel-regularized interval-factor term of the marked
  word/deepest-comb/bare Fourier current. This gives at least
  734,664,038,400,000,000 distinct term/address pairs per target vector.
  It does not say that one selected term survives the aggregate sum,
  lies in the bounded visible carrier, has positive phase, or transports
  terminal-component phase. No scalar row is excluded and LRC(14)
  remains open.
source: codex-2026-07-25-two-sided-address-embedding
depends_on:
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2327-two-colour-marked-unit-c3-triangle
related:
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2329-boundary-triple-rerooting-and-transverse-gain-obstruction
script: 04-computation/lrc14_two_sided_address_embedding_thm2331.py
output: 05-knowledge/results/lrc14_two_sided_address_embedding_thm2331.out
script_sha256: f9afd297c5dc002949b48260bfa06dbbb6f46281c0c4c27f76630f1770095f42
output_sha256: cfc79484cfe134b4a37ee3268064f3e264022eee379de2ffa7ad4e528f0fdd5b
hash_basis: working-tree bytes (LF)
---

# THM-2331 -- every target address embeds termwise in the marked current

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2327 gives a genuine nonzero mixed Fourier triangle

```text
(1_(E_Q))_hat(X)
(1_(D_(c_3)))_hat(m c_3)
conjugate((1_E)_hat(Y)) !=0,

Y=X+m c_3,                    gcd(m,91)=1.           (1)
```

THM-2325 independently gives an enormous bank of exact all-`91`-unit
relations in every prescribed target-vector fibre. The missing statement
was whether one of those addresses can occur inside the interval-factor
expansion of (1), rather than merely sharing its finite-field label.

In fact every address in the bank occurs termwise. The operation is a
two-sided harmonic lift:

```text
relation address r
  + choose a left rectangle harmonic u
  + put v=u+m e_3-r
  + force u and v off every septimal Fourier zero.  (2)
```

Only one scalar equation remains. Two coordinates of the speed word are
nonzero modulo seven, while each endpoint condition forbids at most two
residues. This leaves three completions, uniformly.

The conclusion is deliberately **termwise and Abel-regularized**. It does
not assert that the selected term is a nonzero sub-sum after other terms
are collected, or that it survives cancellation with a controlled sign.

## 1. The two-sided finite-field lift

Let

```text
w=(w_0,...,w_8) in Z^9,
gcd(w_0,...,w_8)=1,

S=sum_i |w_i|,
supp_7(w)={i:w_i!=0 mod 7},
s=|supp_7(w)|>=2.                                  (3)
```

Fix integers `X,m`, the coordinate labelled by the deepest speed `w_3`,
and an exact relation

```text
r in Z^9,                  r.w=0.                  (4)
```

Put

```text
d=m e_3-r,
Y=X+m w_3.                                         (5)
```

For each coordinate define

```text
A_i=F_7 minus {0,-d_i}.                             (6)
```

If `d_i=0`, the two deleted residues coincide, so `|A_i|=6`; otherwise
`|A_i|=5`. A vector `u modulo 7` lies in the product of the `A_i` exactly
when

```text
u_i!=0 mod 7,
v_i:=u_i+d_i!=0 mod 7,              every i.       (7)
```

Choose two distinct indices `a,b` in `supp_7(w)`. Fix the other seven
coordinates arbitrarily in their allowed sets. There are at least `5^7`
such choices. For every fixed choice, the remaining equation

```text
w_a u_a+w_b u_b
 =X-sum_(i notin {a,b})w_i u_i             mod 7   (8)
```

has at least three solutions in `A_a x A_b`. Indeed, every `u_a in A_a`
determines a unique `u_b`; this map is injective, and at most

```text
7-|A_b|<=2
```

of the at least five choices of `u_a` land outside `A_b`. Therefore

```text
# {u mod 7:
     u.w=X,
     u_i!=0,
     u_i+d_i!=0 for all i}
 >=3*5^7
 =234,375.                                         (9)
```

The constant is sharp for the abstract lemma. Take only two supported
speed coordinates, both equal to one, take every `d_i=1`, and choose the
right-hand side on (8) so that the two five-element sets have only three
valid pairs. The exact companion exhausts all two-coordinate data and
finds this equality case.

The support threshold is also exact. With only one supported coordinate,
the scalar equation fixes that coordinate and can force either forbidden
residue. In THM-2325's stronger all-unit relation setting, `s=1` already
makes the address bank empty.

## 2. Exact Bezout lift and height

Choose a Bezout vector

```text
z.w=1,                       B(w)=||z||_infinity.   (10)
```

Take any solution from (9), represented in `[-3,3]^9` by `u_0`. Since
`u_0.w=X mod 7`, the integer

```text
h=(X-u_0.w)/7                                      (11)
```

is defined. Set

```text
u=u_0+7h z,
v=u+d.                                             (12)
```

Then

```text
u.w=X,
v.w=X-r.w+m w_3=Y,                                 (13)

u+m e_3-v=r.                                       (14)
```

The correction in (12) is divisible by seven, so all eighteen nonzero
conditions in (7) survive. Distinct residue choices in (9) give distinct
exact pairs. The explicit height invoice is

```text
||u||_infinity
 <=3+B(w)(|X|+3S),

||v||_infinity
 <=3+B(w)(|X|+3S)+||r||_infinity+|m|.              (15)
```

No positivity, Fourier analysis, or target-plane structure entered this
lift. It is a sharp affine avoidance statement over `F_7`, followed by
one exact integral correction.

## 3. The canonical rectangle has exactly the needed support

On a strict shallow-owner branch, THM-2302 writes the exclusive source
`E=E_j`, up to null endpoints, as a nine-factor circular rectangle. Eight
factors are the centered danger interval

```text
D={t:||t||<1/14}
```

or its complement: five unit-speed complements, the selected shallow
danger factor, and the two other blocker complements. Their coefficients
are

```text
d_hat(0)=1/7,
d_hat(n)=sin(pi*n/7)/(pi*n),             n!=0,

(1-d)_hat(0)=6/7,
(1-d)_hat(n)=-d_hat(n),                  n!=0.      (16)
```

The ninth factor is the wider guard window

```text
C={t:||t||>1/7}.
```

Its coefficients are

```text
c_hat(0)=5/7,
c_hat(n)=-sin(2*pi*n/7)/(pi*n),          n!=0.      (16a)
```

Since multiplication by two permutes the nonzero residues modulo seven,
both (16) and (16a) are nonzero exactly at

```text
n=0 or 7 does not divide n.                         (17)
```

The lift in Section 1 is stronger than necessary: all coordinates of both
`u` and `v` are seven-units. Hence the rectangle monomials

```text
prod_i I_i_hat(u_i),
prod_i I_i_hat(v_i)                                 (18)
```

are nonzero, and (13) places them respectively in the ordinary
frequencies `X` and `Y`.

Write

```text
W=T^(-(lambda_j+1))Q,
E_Q=E intersection W.                              (19)
```

The canonical word `Q` has positive measure, so

```text
(1_W)_hat(0)=measure(W)=measure(Q)>0.               (20)
```

Consequently the `u` monomial in (18), together with the zero Fourier
mode of the extra word factor, is a nonzero term in the factor expansion
of `(1_(E_Q))_hat(X)`. The `v` monomial is a nonzero term in
`(1_E)_hat(Y)`.

Finally, THM-2327 gives `7` not dividing `m`, and hence

```text
(1_(D_(c_3)))_hat(m c_3)
 =sin(pi*m/7)/(pi*m)!=0.                           (21)
```

Equations (14), (18), (20), and (21) show that the selected term in the
expanded mixed current (1) is nonzero and carries the exact relation
address `r`.

## 4. Why Abel regularization is the correct statement

The atomic interval Fourier series are not absolutely summable. For
`0<rho<1`, multiply every atomic coefficient at index `n` by

```text
rho^|n|.                                           (22)
```

The resulting factor expansions are absolutely convergent. Their products
converge in `L^1` to the exact rectangle and word products as `rho` tends
to one from below, so their Fourier coefficients Abel-converge to the
three exact coefficients in (1).

For every fixed `rho`, the term selected above equals

```text
rho^(sum_i(|u_i|+|v_i|)+|m|)
  *measure(W)
  *(sin(pi*m/7)/(pi*m))
  *prod_i I_i_hat(u_i)
  *conjugate(prod_i I_i_hat(v_i)),                 (23)
```

and is nonzero. Its address identity (14) is independent of `rho`.
Therefore `r` genuinely occurs as a nonzero term of an absolutely
convergent regularization of the exact current.

This formulation does **not** say:

- that the term (23) is a nonzero aggregate after collecting other terms;
- that it has positive real part;
- that its contribution survives as `rho` tends to one without
  cancellation;
- that its endpoint harmonics `u,v` have bounded visible height.

Those are different quantifiers. The theorem closes term/address
participation, not cancellation or visibility.

## 5. Every prescribed target-vector bank embeds

Now assume the scalar and owner-packet hypotheses of THM-2325. For every
nonzero target vector

```text
q in K_13/L_13,
```

that theorem supplies at least

```text
R_q=3,134,566,563,840                              (24)
```

distinct exact relations `r` such that

```text
r mod 13 lies in q+L_13,
every r_i is a unit modulo 91,
||r||_infinity<=45(1+S B(w)).                      (25)
```

Apply Sections 1--4 to the marked triangle (1). Each `r` has at least
`234,375` distinct two-sided residue lifts. Since distinct `r` or distinct
`u modulo 7` give distinct term/address pairs, every prescribed vector
fibre contains at least

```text
234,375 R_q
 =734,664,038,400,000,000                          (26)
```

nonzero Abel-regularized term/address pairs. Every projective direction
contains twelve vector fibres and hence at least

```text
8,815,968,460,800,000,000                          (27)
```

such pairs.

Thus the complete target-vector label is not merely copied by a relation
which lives elsewhere: that very relation can be embedded in the marked
mixed current's interval-factor expansion. In particular one may choose
the boundary vector selected by rerooting the THM-2327 triangle, or any
other nonzero vector selected by an independent target-polarized current.

The remaining proof object is narrower:

```text
termwise address participation              PROVED here;
cancellation-free address sub-sum            OPEN;
bounded visible/Jackson membership           OPEN;
terminal-component phase transport           OPEN;
scalar-row exclusion                         OPEN.   (28)
```

No scalar profile is excluded. The exact ledger remains `165`, and
LRC(14) remains open.

## 6. Exact companion

The companion exhausts every two-coordinate septimal speed, displacement,
and right-hand side, verifies that the minimum completion count is exactly
three, checks the separate danger/safe and guard support laws, the general
`3*5^7` count and both term-bank constants, and constructs an exact
positive control on THM-2325's hostile scalar word.
For that control it verifies the target-axis relation, the exact
two-sided lift, all eighteen septimal nonvanishing predicates, the height
invoice, and the Abel-term support predicate. Every load-bearing check
raises explicitly under ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_two_sided_address_embedding_thm2331.py
python3 -O 04-computation/lrc14_two_sided_address_embedding_thm2331.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_two_sided_address_embedding_thm2331.out
```

byte-for-byte after LF normalization.
