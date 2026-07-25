---
id: THM-2331
title: "Two-sided septimal address embedding in a marked current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Let
  w be a primitive nine-coordinate scalar word with at least two
  coordinates nonzero modulo seven. Given an exact relation r with
  r.w=0, any integer frequency triangle X+m c_3=Y, and any requirement
  that both endpoint rectangle harmonics avoid the septimal zero set,
  there are at least 3*5^7=234,375 residue choices u modulo seven such
  that u.w=X and v=u+m e_(c_3)-r has v.w=Y, with every coordinate of u
  and v nonzero modulo seven. Here e_(c_3) is the coordinate vector
  labelled by c_3, independent of its numeric position in w. An exact
  Bezout lift preserves these conditions.
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
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2327-two-colour-marked-unit-c3-triangle
related:
  - THM-2301-essential-affine-arrangement-and-visible-rank-six-address-bank
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2329-boundary-triple-rerooting-and-transverse-gain-obstruction
script: 04-computation/lrc14_two_sided_address_embedding_thm2331.py
output: 05-knowledge/results/lrc14_two_sided_address_embedding_thm2331.out
script_sha256: 63a1d5ad694eaffcf6afb2a43501f04b41bb7d1af0cf0c07022c14bc97f75561
output_sha256: 72c6d0d2d227e46f6800ff33dab0a671eacc2e26f64dc00cb8b5fbbceb6275b7
hash_basis: working-tree bytes (LF)
---

# THM-2331 -- every target address embeds termwise in the marked current

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

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
  + put v=u+m e_(c_3)-r
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

Fix integers `X,m`, let `i_3` be the coordinate index labelled by the
deepest speed

```text
w_(i_3)=c_3,
```

and fix an exact relation

```text
r in Z^9,                  r.w=0.                  (4)
```

Put

```text
d=m e_(i_3)-r,
Y=X+m c_3.                                         (5)
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
v.w=X-r.w+m c_3=Y,                                 (13)

u+m e_(i_3)-v=r.                                   (14)
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

For `ell in {1,2}`, put

```text
J_ell={t:||t||<ell/14}.
```

Thus `J_1` is a danger arc of length `1/7`, while `J_2` is the centered
arc of length `2/7` whose complement is the guard. With the Fourier
convention used in the LRC canon,

```text
(1_(J_ell))_hat(0)=ell/7,
(1_(J_ell))_hat(n)=sin(pi*ell*n/7)/(pi*n),    n!=0,

(1-1_(J_ell))_hat(0)=1-ell/7,
(1-1_(J_ell))_hat(n)
 =-(1_(J_ell))_hat(n),                         n!=0. (16)
```

Since both `1` and `2` are units modulo seven, every one of these atomic
factor coefficients is nonzero exactly at

```text
n=0 or 7 does not divide n.                         (17)
```

This includes the wider guard: its zero coefficient is `5/7`, not
`6/7`.

Write the labelled scalar word as

```text
w=(H,q_1,...,q_5,c_1,c_2,c_3)=(w_0,...,w_8).
```

On a strict shallow-owner branch, the exact set differences in THM-2302
write `E=E_j`, up to null endpoints, as

```text
1_E(t)=prod_(i=0)^8 chi_i(w_i t),                  (18)
```

where the guard factor `chi_0` is the complement of `J_2`, and each
remaining factor is `J_1` or its complement. In particular (18) retains
the minus signs from every exclusive-owner complement; no inclusion-
exclusion term has been dropped.

Likewise, every pure or double terminal word `Q` in THM-2305 has an exact
nine-factor presentation

```text
1_Q(t)=prod_(i=0)^8 psi_i(w_i t),                  (19)
```

with the same guard factor and with the appropriate danger/complement
choice at each of the other eight coordinates. Every zero coefficient
`(psi_i)_hat(0)` is one of

```text
1/7, 5/7, 6/7,
```

and is strictly positive.

Set

```text
k=lambda_j+1,               R=13^k.
```

The transported word and the marked source therefore have the fully
atomic presentations

```text
1_W(t)=prod_(i=0)^8 psi_i(R w_i t),

1_(E_Q)(t)
 =prod_(i=0)^8 chi_i(w_i t)
  prod_(i=0)^8 psi_i(R w_i t).                     (20)
```

Choose the left modes to be `u` from Section 2 and choose all nine
transported-word modes

```text
beta_i=0.                                          (21)
```

Their total ordinary frequency is

```text
sum_i (u_i+R beta_i)w_i=u.w=X,
```

and their coefficient is

```text
prod_i (chi_i)_hat(u_i)(psi_i)_hat(0).             (22)
```

All coordinates of `u` are seven-units, so (17) and positivity of the
nine zero modes make (22) nonzero. Similarly, the bare-source modes `v`
give the nonzero coefficient

```text
prod_i (chi_i)_hat(v_i)                            (23)
```

at ordinary frequency `Y`.

Finally THM-2327 gives `7` not dividing `m`. If `d=1_(J_1)`, the deepest
comb term is

```text
(1_(D_(c_3)))_hat(m c_3)
 =d_hat(m)
 =sin(pi*m/7)/(pi*m)!=0.                           (24)
```

The full marked-current term is the product of (22), (24), and the
conjugate of (23). It is nonzero, and its labelled coordinate address is

```text
u+R beta+m e_(i_3)-v
 =u+m e_(i_3)-v
 =r.                                               (25)
```

Thus the claim uses `9+9+9=27` atomic interval factors, plus the deepest
comb. It does **not** replace the terminal word by its aggregate zero
coefficient: `(1_W)_hat(0)=measure(W)` may already contain cancellation
among its atomic terms.

## 4. Why Abel regularization is the correct statement

The atomic interval Fourier series are not absolutely summable. Apply
Poisson/Abel regularization separately to the eighteen endpoint factors
in (20), an independent nine-factor copy of (18) for the bare endpoint,
and the deepest comb: `28` factors in the current term. For `0<rho<1`,
multiply each factor's coefficient at base index `n` by

```text
rho^|n|.                                           (26)
```

Every regularized factor has an absolutely convergent Fourier series,
takes values in `[0,1]`, and converges in `L^1` to its original interval
or complement indicator. For any finite family of `[0,1]`-valued factors,
the bounded product telescope gives

```text
||prod_j f_(j,rho)-prod_j f_j||_1
 <=sum_j ||f_(j,rho)-f_j||_1.                      (27)
```

Hence the eighteen-factor marked product converges in `L^1` to `1_(E_Q)`
and the nine-factor bare product converges in `L^1` to `1_E`. Their
Fourier coefficients converge to the exact endpoint coefficients in (1);
the same one-factor statement applies to the deepest comb.

For every fixed `rho`, the fully atomic term selected above equals

```text
rho^(sum_i(|u_i|+|v_i|)+|m|)
  *prod_i (psi_i)_hat(0)
  *(sin(pi*m/7)/(pi*m))
  *prod_i (chi_i)_hat(u_i)
  *conjugate(prod_i (chi_i)_hat(v_i)),             (28)
```

because the nine chosen word modes in (21) contribute no Abel exponent
and have strictly positive coefficients. Its address identity (25) is
independent of `rho`. Therefore `r` genuinely occurs as a nonzero,
fully atomic interval-factor term of an absolutely convergent
regularization of the exact current.

This formulation does **not** say:

- that the term (28) is a nonzero aggregate after collecting other terms
  with the same relation address;
- that it has positive real part;
- that its nonzero termwise limit survives cancellation in the
  same-address grouped sum, either before or after `rho` tends to one;
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
R_q=3,134,566,563,840                              (29)
```

distinct exact relations `r` such that

```text
r mod 13 lies in q+L_13,
every r_i is a unit modulo 91,
||r||_infinity<=45(1+S B(w)).                      (30)
```

Apply Sections 1--4 to the marked triangle (1). Each `r` has at least
`234,375` distinct two-sided residue lifts. Since distinct `r` or distinct
`u modulo 7` give distinct term/address pairs, every prescribed vector
fibre contains at least

```text
234,375 R_q
 =734,664,038,400,000,000                          (31)
```

nonzero Abel-regularized term/address pairs. Every projective direction
contains twelve vector fibres and hence at least

```text
8,815,968,460,800,000,000                          (32)
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
function-role/target polarization             OPEN;
terminal-component phase transport           OPEN;
scalar-row exclusion                         OPEN.   (33)
```

No scalar profile is excluded. The exact ledger remains `165`, and
LRC(14) remains open.

## 6. Exact companion

The companion exhausts every two-coordinate septimal speed, displacement,
and right-hand side, verifies that the minimum completion count is exactly
three, checks the separate danger/safe and guard support laws, the
`18+9+1=28` factor ledger, both pure and double terminal-word zero-mode
products, the general `3*5^7` count and both term-bank constants, and
constructs an exact positive control on THM-2325's hostile scalar word.
For that control it verifies the target-axis relation, the exact
two-sided lift, all eighteen septimal nonvanishing predicates, the height
invoice, and the Abel-term support predicate. It also records the exact
support-one boundary and a two-term cancellation hostile showing why
termwise participation cannot be promoted to grouped-address survival.
Every load-bearing check raises explicitly under ordinary and optimized
Python.

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

The independent audit rederived the sharp two-pivot floor and equality
case, the labelled address identity, the exact Bezout lift and height
invoice, and both atomic support laws. It checked the exact presentations
of the pure and double THM-2305 words, the `18+9+1` factor count, the
positive transported zero-mode products, and the bounded `L^1` product
telescope. It reproduced the ordinary, optimized, and stored transcripts
byte-for-byte and confirmed that termwise limiting nonvanishing does not
imply same-address grouped survival, target polarization, visibility, or
terminal-component phase transport.
