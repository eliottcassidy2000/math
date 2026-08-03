---
id: THM-3253
title: "Positive owner-mass Newton cyclicity and maximal common Heisenberg module"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  For every integer dilation in THM-3246's (3,5;1,2) lane, compactify its
  168 strictly positive cell masses by zero on the THM-3234 Singer plane.
  The resulting 13 by 13 matrix is nonsingular in all 672 scalar,
  Frobenius and reflection Singer gauges.  Exact degree-26 Newton
  certificates prove this for every dilation, including two phase classes
  with sharp finite heads and shifted positive tails.  Placing the mass
  matrix on one central slice gives a nonnegative packet whose H_13 orbit
  spans the sharp 2,041-dimensional common submodule from THM-3250.  The
  plane relocation remains abstract and is not a canonical endpoint current.
source: root/2026-08-03
audit: >
  The assertion-independent companion pins THM-3234, promoted THM-3246 and
  promoted THM-3250; rebuilds all 168 numerator quadratics and their
  positivity; verifies the deterministic Singer, scalar, anti-diagonal,
  Frobenius and reflection reductions; computes exact Bareiss determinants;
  and checks every finite head and all 27 Newton coefficients in each tail.
  It additionally tests all 672 restricted gauges at four independent
  dilations and verifies the charged/neutral dimension arithmetic.  Normal,
  optimized and stored transcript replay and the LF hashes are required.
  Independent hostile audit is pending.
depends_on:
  - THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate
  - THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word
  - THM-3250-charged-heisenberg-blowup-address-intertwiner-and-pointed-multiplicity-gate
related:
  - THM-3252-singer-compactified-owner-hodge-word-universal-charged-cyclicity
  - THM-3247-heisenberg-central-fourier-decomposition-and-canonical-current-cyclicity
script: 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
output: 05-knowledge/results/lrc_positive_owner_mass_newton_cyclicity_thm3253.out
script_sha256: e90c5cee37b974ee9f1a9607df119b04f0dd7fc119e9dde42305ac5bc09be7fe
output_sha256: 1cb39817ae45f1e17a2c1c3d6166b7e1daf9ba1bd3cdfa425b5eacdc6cd50b4b
hash_basis: LF-normalized bytes
---

# THM-3253 -- positive owner-mass Newton cyclicity and maximal common Heisenberg module

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-3252 shows that the signed second-corrector word clears THM-3250's
charged determinant gate.  The determinant is not merely an asymptotic
shadow.  The actual positive overlap masses from THM-3246 are already cyclic
at every integer dilation in a natural 672-gauge family.

The resulting central-slice packet is nonnegative.  It does not give a
positive intertwiner to the blowup carrier: instead, its orbit realizes
exactly the largest source submodule which can occur on both sides of the
neutral-character mismatch.

## 1. Positive all-dilation owner masses

In THM-3246's lane

```text
(P,Q;e,f)=(3,5;1,2),
D_g=(504g-1)(840g-2),                                  (1)
```

write

```text
N_(g,j)=D_g I_(g,j),             0<=j<168.             (2)
```

For every integer `g>=1`, THM-3246 gives

| cells `j` | `N_(g,j)` |
|---|---|
| `0<=j<=5` or `162<=j<=167` | `12096g^2-1032g+2` |
| `6<=j<=23` or `144<=j<=161` | `12096g^2-24g` |
| `24<=j<=71` | `(16044-168j)g^2+48g` |
| `72<=j<=95` | `4032g^2+96g` |
| `96<=j<=143` | `(168j-12012)g^2+48g` |

Every displayed quadratic is strictly positive for `g>=1`: its value at
one and its derivative on `[1,infinity)` are positive.  Also

```text
N_(g,167-j)=N_(g,j).                                   (3)
```

Since `D_g>0`, the 168 masses `I_(g,j)` are strictly positive.

## 2. Restricted Singer gauges

Use THM-3234's deterministic plane

```text
F_169=F_13[u]/(u^2-2),          alpha=1+2u.             (4)
```

For

```text
a in {1,13,155,167}={+/-1,+/-13} mod 168,
b in Z/168Z,                                              (5)
```

define a mass matrix on `F_13^2` by

```text
A_g^(a,b)(0,0)=0,
A_g^(a,b)(alpha^(b+aj))=I_(g,j).                        (6)
```

Thus `(6)` is nonnegative, is strictly positive on the punctured plane, and
has only the abstract completion point as a zero.  The theorem is

```text
det A_g^(a,b) != 0                                      (7)
```

for every integer `g>=1` and all

```text
4*168=672                                                (8)
```

gauges in `(5)`.

This is deliberately a restricted family, not all 8,064 primitive Singer
generator gauges from THM-3252.

## 3. Phase and generator reductions

Let `B_g^b` be the cleared integer matrix obtained from `(6)` with `a=1`
and `I_(g,j)` replaced by `N_(g,j)`.  Put

```text
F_b(g)=det B_g^b.                                       (9)
```

Three exact symmetries reduce `(7)` to seven phase classes.

First,

```text
alpha^14=6 in F_13^*.                                   (10)
```

Changing `b` by 14 therefore applies the same scalar permutation to rows
and columns.  Its two determinant signs square, so

```text
F_(b+14)(g)=F_b(g).                                     (11)
```

Second,

```text
alpha^7=9u,
(x,y) |-> alpha^7(x+yu)=(5y,9x).                       (12)
```

This is transpose followed by multiplication permutations on the two axes.
Their signs are the Legendre symbols of 5 and 9, namely `-1` and `+1`.
Hence

```text
F_(b+7)(g)=-F_b(g).                                    (13)
```

It remains to handle `b=0,...,6`.

Finally, Frobenius is

```text
(x,y)->(x,-y),                  alpha^j->alpha^(13j).   (14)
```

It is a column permutation, so nonvanishing for multiplier 1 gives
nonvanishing for multiplier 13.  Reflection `(3)` gives, for either
`a=1` or `13`,

```text
A_g^(-a,b)=A_g^(a,b+a)                                  (15)
```

after relabelling `j` by `167-j`.  Equations `(11)--(15)` prove that the
seven phase certificates below imply all 672 cases.

## 4. Exact Newton certificates

Every entry of `B_g^b` is quadratic in `g`, so

```text
degree F_b <= 26.                                       (16)
```

For a polynomial of degree at most 26 and an integer base `r`, Newton's
identity is

```text
F_b(r+n)=sum_(k=0)^26 Delta^k F_b(r) binom(n,k)          (17)
```

for every integer `n>=0`.  The binomial coefficients are nonnegative.

Exact fraction-free determinants give:

1. For every `b=0,1,2,3,4`, all 27 numbers

   ```text
   Delta^k F_b(1),       0<=k<=26,                      (18)
   ```

   are strictly negative.  Their joint digest is

   ```text
   1aeb6d4070908447584ff0fee52c2ab7e7f0d2287f766bf69b962ef6f2815e16. (19)
   ```

   Hence `F_b(g)<0` for every integer `g>=1`.

2. For `b=5`, direct evaluation gives

   ```text
   F_5(g)<0,                 1<=g<=17.                  (20)
   ```

   All 27 coefficients `Delta^k F_5(18)` are strictly positive, so

   ```text
   F_5(g)>0,                 g>=18.                     (21)
   ```

   The exact head and tail digests are

   ```text
   b228dcdcca8a65cc11b61894e5cf07921238596ff03680873d88fe427d444274,
   f7a36db8a325e9bc576343143455035dab7e8bc1b16341159406c7427b0d52af. (22)
   ```

3. For `b=6`,

   ```text
   F_6(1)<0,                                               (23)
   ```

   while all 27 coefficients `Delta^k F_6(2)` are strictly positive.
   Their digest is

   ```text
   5fb6da1e11014bec9e2e97daa588f982ff318d6a437ebfc08da15141c2da40b4. (24)
   ```

   Thus `F_6(g)>0` for every `g>=2`.

Equations `(18)--(24)` prove `F_b(g)!=0` for the seven representatives;
Section 3 proves `(7)`.  Since

```text
det A_g^(a,b)=det B_g^(a,b)/D_g^13,                     (25)
```

the positive rational mass matrices have the same nonvanishing.

The sign changes between `g=17,18` in phase 5 and between `g=1,2` in phase
6 explain why a single base-one one-sign Newton certificate would be false.
The shifted tails are load-bearing, not cosmetic.

## 5. A nonnegative packet attaining the common-module ceiling

Fix any case in `(7)` and abbreviate its mass matrix by `A`.  On THM-3250's
exact-address carrier define the rational packet

```text
W_A=sum_(s,t in F_13) A_(s,t)[s,t,0].                  (26)
```

It is nonnegative.  In the central Fourier basis,

```text
[s,t,0]=(1/13)sum_(kappa in F_13) E_(s,t)^kappa.        (27)
```

Hence every charged block of `(26)` has coefficient matrix `A/13`.
By `(7)` and THM-3250 it is cyclic and spans all 169 dimensions.  The twelve
central idempotents lie in `K[H_13]`, so all block projections belong to the
orbit span.  The charged blocks therefore contribute

```text
12*169=2028.                                            (28)
```

In the neutral block `H_13` translates only the `s` coordinate and leaves
`t` as multiplicity.  For a Fourier row vector `ell_a` in `s`, invertibility
of `A` implies

```text
ell_a A != 0                         for every a.        (29)
```

Thus `(26)` has a nonzero component in each of the thirteen common neutral
characters `chi_(a,0)`.  An orbit in a one-dimensional isotypic character
contributes only one dimension, so its neutral orbit span is exactly 13.
Combining `(28)` and `(29)`,

```text
dim_K span(H_13.W_A)=2028+13=2041.                     (30)
```

THM-3250 proves that 2041 is the maximum rank of any equivariant map from
the exact-address carrier to the regular nonvertical blowup carrier.  The
module in `(30)` has all charged blocks plus one copy of every common neutral
character, so it is precisely a maximal common submodule.  The positive
packet uses the neutral overlap; it does not contradict or eliminate the
remaining 156-dimensional neutral mismatch.

## 6. Scope

The entries in `(6)` are genuine positive THM-3246 cell masses, but placing
them at Singer-plane points and on the central slice `delta=0` is an abstract
relocation.  No physical LRC owner-to-plane map, canonical endpoint packet,
Boolean observable, positive equivariant map, Markov clutch, or compatibility
with the full affine/Singer action is constructed.  In particular `(26)` is
not the canonical current reserved in THM-3247.

The theorem treats one ordered `(3,5;1,2)` lane and 672 restricted gauges.
It proves no row exclusion, arbitrary-radial NC2 theorem, Gaussian Moment
Conjecture, or `LRC(14)` decrement.  Its exact gain is that positivity of the
owner weights and maximal charged/neutral cyclic rank coexist at every
dilation.  The remaining obstruction is realization and equivariant
transport, not matrix rank.

## 7. Exact companion

Run

```text
python 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
python -O 04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py
```

and compare LF-normalized bytes with the declared output.  The companion
uses exact integer and finite-field arithmetic only, with no floating point,
randomness, discovery cache, or optimization-sensitive assertions.

QED, pending independent hostile audit.
