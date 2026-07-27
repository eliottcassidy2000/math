---
id: THM-2527
title: "Owner-weighted all-mode odd bank and Boolean cut coordinate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A rational positive
  BV owner placed after an active base-thirteen collision need only make the
  total weighted replica defect positive: rationality and Phi_13
  irreducibility then force all twelve weighted root colours positive at
  once.  Applying the guard-selected A_tau H_tau transform on the same
  ancestry integrand gives a nonzero, lossless, reflection-odd bank with
  explicit BV thresholds and norm floors.  At the live Boolean depth-one
  collision, the fixed coordinate -4 tau is pointwise nonnegative on every
  root fibre.  Its cyclic cut score is zero exactly for constant masks and
  is at least 13/42 of the fibre drift, sharply.  A Boolean owner therefore
  gives a finite positive Boolean layer factorization of that coordinate.
  This couples the late owner to the odd bank, but the factor observes the
  full thirteen-root mask rather than one semantic source--arrival pair; no
  owner loop, row exclusion, or LRC(14) proof is claimed.
source: codex-2026-07-27-owner-weighted-odd-bank
depends_on:
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2525-unit-guard-collision-floor-and-word-owner-cross-cospan-collapse
  - THM-2526-affine-skew-orientation-gauge-boundary
related:
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2524-translated-chi7-hamilton-polarization-inversion
script: 04-computation/lrc14_owner_weighted_odd_bank_thm2527.py
output: 05-knowledge/results/lrc14_owner_weighted_odd_bank_thm2527.out
script_sha256: c7ec707e6912e99ebf076ce7765ead131daf369923e508cd6fa0c7b5cdd5ae49
output_sha256: 9fccaf4e7e527edb8932f67440992b2a36cbe25db21588e6bfd29b0c849f86d6
hash_basis: working-tree bytes (LF)
---

# THM-2527 -- a late owner survives the guard-oriented collision transform

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2522 separated the intrinsic collision clock from a freely delayed
owner, while THM-2526 supplied a lossless odd transform on the retained
guard sheet.  The missing compatibility is exact.  The late owner is
constant on every root translation used by the collision table, so it can
be kept inside the same integral while the odd transform is applied.

There are two distinct conclusions.

```text
rational owner weighting + positive total defect
  -> all twelve Fourier colours positive
  -> a lossless all-mode odd bank;

Boolean live root mask
  -> one fixed odd coordinate is pointwise nonnegative
  -> a finite positive Boolean layer factorization.                    (1)
```

The second statement is stronger than merely choosing a positive entry of
a signed vector after integration.  It is still a thirteen-predecessor
mask predicate, not a single ordered source--arrival cospan.

## 1. The total weighted defect is the right mixing target

Use THM-2522's Perron notation.  Let `F` be a rational step function with

```text
0<=F<=B,                  V=Var(F),

f=F_m=P^mF,
D_m=(I-E)f,
d_m=||D_m||_2^2>0.                                             (2)
```

Let the future owner, or the entire future owner--word block, be a rational
BV step function

```text
0<=G<=1,             rho=integral_T G>0,             W=Var(G). (3)
```

Put it `R` steps after the fixed collision and write

```text
K=13^(R+1),
w(y)=G(Ky).                                                   (4)
```

The owner-weighted collision profile of THM-2522 equations (48)--(50) is

```text
B^[m,R]_u
 =integral_T w(y)f(y)f(y+u/13)dy,               u in F_13,    (5)

q_a:=Bhat^[m,R](a)
 =integral_T w(y)|D_(m,a)(y)|^2dy>=0,
                                      a in F_13^*.             (6)
```

Let

```text
Bbar=1/13 sum_u B^[m,R]_u,
b_u=B^[m,R]_u-Bbar.                                           (7)
```

Since `K` is divisible by `13`,

```text
w(y+u/13)=w(y).                                               (8)
```

Changing variables by a root translation in (5) proves

```text
B^[m,R]_(-u)=B^[m,R]_u,
b_(-u)=b_u.                                                   (9)
```

The same invariance orthogonalizes the twelve character pieces even with
the owner weight present.  If `a!=c`, translate by `v/13` in
`integral w D_(m,a) conjugate(D_(m,c))`; the integral is multiplied by a
nontrivial thirteenth root and hence vanishes.  Using
`D_m=sum_(a!=0)D_(m,a)` gives the exact total identity

```text
b_0
 =B^[m,R]_0-Bbar
 =sum_(a!=0)q_a
 =integral_T w(y)|D_m(y)|^2dy.                              (10)
```

This is the owner-weighted replica defect, not merely a sum of unrelated
lower bounds.

Perron contraction and root averaging give

```text
||D_m||_infinity<=B,
Var(D_m)<=2V/13^m,
Var(|D_m|^2)<=4BV/13^m.                                     (11)
```

Apply the periodic two-BV covariance estimate with frequency `K`:

```text
|b_0-rho d_m|
 <= BVW/(3*13^m*K)
 =:epsilon_(m,R).                                           (12)
```

In particular,

```text
b_0>=rho d_m-epsilon_(m,R).                                 (13)
```

Thus the explicit total-defect threshold

```text
K>BVW/(3 rho 13^m d_m)                                     (14)
```

makes `b_0>0`, while the doubled threshold

```text
K>2BVW/(3 rho 13^m d_m)                                    (15)
```

gives

```text
b_0>=rho d_m/2.                                             (16)
```

Equivalently, the first admissible delay in (14) is

```text
R_all=min{R>=0:13^(R+1)>BVW/(3 rho 13^m d_m)}.              (17)
```

No asymptotic choice is hidden in this definition.

## 2. Rationality turns total positivity into twelve-colour positivity

Every entry of (5) is rational.  Suppose one nontrivial normalized Fourier
coefficient vanished.  Then

```text
sum_(u=0)^12 B^[m,R]_u X^u                                  (18)
```

would vanish at a primitive thirteenth root.  Irreducibility of
`Phi_13(X)=1+X+...+X^12` over `Q` forces all thirteen coefficients in
(18) to be equal.  That would give `b_0=0`, contrary to (14).  Together
with the nonnegativity in (6), this proves

```text
b_0>0 and B^[m,R] rational
  -> q_a>0                         for every a!=0.            (19)
```

This is a genuine strengthening of applying THM-2522 equation (53) to each
colour separately.  The old route paid for
`e_m=min_a||D_(m,a)||_2^2`, which may be tiny.  Equations (10), (12), and
(19) pay only for the total energy `d_m` and nevertheless recover all
twelve modes.

At THM-2522's endpoint horizon, equation (46) there gives

```text
d_m=sum_(a!=0)||D_(m,a)||_2^2
 >=3 S_C/(13^2*13^(2m)*d).                                  (20)
```

Hence a denominator-explicit sufficient form of (14) is

```text
K>13^2*d*13^m*B*V*W/(9 rho S_C),                            (21)
```

and doubling its right side gives the half-drift conclusion (16).

Rationality in (19) is load-bearing.  For a general real profile, positive
total energy does not prevent selected root characters from vanishing.

## 3. The owner-weighted odd bank lives on the same integrand

Fix a nonzero slope `tau`.  In the live application it will be the
guard-selected slope `tau_H` of THM-2526.  Let `A_tau` be the symmetric
Fano--Hamilton operator of THM-2523 and let

```text
H_tau=sum_(s=1)^6(T_(-2tau s)-T_(2tau s))                   (22)
```

be THM-2526's cyclic Hilbert/tournament operator.  Define

```text
R_tau=13A_tau b,
O_tau=H_tau R_tau=13A_tau H_tau b.                           (23)
```

For each `y`, put

```text
c_y(u)=f(y)f(y+u/13).                                       (24)
```

Because the operators act only on the root-displacement coordinate, (5)
and (23) give the same-ancestry formula

```text
O_tau(t)
 =13 integral_T w(y)(A_tau H_tau c_y)_t dy.                 (25)
```

The future owner has not been averaged away or reattached after the fact:
the identical factor `w(y)` occurs in every signed term of the odd bank.

Let `zeta=zeta_13`.  Equations (9), (19), and the Fourier multipliers of
THM-2523/2526 give

```text
Ohat_tau(a)
 =13 lambda_(tau,a)
      [(zeta^(a tau)-1)/(zeta^(a tau)+1)]q_a,
                                      a!=0.                  (26)
```

Every factor on the right is nonzero after (14).  Moreover reflection
`Jx(t)=x(-t)` satisfies

```text
JO_tau=-O_tau.                                               (27)
```

Thus (14) produces a nonzero reflection-odd bank with all twelve primitive
Fourier modes.  The transform is lossless because both `A_tau` and `H_tau`
are invertible on the augmentation module.

There are useful norm floors before using Booleanity.  At `tau=1`, put

```text
ell=1/845*(0,-17,-37,-37,-33,-15,1,-1,15,33,37,37,17).      (28)
```

For `F_1=A_1H_1`, exact multiplication gives

```text
ell^T(13F_1)=e_0^T-(1/13)1^T.                               (29)
```

Since `b` is centred, (23) and (29) imply

```text
b_0=ell dot O_1.                                             (30)
```

Multiplicative relabelling gives a row `ell_tau` with the same norms for
every nonzero `tau`.  The exact values

```text
||ell_tau||_1=56/169,
||ell_tau||_2^2=8684/845^2                                  (31)
```

therefore yield, with the unnormalized Euclidean norm,

```text
max_t |O_tau(t)|
 >=169/56 (rho d_m-epsilon_(m,R)),                           (32)

||O_tau||_(ell^2)
 >=845/(2 sqrt(2171)) (rho d_m-epsilon_(m,R)).               (33)
```

By oddness, the maximum absolute value in (32) occurs with both signs.
Under (15), the right sides may respectively be replaced by

```text
169 rho d_m/112,
845 rho d_m/(4 sqrt(2171)).                                  (34)
```

## 4. Booleanity selects one fixed positive odd coordinate

Now specialize to the live first collision:

```text
m=0,
F:T->{0,1}.                                                  (35)
```

Disintegrate one root fibre by `z=13y` and put

```text
e_r(z)=F((z+r)/13),
n(z)=sum_r e_r(z),
c_u(z)=sum_r e_r(z)e_(r+u)(z),                               (36)

g_R(z)=G(13^R z).                                            (37)
```

The root invariance (8) gives

```text
B^[0,R]_u=1/13 integral_T g_R(z)c_u(z)dz.                   (38)
```

For a Boolean mask `e in {0,1}^13`, define its guard-oriented cut score by

```text
Psi_tau(e)=(A_tau H_tau c)_(-4tau)
          =-(A_tau H_tau c)_(4tau).                         (39)
```

At slope one this is the explicit integer

```text
Psi(e)
 =7c_0-12c_1+8c_2-6c_3+7c_4-6c_5+2c_6.                    (40)
```

The exact Boolean cut lemma is

```text
0<=Psi(e)<=98,

Psi(e)=0 iff n in {0,13},

42 Psi(e)>=n(13-n).                                         (41)
```

The last constant is sharp.  Equality away from constant masks occurs for
exactly `52` masks, all with `n=6` or `n=7`.  They form four translation
necklaces, represented in root order by

```text
1111000110000,
1111001110000,
1100011110000,
1110011110000.                                               (42)
```

Their alternating run lengths are rotations/reflections of
`(4,3,2,4)`, `(4,2,3,4)`, `(2,3,4,4)`, and `(3,2,4,4)`;
complement pairs the first with the fourth and the second with the third.

There is a useful cut-polytope description.  Put

```text
delta_d(e)=sum_r(e_r-e_(r+d))^2.                             (43)
```

Then

```text
2Psi
 =12delta_1-8delta_2+6delta_3-7delta_4+6delta_5-2delta_6.    (44)
```

Thus (41) is a cyclic cut inequality, not positive semidefiniteness on the
whole real root module.  The exact referee checks every one of the `8192`
vertices of this cut polytope, including the zero classification and the
four equality necklaces.  Once strict positivity and integrality are
known, the relative constant has a short mechanism:

```text
Psi>=1 on a mixed mask,
n(13-n)<=42,

therefore Psi>=n(13-n)/42.                                  (45)
```

Equality requires `Psi=1` and `n in {6,7}`, explaining the boundary in
(42).  Booleanity is essential: the rational `[0,1]` vector

```text
(1/2,0,0,0,0,1,1,1,1,0,0,0,0)                              (46)
```

has the same quadratic score `Psi=-1/4`.

Multiplicative relabelling proves (39)--(41) for every nonzero `tau`.
Combining (23), (38), and (39) cancels the branch factor `13` and gives the
fixed positive coordinate

```text
O_tau(-4tau)
 =integral_T g_R(z)Psi_tau(e(z))dz
 >=0.                                                        (47)
```

It vanishes exactly when the owner weight is supported almost everywhere
on constant root fibres.  Since

```text
b_0=1/169 integral_T g_R(z)n(z)(13-n(z))dz,                  (48)
```

the sharp pointwise cut inequality gives

```text
O_tau(-4tau)>=169/42 b_0.                                   (49)
```

This improves the adaptive general coordinate floor (32) on the Boolean
depth-one sheet, and names the positive coordinate before integration.

## 5. The positive coordinate has a finite Boolean lift

Suppose now that `G` is itself a Boolean owner or owner--word block.  Since
`Psi` is integer-valued between zero and `98`, pointwise

```text
g_R(z)Psi_tau(e(z))
 =sum_(j=1)^98
    1_{g_R(z)=1 and Psi_tau(e(z))>=j}.                       (50)
```

Therefore

```text
O_tau(-4tau)
 =sum_(j=1)^98 measure{
      z:g_R(z)=1 and Psi_tau(e(z))>=j}.                      (51)
```

This is a positive Boolean factorization on the same ancestry fibre.  The
first layer is simply

```text
{late owner occurs and the collision root mask is nonconstant}.       (52)
```

Whenever (14) holds, at least one layer in (51), and in fact the first
layer, has positive measure.

The factorization does **not** say that one original pair term in the
signed formula for `A_tau H_tau` equals the whole coordinate.  It observes
all thirteen predecessor bits and repackages their cyclic cut score.  It
therefore supplies a Boolean collision-mask witness coupled to the owner,
but not yet a distinguished old-to-new, source-to-arrival, or
blocker-to-owner edge.

## 6. Uniform application to all `165` live rows

For every live THM-2349 owner--word event, THM-2522 proves `m=0`, while
THM-2525 gives, with `mu=integral F`,

```text
d_0=||D_0||_2^2>=3mu/13.                                    (53)
```

Here `B=1`.  Choose the lawful guard slope `tau=tau_H` from THM-2526.  A
row-explicit sufficient all-mode threshold is therefore

```text
K=13^(R+1)>13 V W/(9 rho mu).                               (54)
```

At every such delay, the same owner-weighted ancestry integrand produces
all of the following simultaneously:

```text
q_a>0 for every a!=0;

O_(tau_H) reflection-odd and nonzero on all twelve modes;

O_(tau_H)(-4tau_H)>0;

the positive Boolean layer factorization (51).               (55)
```

More quantitatively, equations (12), (49), and (53) give

```text
O_(tau_H)(-4tau_H)
 >=13 rho mu/14-169 V W/(126K).                              (56)
```

Under the doubled explicit threshold

```text
K>26 V W/(9 rho mu),                                        (57)
```

this becomes the uniform relative floor

```text
O_(tau_H)(-4tau_H)>=13 rho mu/28.                            (58)
```

The old shallow and deep source labels remain retained because the
collision is still the same THM-2522 depth-one collision and the owner was
only moved later.  No rebase or new ancestry address occurs.

## 7. Exact boundary of the advance

This theorem closes the algebraic owner-coupling seam left open in
THM-2526 and does slightly more: on the live Boolean sheet it turns one
fixed odd coordinate into a positive finite sum of Boolean mask events.
Three boundaries remain load-bearing.

1. The direction `-4tau_H` is transported from the physical guard sheet.
   Averaging the two converse guard conventions still kills the odd bank.
2. The Boolean factor in (51) is a predicate of the full root mask.  It is
   not one ordered predecessor pair and does not by itself emit a semantic
   owner loop.
3. Positivity of this mask event is another expression of nonreplication.
   No theorem here identifies it with a forbidden scalar-cover
   configuration or removes one of the `165` rows.

Thus the remaining target is no longer “make the late owner appear in an
odd Boolean object.”  It is to map one positive layer of (51), with the
retained source/deep sidecars, into a scalar-cover contradiction or a
genuine ordered owner loop.

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_owner_weighted_odd_bank_thm2527.py
python3 -O 04-computation/lrc14_owner_weighted_odd_bank_thm2527.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_owner_weighted_odd_bank_thm2527.out
```

byte-for-byte.  The referee uses only integer and rational arithmetic.  It
checks:

- the exact `A_tau H_tau` matrices, rank, covariance, and inverse row;
- every Boolean root mask, the cyclic cut formula, strict/zero boundary,
  sharp `13/42` floor, and all four equality necklaces;
- the `5434` nonempty masks having a cyclic three-zero run, including all
  `52` sharp equality masks;
- the non-Boolean hostile control (46);
- the positive layer identity and a rational weighted-disintegration
  control on the same integrand.

The BV estimates and cyclotomic step are symbolic proofs in Sections 1--2;
the finite cut lemma is the exhaustively verified part.  An independent
audit reproduced the BV constants, cyclotomic implication, endpoint
substitution, inverse-row norms, live `13/28` floor, and the byte-identical
normal/optimized/stored referee output.
