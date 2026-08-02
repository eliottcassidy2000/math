---
id: THM-3120
title: "Pole-prefix Newton flag and 115-support all-degree row positivity"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  For each of
  the 115 integer supports 1<=a<=10, a<b<=min(3a+4,21), both THM-3110 row
  generating functions admit an exact reduced positive pole-prefix flag.
  All 8,241 ordinary/confluent Newton coefficients are strictly positive,
  proving the complete-homogeneous row response in every degree n>=5 for
  those supports.  Arbitrary supports and nonrow histograms remain OPEN.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
related:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
script: 04-computation/gmc_product_gamma_row_pole_flag_thm3120.py
output: 05-knowledge/results/gmc_product_gamma_row_pole_flag_thm3120.out
script_sha256: 93828710e859ef2215ee164708b835d339594d241dcd70a090050ad3187901dd
output_sha256: af259dde1f130393d44792547104d535626bdadbe49fdc574ea41ada591248aa
hash_basis: LF-normalized bytes
---

# THM-3120 -- pole-prefix Newton flag and 115-support all-degree row positivity

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-3110 reduces each anchored product-Gamma width-three invariant to an
exact signed bank of residual root alphabets.  Its finite Schur certificate
reaches total degree twelve, and the row-ray scout remains positive far past
that range.  The row ray in fact has a finite exact positivity certificate:
reverse its reduced numerator, interpolate it on the descending pole flag,
and retain the resulting Newton coefficients.

The certificate below proves all degrees, but currently only on a stated
finite support universe.  Its value is not the size of that universe.  It
replaces an unbounded coefficient check by one finite Newton flag for each
support and exposes a precise universal target for arbitrary `a,b`.

## 1. The reduced row rational function

Fix one of the two THM-3110 banks `i in {1,2}` and an integer support

```text
0<a<b.                                                        (1)
```

For a bank row `R`, let `S_R` be its residual root multiset after the exact
THM-3110 chamber-common-root deletion, and let `epsilon_R` be its signed
integer coefficient.  Zeros in `S_R` may be retained or discarded.  Since

```text
sum_(n>=0) h_n(S_R)t^n=product_(r in S_R)(1-rt)^(-1),         (2)
```

the complete-homogeneous row response has generating function

```text
F_i(t)=sum_(n>=0) Phi_i(h_n)t^n
      =sum_R epsilon_R product_(r in S_R)(1-rt)^(-1).         (3)
```

Let

```text
m(r)=max_R mult_r(S_R),
D_raw(t)=product_(r>0)(1-rt)^m(r).                            (4)
```

Clear `(3)` over `D_raw` and cancel every common linear factor between its
numerator and denominator.  Write the resulting coprime fraction as

```text
F_i(t)=N(t)/D(t),                  D(0)=1.                    (5)
```

Every operation in `(3)--(5)` takes place in `Z[t]`.  No numerical root
finding or rational-function simplification is used.

In all cases certified below,

```text
N(t)=t^5P(t),                 d=deg P,
D(t)=product_(ell=1)^E(1-r_ell t),
r_1>=r_2>=...>=r_E>0,        E>=d+1.                         (6)
```

Repeated poles are repeated in the displayed list.

## 2. The pole-prefix flag

Put

```text
q_0(t)=1,
q_j(t)=product_(ell=1)^j(1-r_ell t).                         (7)
```

The polynomials

```text
B_k(t)=t^k q_(d-k)(t),                0<=k<=d,                (8)
```

form a triangular basis of the polynomials of degree at most `d`: `B_k`
has valuation `k` and coefficient one at `t^k`.  Hence there are unique
numbers `c_0,...,c_d` such that

```text
P(t)=sum_(k=0)^d c_k t^k q_(d-k)(t).                         (9)
```

The exact certificate is simply

```text
c_k>0,                         0<=k<=d.                      (10)
```

Dividing `(9)` by `D` gives the promised positive rational decomposition:

```text
F_i(t)=t^5 sum_(k=0)^d
        c_k t^k / [D(t)/q_(d-k)(t)].                         (11)
```

Every inverse factor on the right has nonnegative integer coefficients.
For `n=5+s` with `0<=s<=d`, the `k=s` summand contributes its positive
constant term.  For `s>d`, the `k=d` summand contributes a positive
coefficient of `1/D`; this is strictly positive because `E>=d+1` and every
remaining pole is positive.  Therefore `(10)` proves

```text
Phi_i(h_n)>0                    for every n>=5.               (12)
```

This is an all-degree argument.  It does not extrapolate from a finite
coefficient prefix.

## 3. Reciprocal Newton and confluent divided differences

The flag has a second, interpolation-theoretic form.  Reverse the numerator:

```text
Q(x)=x^dP(1/x).                                               (13)
```

Equation `(9)` is equivalent term by term to

```text
Q(x)=sum_(j=0)^d c_(d-j)
       product_(ell=1)^j(x-r_ell).                            (14)
```

Consequently

```text
c_(d-j)=Q[r_1,...,r_(j+1)],             0<=j<=d,              (15)
```

where a repeated-node divided difference is interpreted confluently:

```text
Q[r,...,r]  (j+1 copies) = Q^(j)(r)/j!.                       (16)
```

Thus the pole-prefix coefficients are not fitted positive weights.  They
are the unique ordinary/confluent Newton coefficients of the reciprocal
numerator on the first `d+1` nodes of the actual reduced pole multiset.
The companion independently computes both triangular expansions and checks
their equality exactly.

## 4. Minimal recurrence

Write

```text
D(t)=sum_(j=0)^E d_jt^j
    =sum_(j=0)^E (-1)^j e_j(r_1,...,r_E)t^j,   d_0=1,         (17)
a_n=Phi_i(h_n).
```

Then

```text
sum_(j=0)^E d_j a_(n-j)=0                 for n>5+d.          (18)
```

Because `(5)` is coprime, `(18)` has minimal constant-coefficient order
`E`.  This explains why a direct recurrence induction looked expensive:
even in the modest support bank below, `E` reaches `132` for `I_1` and
`133` for `I_2`.  The Newton flag proves positivity without asking the
alternating recurrence itself to preserve a positive cone.

## 5. Exact 115-support theorem

Let

```text
U={(a,b): 1<=a<=10, a<b<=min(3a+4,21)}.                      (19)
```

This set has `115` supports.  On both banks and every support in `U`, exact
integer arithmetic proves `(6)` and `(10)`.  The complete census is

```text
bank   cases  flag coefficients  degree d  recurrence E  min(E-d)
I1      115          3953          2..67       11..132         9
I2      115          4288          2..68       11..133         9.            (20)
```

All `8,241` flag coefficients are positive; none is zero or negative.  The
smallest coefficients are

```text
I1: c_0=36 at (a,b)=(1,2),
I2: c_0=32 at (a,b)=(1,2).                                  (21)
```

Combining `(11)--(12)` proves, for every `(a,b) in U`,

```text
Phi_1(h_n)>0 and Phi_2(h_n)>0             for every n>=5.     (22)
```

Since the THM-3110 forced divisors and distinguished-row normalizations are
positive, the corresponding normalized one-row responses are positive as
well.  No assertion of normalized row monotonicity is needed or made.

## 6. Smallest exact control

At `(a,b)=(1,2)`, both reduced denominators are

```text
D=(1-t)^4(1-2t)^3(1-3t)^2(1-4t)(1-5t),
q_1=1-5t,                 q_2=(1-5t)(1-4t).                  (23)
```

For `I_1`,

```text
P_1=36(1-3t-2t^2)
   =36q_2+216tq_1+288t^2.                                  (24)
```

For `I_2`,

```text
P_2=32+24t+40t^2
   =32q_2+312tq_1+960t^2.                                  (25)
```

The raw numerator in `(24)` already has two negative monomial
coefficients.  Positivity resides in the pole-adapted flag, not in the
ordinary monomial basis.

## 7. The exact cancellation boundary

It is unsafe to identify `r_1` with the largest pole of the uncleared
envelope.  There is a genuine cancellation family.  Put

```text
a=2k,              b=3k,              R=4k=2a.               (26)
```

The raw envelope has multiplicity three at `R`.  Only two `I_1` atoms can
contribute to the numerator at `t=1/R`; their bank coefficients are `2,-4`
and their local constants, normalized by the first, are

```text
1, (-1)^k/2.                                                (27)
```

Only three `I_2` atoms contribute; their coefficients and normalized local
constants are

```text
(1,-4,4),                  (1,(-1)^k/2,1/4).                 (28)
```

In either bank the normalized local residue is therefore

```text
2-2(-1)^k.                                                   (29)
```

The factor `(1-4kt)` cancels exactly when `k` is even.  Formula
`(27)--(29)` follows by grouping the missing root intervals in the two or
three displayed atoms; the companion checks every factor as an exact
rational product for `1<=k<=10`.

Within the finite universe `(19)`, the only cancellations are one copy of

```text
(1-8t)  at (a,b)=(4,6),
(1-16t) at (a,b)=(8,12),                                    (30)
```

in each bank.  The flag always uses the poles remaining *after* these
cancellations.

## 8. Pole deletion commutes algebraically, not positively

Multiplication by `1-Mt` is plethystic subtraction of one virtual letter:

```text
H_X(t) product (1-Mt)=H_(X-M)(t).                            (31)
```

This gives a clean comparison with labelled occupancy deletion.  If
`A` denotes the factorial-normalized labelled-deletion operator satisfying

```text
A Mtilde_(N+1)[X]=p_1[X]Mtilde_N[X],                         (32)
```

and `P_M` substitutes `X-M` for `X`, then naturality gives the strict
commutator identity

```text
A P_M Mtilde_(N+1)[X]
 =P_M A Mtilde_(N+1)[X]
 =(p_1[X]-M)Mtilde_N[X-M].                                  (33)
```

Every prefix `q_j` is a composition of such virtual-letter subtractions,
so it commutes in the same sense.  This identity is termwise and precedes
the signed bank summation.

It is not a positive deletion theorem.  At `(a,b)=(1,2)`, the top pole
`M=5` is absent from `21/24` `I_1` alphabets and `21/25` `I_2` alphabets.
Thus `(31)` cannot mean deletion of a letter physically present in each
atom.  Virtual alphabets are signed, and `(33)` supplies no Hasse-positive
current or row-to-nonrow transport.  This is exactly the distinction needed
when comparing the pole flag with the separately reserved labelled-deletion
program.

## 9. Sharp hostiles to stronger mechanisms

Three tempting strengthenings fail already at `(a,b)=(1,2)`.

1. Ordinary numerator coefficient positivity fails by `(24)`.
2. Sequential coefficientwise removal of the two largest poles fails:

   ```text
   [t^17](1-5t)(1-4t)F_1(t)=-5,901,696.                     (34)
   ```

   The full positive sum `(11)` survives; its individual prefixes do not
   form a nested coefficientwise-positive chain.
3. Giving independent cycle variables `z_j` to the power sums and expanding

   ```text
   sum_R epsilon_R exp(sum_(j>=1)p_j(S_R)z_j)                (35)
   ```

   does not produce a positive cycle cone.  In both banks the coefficient
   of `z_1z_2z_3` is exactly `-84`.

The mechanism is therefore specifically a one-variable ordered pole flag.
It neither termwise orients the signed atoms nor positively refines every
cycle coordinate.

## 10. Exact scope and next theorem target

This candidate proves:

1. the general implication “positive reduced pole flag implies all-degree
   row positivity”;
2. the exact Newton/divided-difference and suffix-denominator forms of that
   flag;
3. minimal row recurrences and the cancellation boundary; and
4. positivity of both complete-homogeneous row responses in every degree
   `n>=5` on all `115` supports in `(19)`.

It does **not** prove that the flag coefficients are positive for arbitrary
integers `0<a<b`.  It also does not prove that the row response is minimal
among partitions, that the normalized row response is increasing, that a
positive Hasse transport exists, or that any nonrow or multi-active-layer
histogram is positive.  In particular it does not close the residual
THM-3110 bank, arbitrary anchored product-Gamma width three, SFC, GMC(2),
NC2, LRC(14), JC(2), or DC(2).

The sharp universal next question is finite and explicit:

```text
For every 0<a<b, are all confluent divided differences
Q_i[r_1,...,r_(j+1)] strictly positive?                       (36)
```

A chamber formula or variation-diminishing theorem for `(36)` would extend
`(22)` from the exact support bank to every anchored triple without any
additional degree induction.

## 11. Exact companion

Reproduce the certificate with

```text
python 04-computation/gmc_product_gamma_row_pole_flag_thm3120.py
python -O 04-computation/gmc_product_gamma_row_pole_flag_thm3120.py
```

The companion reconstructs the THM-3110 banks, clears and reduces every
rational function over `Z[t]`, checks coprimality and minimal recurrence,
recovers all `8,241` flag coefficients twice (triangular pole basis and
confluent Newton differences), verifies the positive rational decomposition,
replays the exact cancellation products and hostile controls, and pins the
complete census.  Both executions must byte-match the stored output.

**QED (candidate pending independent audit).**
