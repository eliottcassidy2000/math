---
id: THM-3432
title: "Order-two fixed/half parity transplant"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  An elementary all-modulus
  parity decomposition identifies a
  fixed-zero cover on 2R containing the order-two block with a half-twist
  cover on R.  Joint quotient period is preserved, but augmented
  primitiveness needs one parity breaker when R is odd.  No mixed rank-seven
  classification or LRC(14) consequence is claimed.
source: root-lrc-order-two-parity-transplant-2026-08-15
audit: independent strict-endpoint, parity-chart, rank, period/gcd, normalized-order, affine-unit, multiplicity-loss, R15, Q474, solver, hash, normal/-O, security, and documentation audit CLEAN
depends_on: []
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
script: 04-computation/lrc_order2_fixed_half_parity_transplant_thm3432.py
output: 05-knowledge/results/lrc_order2_fixed_half_parity_transplant_thm3432.out
script_sha256: 5be8316df17232cb797c830207190b793faf2088a9d9c172a509f6ef986207e9
output_sha256: f4f4ccfccc0ff82a7920abb9d69245539c0eb498a383fe4895a5c08f01088ca9
semantic_sha256: 509ab5d7e01e3ffad5ba93bbfdd9081c84697db9457213658d750bbf40c681a9
hash_basis: LF-normalized bytes
---

# THM-3432 -- order-two fixed/half parity transplant

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Write `||x||` for distance to the nearest integer.  Fix an integer `R>=2`
and a strict arc radius

```text
0<delta<=1/2.                                               (1)
```

For a residue `r modulo 2R`, define the fixed-zero and half-twist blocks

```text
Z_n^delta(r)={ell in Z/nZ: ||r ell/n||<delta},
H_R^delta(r)={j in Z/RZ: ||r(2j+1)/(2R)||<delta}.          (2)
```

Let `e(j)=2j` and `o(j)=2j+1` map `Z/RZ` onto the even and odd sheets of
`Z/(2R)Z`.  Then, for every `r`,

```text
e^(-1) Z_(2R)^delta(r)=Z_R^delta(r),
o^(-1) Z_(2R)^delta(r)=H_R^delta(r),                       (3)

Z_(2R)^delta(r)
 =e(Z_R^delta(r)) disjoint_union o(H_R^delta(r)),
Z_(2R)^delta(R)=e(Z/RZ).                                  (4)
```

Consequently, for every finite labelled family `A` of distinct transverse
residues modulo `2R` (`R` does not divide any `r in A`),

```text
{R} union A covers Z/(2R)Z at fixed zero
iff
A covers Z/RZ at half twist.                               (5)
```

If `rho_H(R,delta)` is the minimum size of such a half family and
`rho_(Z,2)(2R,delta)` is the minimum fixed-zero size **subject to containing
the order-two residue `R`**, then

```text
rho_(Z,2)(2R,delta)=1+rho_H(R,delta).                       (6)
```

The equality includes the value infinity.  It is not a statement about the
unconstrained fixed-zero minimum.

At the LRC radius `delta=1/14`, `(3)--(6)` are the promised fixed/half parity
transplant.

## 2. Proof of the parity decomposition

On an even sheet,

```text
||r e(j)/(2R)||=||rj/R||,                                  (7)
```

and on an odd sheet the left side is exactly the defining expression for
`H_R^delta(r)`.  This proves `(3)` and the first line of `(4)` without any
endpoint relaxation.

For `r=R`, the phase on sheet `ell` is `ell/2`.  It has distance zero for
even `ell` and distance `1/2` for odd `ell`.  Because the inequality is
strict and `delta<=1/2`, its block is exactly the even sheets.  This proves
the second line of `(4)`.

The residue `R` therefore covers the whole even chart, while the remaining
odd chart is covered exactly when the half blocks of `A` cover.  This proves
`(5)`.  Deleting or adjoining the one labelled residue `R` is inverse on the
two cover families, so minima differ by exactly one, proving `(6)`.  **QED.**

## 3. Multiplicity, XOR, and the discarded zero chart

The transplant is stronger than an OR identity but weaker than a partition
identity.  For a family `A`, let `mu_Z`, `mu_0`, and `mu_H` be the block
multiplicities of `{R} union A`, the downstairs zero family `A`, and the
downstairs half family `A`.  Equation `(3)` gives pointwise

```text
mu_Z(e(j))=1+mu_0(j),
mu_Z(o(j))=mu_H(j).                                       (8)
```

Thus the half chart retains every odd-sheet overlap, but it forgets all
even-sheet overlaps contributed by `A`.  In particular, half OR coverage is
equivalent to upstairs OR coverage in `(5)`, while XOR or exact-partition
claims additionally require the zero-chart multiplicities.  This is the
precise information loss in the Q474 order-two reduction; an odd-sheet
certificate alone is not a full multiplicity certificate upstairs.

## 4. Joint period and the augmented parity breaker

Put

```text
h=gcd(R,r:r in A),
m_R(r)=R/gcd(R,r).                                        (9)
```

The elementary lcm/gcd identity gives

```text
lcm_(r in A) m_R(r)=R/h,
lcm(ord_(2R)(R),ord_(2R)(r):r in A)=2R/h.                (10)
```

Hence the half family has joint quotient period `R` if and only if the
upstairs family has joint additive period `2R`.  This primitive-period
predicate is preserved exactly.

The augmented gcd gate is subtler.  Its two values are

```text
g_fixed=gcd(2R,R,r:r in A)=gcd(R,r:r in A)=h,
g_half =gcd(2R,r:r in A).                                 (11)
```

Assume the equivalent joint-period conditions, so `h=1`.  If `R` is even,
some selected `r` must be odd and `g_half=1`; fixed and half augmented gates
agree.  If `R` is odd, however,

```text
g_half=2 iff every r in A is even,
g_half=1 otherwise.                                      (12)
```

Upstairs, the odd order-two residue `R` itself breaks the prime two; after
deleting it, that breaker can disappear.  Therefore augmented primitiveness
is preserved for even `R`, while odd `R` requires the extra sidecar “some
half residue is odd.”

The smallest useful hostile is `R=15`:

```text
A=(2,4,6,8,10,14)                                        (13)
```

covers the half layer and has joint quotient period `15`, but
`gcd(30,A)=2`.  The transplanted fixed family

```text
(15,2,4,6,8,10,14)                                       (14)
```

covers `Z/30Z` and has augmented gcd one.  Thus replacing `(11)` by a claim
that the augmented breaker is automatically preserved would be false.

## 5. The order-doubling sidecar

For one residue put

```text
d_r=gcd(R,r),        m=R/d_r,        s=r/d_r,
gcd(m,s)=1.                                             (15)
```

Then

```text
ord_(2R)(r)=2m/gcd(2,s)
 =m   if s is even,
 =2m  if s is odd.                                      (16)
```

In the first branch `m` is necessarily odd.  Thus the downstairs quotient
order `m` does not determine the upstairs additive order: the parity of the
normalized coefficient `s` is the missing bit.  It is also the branch bit
behind the two formulas for maximal half-block size in THM-3416.

## 6. Unit and reflection covariance

Let `u` be a unit modulo `2R`.  It is odd, fixes the order-two residue
`uR=R`, and acts on the odd chart by the affine permutation

```text
phi_u(j)=u j+(u-1)/2  (mod R),
u(2j+1)=2phi_u(j)+1.                                    (17)
```

Accordingly,

```text
Z_(2R)^delta(ur)=u^(-1) Z_(2R)^delta(r),
H_R^delta(ur)=phi_u^(-1) H_R^delta(r).                  (18)
```

The transplant therefore respects the full common-unit orbit, but the
downstairs action is affine rather than purely multiplicative.  The unit
`u=-1` gives the half reflection `j -> -j-1`.  Coefficient sign itself leaves
each block unchanged because `||-x||=||x||`.

## 7. Exact boundaries and the Q474 application

Direct independent minimum-cover searches at `delta=1/14` give:

| `R` | `rho_H(R)` | fixed minimum containing order two |
|---:|---:|---:|
| 8, 9 | 4 | 5 |
| 10, 12 | 5 | 6 |
| 11, 15, 23, 25 | 6 | 7 |
| 13, 29 | 7 | 8 |

These are controls, not bounded evidence for `(6)`, whose proof is the set
identity `(3)`.  The named half atoms from THM-3415/3416/3421 transplant
verbatim by adjoining `R`.

The FINITE-EXACT
[`Q=474` certificate](../../04-computation/fixed_zero_q474_primitive_cap7_exact.py)
uses precisely the `R=237` instance.  Its order-two fixed block pays all even
sheets, and the odd restrictions of the other six blocks form a half-twist
six-cover problem on `237`.  The independent half-`237` fresh-gain deficit in
that artifact therefore excludes the entire order-two branch upstairs.
THM-3432 supplies the all-modulus mechanism; it does not promote the Q474
finite decision to an all-q exclusion.

## 8. Connection and loss ledger

| field | exact content |
|---|---|
| source | labelled fixed-zero family on `2R` containing residue `R` |
| target | paired zero/half families on `R` |
| map | restrict sheets along `e(j)=2j`, `o(j)=2j+1` |
| preserved | labels, strict endpoints, odd-chart OR/multiplicity, joint period, common-unit orbit |
| destroyed by half projection | even overlaps, fixed XOR/partition data, and possibly the prime-two augmented breaker |
| restoring sidecar | downstairs zero multiplicities plus normalized-coefficient parity |
| cheapest hostile | all-even `R=15` family `(13)--(14)` |

This is a parity decomposition, not a tournament: there is no intrinsic
orientation between owners.

## 9. Exact companion and scope

The standard-library companion verifies `(3)--(4)`, `(8)`, `(12)`, `(16)`,
and `(18)` for five strict arc widths through `R=96` (`5,990,700` sheet
cells), then independently searches the ten boundary ranks in Section 7.
Ordinary and optimized runs byte-match the stored output:

```bash
python3 -B 04-computation/lrc_order2_fixed_half_parity_transplant_thm3432.py
python3 -B -O 04-computation/lrc_order2_fixed_half_parity_transplant_thm3432.py
```

Nothing here says every fixed-zero cover contains the order-two owner, that
the unconstrained fixed minimum equals a half minimum plus one, that half
projection preserves XOR or augmented primitiveness without its sidecars,
that mixed quotient-order rank seven is classified, that an arbitrary common
physical time is one of the two THM-3405 gauge layers, or that `LRC(14)` is
proved.
