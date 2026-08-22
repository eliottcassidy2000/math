---
id: THM-3425
title: "Half-twist rank-six primitive breaker-profile closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Outside quotient degrees
  divisible by 8, 9, 10, or 12, a half-twist cover by at most six blocks has
  full quotient-order joint period exactly at Q=11,15,22,23,25.  Consequently
  the primitive augmented half-twist cap-six support is exactly the multiples
  of 8,9,10,12 together with those five isolated degrees.  In the same
  lower-base-free sector, primitive fixed-zero cap-six support is exactly
  Q=15.  This is not a standalone all-Q fixed-zero classification and gives
  no arbitrary-time or LRC(14) conclusion.
source: codex2-major-frontiers-2026-08-15
audit: independent proof reconstruction CLEAN after family-level parity and candidate-type wording repairs; Q27 positive and Q33/Q46 hostile checks; normal/optimized/stored-output replay
depends_on:
  - THM-3414-fixed-zero-six-owner-base-classification
  - THM-3416-zero-mode-cochain-global-rank-six-support
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures
script: 04-computation/lrc_half_twist_rank6_joint_period_probe_20260815.py
output: 05-knowledge/results/lrc_half_twist_rank6_joint_period_probe_20260815.out
script_sha256: df87baa9e752b3bda4f55fe29f7bda0219b2fc7978ffc47c3cd37c0f4319cf75
output_sha256: c643d6ddd1c3fe03daa06850cf8aae91287eab90229bcf2c95ffa0469a8dee79
semantic_sha256: 4ac644546e4f81631b5c12a404779397c5c63da67f3e4dd5c1b88da1fa8beda1
hash_basis: LF-normalized bytes
---

# THM-3425 -- half-twist rank-six primitive breaker-profile closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and scope

For a half-twist block on `Q` sheets with residue `r`, put

```text
m(r)=Q/gcd(Q,r).                                         (1)
```

This is its quotient order.  For a selected family `R`, define its joint
period

```text
L(R)=lcm_(r in R) m(r).                                  (2)
```

Retain the literal strict masks `B^1_(Q,r)` from THM-3416.  If none of
`8,9,10,12` divides `Q`, then

```text
|R|<=6,       union_(r in R) B^1_(Q,r)=Z/QZ,
L(R)=Q

iff Q in {11,15,22,23,25}.                              (3)
```

Let `r^prim_(1/2)(Q)` be the minimum half-twist augmented cover rank, where
the positive-lift breaker is

```text
gcd(2Q,r:r in R)=1.                                     (4)
```

Then the all-`Q` consequence is

```text
r^prim_(1/2)(Q)<=6
iff 8|Q or 9|Q or 10|Q or 12|Q
    or Q in {11,15,22,23,25}.                           (5)
```

Finally, in the same lower-base-free sector,

```text
r^prim_0(Q)<=6 iff Q=15.                                (6)
```

Formula `(6)` is not promoted to a standalone all-`Q` fixed-zero formula.
The theorem concerns the primitive Boolean layers of zero complete mode
cochains.  It does not classify arbitrary physical common times or close an
LRC(14) row.

## 2. The breaker is joint period plus one parity bit

Prime by prime,

```text
v_p(m(r))=v_p(Q)-min(v_p(Q),v_p(r)).                    (7)
```

Therefore

```text
L(R)=Q/gcd(Q,r:r in R),                                 (8)
```

and in particular

```text
gcd(2Q,R)=1
iff L(R)=Q and at least one residue in R is odd.         (9)
```

When `Q` is even, `L(R)=Q` already forces an odd residue.  When `Q` is odd,
the parity sidecar is genuine.  At `Q=15`, the all-even family

```text
(2,4,6,8,10,14)                                        (10)
```

covers all sheets and has joint period 15, but its augmented gcd is two.
Thus `(3)` is stronger than needed for the negative implication in `(5)`,
while existence in the five positive degrees still requires an explicit odd
breaker.

## 3. Positive families and lower-base pullback repair

The lower-base-free atoms are:

| `Q` | half-twist residues | quotient-order profile |
|---:|---|---|
| 11 | `(1,2,3,5,7,9)` | `11^6` |
| 15 | `(1,4,6,7,8,10)` | `(3,5,15^4)` |
| 22 | `(1,4,9,10,13,21)` | `(11^2,22^4)` |
| 23 | `(1,4,5,7,9,11)` | `23^6` |
| 25 | `(1,9,10,11,19,21)` | `(5,25^5)` |

Every row contains residue one, covers literally, has joint period `Q`, and
passes `(4)`.

For the four lower bases use the primitive atoms

```text
8:  (1,3,5,7),             9:  (1,5,6,7),
10: (1,3,4,7,9),           12: (1,5,7,8,11).           (11)
```

At `Q=ka`, scaling the residues in the `a`-atom by `k` pulls back the sheet
cover but leaves joint period only `a`.  If `k>1`, adjoining residue one
restores joint period `Q` and augmented gcd one.  This uses five owners above
the rank-four atoms and six above the rank-five atoms.  It proves the positive
direction of `(5)` for every multiple, not merely at the four bases.

The useful correction control is

```text
Q=27: (1,3,15,18,21),      orders (3,9,9,9,27).        (12)
```

It is a primitive five-cover.  The fact that `Q=27` does not occur in a list
of *exact rank-six* primitive instances is not evidence against cap-six
support.

## 4. Capacity leaves eight high-density orders

Let `h(m)` be THM-3416's exact maximum half-twist block size on an order-`m`
quotient.  Its formula and bound are

```text
h(m)<= (m+6)/7.                                         (13)
```

Six-cover mass, followed by exclusion of orders carrying one of
`8,9,10,12`, leaves

```text
{3,5,11,15,17,22,23,29}.                               (14)
```

The order-three and order-five anchored complement classifiers inherited
from THM-3416 sharpen as follows in the present sector:

- an order-three anchor needs `2/15` from each remaining slot; equality is
  possible only at orders 5 and 15, so full joint period forces `Q=15`;
- with no order three, an order-five anchor needs `4/25`; equality is possible
  only at order 25, so full joint period forces `Q=25`.

It remains to remove maximal orders 15, 17, and 29.  For a maximal order-`a`
anchor, let `c(a,m)` be the largest fraction of an order-`m` block lying in
the anchor complement after pullback to `lcm(a,m)`.  The five-slot quotas are

```text
theta_15=4/25,       theta_17=14/85,       theta_29=24/145. (15)
```

Exact finite evaluation plus `(13)` gives:

- for `a=15`, only order 3 reaches the quota;
- after orders 3 and 5 are removed, an order-17 anchor can reach equality
  only with a maximal order-15 block;
- the same is true for an order-29 anchor.

The last two equalities are CRT products:

```text
c(17,15)=14/85,             c(29,15)=24/145.           (16)
```

All five companions would have to attain equality.  Thus 17 and 29 reduce to
the order-15 branch, which in turn forces the already handled order-three
branch.  Every remaining noncore block has density at most `7/43`.

## 5. Weighted reflection cores at 11 and 23

The surviving dense orders are `11,22,23`.  They require weights that retain
the reflection-fixed sheet.

### 5.1 The 11/22 core

On `Z/11Z`, put

```text
w_11(5)=1,             w_11(x)=1/2 otherwise,
sum w_11=6.                                               (17)
```

For an order-`m` quotient block `B`, pull it to `lcm(11,m)` and normalize its
weight to a score `sigma_11(m,B)`.  If `11` does not divide `m`, CRT gives

```text
sigma_11=6|B|/m.                                        (18)
```

If `11|m` and `N_11(B)` counts block points congruent to the fixed sheet, then

```text
sigma_11=11(|B|+N_11(B))/(2m),
N_11(B)<=ceil(|B|/11).                                  (19)
```

The exact small scores are

```text
m:             11, 22, 33, 44, 55, 66
max sigma_11:   1,  1,  1, 3/4, 4/5, 5/6.             (20)
```

For every multiple `m>=77`, `(13)` and `(19)` give

```text
11(h+ceil(h/11)) <= (12m+142)/7 < 2m,                  (21)
```

so the score is strictly below one.  After the earlier exceptional branches
and the 23-core are removed, equality occurs only at orders `11,22,33`.
A six-cover must spend exactly one unit in every slot, hence all selected
orders lie in that set and their lcm divides 66.

The two nonpositive residual moduli are exact finite gates:

```text
Q=33: 26 candidate types, 230230 six-subsets,
      one literal cover, joint period 11;

Q=66: 36 candidate types, 1947792 six-subsets,
      thirteen literal covers, joint period 11 or 22.   (22)
```

The word “candidate” matters: only 16 and 26 of those respective supersets
have score exactly one.  Enumerating the larger supersets is harmless and
makes `(22)` a stronger hostile check.  Neither modulus has a cover of joint
period `Q`; only 11 and 22 survive.

### 5.2 The 23 core

On `Z/23Z`, put

```text
w_23(11)=1/2,          w_23(x)=1/4 otherwise,
sum w_23=6.                                               (23)
```

The corresponding score is `6|B|/m` when `23` does not divide `m`, and

```text
sigma_23=23(|B|+N_23(B))/(4m),
N_23(B)<=ceil(|B|/23)                                   (24)
```

when `23|m`.  Its exact small table is

```text
m:             23, 46, 69
max sigma_23:   1, 3/4, 5/6.                            (25)
```

For every remaining multiple from 92 onward,

```text
23(h+ceil(h/23)) <= (24m+298)/7 < 4m.                  (26)
```

Thus equality occurs only at order 23, forcing joint period 23 in the pure
core.

## 6. Mixed cores leave a positive CRT cylinder

Suppose `a>=1` blocks from the 11/22 core and `b>=1` order-23 blocks occur.
The first missed set is typed modulo 11 when `Q` is odd and modulo 22 when
`Q` is even; this parity-sensitive choice is required for an equal-fibre CRT
projection.  In either case the two coprime core coordinates leave a missed
cylinder of density at least

```text
((11-2a)/11) ((23-4b)/23).                              (27)
```

The other `6-a-b` blocks cover at most `7/43` each.  Subtracting that total
from `(27)` gives

```text
(344ab-207a-121b+253)/10879.                            (28)
```

For positive `a,b` with `a+b<=6`, `(28)` is positive; its minimum is

```text
269/10879                                                (29)
```

at `a=b=1`.  Hence mixed cores cannot cover.  Sections 4--6 prove the
negative direction of `(3)`.

## 7. Fixed-zero lower-base-free corollary

THM-3414 forces a fixed-zero cap-six cover to have one of the bases
`15,16,18,20,24`.  Under the hypothesis in `(3)`, only 15 can divide `Q`.
Moreover `Q` is odd: an even multiple of 15 would also be a multiple of 10.

For odd `Q`, send every fixed-zero residue `u` to the half-twist residue
`2u`.  The sheet bijection

```text
x=2ell+1 mod Q                                           (30)
```

identifies the masks, and multiplication by two preserves every quotient
order.  A primitive fixed-zero family therefore yields a half-twist cover of
joint period `Q`.  Formula `(3)` and `15|Q` force `Q=15`.  The positive
fixed-zero witness is

```text
(1,2,3,4,5,7).                                          (31)
```

This proves `(6)` without asserting an all-`Q` fixed-zero primitive formula.

## 8. Exact companion and validity boundary

Run

```bash
python3 04-computation/lrc_half_twist_rank6_joint_period_probe_20260815.py
python3 -O 04-computation/lrc_half_twist_rank6_joint_period_probe_20260815.py
```

The companion freezes every competitive complement value for anchors
15, 17, and 29; the full 11/23 score classification and analytic tails; all
2,178,022 candidate six-subsets in `(22)`; all fifteen mixed-core gaps; the
five positive atoms; 160 scaled lower-base breaker controls; the primitive
`Q=27` correction; the `Q=33,46` period-loss hostiles; and the odd-`Q`
all-even parity hostile `(10)`.  Normal and optimized transcripts are
identical and match the stored output.

The independent audit reconstructed the proof and computations, caught the
family-level distinction in `(9)`, corrected “equality types” to “candidate
types,” and checked the parity-sensitive modulo-11/modulo-22 routing in the
mixed core.  No claim in this theorem lowers the LRC(14) ledger.
