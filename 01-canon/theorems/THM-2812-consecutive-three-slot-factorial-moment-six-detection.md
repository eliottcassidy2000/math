---
id: THM-2812
title: "Consecutive three-slot factorial moment and Gaussian moment-six detection"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT HOSTILE AUDIT.  For L(s^k)=k!, every nonzero
  H=s^n(x+ys+zs^2), n>=0, has a nonzero power moment among j=1,2,3.
  For n>=1 its polynomial two-charge lift has a nonzero Gaussian moment
  of order at most six, and an exact full-support witness first dies at
  order six.  This proves the sharp HYP-8765 cutoff on consecutive
  three-slot two-charge envelopes only; arbitrary three-slot supports,
  arbitrary Wick channels, and the full effective conjecture remain open.
source: root/consecutive-three-slot-factorial-detection-2026-07-28
depends_on: []
related:
  - THM-1790-the-emp-floor-detection-depth-at-least-degree-plus-one
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go
  - HYP-8765-gmc2-radial-channel-return-tower
script: 04-computation/gmc_consecutive_three_slot_factorial_detection_thm2812.py
output: 05-knowledge/results/gmc_consecutive_three_slot_factorial_detection_thm2812.out
script_sha256: f5452eb57301afe157d2fcbddb84d219d4d516aa1e62fccd418da731ac121f7f
output_sha256: efe19b781bea32cc75a472316421c3f9b7d3dd9cbe74b873410629acfdee36df
hash_basis: LF-normalized bytes
---

# THM-2812 -- consecutive three-slot radial support closes at moment three

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT HOSTILE AUDIT.**

THM-2173 proves the sharp lower side: every prescribed three-dimensional
radial envelope contains a full-support polynomial whose first two factorial
moments vanish.  It leaves the matching third-moment upper bound as a
three-term Strong Factorial problem.

This theorem closes that upper bound on every translated consecutive
three-slot envelope.  The proof is a division-free elimination after the
linear moment is removed.  Its decisive feature is not support counting:
the second-moment conic has two nonreal points, while the third moment reduces
to a nonzero real linear function on that conic.

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^k)=k!.           (1)
```

For every integer `n>=0` and every `(x,y,z) in C^3`, put

```text
H=s^n(x+y s+z s^2).                                    (2)
```

Then

```text
L(H)=L(H^2)=L(H^3)=0             implies             H=0. (3)
```

Equivalently, every nonzero `(2)` is detected by one of its first three
factorial power moments.

## 2. Coordinate-boundary supports

A one-term polynomial has a nonzero first moment.  Suppose instead that

```text
H=A s^u+B s^v,                    0<=u<v,   AB!=0.       (4)
```

If `L(H)=0`, then

```text
B/A=-u!/v! in R_<0.                                      (5)
```

Thus `H=A R` for a nonzero real polynomial `R`.  The integral realization

```text
L(f)=integral_0^infinity f(s)e^(-s) ds                  (6)
```

gives

```text
L(H^2)
 =A^2 integral_0^infinity R(s)^2e^(-s) ds !=0.          (7)
```

Hence every one- or two-coordinate boundary of `(2)` already satisfies the
theorem.  It remains to treat exact three-term support, where `xyz!=0`.

## 3. Eliminate the first moment

For `j=1,2,3`, normalize

```text
F_j=L(H^j)/(jn)!.                                       (8)
```

Because scaling `(x,y,z)` does not affect common projective nullity, set

```text
x=1,                         t=y/x.                     (9)
```

The first normalized moment is

```text
F_1=x+(n+1)y+(n+1)(n+2)z.                              (10)
```

Thus `F_1=0` gives

```text
z/x=-(1+(n+1)t)/((n+1)(n+2)).                          (11)
```

After substituting `(11)`, the second moment is

```text
F_2=Q_n(t)/(n+1),                                      (12)

Q_n(t)
 =2(n+1)^2(2n+1)t^2
  +6(n+1)(2n+1)t
  +(9n+5).                                             (13)
```

Its discriminant is

```text
disc(Q_n)=-4(n+1)^2(2n+1)<0.                           (14)
```

Therefore both roots of `Q_n` are nonreal.

## 4. The cubic has a real linear remainder

Substitute `(11)` into `F_3` and divide the resulting cubic in `t` by
`Q_n`.  Exact Euclidean division gives the remainder

```text
-(A_n t+B_n)/((n+1)^2(2n+1)),                          (15)

A_n=36n^4+57n^3+15n^2-9n-3,
B_n=52n^3+49n^2+8n-1.                                 (16)
```

At `n=0`,

```text
(A_0,B_0)=(-3,-1).                                     (17)
```

For `n>=1`, both entries in `(16)` are strictly positive.  Indeed

```text
A_n=(36n^4-9n)+(57n^3+15n^2-3)>0,
B_n>=52+49+8-1>0.                                     (18)
```

Thus `(15)` is a nonzero real linear polynomial for every integer `n>=0`;
its only root is real.  It cannot vanish at either nonreal root of `Q_n`.
Consequently `F_1,F_2,F_3` have no common exact-three-support projective
zero.  Together with Section 2, this proves `(3)`.

## 5. The bound is sharp

At `n=1`, take

```text
t=-3/4+i sqrt(3)/12,
u= 1/12-i sqrt(3)/36,
H_*=s(1+t s+u s^2).                                   (19)
```

All three coefficients are nonzero, and direct exact evaluation gives

```text
(L(H_*),L(H_*^2),L(H_*^3))
 =(0,0,-18-4i sqrt(3)).                                (20)
```

Thus three powers are necessary uniformly, not merely sufficient.

## 6. Two-charge Gaussian lift

Assume `n>=1`, put `s=ZW`, and define the polynomial

```text
h=s^(n-1)(x+y s+z s^2),
P=W+Z h(s).                                             (21)
```

The two summands have charges `-1` and `+1`.  Every odd moment vanishes.
In degree `2j`, charge balance forces exactly `j` copies of each summand, so

```text
E[P^(2j)]=binom(2j,j)L((s h)^j)
          =binom(2j,j)L(H^j).                           (22)
```

Equations `(3)` and `(22)` prove that every nonzero polynomial `(21)` has a
nonzero Gaussian moment of order at most six.

For the witness `(19)`, moments one through five vanish and

```text
E[P^6]=20L(H_*^3)=-360-80i sqrt(3) !=0.                (23)
```

Hence six is the exact uniform bound on this family.  In HYP-8765 notation,
`P` has `k=4` monomials and primitive return `R=2`, so `(23)` proves the
conjectured sharp cutoff

```text
(k-1)R=6                                               (24)
```

on every consecutive three-slot two-charge envelope.

## 7. Information ledger and boundary

| item | exact content |
|---|---|
| source | `s^n span{1,s,s^2}`, arbitrary complex coefficients |
| map | eliminate `L(H)`, then reduce `L(H^3)` modulo the second-moment conic |
| preserved | coefficient phase and exact factorial weights |
| decisive invariant | nonreal conic roots versus nonzero real linear remainder |
| sharp hostile | `(19)--(20)`, whose first two moments vanish |
| Gaussian consequence | two-charge detection by moment six, sharply |
| first unproved extension | nonconsecutive three-term support |

This theorem does not prove the Strong Factorial Conjecture for arbitrary
three-term supports.  It does not separate arbitrary Wick channels inside a
scalar moment, prove general HYP-8765, or replace THM-2022's already complete
proof of NC2/GMC2.

## 8. Exact companion

Run

```bash
python 04-computation/gmc_consecutive_three_slot_factorial_detection_thm2812.py
python -O 04-computation/gmc_consecutive_three_slot_factorial_detection_thm2812.py
```

Both executions byte-match

```text
05-knowledge/results/gmc_consecutive_three_slot_factorial_detection_thm2812.out.
```

The companion derives `(10)--(16)` symbolically, checks the integer positivity
boundary through `n=200` as a hostile control in addition to the proof
`(17)--(18)`, and verifies the sharp factorial and Gaussian witness exactly
in `Q(sqrt(-3))`.  It has explicit exception gates, no truth-bearing Python
assertions, no floating point, and no scratch dependency.

**Awaiting independent hostile audit; not QED.**
