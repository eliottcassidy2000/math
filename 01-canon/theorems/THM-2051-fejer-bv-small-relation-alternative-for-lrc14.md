---
id: THM-2051
title: Fejer--BV whole-product approximation gives a small-relation-or-positive-BONF5 alternative for LRC(14)
status: PROVED. If a thirteen-speed row has no exact relation of support two through five and coefficient height at most 2^20, its continuous quintic Bonferroni functional is positive, hence it has a positive-measure strict lonely set. The proof keeps each whole signed centered product, annihilates its finite Fejer part, and controls the remainder in L1. It bypasses the THM-946 absolute strip/slab estimates for this dichotomy, but does not classify the small-relation branch or prove LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on:
  - THM-604
  - THM-699
  - THM-935
  - THM-1092
related:
  - THM-940
  - THM-946
  - THM-1645
  - THM-2050
  - HYP-8841
---

# THM-2051 -- Fejer--BV small-relation alternative

Let `T=R/Z` carry Haar probability measure, and put

```text
D={t: ||t||<1/14},   h=1_D,   a=1/7,   u=h-a.
```

Let `v_1,...,v_13` be distinct positive integer speeds. Define

```text
S_k = sum_(|A|=k) integral_T product_(i in A) h(v_i t) dt,
BONF5(v) = sum_(k=0)^5 (-1)^k S_k.
```

## 1. The theorem

**Theorem.** Suppose there is no integer relation

```text
sum_{i in A} k_i v_i=0,
2<=|A|<=5,             0<|k_i|<=H=2^20.              (1)
```

Then

```text
BONF5(v)>0.                                               (2)
```

Consequently the weak safe set has positive measure, and after discarding its
finitely many wall points there is a positive-measure set of phases satisfying
`||v_i t||>1/14` for every `i`.

Equivalently, every thirteen-speed row obeys the rigorous alternative

```text
positive-measure strict loneliness
    OR
an exact support-2-to-5 relation of coefficient height at most 2^20.   (3)
```

## 2. Exact centered expansion

For a nonempty subset `A`, set

```text
delta_A = integral_T product_(i in A) u(v_i t) dt.
```

Expanding `h=a+u` and grouping by the exact centered support gives

```text
BONF5
 = B_eq
   +(24/343) sum_(|A|=2) delta_A
   -(24/49)  sum_(|A|=3) delta_A
   -(2/7)    sum_(|A|=4) delta_A
   -         sum_(|A|=5) delta_A,                    (4)

B_eq=sum_(k=0)^5 (-1)^k C(13,k)/7^k=2052/16807.     (5)
```

Indeed, the coefficient of a centered support of size `s` is

```text
(-1)^s E_s,
E_s=sum_(j=0)^(5-s) (-1)^j C(13-s,j)/7^j,
```

and `(E_2,E_3,E_4,E_5)=(24/343,24/49,-2/7,1)`. The
size-one centered moments vanish. This is THM-935's *continuous* centered
identity. THM-940 is a useful discrete analogue, but is not used as an
identification between the two settings.

## 3. A quantitative Fejer--BV lemma

Let `F_H` be the normalized Fejer kernel of order `H`, put `N=H+1`, and set

```text
p_H=F_H*u.
```

Then `p_H` is a mean-zero trigonometric polynomial with Fourier support
`0<|n|<=H`; explicitly its multiplier is
`max(1-|n|/(H+1),0)`. Positivity of `F_H` gives

```text
||u||_infinity, ||p_H||_infinity <= 6/7.              (6)
```

Moreover

```text
||u-p_H||_1 <= epsilon_H
 := (1+2 log(H+1))/(2(H+1)).                           (7)
```

*Proof.* Translation of the interval `D` gives

```text
||u-u(.-y)||_1=||h-h(.-y)||_1<=2||y||.
```

For `y` represented in `[-1/2,1/2]`, the elementary bounds
`|sin(pi N y)|<=1` and `|sin(pi y)|>=2|y|` give

```text
F_H(y)<=min(N,1/(4N|y|^2)).                            (8)
```

Therefore

```text
||u-F_H*u||_1
 <=2 integral_T ||y||F_H(y)dy
 <=2[2 integral_0^(1/(2N)) yN dy
       +2 integral_(1/(2N))^(1/2) y/(4Ny^2) dy]
 =(1+2 log N)/(2N).
```

This proves (7). The estimate is the whole-object BV contraction principle of
THM-699 applied to a Fejer approximate identity; it never prices individual
relation-lattice atoms. QED.

## 4. Finite annihilation and the signed tail bound

Fix `A` with `2<=s=|A|<=5`. Expanding the finite Fourier polynomial gives

```text
integral_T product_(i in A) p_H(v_i t)dt
 = sum_(0<|n_i|<=H, sum_i n_i v_i=0)
       product_i p_H_hat(n_i).                         (9)
```

Hypothesis (1) makes the sum empty. Hence its value is zero. Multiplication by
each nonzero integer `v_i` preserves Haar measure, so product telescoping and
(6)--(7) give

```text
|delta_A|
 = |integral product_i u(v_i t)-integral product_i p_H(v_i t)|
 <= s(6/7)^(s-1) epsilon_H.                            (10)
```

Taking absolute values only *after* each signed product has been summed, (4)
and (10) yield

```text
|BONF5-B_eq| <= K epsilon_H,                            (11)

K=sum_(s=2)^5 C(13,s)|E_s|s(6/7)^(s-1)
 =1477008/343.                                         (12)
```

This is the step unavailable to the termwise absolute `T_4/T_5` program: all
high-frequency relations are already packaged inside one L1 error.

## 5. Exact numerical closure at `H=2^20`

Now `N=2^20+1<2^21` and `log 2<1`, so

```text
epsilon_H < 43/(2N).
```

Consequently

```text
K epsilon_H
 < 31755672/359661911
 < 2052/16807=B_eq,                                    (13)

B_eq-31755672/359661911
 =595652076/17623433639>0.                             (14)
```

Equations (11)--(14) prove (2). Finally, the odd Bonferroni inequality gives

```text
measure{t: no bad interval is active} >= BONF5(v)>0,
```

which proves the theorem. QED.

## 6. What this changes, and what it does not

THM-946's punctured-line, strip, and slab estimates remain open as estimates
for the *absolute relation masses*. They are no longer required to prove the
coarse universal alternative (3). The physical-space telescope preserves the
cancellation that those estimates discarded.

The result gives HYP-8841 a finite circuit sidecar: a row either terminates at
a strict positive-measure exit, or it lands on one of finitely many relation
templates `(A,(k_i))` with `|A|<=5` and `|k_i|<=2^20`. It does not show that a
move inside such a hyperplane decreases a Noetherian height, and it does not
classify the structured branch. The tight AP is deliberately not in the
dissociated branch: it has many tiny support-three relations.

The theorem also sharpens the CT bridge of THM-1645/HYP-8840. Constant-term
information can prove a *strict* exit when its low-frequency relation packet
is empty; THM-2047/2050 remain decisive warnings that volume and local germs
cannot recognize the measure-zero AP boundary. Thus this result attacks the
strict branch, not Wall A itself.
