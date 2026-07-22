---
id: THM-2051
title: Fejer--BV whole-product approximation gives a genuine-higher-relation-or-positive-BONF5 alternative for LRC(14)
status: PROVED. The basic form says that no exact relation of support two through five and coefficient height at most 2^20 forces positive continuous BONF5. The stronger pair-exact form uses THM-965's sharp covariance floor and requires dissociation only for supports three through five, at height 2^21. Hence every thirteen-speed row either has a positive-measure strict lonely set or a genuine 3-, 4-, or 5-speed relation of bounded height. The proof bypasses the THM-946 absolute strip/slab estimates, but does not classify the higher-relation branch or prove LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on:
  - THM-604
  - THM-699
  - THM-935
  - THM-1092
related:
  - THM-940
  - THM-946
  - THM-965
  - THM-1645
  - THM-2050
  - HYP-8841
script: 04-computation/lrc14_fejer_bv_higher_circuit_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_fejer_bv_higher_circuit_referee_codex_20260721.out
script_sha256: 9da37dfd4bc3e87538202080cd5d46f15c738c1086a56ee1304c9309bbd7e70b
output_sha256: 55d38afb379280ab51f61ef3675d5e5fafa2b34af229998d735466787cad3b62
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

**Theorem A.** Suppose there is no integer relation

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

## 6. Pair-exact sharpening: only genuine higher relations remain

The pair term in (4) can be paid directly instead of approximated. This
removes the support-two hypothesis, at the cost of moving the higher-support
Fejer horizon from `2^20` to `2^21`.

**Theorem B (stronger alternative).** Suppose there is no relation

```text
sum_(i in A) k_i v_i=0,
3<=|A|<=5,             0<|k_i|<=H_1=2^21.             (15)
```

Then `BONF5(v)>0`. Consequently every thirteen-speed row satisfies

```text
positive-measure strict loneliness
    OR
an exact support-3-to-5 relation of coefficient height at most 2^21.   (16)
```

*Proof.* For distinct speeds `a,b`, divide by `g=gcd(a,b)` and write the
coprime reduced pair as `A<B`. THM-965 gives its centered pair moment exactly:

```text
delta_(A,B)
 =[fold_14(A+B)-fold_14(B-A)]/(196AB),
fold_14(r)=(r mod 14)(14-(r mod 14)).                  (17)
```

Since `0<=fold_14<=49`, if `AB>=27` then

```text
delta_(A,B)>=-1/(4AB)>=-1/108>-6/637.                 (18)
```

The remaining coprime pairs have `AB<=26`, hence `A<=4`. Direct substitution
in (17) gives the exact minima at each possible `A`:

| `A` | minimizing `B` | minimum `delta_(A,B)` |
|---:|---:|---:|
| 1 | 13 | `-6/637` |
| 2 | 11 | `-4/539` |
| 3 | 8 | `-1/392` |
| 4 | 5 | `2/245` |

Thus every pair moment is at least `-6/637`, with equality exactly at reduced
ratio `1:13`. Since the pair coefficient in (4) is positive,

```text
(24/343) sum_(|A|=2) delta_A
 >=(24/343) C(13,2)(-6/637)
 =-864/16807.                                          (19)
```

After paying (19), the equilibrium reserve is

```text
2052/16807-864/16807=1188/16807.                       (20)
```

Apply the proof of Sections 3--4 only to supports `s=3,4,5`. Their telescope
constant is

```text
K_high=sum_(s=3)^5 C(13,s)|E_s|s(6/7)^(s-1)
      =10316592/2401.                                  (21)
```

At `H_1=2^21`, put `N_1=H_1+1<2^22`. Since `log 2<1`, (7) gives

```text
epsilon_(H_1)<45/(2N_1),
K_high epsilon_(H_1)<25791480/559473817.               (22)
```

The exact remaining margin is

```text
1188/16807-25791480/559473817
 =96283836/3916316719>0.                               (23)
```

This proves Theorem B. The finite pair core, all coefficients, and both
rational margins are independently checked by the stored optimization-safe
referee. QED.

## 7. Geometry of the surviving relation branch

The relations in (16) are height-filtered signed hyperedges, not ordinary
matroid circuits: thirteen scalar speeds have rational rank one, so their
ordinary linear-matroid circuits are already pairs. Two elementary facts give
the useful structure.

### 7.1 Balanced versus anchored relations

For a surviving relation on `x_1<...<x_s`, put `sigma=sum_i k_i`. Then

```text
sigma x_1+sum_(i=2)^s k_i(x_i-x_1)=0.                 (24)
```

If `sigma=0`, the relation depends only on gaps and survives every common
translation. THM-1225 warns that this translation-blind branch cannot by
itself characterize LRC extremality; it needs the divisor/owner-height
sidecar. If `sigma!=0`, (24) anchors absolute position and gives

```text
x_1 <= H_1(s-1)(x_s-x_1)/|sigma|.                     (25)
```

Thus the unbalanced branch carries an explicit magnitude-to-span constraint,
while the balanced branch is the AP/Freiman-facing one.

### 7.2 Twelve independent relations force a finite speed box

Suppose the same primitive speed vector `v in Z_{>0}^13` satisfies twelve
linearly independent relation rows. If each row has support at most five and
coefficient height at most `H`, place them in a `12 x 13` matrix `K`. The
signed `12 x 12` maximal minors form a nonzero integer generator of `ker K`.
Primitivity divides out their common gcd, and Hadamard's inequality gives

```text
max_i v_i <= (sqrt(5)H)^12=5^6 H^12.                  (26)
```

Therefore twelve independent harvested relations reduce LRC to a finite
(albeit enormous) exact enumeration. The missing termination supplier is now
precise: an unresolved move must produce a *new independent* relation. Merely
rediscovering the same hyperedge gives no descent, and coefficient elimination
may leave the height box.

## 8. What this changes, and what it does not

THM-946's punctured-line, strip, and slab estimates remain open as estimates
for the *absolute relation masses*. They are no longer required to prove the
coarse universal alternative (3). The physical-space telescope preserves the
cancellation that those estimates discarded.

Theorem B gives HYP-8841 a finite higher-relation sidecar: a row either
terminates at a strict positive-measure exit, or it lands on one of finitely
many relation templates `(A,(k_i))` with `3<=|A|<=5` and
`|k_i|<=2^21`. It does not show that a move inside such a hyperplane decreases
a Noetherian height, and it does not classify the structured branch. The tight
AP is deliberately not in the dissociated branch: it has many tiny
support-three relations.

The theorem also sharpens the CT bridge of THM-1645/HYP-8840. Constant-term
information can prove a *strict* exit when its low-frequency relation packet
is empty; THM-2047/2050 remain decisive warnings that volume and local germs
cannot recognize the measure-zero AP boundary. Thus this result attacks the
strict branch, not Wall A itself.
