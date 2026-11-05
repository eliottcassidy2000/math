---
id: THM-2048
title: The LRC peel discrepancy has an exact fiber-occupancy decomposition and a quantized variance tax
status: PROVED. This strengthens the necessary THM-731/732 obstruction for a zero-measure cover; it is not a proof of LRC(14).
source: codex-2026-07-21-LRC-unrelated-transfer
depends_on:
  - THM-731
  - THM-732
related:
  - HYP-8815
  - THM-2047
  - THM-2000
---

# THM-2048 -- fiber quantization inside the peel discrepancy

## 0. Abstract finite-orbit lemma

The mechanism is not specific to the number fourteen. Let a probability space
carry a measure-preserving cyclic action `T` of order `v`. Let `D` be an
invariant measurable set of measure `rho>0`, let `E` be contained in `D` modulo
null sets, put `a=|E|`, and define the orbit projection

```text
P f = (1/v) sum_{j=0}^{v-1} f composed with T^j,  f=1_E.
```

Then, with `theta=fractional_part(va/rho)`,

```text
||P f-a||_2^2
 = ((1-rho)/rho)a^2 + ||P f-(a/rho)1_D||_2^2
 >= ((1-rho)/rho)a^2 + (rho/v^2)theta(1-theta).             (A)
```

Indeed `P f` is supported on `D`, has mean `a`, and is orthogonal there to the
constant `(a/rho)1_D`; this gives the equality. Under normalized measure on
`D`, the integer-valued random variable `vP f` has mean `va/rho`, so the
two-neighbor integer variance bound gives the inequality. Equality in the
first bound requires uniform integral orbit occupancy.

This abstract form is available in any support, tiling, coding, or cyclic
packet problem only after an invariant carrier `D` and a genuine containment
`E subset D` have been proved. It is a reusable projection lemma, not a claim
that those other problems automatically have the LRC containment.

Fix the threshold `1/14`, a finite speed set, and a peeled speed `v`. Let

```text
E = G'_{~v},              mu = |E|,
D_v = {t: ||v t||<1/14},  |D_v|=1/7,
f = 1_E,
P_v f(t) = (1/v) sum_{j=0}^{v-1} f(t+j/v).
```

THM-731's discrepancy is

```text
disc_v = sum_{m!=0} |f_hat(mv)|^2
       = ||P_v f-mu||_2^2.                                 (1)
```

Suppose the full weak-safe set has measure zero. In particular this necessary
condition holds for every hypothetical LRC(14) counterexample. The peeling
identity gives `|E\D_v|=0`, so `E` is contained in `D_v` modulo null sets.

## 1. Exact orthogonal decomposition

**Theorem.** Under the zero-measure hypothesis,

```text
disc_v = 6 mu^2 + ||P_v f-7mu 1_{D_v}||_2^2.                (2)
```

*Proof.* The set `D_v` is invariant under translation by `1/v`. Since `f`
vanishes almost everywhere off `D_v`, so does `P_v f`. Translation invariance
also gives `integral P_v f=integral f=mu`. Hence

```text
integral_{D_v}(P_v f-7mu)=mu-7mu|D_v|=0.
```

On `D_v`, split `P_v f-mu=(P_v f-7mu)+6mu`; off `D_v` it equals `-mu`.
The displayed mean-zero identity kills the cross term, and

```text
(6/7)mu^2 + (1/7)(6mu)^2 = 6mu^2.
```

This proves (2). ∎

Thus THM-731's old necessary inequality `disc_v>=6mu^2` is a Pythagorean
projection statement: equality means uniform conditional occupancy on the
`v`-fibers inside `D_v`.

## 2. Quantized variance tax

The function `P_v f` takes values in
`{0,1/v,...,1}`. Put

```text
theta_v = fractional_part(7v mu).
```

**Theorem.** Under the same hypothesis,

```text
disc_v >= 6mu^2 + theta_v(1-theta_v)/(7v^2).                (3)
```

Equality in the old bound `disc_v=6mu^2` is possible only if
`7vmu` is integral and almost every `v`-orbit based in `D_v` has the same
number `7vmu` of points in `E`.

*Proof.* Give `D_v` normalized probability measure `7 dt` and set
`X=vP_v f`, an integer-valued random variable in `{0,...,v}`. Its mean is

```text
E_D[X] = 7v integral_{D_v} P_v f = 7vmu.
```

For any integer-valued random variable of mean `a+theta`, `a` integral and
`0<=theta<1`, the pointwise inequality `(X-a)(X-a-1)>=0` gives
`Var(X)>=theta(1-theta)`. Therefore

```text
||P_v f-7mu 1_{D_v}||_2^2
  = Var_D(X)/(7v^2)
  >= theta_v(1-theta_v)/(7v^2).
```

Combine with (2). Equality in the old bound makes this variance zero, so `X`
is almost surely constant; its value is the integral `7vmu`. ∎

## 3. Combination with the Bernoulli edge bound

Let `r_v` be the number of positive-length interval components of `E`, with
null isolated points discarded. THM-732 gives

```text
disc_v <= r_v^2/(3v^2).                                    (4)
```

Combining (3) and (4), every zero-measure packet—and hence every hypothetical
counterexample—must satisfy, for every peel,

```text
6(vmu)^2 + theta_v(1-theta_v)/7 <= r_v^2/3.                 (5)
```

Equivalently, if one peel violates (5), then the full safe set has positive
measure and supplies a strict lonely interval. This is a finite exact
arithmetic test because the threshold endpoints, `mu`, and `r_v` are rational
and computable from the signed wall arrangement.

Indeed, the equality phases `||wt||=1/14` form a finite wall set. Removing
them from a positive-measure weak-safe set leaves a positive-measure strict
safe set, which is open and hence contains an interval.

The new term is small but rigid: it vanishes exactly on integral conditional
fiber occupancy and otherwise charges the distance to the nearest such
occupancy. It turns the qualitative discrepancy direction corrected in
MISTAKE-221 into a quantized obstruction without asserting that the obstruction
alone closes Wall A.

## 4. An exact strict-gain Cover14 certificate

The tax can decide a row on which the old uniform `r^2/(3v^2)` tail bound is
inconclusive. Take

```text
S = {1,8,11,12,14,17,22,26,35,40,54,90,93},   v=93,
B = S\{93}.
```

This row is primitive, and every modulus `q=2,...,14` divides at least one
member of `S`, so it lies in the genuine Cover14 residual rather than the easy
missed-modulus branch. Exact interval arithmetic on the good set of `B` gives

```text
mu = 35517/280280,       r_v = 50,
theta_v = {7*93*mu} = 19801/40040.
```

The old necessary tail inequality is still allowed, with positive slack

```text
r_v^2/3 - 6(93mu)^2 = 3001166951/117835317600 > 0.
```

But the new integer-fiber tax is

```text
theta_v(1-theta_v)/7 = 400752439/11222411200,
```

and it exceeds that slack by

```text
2413467317/235670635200 > 0.
```

Thus (5) fails, so the full good set of `S` cannot have measure zero. This is a
strict gain over the *uniform tail corollary* of THM-732: it proves a lonely
interval using only `(mu,r_v,v)` and the integrality of orbit occupancy, without
evaluating the Bernoulli edge-pair sum.

For an independent audit, the full exact interval computation gives measure

```text
|G'| = 589595/5213208 > 0,
```

while the exact Bernoulli computation gives

```text
disc_93 = 2220515193157/2038315323844800,
6mu^2-disc_93 = 27737954067563/291187903406400 > 0.
```

So the stronger per-row THM-731/732 certificate also closes this example. The
new content is the cheaper quantized-tail decision, not a family inaccessible
to exact Bernoulli arithmetic and not a proof of the general Cover14 case.
