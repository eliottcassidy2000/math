---
id: THM-4435
title: "Four-node metric blindness and universal Hermite precision"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 third research wave
---

# THM-4435 -- Four-node metric blindness and universal Hermite precision

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof, complete minor certificate and independent inverse audits](../../05-knowledge/results/four-node-metric-and-hermite-precision-overnight-hexagon-sep05.md)
are part of this theorem.

## Sharp two-jet precision at arbitrary integer nodes

Let n>=1, let x_1,...,x_n be distinct integers, and let E send an integer
polynomial of degree<2n to its values and first Hasse derivatives at those
nodes. Put F(X)=product_i(X-x_i) and d_i=F'(x_i). Then the largest integer
Smith factor of E is exactly

```text
s_max = lcm_i [ |d_i|^3 / gcd(|d_i|,F''(x_i)) ].
```

F'' is the **ordinary** second derivative. For every prime p, with
v_p(0)=infinity, the sharp uniform precision loss is

```text
L_p=max_i{2v_p(d_i), 3v_p(d_i)-v_p(F''(x_i))}.
```

For every N>=1, data modulo p^(N+L_p) determine every coefficient modulo
p^N, and loss L_p-1 fails uniformly. This is the fixed degree box and full
observer, not a moving-module statement.

The proof constructs every inverse column from the classical Hermite
cardinal polynomials. The leading derivative/value columns attain the
two possible denominator types1/d_i^2 and-F''(x_i)/d_i^3; all other
coefficients are integral combinations of them. The least inverse
denominator is precisely the largest Smith factor. The basis is classical;
no external priority claim is made.

If x_i has a unique closest p-adic neighbour of depth f_i and
S_i=sum_(j!=i)v_p(x_i-x_j), its local contribution is
2S_i+max(0,f_i-v_p(2)). Tied closest neighbours retain the F'' coordinate;
a general metric-only largest-factor law is not asserted.

## The full four-node partition is not metric-only

For every integer e>=0 consider

```text
A_e=2^e*(0,1,2,5),       B_e=2^e*(0,1,3,4).
```

The correspondence from(0,1,2,5) to(1,0,3,4) is a dyadic distance isometry.
Every corresponding three-node Smith restriction is identical by
[THM-4429, arbitrary three-node law](THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md).
Yet at e=3 the full exponent lists are respectively

```text
(0,0,4,7,12,16,19,26),    (0,0,4,7,12,17,18,26).
```

Their determinant valuations and largest exponents coincide. For all e>=0,
the six residual determinantal valuations are

```text
delta1=min(2e,e+1), delta2=min(4e,3e+2),
delta3=min(7e+2,6e+6),
delta4(A)=min(12e+4,11e+6), delta4(B)=min(12e+4,11e+7),
delta5=17e+7, delta6=24e+12.
```

Thus the partitions coincide exactly for e<=2 and differ for every e>=3;
the common largest exponent is7e+5. The missing coordinate is a unit-sensitive
quadratic coefficient of an intermediate minor, not another pair distance.

All923 nonempty residual minors of each base were evaluated by independent
determinant definitions. Their homogeneous scale costs give a finite proof
for every depth, not a depth extrapolation. Full integer Smith/inverse
audits and explicit n=1/vanishing-F'' controls passed. The dyadic
counterexample does not contradict the incoming **ternary** two-plus-two
metric theorem or THM-4429. General higher-jet and full higher-node
partition descriptions remain OPEN.
