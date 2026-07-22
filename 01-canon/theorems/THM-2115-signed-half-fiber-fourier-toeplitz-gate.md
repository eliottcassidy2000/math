---
id: THM-2115
title: "Signed Fourier-Toeplitz gate on the guard half-fiber"
status: >
  PROVED. On p.X=1/2, terminal lifts c_i=a_i p+n_i q become
  binary-shifted combs with shifts e_i=a_i mod 2. If n such combs cover the
  fiber at radius 1/14, their excess multiplicity has mean n/7-1 and an
  explicit positive-semidefinite Toeplitz moment sequence. At rank eight
  every divisor-frequency coefficient has absolute value at most 1/7.
  Violation gives a strict mixed interval. An exact eight-character row which
  passes every THM-2105 clock m=2,...,7 and saturates its half-fiber Hunter
  tree is nevertheless excluded at frequency 84, where the 21/42 harmonic
  packet has coefficient greater than 1/7. This supplies a new higher-order
  route after clocks and pair trees, not a universal rank-eight closure.
source: codex-2026-07-22-LRC-signed-half-fiber-toeplitz
depends_on:
  - THM-641
related:
  - THM-645
  - THM-1185
  - THM-2105
  - THM-2114
  - HYP-2974
  - HYP-3954
---

# THM-2115 -- signed Fourier-Toeplitz gate on the half-fiber

Let `p,q:T^2->T` be independent integer characters, take guard `g=p`, and
write the transverse terminal characters as

```text
c_i=a_i p+n_i q,              a_i in Z, n_i in Z\{0}.  (1)
```

On the guard half-fiber `p.X=1/2`, put

```text
e_i=a_i mod 2 in {0,1},       v_i=|n_i|,
B_i(beta)=1_{||v_i beta+e_i/2||<1/14},
N(beta)=sum_i B_i(beta),      H(beta)=N(beta)-1.        (2)
```

The replacement of `n_i` by its absolute value is harmless because the two
possible centers `0,1/2` are fixed by negation.

If the strict mixed-safe cell is empty, then `N(beta)>=1` away from the finite
set of terminal equality endpoints. Therefore

```text
H(beta)>=0 almost everywhere.                           (3)
```

This endpoint qualification is exact: emptiness of the strict-safe cell gives
a cover by the **closed** danger bands, while (2) uses open bands. They differ
only on finitely many points and hence have the same Fourier integrals.

## 1. Exact signed divisor coefficients

Use the Fourier convention

```text
hat f(r)=integral_T f(beta) exp(-2 pi i r beta) d beta. (4)
```

For the basic arc `b(x)=1_{||x||<1/14}`,

```text
hat b(0)=1/7,
hat b(k)=sin(pi k/7)/(pi k),             k!=0.          (5)
```

Since

```text
B_i(beta)=b(v_i beta+e_i/2),                            (6)
```

pullback of the Fourier series gives, for `r!=0`,

```text
hat H(r)=sum_(i:v_i|r)
  (-1)^(e_i r/v_i)
  sin(pi r/(7v_i))/(pi r/v_i).                         (7)
```

If there are `n` terminals, the constant coefficient is

```text
hat H(0)=mu_n=n/7-1.                                   (8)
```

Equations (7)--(8) are exact. They retain three pieces of information erased
by an unlabelled valuation profile: divisibility by the chosen frequency,
the harmonic quotient `r/v_i`, and the lift parity `e_i`.

## 2. Toeplitz positivity and the scalar gate

For any finite list of integer frequencies `r_0,...,r_s`, define

```text
T_R=(hat H(r_a-r_b))_(0<=a,b<=s).                      (9)
```

Then every cover satisfies

```text
T_R is positive semidefinite.                          (10)
```

Indeed, for arbitrary complex `z_0,...,z_s`,

```text
sum_(a,b) conjugate(z_a) z_b hat H(r_a-r_b)
 =integral_T H(beta)|sum_a z_a exp(2 pi i r_a beta)|^2 d beta
 >=0                                                   (11)
```

by (3). The two-frequency principal minor gives the reusable scalar test

```text
|hat H(r)|<=mu_n                  for every r!=0.       (12)
```

At the critical half-fiber rank `n=8`, this is

```text
|sum_(i:v_i|r) (-1)^(e_i r/v_i)
 sin(pi r/(7v_i))/(pi r/v_i)| <= 1/7.                  (13)
```

Conversely, if one finite Toeplitz section has a negative quadratic form,
then (3) is false on positive measure. Because `H` is integer-valued and
piecewise constant off finitely many endpoints, `N=0` on an open interval.
Every terminal distance is strictly greater than `1/14` in a smaller
interval. Together with `p.X=1/2`, this is a strict mixed escape.

Thus a negative Toeplitz certificate is stronger than a single boundary
witness: it supplies the open cell needed by THM-2097's geodesic finiteness
step.

## 3. An isolated harmonic packet corollary

Formula (13) is especially cheap when few quotient heights divide `r`. For
example, suppose the only terminal heights dividing `4d` are

```text
(v,e)=(d,0),                    (2d,1).                 (14)
```

At frequency `r=4d`, their contribution is

```text
sin(4pi/7)/(4pi)+sin(2pi/7)/(2pi)
 =sin(3pi/7)/(4pi)+sin(2pi/7)/(2pi)>1/7.               (15)
```

The inequality is elementary. On `[0,pi/2]`, concavity gives
`sin x>=2x/pi`, so the left side of (15) is at least

```text
(6/7)/(4pi)+(4/7)/(2pi)=1/(2pi)>1/7.                  (16)
```

The last inequality follows from `pi<7/2`. Therefore the divisor-isolated
dyadic packet (14) is incompatible with a half-fiber cover. Additional
heights dividing `4d` can cancel this coefficient, so the isolation
hypothesis must not be dropped.

## 4. A clock-admissible, tree-equality row killed at frequency 84

Take the guard and eight terminals

```text
p=(1,0),
c_i=(e_i,v_i),
(v_i,e_i)=
 (145,0),(98,1),(76,1),(105,1),
 (71,1),(156,1),(21,0),(42,1).                         (17)
```

The parameter direction `(1,1)` specializes these to the odd guard `1` and
the distinct positive terminal speeds

```text
145,99,77,106,72,157,21,43.                            (18)
```

So (17) is an honest rank-eight coefficient-plane row.

### It passes every small affine clock

For THM-2105's clock `m`, the exact solution sets of

```text
v_i j+m e_i=0 mod 2m                                  (19)
```

contain the following covering subfamilies:

```text
m=2:  98/1->{1,3}, 105/1->{2}, 145/0->{0};
m=3: 105/1->{1,3,5}, 21/0->{0,2,4};
m=4:  76/1->{1,3,5,7}, 98/1->{2,6},
      145/0->{0}, 105/1->{4};
m=5: 145/0->{0,2,4,6,8}, 105/1->{1,3,5,7,9};
m=6:  42/1->{1,3,5,7,9,11}, 21/0->{0,4,8},
      105/1->{2,6,10};
m=7:  21/0->{0,2,4,6,8,10,12},
      105/1->{1,3,5,7,9,11,13}.                       (20)
```

Thus every clock `m=2,...,7` is completely blocked.

### It saturates the half-fiber Hunter tree

Let `w_ij=measure(B_i intersect B_j)`. Applying THM-641's exact anchored
pair law, the lower-triangular matrix of normalized weights `49w_ij`, in the
order (17), is

```text
1
5509/5510  1
1          1          1
2058/2059  1          2695/2698  1
2261/2262  1          245/247    1  5537/5538
1          1          1          0  1          1
1          1          1          0  1          1  0.   (21)
```

Every edge has weight at most `1/49`, and the equality edges contain a
spanning tree. Hence the maximum spanning-tree weight is exactly

```text
tau=7/49=1/7.                                           (22)
```

Hunter's necessary inequality for eight sets of mass `1/7` is precisely
`tau<=8/7-1=1/7`, so the pair-tree test saturates and says nothing.

### The frequency-84 certificate

Among the heights in (17), exactly `21` and `42` divide `84`. Their parities
are the isolated packet (14), with `d=21`. Formula (7) gives

```text
hat H(84)=sin(4pi/7)/(4pi)+sin(2pi/7)/(2pi)>1/7,       (23)
```

contradicting (13). Thus the row has a strict mixed interval despite passing
all six affine clocks and saturating the strongest half-fiber pair-tree
invoice. This is an exact higher-order separation, not a numerical near miss.

## 5. Relation to the live rank-eight carrier

THM-2114's finite-ring needles sample the full two-torus and its
maximum-basis law retains which active label sets are connected in every
maximum tree. THM-2115 instead freezes the guard at its antipode and tests the
binary-shifted one-dimensional fiber by harmonic divisor packets. The routes
are complementary:

```text
finite-ring point exists        -> explicit torsion escape;
half-fiber Hunter tau>1/7       -> pairwise escape;
half-fiber Toeplitz violation   -> open harmonic escape;
all tests pass                  -> retain owners/higher intersections.      (24)
```

The arbitrary shifted Lonely Runner Conjecture is false (as already warned in
HYP-3954). Nothing here asserts a theorem for arbitrary initial positions.
The shifts in (2) are only `0` and `1/2`, forced by the integral lift parity
on the guard half-fiber; that binary restriction is load-bearing.

There is also a structural ceiling. As THM-1185 explains in the unshifted
setting, all finite Toeplitz sections are PSD exactly when the multiplicity
symbol is nonnegative almost everywhere. Toeplitz tests are therefore ideal
for proving an **open** safe interval, but they cannot see a genuinely
boundary-only lonely point. Such equality atoms must return to the finite
clock, endpoint-owner, or maximum-basis lanes.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that, after a clock cover and pair-tree equality,
only spatial triple intersections remain. Equation (7) reveals a different
higher-order carrier: several labels can collide in one divisor frequency
even when their pairwise interval weights look independent.

Candidate tournament vertices were terminals, pair edges, divisor packets,
frequencies, Toeplitz minors, clock residues, and proof obligations. The
faithful vertices here are frequencies `r`, decorated by the signed divisor
packet `{(v_i,e_i,r/v_i):v_i|r}`. Orienting two frequencies by the first one
whose scalar bound is larger, with increasing `r` as the tie Hamiltonian path,
gives a search scheduler; its scores, cycles, SCCs, edge flips, and path count
do not determine PSD of a larger Toeplitz section. The quotient preserves
`cover => H>=0 => Toeplitz PSD` and destroys spatial endpoint ownership, which
must remain a sidecar after harmonic tests pass.

THM-2115 proves the signed harmonic gate and one exact separation row. It does
not prove that every clock-admissible rank-eight row has a violating frequency
or negative Toeplitz minor, classify the PSD equality fibers, or prove
LRC(14). QED.
