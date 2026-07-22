---
id: THM-2085
title: "Explicit height-57 rank-seven relation gate from Selberg box sandwiches"
status: >
  PROVED from the classical Vaaler interval-sandwich lemma and THM-2081.
  If the safe set of seven speeds is contained in one guard, then some
  guard/two-speed triple has a nonzero integer relation of coefficient height
  at most 57. The proof is uniform in all speeds. Degree 57 is the first
  degree at which this particular coordinatewise Selberg/Hunter certificate
  has positive margin; no optimality for the true relation height is claimed.
source: codex-2026-07-22-LRC14-Selberg-relative-Hunter
depends_on:
  - THM-2081
  - THM-2083
related:
  - THM-537
  - THM-685
  - THM-2051
  - THM-2082
  - THM-2086
  - THM-2087
script: 04-computation/lrc14_selberg_rank7_relation_gate_referee_codex_20260722.py
output: 05-knowledge/results/lrc14_selberg_rank7_relation_gate_referee_codex_20260722.out
script_sha256: fba7b10cb67c8929afe3ca512f471b6f3bb3e901f79027cace56514e4396ba69
output_sha256: 2d87d9f0fb4c45f9d629e62b366c932b13a19bc4b4ad98ad03caa7cb1c70b652
hash_basis: working-tree bytes (LF)
external: >
  J. D. Vaaler, Some extremal functions in Fourier analysis,
  Bull. Amer. Math. Soc. (N.S.) 12 (1985), 183--216, Theorem 19.
---

# THM-2085 -- explicit height-57 rank-seven Selberg gate

Use the notation of THM-2081:

```text
D_q={t:||qt||<1/14},
E_h={t:||ht||<1/7},
C_h=(R/Z) minus E_h,
G_Q=intersection_(q in Q) D_q^c.                       (1)
```

For positive integers `u_1,...,u_k`, let

```text
lambda(u_1,...,u_k)
 =min{||m||_infinity:0!=m in Z^k, sum m_r u_r=0}.       (2)
```

For a labelled seven-set `Q={q_1,...,q_7}`, put

```text
Lambda(h,Q)=min_(i<j) lambda(h,q_i,q_j).                (3)
```

## 1. The explicit theorem

For every positive integer `h` and every seven distinct positive integers
`q_1,...,q_7`,

```text
G_Q subset E_h     implies     Lambda(h,Q)<=57.         (4)
```

Thus every rank-seven terminal obstruction in the THM-2073 dyadic tower has
indices `i<j` and integers `a,b,c`, not all zero, such that

```text
a h+b q_i+c q_j=0,
max(|a|,|b|,|c|)<=57.                                  (5)
```

Parity, divisor-completeness, and a height cutoff are not used in (4).

## 2. Cited one-dimensional interval sandwich

We use the periodic interval form of Vaaler's theorem.  If `J` is an interval
in `R/Z` and `H>=1`, there are real trigonometric polynomials

```text
L_(J,H)<=1_J<=U_(J,H)                                  (6)
```

apart from an immaterial endpoint convention, whose Fourier frequencies lie
in `[-H,H]` and whose constant coefficients are

```text
integral L_(J,H)=|J|-epsilon,
integral U_(J,H)=|J|+epsilon,
epsilon=1/(H+1).                                       (7)
```

In particular `U_(J,H)>=0`.  This is the standard Selberg interval pair,
equivalently obtained from Vaaler's approximation/error polynomials.  Only
the order, Fourier support, and constant coefficients in (6)--(7) are used.

Endpoints cause no problem: every quantity below is a Haar measure and all
box boundaries have measure zero on both the ambient torus and the
one-parameter integer orbit.

## 3. A signed tensor minorant for a box

Let `J_1,...,J_k` be circle intervals, write `chi_r=1_(J_r)`, and abbreviate
their Selberg pairs by `L_r,U_r`.  Define

```text
P^-(x_1,...,x_k)
 =product_r U_r(x_r)
  -sum_r (U_r(x_r)-L_r(x_r))*product_(s!=r) U_s(x_s).
                                                               (8)
```

Then

```text
P^-<=product_r chi_r                                  (9)
```

pointwise away from endpoints.  Indeed, the exact product telescope gives

```text
product chi
 =product U
  -sum_r (U_r-chi_r)
       product_(s<r) chi_s product_(s>r) U_s.          (10)
```

Here `0<=chi_s<=U_s` and
`0<=U_r-chi_r<=U_r-L_r`.  Enlarging the nonnegative factors in every
subtracted summand gives (9).

This is the load-bearing distinction from THM-537's real-analyticity wall.
We do **not** ask for a nonnegative trigonometric-polynomial minorant supported
inside an interval.  The functions `L_r`, and consequently `P^-`, may be
negative.  Inequality (9) is instead a signed, pointwise-valid tensor
correction.

If the side lengths are `alpha_1,...,alpha_k`, then (7)--(8) give the ambient
Haar integral

```text
B_H(alpha_1,...,alpha_k)
 =product_r(alpha_r+epsilon)
  -2epsilon*sum_r product_(s!=r)(alpha_s+epsilon).     (11)
```

## 4. Relation-free characters make the bound exact

Let `v=(v_1,...,v_k)` be an integer vector with `lambda(v)>H`. Every Fourier
monomial of (8) has frequency `m in [-H,H]^k`.  On the orbit

```text
t ->(v_1t,...,v_kt),                                   (12)
```

its integral is zero unless `m.v=0`.  The relation-height assumption makes
`m=0` the only survivor.  Hence orbit integration of (8) equals its ambient
torus constant coefficient, and (9) proves

```text
measure{t:v_rt in J_r for every r}
 >=B_H(alpha_1,...,alpha_k).                           (13)
```

This is finite Fourier exactness, not an asymptotic equidistribution claim.

## 5. Degree 57 versus the relative-Hunter deficit

Assume for contradiction that `Lambda(h,Q)>57`.  Take

```text
H=57,                 epsilon=1/58.                    (14)
```

For every `i`, the two-frequency vector `(h,q_i)` has relation height greater
than `57`; otherwise append a zero coefficient to obtain a forbidden triple
relation. Apply (13) to the side lengths `2/7,1/7`.  With

```text
u_2=2/7+epsilon,       u_1=1/7+epsilon,                 (15)
```

we obtain

```text
I_i=measure(E_h intersect D_(q_i))
 >=u_2*u_1-2epsilon*(u_1+u_2)
 =5363/164836.                                          (16)
```

For every `i<j`, apply (13) to `C_h,D_(q_i),D_(q_j)`,
whose side lengths are `5/7,1/7,1/7`.  With

```text
u_0=5/7+epsilon,                                       (17)
```

this gives every restricted edge the lower bound

```text
w_ij=measure(C_h intersect D_(q_i) intersect D_(q_j))
 >=u_0*u_1^2-2epsilon*(u_1^2+2u_0*u_1)
 =655135/66923416.                                     (18)
```

Every spanning tree on seven vertices has six edges, so the maximum
spanning-tree weight from THM-2081 satisfies

```text
tau_h(Q)>=6*(655135/66923416).                          (19)
```

At the same time (16) bounds the THM-2081 deficit by

```text
Delta_h(Q)=2/7-sum_i I_i
 <=2/7-7*(5363/164836).                                (20)
```

Subtracting (20) from (19) gives

```text
tau_h(Q)-Delta_h(Q)
 >=6435/8365427
 >0.                                                   (21)
```

THM-2081 says that containment `G_Q subset E_h` forces
`tau_h(Q)<=Delta_h(Q)`, contradicting (21).  Therefore (4) holds. QED.

The exact referee also checks that the right side of (21), with the same
coordinatewise construction, is nonpositive for every `1<=H<=56`.  Thus 57
is the smallest degree closing **this certificate**.  It is not asserted to
be the smallest possible universal relation height.

## 6. Frontier effect

THM-2083's ineffective compactness conclusion is now explicit.  The surviving
depth-four branch is contained in the finite ledger of primitive coefficient
templates

```text
(a,b,c) in [-57,57]^3 minus {(0,0,0)},                 (22)
```

modulo common gcd and simultaneous sign.  Speeds on a template remain
unbounded, so this is not yet LRC(14).  The next exact target is to combine the
template equation with positivity, odd guard parity, divisor-completeness,
and the mirror-complement owner addresses of THM-2079. THM-2082's translated-
prime grids can discharge templates with sparse noncarrier incidence; the
claimed THM-2086 Fourier cone targets lacunary templates. The residual
projective planes still require endpoint/CRT inequalities or a sharper
relative-Hunter channel argument.

THM-2087 subsequently extracts the graph information implicit in this proof:
the relation-free pair graph must be disconnected, so a complete cut of at
least six height-57 relations exists. It converts that cut into either a
bounded guard ratio or a guard-anchored two-anchor star of height `6498`.

## 7. Assumption challenge and Tournament Analysis

The challenged assumption is that a useful spectral certificate must be a
nonnegative factorwise minorant.  That is false here: a signed box polynomial
can be pointwise valid after its defects are labelled by coordinates.  What
must remain nonnegative are the upper factors in the product telescope, not
the final minorant.

The natural vertices in this step are Fourier-coordinate defects, not runners,
arcs, or wall events.  A tournament on them is all ties: every coordinate is
charged once in (8).  The quotient preserves box membership and destroys the
endpoint order, which is harmless only in the relation-free branch.  The
weighted graphic matroid re-enters after (18), where all runner-pair edges
receive the same certified floor and any six-edge spanning tree suffices. QED.
