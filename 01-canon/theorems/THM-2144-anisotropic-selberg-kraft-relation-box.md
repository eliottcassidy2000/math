---
id: THM-2144
title: "Anisotropic Selberg--Kraft relation boxes"
status: >
  PROVED from THM-2085's signed Selberg tensor construction. A zero-measure
  safe set forces an integer relation in every strictly Kraft-admissible
  anisotropic coefficient box. For thirteen speeds this gives the first
  equal-cap certificate H=29, the certificate-minimal total cap budget 367,
  and the dichotomy rank(W_105)>=2 or a signed subset-sum equality. These are
  certificate statements, not a proof of LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2085
related:
  - THM-537
  - THM-2051
  - THM-2052
  - THM-2054
script: 04-computation/lrc14_anisotropic_selberg_kraft_referee_codex_20260724.py
output: 05-knowledge/results/lrc14_anisotropic_selberg_kraft_referee_codex_20260724.out
script_sha256: 4d62acd3eeb39e6df8e61249887a876c4439f6e40990f23fb433b67fb6724795
output_sha256: 2015ab9bbdab3eecf7d0cf8df4e20583e9e85e31ba73dd7e9e5f6f8f45b07e30
hash_basis: working-tree bytes (LF)
---

# THM-2144 -- anisotropic Selberg--Kraft relation boxes

For positive integers `v_1,...,v_k`, let

```text
S_h(v)={t in R/Z: ||v_i t||_(R/Z)>=h for every i}.      (1)
```

For integer caps `H_i>=1`, write

```text
R_H(v)={m in Z^k: m.v=0 and |m_i|<=H_i for every i}.    (2)
```

The zero vector belongs to `R_H(v)` and is never called a relation below.

## 1. Master anisotropic theorem

Let `0<h<1/2`, put

```text
alpha=1-2h,          epsilon_i=1/(H_i+1),
u_i=alpha+epsilon_i.                                  (3)
```

If `R_H(v)={0}`, then

```text
mu(S_h(v))
 >= product_i u_i
    * (1-sum_i 2epsilon_i/u_i).                        (4)
```

Consequently,

```text
mu(S_h(v))=0
and
sum_i 2epsilon_i/u_i<1                                (5)
```

force a nonzero relation `m.v=0` with `|m_i|<=H_i`.

### Proof

For the circle interval `J=[h,1-h]` of length `alpha`, choose in coordinate
`i` the degree-`H_i` Selberg pair `L_i<=1_J<=U_i` used in THM-2085. Its
constant coefficients are

```text
integral L_i=alpha-epsilon_i,
integral U_i=alpha+epsilon_i=u_i.                       (6)
```

The signed tensor polynomial

```text
P^-(x)
 =product_i U_i(x_i)
  -sum_i (U_i(x_i)-L_i(x_i))*product_(j!=i) U_j(x_j)   (7)
```

satisfies

```text
P^-(x)<=product_i 1_J(x_i)                              (8)
```

away from endpoints. This is exactly THM-2085's product-telescope argument;
the tensor polynomial is signed, so THM-537's nonnegative analytic-minorant
obstruction is not implicated.

The Fourier support of `P^-` lies in the anisotropic box
`product_i[-H_i,H_i]`. On the integer orbit

```text
t |-> (v_1t,...,v_kt),                                  (9)
```

a character `m` has integral zero unless `m.v=0`. Under `R_H(v)={0}`, only
the constant character survives. Integrating (7) therefore gives

```text
integral P^-(v_1t,...,v_kt)
 =product_i u_i-sum_i 2epsilon_i*product_(j!=i)u_j,     (10)
```

which factors as the right side of (4). The preimage of every interval
endpoint is finite because each `v_i` is positive, so endpoint conventions
do not change Haar measure. Equations (8)--(10) prove (4), and (5) follows by
contradiction. QED.

This is finite Fourier exactness. It uses neither asymptotic decorrelation nor
independence of the coordinate events.

## 2. Equal caps and the height-29 LRC(14) gate

For a common cap `H`, write `epsilon=1/(H+1)`. Formula (4) becomes

```text
mu(S_h(v))
 >=(alpha+epsilon)^(k-1)
   * (alpha-(2k-1)epsilon).                            (11)
```

Take `k=n-1` and `H=2n+1`. If no relation has height at most `2n+1`, then
(11) is positive for every

```text
h<5/(4(n+1)).                                           (12)
```

Thus

```text
M(v):=sup_t min_i ||v_i t|| >=5/(4(n+1)).              (13)
```

The endpoint in (13) is obtained by a supremum; the displayed tensor bound
vanishes there. For `n>=4`, the right side is at least `1/n`.

For the thirteen relative speeds of LRC(14), at `h=1/14` and `H=29`,

```text
alpha=6/7, epsilon=1/30,
B_29=(187/210)^12*(1/42)>0.                             (14)
```

At `H=28`, the final factor in (11) is

```text
6/7-25/29=-1/203.                                      (15)
```

Therefore `29` is the first common integer cap accepted by this particular
equal-degree tensor certificate. More strongly, absence of a height-29
relation gives `M(v)>=1/12` by (13). No positive-measure assertion is made at
`h=1/12` itself.

## 3. The Kraft spectrum and total cap budget 367

At the LRC(14) threshold `h=1/14`, condition (5) is exactly

```text
sum_(i=1)^13 14/(6H_i+13)<1.                            (16)
```

This is the anisotropic Selberg--Kraft test. Define

```text
f(H)=14/(6H+13).                                       (17)
```

The sequence is strictly decreasing and discretely convex. If
`A>=B+2`, the balancing exchange

```text
(A,B) |-> (A-1,B+1)                                    (18)
```

strictly decreases `f(A)+f(B)`. Hence a balanced integer profile minimizes
the left side of (16) at fixed total cap.

The adjacent exact profiles are

```text
11 f(28)+2 f(29)=33866/33847>1,   total 366,            (19)
10 f(28)+3 f(29)=33782/33847<1,   total 367.            (20)
```

Thus every thirteen-speed vector with zero-measure safe set has, for any
chosen placement of three cap-29 coordinates and ten cap-28 coordinates, a
nonzero relation in that box. Total cap `367` is minimal among profiles
certified by (16). This is certificate-optimality only; it does not assert
that the true minimal relation box has total cap 367.

## 4. Rank 105 or a signed subset-sum

For a thirteen-speed vector define the bounded relation span

```text
W_105(v)
 =span_Q{m in Z^13: m.v=0 and ||m||_infinity<=105}.     (21)
```

The one-small-coordinate profiles satisfy

```text
f(1)+12f(104)=1730/1729>1,
f(1)+12f(105)=12194/12217<1.                            (22)
```

Suppose `mu(S_(1/14)(v))=0`. For each coordinate `i`, apply (16) with cap
`1` at `i` and cap `105` elsewhere. This supplies a nonzero
`r^(i) in W_105(v)` with

```text
|r^(i)_i|<=1.                                           (23)
```

If `dim_Q W_105(v)=1`, let `m` be the primitive generator of the rank-one
lattice `W_105(v) intersect Z^13`. Then `r^(i)=a_i m` for a nonzero integer
`a_i`. Equation (23) gives `m_i in {-1,0,1}` for every `i`. Since the speeds
are positive and `m.v=0` is nonzero, both signs occur. The disjoint nonempty
sets

```text
A={i:m_i=1},        B={i:m_i=-1}                       (24)
```

therefore obey

```text
sum_(i in A) v_i=sum_(j in B) v_j.                     (25)
```

We have proved the dichotomy

```text
dim_Q W_105(v)>=2
or
v has a disjoint nonempty signed subset-sum equality.  (26)
```

The number `105` is the first uniform complementary cap in the particular
`(1,H,...,H)` Kraft profile. It is not claimed to be the optimal rank-forcing
height.

## 5. Scope and frontier effect

The theorem applies to arbitrary positive integer speeds; distinctness and
primitive normalization are unnecessary. The forced relations may have full
support. In particular:

1. (16) is an anisotropic spectrum of relation boxes, not a sparse-relation
   theorem.
2. (26) is the first unconditional rank-or-subset-sum fork produced by this
   tensor mechanism, but it does not raise the sparse bounded-relation rank in
   THM-2052.
3. Neither branch of (26) bounds the speeds or proves LRC(14).
4. The useful next target is a relative version modulo an already known
   relation lattice: force a new independent relation, not merely another
   vector in the same span.

The companion script checks (14)--(15), the endpoint in (13), the exact
balancing exchanges through the relevant range, (19)--(20), and (22) using
rational arithmetic. Normal and optimized Python produce identical output.
