---
id: THM-2264
title: "Signed Selberg--Kraft relation spectrum"
status: >
  SUPERSEDED / EXACT DUPLICATE OF THM-2144. The anisotropic signed-box
  theorem, height-29 gate, total Kraft budget 367, and height-105
  rank-or-subset-sum fork below were already proved in THM-2144. Section 3's
  table of selected coordinatewise Kraft specializations is a convenient
  corollary only. THM-2164 is strictly stronger at the final rank step:
  distinct zero-safe LRC rows satisfy rank(W_105)>=2 unconditionally.
source: codex-2026-07-25-signed-selberg-kraft
depends_on:
  - THM-2144
related:
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2164-relative-packet-rank-harvesting
  - THM-2085-explicit-height-57-rank-seven-selberg-gate
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
script: 04-computation/lrc14_signed_selberg_kraft_relation_spectrum_referee_codex_20260725.py
output: 05-knowledge/results/lrc14_signed_selberg_kraft_relation_spectrum_referee_codex_20260725.out
script_sha256: 62eb53f3c859274bb2c396d8954d0a2a4c06ef54e2d7aba574760eb5859ae91b
output_sha256: 26ec7df7421a2678f21385de742094076b93897e0f7dde0a3991d42d51e9d0ca
hash_basis: working-tree bytes (LF)
---

> **SUPERSEDED / EXACT DUPLICATE.** Sections 1, 2, and 4 restate
> THM-2144, while Section 3 merely tabulates direct substitutions in
> THM-2144's Kraft inequality. For distinct zero-safe LRC rows, THM-2164
> removes the subset-sum alternative and proves `rank W_105>=2`.
> This file is retained only to preserve correction lineage and the
> convenience table; it is not an independent proved dependency.

# THM-2264 -- signed Selberg--Kraft relation spectrum

Let `v=(v_1,...,v_k)` be a vector of positive integers and put

```text
M(v)=max_(t in R/Z) min_i ||v_i t||.                    (1)
```

For an integer relation `m.v=0`, its anisotropic height box is described by
coordinate bounds `|m_i|<=H_i`.

## 1. The anisotropic signed-box theorem

Fix `0<h<1/2`, put `alpha=1-2h`, and choose integers `H_i>=1`. If

```text
sum_i 2/((H_i+1)alpha+1)<1,                            (2)
```

then exactly one of the following conclusions is needed:

```text
M(v)>=h,                                                (3)
```

or there is a nonzero `m in Z^k` such that

```text
m.v=0,                  |m_i|<=H_i for every i.         (4)
```

More precisely, if (4) does not exist, the closed `h`-safe set has Haar
measure at least

```text
B(alpha;H_1,...,H_k)
 =product_i(alpha+epsilon_i)
    * (1-sum_i 2epsilon_i/(alpha+epsilon_i))>0,
epsilon_i=1/(H_i+1).                                   (5)
```

### Proof

For `J=[h,1-h]`, take in coordinate `i` the degree-`H_i` Vaaler lower and
upper interval polynomials `L_i<=1_J<=U_i`. Their constant coefficients are
`alpha-epsilon_i` and `alpha+epsilon_i`. THM-2085 proves the pointwise signed
tensor inequality

```text
P^-=product_i U_i
    -sum_i (U_i-L_i) product_(j!=i) U_j
   <=product_i 1_J(x_i).                                (6)
```

The proof there is coordinatewise and does not require equal degrees. Its
constant coefficient is exactly (5).

Every Fourier character of `P^-` has coordinate frequencies
`|m_i|<=H_i`. On the orbit

```text
t |-> (v_1t,...,v_kt),                                  (7)
```

the integral of a character is zero unless `m.v=0`. In the absence of (4),
only the constant character survives. Integrating (6) therefore gives (5).
Condition (2) is precisely its positivity condition, proving (3)--(4).
Endpoints are harmless: the orbit meets only finitely many interval endpoints,
and all statements are Haar-measure statements. QED.

## 2. A uniform relation-height theorem

Let `n>=5` be the total number of runners, so `k=n-1`. Take

```text
H_i=2n+1,             epsilon=1/(2n+2).                 (8)
```

For equal side length `alpha`, (5) simplifies to

```text
B=(alpha+epsilon)^(n-2)
     * (alpha-(2n-3)epsilon).                           (9)
```

It is positive whenever

```text
h<5/(4(n+1)).                                          (10)
```

Consequently every positive integer speed vector satisfies

```text
M(v)>=5/(4(n+1))                                       (11)
```

or has a nonzero relation `m.v=0` with

```text
||m||_infinity<=2n+1.                                  (12)
```

The limiting value in (11) follows by applying (5) for every strict
`h<5/(4(n+1))`. Since `5/(4(n+1))>=1/n` for `n>=4`, any LRC counterexample
in this range lies on one of the bounded hyperplanes (12).

For LRC(14), `n=14`, `k=13`, and `H=29`. At `h=1/14`,

```text
B=(187/210)^12*(1/42)>0.                               (13)
```

Moreover (10) is `h<1/12`. Hence

```text
no relation of height <=29  ==>  M(v)>=1/12,            (14)
```

and, in particular,

```text
M(v)<1/14
 ==> exists 0!=m in [-29,29]^13 with m.v=0.             (15)
```

For this equal-degree certificate, `29` is the first successful degree:
at `H=28` the last factor in (9), evaluated at `h=1/14`, is `-1/203`.
No optimality for the true relation height is asserted.

## 3. The Selberg--Kraft coefficient budget

At `h=1/14`, condition (2) becomes the exact anisotropic criterion

```text
sum_(i=1)^13 14/(6H_i+13)<1.                            (16)
```

Thus every hypothetical LRC(14) counterexample has a relation in every
coefficient box satisfying (16). Useful consequences, with all unlisted
coordinates bounded by `29`, are

```text
number of special coordinates     special bound
              1                         21
              2                         25
              3                         26
              4 or 5                    27
              6 through 10              28.             (17)
```

For example, the tightest displayed rows have positive Kraft margins

```text
1x21+12x29: 23/25993,       5x27+8x29: 1/935,
10x28+3x29: 65/33847.                                  (18)
```

The least possible total coordinate budget supplied by (16) is

```text
sum_i H_i=367.                                          (19)
```

Indeed `H |->14/(6H+13)` is decreasing and strictly convex. At fixed integer
sum its total is minimized by a balanced profile. The balanced sum-366 and
sum-367 profiles give, respectively,

```text
11*(14/181)+2*(14/187)=33866/33847>1,
10*(14/181)+3*(14/187)=33782/33847<1.                  (20)
```

This is optimal only for the signed tensor certificate (6).

## 4. First rank harvest at height 105

Let

```text
W_105(v)=span_Q{m in Z^13:m.v=0, ||m||_infinity<=105}. (21)
```

For any chosen coordinate `i`, use `H_i=1` and `H_j=105` for `j!=i`.
The Kraft sum is

```text
14/19+12*(14/643)=12194/12217<1.                        (22)
```

Therefore a hypothetical LRC(14) counterexample has, for every `i`, a
nonzero relation in `W_105(v)` whose `i`-th coefficient has magnitude at
most one.

If `dim W_105(v)=1`, let its integer lattice be generated by a primitive
vector `m`. Every integer relation on that line is a nonzero integer multiple
of `m`. Applying the preceding relation separately at all thirteen
coordinates forces

```text
m_i in {-1,0,1} for every i.                            (23)
```

Because all speeds are positive, `m` has both signs. Thus:

```text
dim W_105(v)>=2,                                        (24)
```

or there are disjoint nonempty index sets `A,B` with

```text
sum_(i in A) v_i=sum_(j in B) v_j.                      (25)
```

This is a genuine low-height addition to THM-2052's high-rank sparse relation
code. It does not say that the height-29 or height-105 rows are independent
of THM-2052's eleven-dimensional code, so it does not by itself close a
rank-eleven star or the rank-twelve finite box.

## 5. Verification and boundary

The dependency-free exact referee checks (9), (13), (16)--(22), the first
failed equal-degree row, and the sharp total-budget transition. Reproduce with

```text
python 04-computation/lrc14_signed_selberg_kraft_relation_spectrum_referee_codex_20260725.py
python -O 04-computation/lrc14_signed_selberg_kraft_relation_spectrum_referee_codex_20260725.py
```

The theorem forces bounded relations, not a finite speed box. A subsequent
rank-harvesting step still needs either a character outside the known relation
space or an exact classification of the low-height hyperplane intersections.
