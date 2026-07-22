---
id: THM-2065
title: "Two-anchor Fejer relations collapse the strict LRC residual to circuit rays"
status: >
  PROVED REDUCTION; NOT LRC(14). THM-2051 forces every thirteen-speed row
  having no positive-measure strict lonely set to carry a support-three-to-five
  relation of height at most 2^20. On a fixed rank-two coefficient template,
  every such relation is either a persistent identity among the coefficient
  rows or cuts out one primitive projective parameter direction. Therefore a
  template without a height-bounded persistent circuit has only finitely many
  strict-null primitive rows. On a fixed-N THM-2058 interval, each nonzero
  pulled relation selects at most one integer M, after which THM-2062's CRT
  wheel filters hereditary primitivity exactly. The persistent-circuit branch
  remains open and is the sharply isolated structured residual.
source: codex-2026-07-21-LRC-circuit-ray-collapse
script: 04-computation/lrc_two_anchor_fejer_circuit_ray_codex_20260721.py
result: 05-knowledge/results/lrc_two_anchor_fejer_circuit_ray_codex_20260721.out
script_sha256: 07c7e8d3a3f04372efe52678d84b8b7b2418c1897bcb15b12860025fe185e5bd
result_sha256: 54a740db647a7c92fa6168e4253f0fc7357ea3491cdeee715679f3aa1c669b07
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2051
  - THM-2052
related:
  - THM-2058
  - THM-2062
  - THM-2064
  - THM-730
  - HYP-8846
---

# THM-2065 -- Two-anchor Fejer circuit-ray collapse

This theorem composes THM-2051's Fejer--BV relation alternative with the
THM-2052/2058 two-anchor atlas. It is a strict-branch terminal, not a proof of
LRC(14).

## 1. The exact pulled relation packet

Let

```text
r_i=(a_i,m_i) in Z^2,          1<=i<=13,
d=(N,M),                       gcd(N,M)=1,
v_i(d)=r_i.d=a_i N+m_i M.                              (1)
```

Fix an admissible parameter cone on which the thirteen `v_i(d)` are distinct
positive integers. Put `H=2^20`, and let `K_H` be the finite set of primitive
coefficient vectors `k in Z^13`, modulo simultaneous sign, such that

```text
3<=|supp(k)|<=5,       0<|k_i|<=H on supp(k).           (2)
```

For `k in K_H`, define its pulled coefficient row

```text
R_k=sum_i k_i r_i=(A_k,B_k).                            (3)
```

Call `k` **persistent** or **internal** if `R_k=0`. Otherwise define its
exceptional projective ray

```text
e(k)=[B_k:-A_k] in P^1(Q).                              (4)
```

Let

```text
L_str(d)={t in R/Z:||v_i(d)t||>1/14 for every i}.       (5)
```

Then every admissible primitive `d` with `measure(L_str(d))=0` satisfies the
exact dichotomy

```text
some k in K_H is persistent,
or
[d]=e(k) for some nonpersistent k in K_H.               (6)
```

In particular, if the template has no persistent member of `K_H`, its whole
strict-null locus is contained in the explicit finite ray packet

```text
E_H(r)={e(k):k in K_H, R_k!=0}.                          (7)
```

This includes every hypothetical strict LRC(14) counterexample in the
template, and also includes tight boundary rows whose strict safe set has
measure zero.

### Proof

THM-2051B says that `measure(L_str(d))=0` forces an integer relation

```text
sum_i k_i v_i(d)=0                                     (8)
```

with (2). Divide by the coefficient gcd and fix its sign, so `k in K_H`.
Linearity gives

```text
sum_i k_i v_i(d)
 =sum_i k_i(r_i.d)
 =(sum_i k_i r_i).d
 =R_k.d.                                                (9)
```

If `R_k=0`, the relation is persistent. If `R_k=(A_k,B_k)!=0`, its rational
kernel is the line spanned by `(B_k,-A_k)`. Writing

```text
g_k=gcd(|A_k|,|B_k|),                                  (10)
```

the only primitive integer points on that line are

```text
+-(B_k/g_k,-A_k/g_k).                                  (11)
```

Thus `[d]=e(k)`, proving (6)--(7). QED.

The number of canonical patterns, and hence the number of distinct rays, is
at most

```text
R_H=(1/2) sum_(s=3)^5 binom(13,s)(2H)^s
   =26103468074956706973500261957369856.                (12)
```

This deliberately crude bound precedes primitive-coefficient normalization,
duplicate slopes, positivity, owner-sector, collision, determinant-gate, and
hereditary-primitivity filters. In an actual template the packet is normally
far smaller.

## 2. Fixed-N collapse and the CRT wheel

Fix `N>0`. A nonpersistent relation `R_k=(A_k,B_k)` can occur only at

```text
M_k(N)=-N A_k/B_k.                                     (13)
```

If `B_k=0`, it never occurs. If `B_k!=0`, it selects at most one integer `M`,
and one exists exactly when `B_k` divides `-N A_k`. Hence define the finite
exceptional set

```text
E_N(r)={-N A_k/B_k in Z:
          k in K_H, R_k!=0, B_k!=0}.                    (14)
```

When no persistent `K_H`-circuit exists, every strict-null point of every
THM-2058 longitudinal interval `[L,U]` lies in

```text
E_N(r) intersect [L,U].                                (15)
```

Now assume the rows generate `Z^2`. In THM-2062's all-rank-two branch, let
`Q_0` and `W_N subset Z/Q_0Z` be the exact hereditary-primitivity wheel. The
hereditarily primitive strict-null candidates are contained losslessly in

```text
{M in E_N(r) intersect [L,U]:
  gcd(N,M)=1 and M mod Q_0 in W_N}.                     (16)
```

Collision walls, positivity, owner cones, and phase packets may then delete
more candidates. In THM-2062's rank-one-deletion branch, intersect (14)
instead with its two affine `+-1` terminals and one-dimensional coprime wheel.

Thus the two mechanisms have complementary jobs:

```text
THM-2051/2065: continuous strict failure -> finitely many exact rays;
THM-2062:      hereditary gcd conditions -> exact prime-labelled wheel. (17)
```

The wheel filters a circuit ray; it cannot create one. Conversely, its
positive density on the ambient plane says nothing about whether one
particular circuit ray survives.

## 3. Why the persistent branch is exact and necessary

The rows

```text
r_1=(1,0),       r_2=(0,1),       r_3=(1,1)             (18)
```

carry the persistent height-one circuit

```text
r_1+r_2-r_3=0.                                         (19)
```

Consequently `v_1(d)+v_2(d)-v_3(d)=0` for every parameter `d`; the relation
does not cut the plane to a ray. THM-2051 is existential, so the presence of
even one persistent bounded circuit makes its conclusion compatible with
every parameter. One must not discard the persistent branch or assume that
another, nonpersistent relation is also supplied.

The branch is finite and exact at the template level. For a subset `A` of
three to five rows, compute the integer kernel of the `2 x |A|` row matrix
and ask whether it meets the all-nonzero box `[-H,H]^A`. For a rank-two
triple `(r_i,r_j,r_l)`, its primitive kernel is obtained by dividing the
signed-minor vector

```text
(det(r_j,r_l), det(r_l,r_i), det(r_i,r_j))              (20)
```

by its coefficient gcd. Equation (20) also explains why merely saying
"the template has a circuit" is meaningless: every rank-two triple has a
rational circuit. The load-bearing condition is the coefficient-height cutoff
`H=2^20`, together with nonzero support and the owner/deck labels.

Inside THM-2052's rank-eleven template kernel, a nonpersistent `k` occurring
on (4) supplies the missing twelfth scalar relation at that specialized row.
A persistent `k` already lies in the rank-eleven template kernel and gives no
rank gain. This is the exact algebraic distinction the atlas must retain.

## 4. Higher-anchor form and scope

For a rank-`r` parameter template `r_i in Z^r`, the same proof replaces (4)
by the rational hyperplane

```text
{[d] in P^(r-1)(Q):R_k.d=0}.                            (21)
```

Thus rank two is special: a nonzero pulled relation becomes a projective
point rather than a positive-dimensional exceptional family. This is a
general finite-relation pullback principle, but THM-2051 supplies the needed
relation alternative only for the LRC(14) strict-measure problem considered
here.

The same dimensional boundary appears in the incoming planar-Jacobian work
THM-2063/HYP-8905: one effective direction is rigid while multiple directions
permit resonance. There is no reduction between those problems. Their shared
mechanism is only that a homogeneous relation has a zero-dimensional
projective kernel in two parameters.

## 5. What remains

THM-2065 closes the circuit-free template branch up to an explicit finite
candidate list. It does **not** say that a relation implies failure; many
strictly lonely rows have small relations. Nor does it discharge a persistent
circuit.

The remaining structured task is now precise:

> Classify saturated rank-two templates carrying a persistent marked
> support-three-to-five circuit of height at most `2^20`, and show that every
> owner/deck interval surviving the THM-2062 wheel either gains a
> nonpersistent bounded relation, enters a settled seam, or has an explicit
> phase-height, Fejer, pair-sum, or Euler certificate.

THM-730 and the older Freiman/additive-energy files classify special highly
coherent relation systems, but no canonical theorem currently handles one
arbitrary persistent marked circuit with these sidecars.

## Assumption challenge and tournament analysis

The useful vertices are signed relation patterns, not runners. The pairwise
observable is equality of their pulled projective kernels; quotienting by it
deduplicates (7). Any orientation can be tied by support size, height, and
lexicographic coefficient word, producing a transitive scheduling tournament.
Its score histogram only orders the finite checks: it destroys the signed
linear form, internal/noninternal distinction, owner sector, and CRT residue.
There is therefore no mathematical tournament certificate here. The faithful
carrier is the signed circuit hypergraph plus its pulled row and atlas
sidecars.

## Computational audit

The companion uses integer arithmetic only. It exhausts primitive projective
kernels in a bounded box, fixed-N uniqueness, the complete height-two packet
on deterministic five-row templates, direct hereditary-wheel intersections,
the persistent-circuit guardrail, and `10,000` signed-minor identities. It
does not re-prove THM-2051's analytic Fejer--BV alternative. Runtime checks
survive ordinary and optimized Python; the frozen output ends in `PASS`.
