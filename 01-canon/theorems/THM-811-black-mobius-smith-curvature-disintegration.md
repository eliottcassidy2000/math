---
id: THM-811
title: The black complement-line flow has an exact Mobius/Smith master polynomial and negative curvature-energy coupling
status: PROVED (general master polynomial, coefficients, moments, and support law) + FINITE-EXACT (node/edge classifications and q-stratified flow through n=7)
source: codex-2026-07-15-S13
depends_on: [THM-785, THM-790, THM-801]
related: [THM-781, THM-796, THM-809, HYP-6880]
verification:
  - 04-computation/black_mobius_curvature_disintegration_codex_S13.py
  - 05-knowledge/results/black_mobius_curvature_disintegration_codex_S13.out
  - 05-knowledge/results/black_mobius_curvature_disintegration_codex_S13.json
---

# THM-811 — exact black curvature disintegration

For the apex-zero endpoint `t` of a complement line put

```text
S=C3(kappa t)-C3(t),
q0=Omega C3(t),          q1=Omega C3(kappa t),
epsilon=e1+e_n-(n-2),    r=n-4,
M=binom(n-1,2).
```

Then the complete four-variable line polynomial is

```text
sum_(e in E_n) u^S z^q0 w^q1 v^epsilon
 =2^(M-2n+5) u^(n-2)
  (1+w v/u)(1+w/(uv))
  (1+v/u+1/(uv)+z w/u^2)^r.                              (1)
```

It is the common source of THM-785's complement flux/defect law and
THM-801's boundary-curvature polynomial.  In particular, at `u=1`,

```text
sum_e z^q0 w^q1 v^epsilon
 =2^(M-2n+5)(1+wv)(1+w/v)(1+v+v^-1+zw)^r.               (2)
```

The black polynomial is (2) minus THM-801's blue `(q0,q1)` polynomial placed
entirely at `v^0`.  Every nonzero-`epsilon` coefficient is therefore black.

## Proof and the four local states

Write `B` and `T` for the one-counts on the bottom and top marked-path legs.
The apex-zero convention leaves two boundary leg bits and `r` paired internal
positions.  At an internal position the four states contribute

```text
00 -> 1,       10 -> v/u,       01 -> v^-1/u,
11 -> z w/u^2.                                         (3)
```

Each selected boundary bit contributes `w/u` and respectively `v` or `v^-1`.
The remaining `M-(2r+2)` bits are free.  Multiplying these independent local
contributions proves (1).  Equivalently,

```text
epsilon=B-T,                   lambda=B+T-(n-2),
S=-lambda=n-2-B-T,
q0=sum_internal B_i T_i,       q1=q0+B_boundary+T_boundary. (4)
```

Thus `(lambda,epsilon)` are first-order leg totals, while `q` is their
positional quadratic overlap.

THM-801's coarse B3 stratum counts already determine the currents:

```text
B=S3_A+S3_AB,                 T=S3_C+S3_BC,
epsilon=S3_A+S3_AB-S3_C-S3_BC,
lambda=S3_A+S3_AB+S3_C+S3_BC-(n-2).                     (5)
```

They do not determine `q` from `n=6` onward.  This precisely separates B3
linear current from Möbius quadratic curvature.

## Coefficients and curvature uncertainty

Let `j=q1-q0 in {0,1,2}` and put

```text
C_0(v)=C_2(v)=1,             C_1(v)=v+v^-1.
```

For a black line,

```text
N_K(k,k+j,e)
 =2^(M-2n+5) binom(r,k)
  [v^e] C_j(v)(1+v+v^-1)^(r-k)
  -1_(e=0) N_B(k,k+j).                                  (6)
```

Before the blue zero-slice is removed,

```text
Var(epsilon | q0=k,j)=2(r-k)/3+1_(j=1),                 (7)
|epsilon| <= r-k+1_(j=1).                               (8)
```

The coefficient row is symmetric and trinomial-unimodal.  Hence high boundary
cyclicity forces narrow transverse Smith curvature; this is a literal support
law, not merely a correlation.

The exact unconditioned moments are

```text
Cov(q0,epsilon^2)=Cov(q1,epsilon^2)=-r/8,
Cov(q1-q0,epsilon^2)=0.                                  (9)
```

Let `a=beta_n` be the blue line fraction, `p=n-3`, and let `delta_n` be zero
for even `n` and `1/4` for odd `n`.  Conditioning on black gives

```text
Cov_K(q_i,epsilon^2)
 = (-(r/8)(1-a)+a(p/2)delta_n)/(1-a)^2 < 0              (10)
```

for every `n>=5`.  The values at `n=5,6,7` are respectively
`-1/18,-4/15,-1480/3969`.

## Node curvature polynomial

For a merged node `x`, define the black reflection-orbit polynomial

```text
H_x(z,v)=1/2 sum_(t in F_x, t black) z^q(t) v^|epsilon(t)|. (11)
```

Reflection fixes `x,q`, reverses signed `epsilon`, and acts freely on black
tilings, so every pre-halving coefficient is even and (11) has an exact orbit
meaning.  The carrier `(C3(x),H_x)` separates every merged node through `n=7`:

```text
n                    4    5    6     7
nodes/cells           3    10   34    272.
```

At `n=7`, THM-801's all-tiling `q` polynomial gives 238 cells; its black-only
version gives 235.  Adjoining `C3` gives 249, black `|epsilon|` gives 263,
joint black `(q,|epsilon|)` gives 266, and adjoining `C3` gives all 272.  The primitive
joint shape together with black mass and `C3` is also complete, so this is not
an accidental raw-scaling label.

## Curvature controls flow but does not identify an edge

On black reflection orbits at `n=7`, projected node pair plus `(q0,q1)` gives
7,206 cells out of 8,064.  Adding `|epsilon|` gives only 7,248.  By contrast,
projected node pair plus the unordered reflection orbit of `(B2,B3)` is exact
on all 8,064 orbits.  The corresponding exact counts for `n=4..7` are

```text
1/1, 12/12, 240/240, 8064/8064.                          (12)
```

All 16 residual literal `(node pair,B2,B3)` collisions at `n=7` are precisely
the 16 reflection pairs.  Every black projected loop is already separated
literally by `(B2,B3)`.  Thus curvature is the right disintegration coordinate
but cannot replace the positional edge address.  Any reflection-invariant
carrier can classify at most a black edge orbit, never a literal black line.

## The finite source-normalized flow signal

Orient a non-tied black boundary line toward larger `C3`.  Raw multiplicity at
`n=7` points from pure black back to mixed, as THM-785 found.  After dividing
separately by the available source black endpoints, however, mixed-to-black
flow is larger in every source-`q` stratum carrying at least one
mixed--pure-black boundary edge at both `n=6` and `n=7`.  The nonempty
`n=6,q=2` source populations carry no such boundary edge and give the vacuous
equality `0=0`; they are not a strict-dominance row.
At `n=7` the outward/reverse rates for `q=0,...,5` are approximately

```text
.497/.172, .163/.105, .0115/.00735,
.263/.130, .625/.309, .722/0.                             (13)
```

This is finite-exact evidence for an all-size outward-dominance conjecture,
not part of the general theorem.  Stratifying only by `|epsilon|` does not
have this property: the `|epsilon|=1` row can reverse and the largest tails
are reverse-only.  Boundary cyclicity, rather than transverse magnitude, is
therefore the sharper control variable for the phase-volume reversal.

## Preservation boundary

The faithful black object here is a reflection orbit of literal complement
lines with both endpoints, multiplicity, and apex orientation.  The master
polynomial preserves longitudinal flow, transverse Smith current, boundary
Möbius curvature, and their exact joint distribution.  It destroys mirror
position and therefore does not recover a literal edge.  None of these data
preserve the LRC metric witness, owner assignment, threshold side, centered
wall schedule, or sheet carry; THM-808's transported owner/root stalk remains
separate until a path-orbit lift relates the two.
