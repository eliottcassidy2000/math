---
id: THM-4451
title: "LRC14 dyadic three-tail strict component-width caps"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. Every strict
  physical clock-two failure component of three distinct positive odd tails
  has width at most 17/693, uniquely sharp at (1,9,11). In the odd-3-unit
  domain the sharp bound is 19/1001, uniquely at (1,11,13). Quotient widths
  are exactly twice these values. The proof is all-height; the residual exact
  boxes contain 123 and 209 triples. This does not close a residual row and
  LRC(14) remains open.
source: strict endpoint topology continuation + clean-room audit, 2026-09-06
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
  - THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry
related:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4450-lrc14-absorbed-label-overlap-hierarchy-and-component-address-decoder
primary_script: 04-computation/lrc14_dyadic_strict_component_width_thm4451.py
primary_output: 05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451.out
primary_script_sha256: c8e1a17b87fd0c7b71c4dcb378db42f6ad94fd9a419d5c33bb8832cbb2a4c238
primary_output_sha256: e80aef3abc01c484fd81b5e72d8767c9ef6d063ccc6a7f8e3fcf165579749a19
auxiliary_script: 04-computation/lrc14_dyadic_strict_component_topology_thm4451.py
auxiliary_script_sha256: e6f3dd3933076c05d452211ab6d4d4ef4601746e0735a19e9e5564c48c98a951
independent_script: 04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
independent_output: 05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451_independent.out
independent_script_sha256: 7c8ebcbc41beee193c422042a7f09b170c7ecb231325628208b7efb02aa2b92d
independent_output_sha256: 772f293b16edefb6750128327bdb7919df0f287dc3c82d5aa3e9e8883771fa9d
report: 05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451.md
audit: 05-knowledge/results/lrc14_dyadic_strict_component_width_thm4451_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4451 -- LRC14 dyadic three-tail strict component-width caps

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.** This theorem
concerns actual connected components of a strict open set, not a.e.-merged
interval representatives. That distinction repairs a false preliminary
candidate at `(1,7,13)`. The all-height reduction is analytic; exact
computation is used only on the explicitly bounded residual boxes.

For a positive integer `n`, define
\[
 D_n=\{x\in\mathbb R/\mathbb Z:\|nx\|<1/14\},\qquad
 \tau(x)=x+1/2.                                             \tag{1}
\]
For three distinct positive odd tails `T={a,b,c}`, put
\[
 A_T=\bigcup_{n\in T}D_n,\qquad F_T=A_T\cap\tau(A_T),       \tag{2}
\]
and let `lambda(T)` be the supremum of the lengths of the actual open
connected components of `F_T`.

## 1. Sharp component theorems

For every such odd triple,
\[
 \boxed{\lambda(T)\le {17\over693}},                        \tag{3}
\]
with equality uniquely for the individual triple
\[
 T=(1,9,11).                                                 \tag{4}
\]
If no member of `T` is divisible by three, then
\[
 \boxed{\lambda(T)\le {19\over1001}},                       \tag{5}
\]
with equality uniquely for
\[
 T=(1,11,13).                                                \tag{6}
\]
These are physical-`x` widths. In the quotient coordinate `y=2x`, antipodal
components are identified and all widths double, so the sharp quotient caps
are
\[
 {34\over693}\quad\hbox{and}\quad {38\over1001}.           \tag{7}
\]
Unlike Haar mass, component width is not invariant under common odd
dilation. Thus (4) and (6) are individual equality triples, not dilation
rays.

## 2. Two-lane capacity

Put `C_n=D_n union tau(D_n)`. For odd `n`, this is a comb of period `1/(2n)`
and tooth width `1/(7n)`. If `J` is an interval of length `L`, write
`2nL=m+r`, where `m=floor(2nL)` and `0<=r<1`. Exact one-period packing gives
\[
 \operatorname{cap}_n(L):=\sup_J\mu(J\cap C_n)
 ={2m/7+\min(r,2/7)\over2n}.                                \tag{8}
\]
On an interval contained in `F_T`, a first-sheet owner and a second-sheet
owner exist at every point. An odd tail cannot own both sheets, so at least
two distinct lane indicators are active. Hence
\[
 2L\le\sum_{n\in T}\operatorname{cap}_n(L).                 \tag{9}
\]
Writing
\[
 s_n(L)=\operatorname{cap}_n(L)-{2L\over7},
\]
equation (8) gives
\[
 0\le s_n(L)\le {5\over49n},\qquad
 \sum_{n\in T}s_n(L)\ge {8L\over7}.                       \tag{10}
\]

Sort `a<b<c` and put `Delta=8L/7`. Successive use of (10) yields
\[
 a\le {15\over49\Delta},
\quad b\le {10\over49(\Delta-s_a)},
\quad c\le {5\over49(\Delta-s_a-s_b)},                    \tag{11}
\]
whenever the displayed residual denominator is positive. At `L=17/693`,
this leaves smallest tails `1,3,5,7,9` and finitely many bounded families;
at `L=19/1001` in the odd-3-unit domain it leaves `1,5,7,11,13` and finitely
many bounded families.

## 3. Closing the nominally unbounded pairs

Let
\[
 P_{ab}=(D_a\cup D_b)\cap\tau(D_a\cup D_b),                 \tag{12}
\]
and let `P*_{ab}` be the explicit a.e. superset obtained by filling only
isolated deleted endpoints and then merging. Split (2) by its two owners.
The exact owner identity is
\[
 F_{\{a,b,c\}}=P_{ab}\cup\Sigma_{ac}\cup\Sigma_{bc},
 \qquad \Sigma_{ac}\cup\Sigma_{bc}\subset C_c.            \tag{13}
\]
If `P*_{ab}` is empty, (13) gives `F subset C_c`, so every component has
width at most `ell=1/(7c)`. Otherwise let `w` be its maximum component width
and `g` its least positive circular gap. When `ell<g`, no `C_c` tooth meets
two carrier components. The incidence graph is a star forest and every
component has width at most
\[
 \max(\ell,w+2\ell).                                        \tag{14}
\]

The only capacity-unbounded all-odd pairs are
\[
 (1,3),(1,5),(1,7),(3,5),(3,7),(5,7),                      \tag{15}
\]
with exact `(w,g)` values
\[
 (0,1),(0,1),(1/98,6/49),(1/210,5/42),
 (1/98,19/98),(1/98,1/14).                                 \tag{16}
\]
For odd 3-units only `(1,5),(1,7),(5,7)` remain. Equations (13)--(16)
strictly close the infinite tails after finite prefixes. In the bounded
all-odd branch, the formal `(7,11)` cutoff is `165/14>11`, but it contains no
admissible next odd `c` because `13>165/14`.

## 4. Exact residual boxes and strict topology

The reductions leave exactly 123 all-odd triples and 209 odd-3-unit triples.
Two exact representations evaluate every row:

1. union the six oriented owner-cross interval families, joining positive
   overlaps but never a shared deleted endpoint;
2. partition at every defining wall, test strict membership on cells and
   walls, and join adjacent live cells only through a live wall.

They agree on every row. The maxima, unique leaders, and runners-up are
\[
\begin{array}{c|c|c}
\text{domain}&\text{maximum and leader}&\text{runner-up}\\ \hline
\text{all odd}&17/693\text{ at }(1,9,11)&1/42\\
\text{odd 3-units}&19/1001\text{ at }(1,11,13)&29/1547.
\end{array}                                                  \tag{17}
\]
The first equality triple has four components of width `17/693` and four of
width `13/1386`. The second has four each of widths `19/1001`, `17/2002`,
and `3/2002`.

The endpoint hostile explains why strict topology is essential. At
`T=(1,7,13)`, the a.e.-merged interval of width `1/49` crosses `x=1/14`.
At that point no first-sheet owner is strict, so it is absent. The actual
multiset is four components of width `1/91` and eight of width `1/98`.
Thus `1/49 at (1,7,13)` is **REFUTED as a strict component value**, though
the a.e. merge remains legitimate for Haar measure.

The doubling map identifies antipodal component pairs. Direct exact
recomputation also verifies the factor two in (7) on all 332 residual rows.

## 5. LRC14 component-entry consequence

Whenever a clock-two reduction says that a closed connected quotient body
component `K` must lie inside one strict three-tail failure component, (7)
gives the sufficient escape tests
\[
 |K|\ge {34\over693}\quad\text{for arbitrary odd tails},
 \qquad
 |K|\ge {38\over1001}\quad\text{for odd 3-unit tails}.      \tag{18}
\]
Equality in (18) also forces escape: a compact interval cannot occupy the
entire equal-width open failure component. These are conditional component
entry gates, not universal lower bounds on a ten-body component. They do not
close the residual dyadic row.

## 6. Verification and scope

The primary program performs the capacity reduction, the exact owner-union
calculation, and the wall-partition calculation. The independent program is
self-contained and rebuilds the strict walls, both finite boxes, every
carrier value, equality topology, endpoint hostile, and quotient doubling.
Both programs and their optimized runs are line-identical to their frozen
outputs and end in `PASS`.

Oddness is essential in (9): for even `n`, `D_n=tau(D_n)` and one label can
own both sheets. The height-99 discovery scan is provenance, not a theorem
dependency. No finite-height extrapolation, measure-to-connectivity
substitution, or tournament completion is used. LRC(14) remains **OPEN**.
