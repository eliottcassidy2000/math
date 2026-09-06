---
id: THM-4452
title: "LRC14 clock-four one-v2 strict component-width caps"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. For the strict
  four-sheet failure set of tails (2r,a,b), with r,a,b odd and a distinct
  from b, every physical component has width at most
  min(1/(14r),1/(7 max(a,b))). The sharp universal caps are 1/98 for all
  odd labels and 1/110 for odd 3-units; quotient widths are four times these
  values. This gives original-body component entry gates 2/49 and 2/55 but
  no universal body floor. Clock four and LRC(14) remain open.
source: colored four-sheet parity decomposition + clean-room audit, 2026-09-06
depends_on:
  - THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction
related:
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
  - THM-4450-lrc14-absorbed-label-overlap-hierarchy-and-component-address-decoder
  - THM-4451-lrc14-dyadic-three-tail-strict-component-width-caps
primary_script: 04-computation/lrc14_clock_four_one_v2_component_thm4452.py
primary_output: 05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452.out
primary_script_sha256: 03432033c64bdd46f50788f01dbe8545e89e31fccac3b2c2bac23900418548f3
primary_output_sha256: 3f11751c9bff00ecaf39184f811a96d64b700fc14d5b12ad9b6a262bbaa5d8d6
auxiliary_script: 04-computation/lrc14_clock_four_one_v2_component_thm4452_geometry.py
auxiliary_script_sha256: ca1497c66ec0a7dbfc8c9103409e3cddba283abd9a6149088a9b8f6c914a4d65
independent_script: 04-computation/lrc14_clock_four_one_v2_component_thm4452_independent.py
independent_output: 05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452_independent.out
independent_script_sha256: 6ead91cfa71c5c31fd212d8fb12cf67fc834e968f74c08c72fe97e163b75197e
independent_output_sha256: c9e830740e0bf0123244bbf92c2b4bf8fc77bc6213a6195974843a70dd6f3f01
report: 05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452.md
audit: 05-knowledge/results/lrc14_clock_four_one_v2_component_thm4452_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4452 -- LRC14 clock-four one-v2 strict component-width caps

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The proof keeps
the four sheet colors. The corresponding uncolored scalar capacity inequality
is too weak even at the sharp odd-3-unit boundary. The infinite reduction is
analytic, and exact computation is confined to 84 and 30 stated rows.

Put
\[
 D_n=\{x\in\mathbb R/\mathbb Z:\|nx\|<1/14\}.              \tag{1}
\]
Let `r,a,b` be positive odd integers with `a` distinct from `b`, and define
the physical four-sheet tail failure set
\[
 P_4(r;a,b)=\{x:\text{each }x+j/4, 0\le j<4,
       \text{ lies in }D_{2r}\cup D_a\cup D_b\}.            \tag{2}
\]
Write `lambda_4(r;a,b)` for the longest actual open connected-component
width. For a set `E`, use `E-s={x:x+s in E}`.

## 1. Exact colored decomposition

An odd tail kills at most one quarter-lift. The even tail `2r` kills either
the pair of sheets `0,2`, the pair `1,3`, or neither. Since the total sheet
capacity is exactly `2+1+1=4`, covering all four sheets forces saturation:
`2r` owns one parity pair and `a,b` own the other two sheets separately.

Let
\[
 P_{ab}=(D_a\cup D_b)\cap((D_a\cup D_b)-1/2).              \tag{3}
\]
The saturation statement is the pointwise identity, including all strict
walls,
\[
 \boxed{P_4(r;a,b)=
 [D_{2r}\cap(P_{ab}-1/4)]\ \cup\
 [(D_{2r}-1/4)\cap P_{ab}].}                               \tag{4}
\]
The two displayed sets are disjoint. Moreover
\[
 P_{ab}=[D_a\cap(D_b-1/2)]\cup[D_b\cap(D_a-1/2)],          \tag{5}
\]
again disjoint, because an odd label cannot kill antipodal sheets.
Consequently each component in (4) lies inside one `D_(2r)` tooth and one
tooth belonging to the larger of `a,b`. Therefore
\[
 \boxed{\lambda_4(r;a,b)\le
 \min\left({1\over14r},{1\over7\max(a,b)}\right).}         \tag{6}
\]

This is a componentwise colored bound, not a Haar inequality. If the colors
are erased, the scalar three-comb capacity condition at length `1/110` is
already an equality for `(r,a,b)=(1,1,5)`, although the actual carrier is
empty. Thus the uncolored condition cannot prove the sharp theorem.

## 2. Sharp all-height constants

For equality or violation at `1/98`, (6) restricts
\[
 r\in\{1,3,5,7\},\qquad
 \{a,b\}\subset\{1,3,5,7,9,11,13\}.                       \tag{7}
\]
There are 84 unordered rows. Exact strict-wall evaluation gives
\[
 \boxed{\lambda_4(r;a,b)\le1/98},                           \tag{8}
\]
with equality, up to swapping `a,b`, precisely at
\[
 (r;\{a,b\})=(1;\{7,9\}),\ (3;\{7,13\}),\ (5;\{3,7\}).   \tag{9}
\]
The runner-up is `1/105`.

If `r,a,b` are also 3-units, equality or violation at `1/110` is restricted
to
\[
 r\in\{1,5,7\},\qquad
 \{a,b\}\subset\{1,5,7,11,13\}.                            \tag{10}
\]
The 30 unordered rows give
\[
 \boxed{\lambda_4(r;a,b)\le1/110},                         \tag{11}
\]
with equality uniquely at
\[
 (r;\{a,b\})=(5;\{1,11\}).                                \tag{12}
\]
The runner-up is `17/2002`. In fact the only nonempty rows in this box are
\[
\begin{array}{c|ccccc}
(r;\{a,b\})&(1;\{11,13\})&(5;\{1,11\})&(5;\{1,13\})
 &(7;\{1,11\})&(7;\{1,13\})\\ \hline
\lambda_4&17/2002&1/110&1/910&1/539&5/637.
\end{array}                                                 \tag{13}
\]
Every equality row in (9) and (12) has eight physical components of the
equality width.

## 3. Quotient normalization and LRC14 entry

The set (2) is invariant under quarter-turn translation. Its components are
strictly shorter than `1/4`, so the map `y=4x` identifies each four-component
orbit and multiplies component widths by exactly four. The sharp quotient
caps are therefore
\[
 {2\over49}\quad\text{and}\quad {2\over55}.                \tag{14}
\]

For a ten-body `C`, put
\[
 G_C=\{y:\|cy\|\ge1/14\text{ for every }c\in C\}.          \tag{15}
\]
Strict failure of the original row
\[
 4C\cup\{2r,a,b\}                                          \tag{16}
\]
is equivalent to containment of `G_C` in the quotient of (2). If `L_C` is
the longest positive-length connected component of `G_C`, (14) gives the
inclusive sufficient escape gates
\[
 \boxed{L_C\ge2/49\quad\text{(all odd)},\qquad
        L_C\ge2/55\quad\text{(odd 3-units)}.}              \tag{17}
\]
Inclusivity is valid: a closed interval cannot fit inside an open interval
of the same length. The live one-`v_2` clock-four residual uses the second
line of (17).

If `kappa_+(C)` is the number of positive-length components of `G_C`, then
\[
 \mu(G_C)\ge {2\over49}\kappa_+(C)
 \quad\text{or}\quad
 \mu(G_C)\ge {2\over55}\kappa_+(C)                         \tag{18}
\]
is sufficient in the respective domain. These are body component/count
gates, not universal lower bounds on `L_C` or `mu(G_C)`.

For an odd common dilation `d`, quarter sheets are permuted and
`lambda_4(dr;da,db)=lambda_4(r;a,b)/d`. Equality in the universal caps is
therefore attained only by the individual triples (9),(12), not their full
dilation rays.

## 4. Verification and scope

The primary program checks (4) against the original four-sheet predicate on
every defining wall and complementary cell and evaluates both finite boxes.
The independent program imports no author geometry: it reconstructs all
walls, the equality topology, the runner-ups, and quotient normalization.
Normal and optimized runs are line-identical to the frozen outputs and end
in `PASS`.

The distinction from THM-4451 is load-bearing: the even tail can own two
quarter sheets, so the clock-two two-distinct-owner inequality cannot simply
be reused. The colored saturation in (4) repairs that type mismatch. No
universal body component floor is proved, clock four remains open, and
LRC(14) remains **OPEN**.
