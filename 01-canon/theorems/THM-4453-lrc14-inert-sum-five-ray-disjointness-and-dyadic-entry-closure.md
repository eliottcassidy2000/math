---
id: THM-4453
title: "LRC14 inert-sum five-ray disjointness and dyadic entry closure"
status: >
  PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY LOGIC-AUDITED. For a
  doubled body with safe mass at least 4/91, a distinct odd-3-unit tail pair
  whose primitive sum satisfies THM-3818's inert-prime/exponent condition
  cannot fail: the five high-overlap rays are disjoint from the 5,855-ray
  inert atlas. This yields inclusive original-body gates 8/91 for the
  one-v2 clock-four form and 20/117 for the clock-two one-even form, within
  the inert-pair class. It does not close arbitrary pairs or LRC(14).
source: THM-3818/4153/4450 composition + incoming unit-entry audit, 2026-09-06
depends_on:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
  - THM-4450-lrc14-absorbed-label-overlap-hierarchy-and-component-address-decoder
related:
  - THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry
  - THM-4452-lrc14-clock-four-one-v2-strict-component-width-caps
primary_script: 04-computation/lrc14_inert_sum_five_ray_disjointness_thm4453.py
primary_output: 05-knowledge/results/lrc14_inert_sum_five_ray_disjointness_thm4453.out
primary_script_sha256: e456fb062a4e13ed29b7c68c5b926107d4524f6903a0ac1c70b95c6755bfa74d
primary_output_sha256: da67b3db4d130994ff7eb643b2ee595ea23950d630a3d3e3b219502800314c8d
report: 05-knowledge/results/lrc14_inert_sum_five_ray_disjointness_thm4453.md
audit: 05-knowledge/results/lrc14_inert_sum_five_ray_disjointness_thm4453_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4453 -- LRC14 inert-sum five-ray disjointness and dyadic entry closure

**PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY LOGIC-AUDITED.** This
theorem composes a Haar ratio filter with an arithmetic decoder atlas. It
closes the overlap of those two branches; it does not assert that every
LRC(14) residual pair is inert-admissible.

For a finite positive body `H`, let
\[
 G_H=\{y\in\mathbb R/\mathbb Z:\|hy\|\ge1/14
       \text{ for every }h\in H\}.                         \tag{1}
\]
Let `a,b` be distinct positive integers coprime to six. Write
\[
 s=\gcd(a,b),\qquad a=sp,quad b=sq,quad
 1\le p<q,quad\gcd(p,q)=1.                                \tag{2}
\]
Assume every prime divisor `ell` of `p+q` satisfies
\[
 \ell\equiv2\pmod3,qquad v_\ell(p+q)\le2.                 \tag{3}
\]

## 1. Inert-pair entry theorem

If
\[
 \mu(G_H)\ge {4\over91},                                  \tag{4}
\]
then the row
\[
 2H\cup\{a,b\}                                             \tag{5}
\]
has a common phase of clearance at least `1/14`.

Indeed, under `y=2x`, strict failure of (5) puts the compact set `G_H` inside
the proper open pulled two-tail cross-comb for `(p,q)`. Multiplication by the
odd common scale `s` preserves Haar measure. Compact-to-proper-open
containment makes the mass inequality strict, so failure and (4) force the
primitive cross-comb mass to exceed `4/91`.

THM-4153's exact all-height ratio classification, restricted to 3-units,
then leaves precisely
\[
 (1,11),(1,23),(5,11),(1,37),(1,25).                       \tag{6}
\]
Their sums and obstructions to (3) are
\[
\begin{array}{c|c|c}
(p,q)&p+q&\text{obstruction}\\ \hline
(1,11)&12=2^2\cdot3&3\not\equiv2\pmod3\\
(1,23)&24=2^3\cdot3&3\text{ occurs and }v_2=3\\
(5,11)&16=2^4&v_2=4\\
(1,37)&38=2\cdot19&19\equiv1\pmod3\\
(1,25)&26=2\cdot13&13\equiv1\pmod3.
\end{array}                                                 \tag{7}
\]
Thus (3) and (6) are disjoint, a contradiction. The boundary in (4) is
inclusive because failure forced a strict mass inequality.

The coprimality-to-six hypothesis is essential. In the larger all-odd lane,
`(1,9)` survives the filter and its sum `10=2*5` satisfies (3).

## 2. Actual THM-3818 overlap

THM-3818's primitive pair atlas consists of exactly 5,855 ratios with
`p+q<=356` satisfying (3). Hence an actual `11+2` decoder equality row in the
odd-3-unit residual branch is closed by (4) whenever its eleven-coordinate
component is the physical `2H` and its actual two-vertex decoder component is
`{a,b}`. The pair must be that labelled decoder component; the existence of
an unrelated inert-admissible pair elsewhere in the row is insufficient.

At this threshold no normalized-unit, `W=V_dec`, or finite-box hypothesis is
needed. The incoming unit-core theorem remains complementary below (4); it
does not supply a universal body-component lower bound.

## 3. Original-body dyadic gates

For a ten-body `C` and odd 3-units `r,a,b`, consider first
\[
 4C\cup\{2r,a,b\}=2H\cup\{a,b\},\qquad H=2C\cup\{r\}.      \tag{8}
\]
If the reduced pair `(a,b)` satisfies (3), THM-4450 gives
\[
 \mu(G_H)\ge\max\bigl(\mu(G_C)/2,\mu(G_C)-1/8\bigr).       \tag{9}
\]
Therefore the inclusive conditional gate
\[
 \boxed{\mu(G_C)\ge8/91}                                  \tag{10}
\]
forces (4) and closes (8). The alternative term in (9) alone would require
`123/728`. For arbitrary odd-3-unit pairs, (10) only localizes to (6); the
inert condition is what converts localization into safety.

For the clock-two one-even form
\[
 2C\cup\{2r,a,b\}=2H\cup\{a,b\},\qquad H=C\cup\{r\},      \tag{11}
\]
THM-4450 gives
\[
 \mu(G_H)\ge\mu(G_C)-8/63.                                \tag{12}
\]
Thus, again within the inert-pair class,
\[
 \boxed{\mu(G_C)\ge {4\over91}+{8\over63}={20\over117}}  \tag{13}
\]
is an inclusive sufficient original-body gate. It improves `124/693` only
under the added arithmetic hypothesis (3).

There is no corresponding q=2 zero-even absorption: `2C union {a,b,c}` is a
natural `10+3` split, and THM-4451 remains its correctly typed component tool.

## 4. Incoming unit condition and scope

In (8), let `d=gcd(H)`. Since `r` is odd, the normalized body `H/d` contains
the literal coordinate one exactly when
\[
 r\mid c\quad\text{for every }c\in C.                      \tag{14}
\]
This is the cheap test for the incoming actual-unit-core theorem below the
mass gate. Its safe component has scale-dependent width `1/(42 max H)`, not
a constant floor that can be substituted into THM-4451 or THM-4452.

The primary audit independently enumerates all 5,855 atlas ratios, verifies
the five exclusions and the `(1,9)` hostile, checks (9)--(14), tests 45,045
unit/divisibility controls, and verifies the even-tail Boolean identity.
Normal and optimized runs are line-identical to the frozen output and end in
`PASS`. An independent quantifier audit replayed the inherited statements and
confirmed the labelled-component scope. LRC(14) remains **OPEN**.
