---
id: THM-4450
title: "LRC14 absorbed-label overlap hierarchy and component-address decoder"
status: >
  PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY AUDITED. The overlap of one absorbed label with
  a ten-body has a sharp four-tier tenth-order lower bound according to
  gcd(r,6), improving every clock-two one-even deletion loss. For clock four
  a mixed-radius order statistic gives a global hybrid floor. After inherited
  five-ray localization, strict failure is exactly a finite complete-component
  address condition; endpoint-only and single-component criteria are false.
  No residual ray is closed and LRC(14) remains open.
source: dyadic residual overlap continuation + two clean-room audits, 2026-09-06
depends_on:
  - THM-739-pairwise-coprime-bad-overlap-exact-bernoulli-closed-form
  - LEM-042-pair-overlap-law
  - THM-2060-crt-tail-coset-saturation
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
  - THM-4449-lrc14-dyadic-seventh-rounding-energy-and-residual-haar-entry
  - LRCUpTo13
related:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction
primary_script: 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450.py
primary_output: 05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450.out
primary_script_sha256: 28e95f61464adcdf127d5b1d3e3d56949f871e644835e677c56b6326ab9e58fa
primary_output_sha256: a0555be8fe4db9122cc2267f4a14b2de9191eee00ccdd4fa49b32f29b3e995db
independent_script: 04-computation/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.py
independent_output: 05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_cleanroom.out
independent_script_sha256: 832d3012823d374046964affd3028275bd125e259e4f8183ed769255856079e3
independent_output_sha256: 40bd7e00e62e0e2a87c591116e217d95a8e2087373202494f1bc6a4fc1743f98
report: 05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450.md
audit: 05-knowledge/results/lrc14_absorbed_label_overlap_hierarchy_thm4450_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4450 -- LRC14 absorbed-label overlap hierarchy and component-address decoder

**PROVED ELEMENTARY + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** The Bernoulli overlap identity below
is recovered from THM-739 (the repository has a legacy THM-739 ID collision,
so the full slug is essential) and LEM-042. New here are its oriented
denominator classification by `gcd(r,6)`, the sharp tenth order statistic,
the structured-body consequences, and the exact complete-component decoder.
The five primitive dyadic rays and their widths are inherited from THM-4153
and THM-4449, not reclaimed as new. No literature-priority claim is made.

Write
\[
 D_s=\{y\in\mathbb R/\mathbb Z:\|sy\|<1/14\},\qquad
 G_A=(\mathbb R/\mathbb Z)\setminus\bigcup_{a\in A}D_a.       \tag{1}
\]

## 1. Sharp tenth-order overlap hierarchy

For distinct positive integers `c,r`, reduce the **oriented** ratio
\[
 {c\over r}={p\over q},\qquad \gcd(p,q)=1.                    \tag{2}
\]
Then
\[
 \mu(D_c\cap D_r)={1\over49}
 +{B_2(\{(q-p)/14\})-B_2(\{(q+p)/14\})\over pq}.             \tag{3}
\]
The measure is symmetric in `c,r`; the orientation is retained because
`q` divides `r`, so the reduced denominator remembers arithmetic of the
absorbed label.

The radius-`1/14` tooth has Fourier coefficients `1/7` at zero and
`sin(pi k/7)/(pi k)` otherwise. Orthogonality at frequencies `p,q`,
product-to-sum, and the Bernoulli cosine series give (3), recovering the
cited overlap law. Since the oscillation of `B_2` is `1/4`,
\[
 \mu(D_c\cap D_r)\ge {1\over49}-{1\over4pq}.                 \tag{4}
\]

Let `C` contain ten distinct positive labels and suppose `r` is not in `C`.
Then
\[
 \boxed{\max_{c\in C}\mu(D_c\cap D_r)\ge\lambda_g,\qquad
        g=\gcd(r,6),}                                          \tag{5}
\]
where
\[
\begin{array}{c|cccc}
g&1&2&3&6\\ \hline
\lambda_g&1/63&1/70&1/70&1/77.
\end{array}                                                    \tag{6}
\]
Each constant is sharp uniformly over its residue class, already for a
primitive ten-set `C`.

### Exhaustive all-height proof

Because `q` divides `r`, the allowed denominator conditions are
\[
\begin{array}{c|c}
g&\text{allowed }q\\ \hline
1&\gcd(q,6)=1\\
2&3\nmid q\\
3&2\nmid q\\
6&\text{unrestricted}.
\end{array}                                                    \tag{7}
\]
Equation (4) makes the exact searches finite: products beyond `55`, `40`,
`40`, and `33`, respectively, are strictly above the target. Exact
substitution below those cutoffs gives the complete subcritical oriented
ratio lists
\[
\begin{array}{c|l}
g=1&1/13,1/11,2/11,3/11,10,11,12,13\\
g=2&1/13,1/11,2/11,3/11,11/2,11,12,13\\
g=3&1/13,1/11,2/11,3/11,11/3,11,12,13\\
g=6&1/13,1/12,12,13.
\end{array}                                                    \tag{8}
\]
There are fewer than ten in every row. The ten ratios `c/r` are distinct, so
at least one reaches the target, proving (5). The equality ratios are
\[
\begin{array}{c|l}
g=1&9/5,9,27\\
g=2&1/10,3/10,10\\
g=3&10/3,10\\
g=6&1/11,2/11,3/11,11/3,11/2,11.
\end{array}                                                    \tag{9}
\]
The audit supplies explicit primitive ten-sets containing all subcritical
ratios and enough boundary ratios so their maximum is exactly `lambda_g`.
Sharpness here means sharpness of the tenth pair correlation; it does not
assert equality in the union bound used next.

## 2. Clock-two absorbed-body floors

Let `H=C union {r}`. Exact inclusion-exclusion gives
\[
 \mu(G_H)=\mu(G_C)-{1\over7}
      +\mu\!\left(D_r\cap\bigcup_{c\in C}D_c\right).       \tag{10}
\]
Using (5),
\[
\boxed{
\begin{array}{c|ccc}
\gcd(r,6)&1&2\text{ or }3&6\\ \hline
\mu(G_H)-\mu(G_C)&\ge-8/63&\ge-9/70&\ge-10/77.
\end{array}}                                                  \tag{11}
\]
These improve the old uniform subtraction `-1/7` for every absorbed label.

Failure of `2H union {a,b}` places compact `G_H` inside the proper open
two-tail cross-comb. Combining (11) with the general odd-pair cap `4/63`,
the sufficient original-body entries are
\[
\begin{array}{c|ccc}
\gcd(r,6)&1&2\text{ or }3&6\\ \hline
\mu(G_C)\text{ entry}&4/21&121/630&134/693.
\end{array}                                                    \tag{12}
\]
For an odd-3-unit pair, cap `4/77` sharpens these to
\[
 124/693,\qquad139/770,\qquad2/11.                           \tag{13}
\]
At the inherited five-ray localization level `4/91`, they become
\[
 20/117,\qquad157/910,\qquad174/1001.                        \tag{14}
\]
The live residual absorbed label is an odd 3-unit, so `124/693` and `20/117`
replace `15/77` and `17/91`. These are sufficient gates; no universal
ten-body floor is asserted.

## 3. A clock-four hybrid floor

Now let `r` be an odd 3-unit and `H=2C union {r}`. Disintegrate over `y=2x`
and put
\[
 E_r=\{y:\|ry\|<1/7\}.
\]
Both lifts are safe for `2C` when `y` lies in `G_C`, and at most one is
`r`-dangerous. Hence the exact fibre identity is
\[
 \mu(G_H)=\mu(G_C)-{1\over2}\mu(G_C\cap E_r).              \tag{15}
\]
The recovered mixed-radius Bernoulli formula is
\[
 \mu(D_c\cap E_r)={2\over49}
 +{B_2(\{(q-2p)/14\})-B_2(\{(q+2p)/14\})\over pq}.         \tag{16}
\]
Its tenth order statistic is `1/28`: only seven ratios are subcritical,
\[
 1/25,1/13,1/11,2/11,5,6,13,                                \tag{17}
\]
and equality occurs at `4/5,4,12,20`. A primitive equality control is in the
audit. Consequently
\[
 \boxed{\mu(G_{2C\cup\{r\}})\ge
     \max\bigl(\mu(G_C)/2,\ \mu(G_C)-1/8\bigr).}          \tag{18}
\]
In the live residual range `mu(G_C)<1/4`, the half-mass term dominates, so
this global improvement does not lower the clock-four entry gate.

## 4. Exact complete-component decoder

Assume inherited localization has left a primitive odd-3-unit ratio
\[
 (p,q)\in\{(1,11),(1,23),(5,11),(1,37),(1,25)\}.            \tag{19}
\]
Let `(A_k,B_k)` be every primitive open quotient-`y` cross-comb component. At common
scale `t`, the complete pulled cells are
\[
 U_{n,k,t}=((n+A_k)/t,(n+B_k)/t),\qquad0\le n<t.             \tag{20}
\]
Let `J_1,...,J_m` be **all** closed connected components of `G_H`, retaining
zero-length isolated safe points. Strict failure at `(p,q),t` is equivalent
to
\[
 \boxed{\forall i\ \exists(n,k):\ J_i\subset U_{n,k,t}.} \tag{21}
\]
Indeed, (20) are exactly the disjoint open components of the pulled
cross-comb. A connected closed subset of their union lies in one component,
and conversely (21) covers all of `G_H`. Open containment retains endpoint
margins that an a.e. interval merge destroys. The physical-`x` pullback has
twice as many cells, each of half the quotient width.

A safe pulled-cell endpoint immediately refutes (21), but the converse is
false. Four exact structured bodies, two of each clock provenance, give
twenty ray checks in which every ray has a trapped component somewhere
despite a direct full-row witness. Across them, endpoint-only and
single-component tests both fail. They are method hostiles, not LRC
counterexamples.

If `L_H` is the longest positive-length component, inherited THM-4153 gives
`tL_H<beta_(p,q)`. Cited LRC through eleven speeds gives
`L_H>=1/(42 max H)`, so a fixed body has explicit finite scale caps
\[
\begin{array}{c|ccccc}
(p,q)&(1,11)&(1,23)&(5,11)&(1,37)&(1,25)\\ \hline
t/(\max H)&<12/11&<12/23&<54/55&<12/37&<12/25.
\end{array}                                                    \tag{22}
\]
The exact remaining state is the full address vector in (21), together with
THM-2061's parity and determinant sidecars. No uniform body bound or complete
ray escape is proved; the residual clocks and LRC(14) remain **OPEN**.

## 5. Verification

The primary audit checks (3) and (16) against literal wall integration on
3,309 coprime pairs and tests twenty component banks. A clean-room referee
independently rebuilds both overlap atlases, all five quotient cell lists,
scale pullbacks, and every hostile address vector; an orthogonal arithmetic
audit separately checks the transfers. All checks remain active under
optimized Python. Reproduction commands and raw-LF hashes are in the result
note and audit.
