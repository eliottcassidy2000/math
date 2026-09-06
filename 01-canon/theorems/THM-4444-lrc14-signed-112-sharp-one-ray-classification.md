---
id: THM-4444
title: "LRC14 signed (1,1,2) sharp one-ray classification"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
  AUDITED. Every primitive sorted distinct positive ternary-unit triple with
  a signed coefficient-magnitude (1,1,2) relation has min_i E_i equal to its
  physical failure mass and at most 11/140, sharply. Exactly (2,11,20) is
  above 6/77 and exactly (1,5,11) equals 6/77. Entry, the additive family,
  synchronization, and LRC(14) remain open.
source: root low-circuit continuation + independent literal referee, 2026-09-06
depends_on:
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
related:
  - THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits
  - THM-4441-lrc14-signed-122-sharp-ray-closure
primary_script: 04-computation/lrc14_signed_112_sharp_classification_thm4444.py
primary_output: 05-knowledge/results/lrc14_signed_112_sharp_classification_thm4444.out
primary_script_sha256: 304babb2c2a19f5f8bdfb617f9284808ac8676e0a5edc7ea2f0583a3bc47779f
primary_output_sha256: 15323d9c04863b716d88707ce98e9a4f69d8255c78c2c98ffd3c04df4f598758
independent_script: 04-computation/lrc14_signed_112_sharp_classification_thm4444_independent.py
independent_output: 05-knowledge/results/lrc14_signed_112_sharp_classification_thm4444_independent.out
independent_script_sha256: 8efca1b7d107c7317e9e731048513b9fa41fb5da83873b2b5989328ebb8eceb3
independent_output_sha256: 7c6962e5eecd476af8f378ceaf8458e1711e2d0b9aca878aa31f356c6c235074
report: 05-knowledge/results/lrc14_signed_112_sharp_classification_thm4444.md
audit: 05-knowledge/results/lrc14_signed_112_sharp_classification_thm4444_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4444 -- LRC14 signed (1,1,2) sharp one-ray classification

**PROVED ELEMENTARY RELATIVE TO THM-4414 + FINITE-EXACT + INDEPENDENTLY
AUDITED.** This classifies the second residual circuit from THM-4437.
It is a local theorem; `LRC(14)` remains **OPEN**.

## 1. Sharp statement

Let \(w=(a,b,c)\) be primitive, with
\[
 1\le a<b<c,\qquad 3\nmid abc,
\]
and suppose it has a primitive signed relation with coefficient magnitudes
\((1,1,2)\). Use THM-4414's complete carrier set \(\Lambda(w)\), carrier
terms \(e_i(C)\), network projections \(E_i=\sum_Ce_i(C)\), and physical
failure mass \(\mu(F_w)=\sum_C\min_i e_i(C)\).

If \(r(w)\) is the coordinate carrying coefficient two, then exactly
\[
             \min_iE_i(w)=\mu(F_w)=E_{r(w)}(w).                    \tag{1}
\]
The sharp all-height bound is
\[
             \min_iE_i(w)=\mu(F_w)\le {11\over140},               \tag{2}
\]
with equality only at \(w=(2,11,20)\). Relative to the scale-three consumer
threshold,
\[
\begin{aligned}
 \min_iE_i(w)>6/77&\quad\Longleftrightarrow\quad w=(2,11,20),\\
 \min_iE_i(w)=6/77&\quad\Longleftrightarrow\quad w=(1,5,11).
\end{aligned}                                                     \tag{3}
\]

The three disjoint sorted sign cones are:

| family | relation / primitive ray \(u\) | selected coordinate | sharp value | equality |
|---|---|---:|---:|---|
| F1 | \(c=2a+b\), \((2,1,-1)\) | \(a\) | \(58/833\) | \((5,7,17)\) |
| F2 | \(c=2b-a\), \((1,-2,1)\) | \(b\) | \(11/140\) | \((2,11,20)\) |
| F3 | \(c=a+2b\), \((1,2,-1)\) | \(b\) | \(6/77\) | \((1,5,11)\) |

## 2. Sign classification and complete carrier ray

The twelve signed placements modulo global reversal reduce, under positivity
and \(a<b<c\), to
\[
 2a+b=c,\qquad a+c=2b,\qquad a+2b=c.
\]
Their arithmetic parameterizations impose only coprimality and the necessary
nonzero residue classes modulo three. Pairwise cone intersections force
\(b=3a\), \(a=b\), or \(a=0\), so no intersection survives in the typed
universe.

Let \(L=\ker_{\mathbb Z}(w)\), \(r=3/14\), and let \(u\) be the corresponding
primitive relation. THM-4414 represents every carrier as
\[
 C=-w\times e,\qquad |e_i|<r.
\]
Since \(w\) is primitive, \(w\times(-):\mathbb Z^3\to L\) is onto. Choose
\(\eta\in\mathbb Z^3\) with \(w\times\eta=u\), and put
\(\phi(C)=\eta\cdot C\). Its kernel on \(L\) is \(\mathbb Zu\), while
\[
             |\phi(C)|=|u\cdot e|<r\|u\|_1={6\over7}<1.
\]
The integral value \(\phi(C)\) is therefore zero. Thus every carrier lies on
the one ray, and conversely every permitted multiple is live:
\[
 \Lambda(w)=\{ku:k\in\mathbb Z,\ 3\nmid k,
 0<|k|<\min_i{3(w_j+w_k)\over14|u_i|}\}.                            \tag{4}
\]
No live address lies on a roof: a roof equality has right side divisible by
three and left side nonzero modulo three.

## 3. Pointwise selector and universal bulk

Normalize by \(z_i=w_i/c\). On the positive ray define
\[
 f_i(x)=\min\!\left({3\over7},
 {\,{3\over14}(z_j+z_k)-|u_i|x\over z_jz_k}\right),
\]
with zero beyond its strict support. Direct comparison of the capped affine
pieces in each cone shows that the coefficient-two profile is pointwise
least. This proves (1) before summation.

Exact cap-and-trapezoid integration gives the same selected and physical
positive-half integral in all three cones:
\[
                 \int f_{r(w)}=\int\min_i f_i
                 =r^2={9\over196}.                                \tag{5}
\]
The other coordinate integrals retain the cone parameter, but (5) loses
neither the physical overlap nor the selected owner.

For \(R_<(T)=\#\{1\le k<T:3\nmid k\}\), exact three-term blocks give
\[
 -{2\over3}\le R_<(T)-{2T\over3}<{2\over3}.
\]
Layer cake applied to the profiles yields the two-sided deleted-third band
\[
 {3\over49}-{4\over7c}\le
 \min_iE_i(w)=\mu(F_w)
 <{3\over49}+{4\over7c}.                                          \tag{6}
\]
In particular every row of height \(c\ge35\) is strictly below \(6/77\).
Thus the hostile and boundary packets in (3) are finite residue effects, not
a persistent continuum obstruction.

## 4. Exact head and independent replay

The upper side of (6) beats the F1/F2/F3 sharp leaders from heights
\(68,33,35\), respectively. A raw head through \(c=68\) has 361 rows, split
\(121/179/61\), and proves all constants in the table. The packets are
\[
\begin{array}{c|c}
(5,7,17)&(58/833,\ 12/119,\ 346/4165)\\
(2,11,20)&(131/1540,\ 11/140,\ 3/35)\\
(1,5,11)&(6/77,\ 6/77,\ 6/77).
\end{array}
\]
A diagnostic ray extension through \(c=611\) has 28,438 rows, split
\(9,482/14,214/4,742\), with the same leaders and exactly the loci (3).

The independent program imports no candidate code. It pins the earlier
literal six-sheet engine, independently constructs all twelve signed
directions, and checks all 28,438 rows against raw carrier enumeration and
a generic breakpoint integrator. It reproduces the ray, selector identity,
continuum integral, family leaders, and the complete \(6/77\) locus in
7,484,368 literal-engine checks plus 199,075 new gates.

## 5. Consequence and scope

Together, THM-4437, THM-4441, and this theorem show that outside the additive
\((1,1,1)\) family the only strict local \(6/77\) failure is
\((2,11,20)\), and the only boundary packet is \((1,5,11)\). This is the
strongest scalar local classification; it does not show that an incoming body
avoids either physical failure set.

The additive family, arbitrary chart entry, synchronization, and LRC(14)
remain open. The bounded-body theorem THM-4442 closes these packets only
after entry into its stated ten-body chart.

## 6. Reproduction

```powershell
python -B 04-computation/lrc14_signed_112_sharp_classification_thm4444.py
python -B -O 04-computation/lrc14_signed_112_sharp_classification_thm4444.py
python -B 04-computation/lrc14_signed_112_sharp_classification_thm4444_independent.py
python -B -O 04-computation/lrc14_signed_112_sharp_classification_thm4444_independent.py
```

All theorem decisions are exact integer or rational comparisons.
