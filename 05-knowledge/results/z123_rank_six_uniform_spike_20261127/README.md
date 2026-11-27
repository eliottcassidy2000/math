# A rank-six counterexample from the z123 flag calculation

The linked [z123 program](https://github.com/guoriyue/z123/blob/main/z123.py) tests a rank-seven matroid obtained by truncating the direct sum of **two** rank-four Gabidulin matroids. This note replaces one of those expensive summands by a rank-three uniform matroid. After one truncation the rank is **six**, and a violation survives even when the second summand is just three free elements. This is a new family relative to the linked construction; no priority or globally minimum rank is asserted.

## 1. Statement and scope

Write \(G_q=G_{4,q}\) for the rank-four Gabidulin matroid of [Bérczi–Coron–Eur–Li, §4.2, Definition 4.12, Proposition 4.13 and Corollary 4.14](https://arxiv.org/html/2610.17135#S4.SS2), with \(n=6,\ m=q^3,\ b=q^9\) in their notation. For prime powers \(q\) and integers \(N\ge3\) define
\[
M_{q,N}=\operatorname{Trunc}(G_q\oplus U_{3,N}),\qquad
Z_{M_{q,N}}(t)=1+z_1t+z_2t^2+z_3t^3+z_2t^4+z_1t^5+t^6.
\]
The rank decreases from seven to six. All constructions use actual matroids; none assumes an arbitrary flag vector is realizable. The results concern the ordinary **\(Z\)-polynomial**, not the Chow polynomial, the inverse \(Z\)-polynomial, or the coefficient sequence counting independent sets.

**Results.**
- At \((q,N)=(181,76)\) the index-two gap is
  \[
  z_2^2-z_1z_3=-3047994723191979213261875140252924091840840461600135485.
  \]
  In fact this index fails for every prime power \(q\ge181\) and every \(N\ge76\). For all field sizes \(2\le q\le179\) it passes at this index for every \(N\ge3\). At \(q=181\), it passes for \(3\le N<76\). These checks do not assert that a pass at one index implies log concavity.
- With the second summand **only three free elements**, the gap is negative for every prime power \(q\ge239\). At \(q=239\) the three exact coefficients are \(218938492224275943517703278,\ 2999379799972648156166575529,\ 41102872066332189025600618271\). This reduction is in rank and in the complexity of the second summand, not in the size of the whole ground set relative to z123.
- For each of these two fixed-\(N\) families throughout the stated failing range, the polynomial has exactly **two real roots**, both negative, and four nonreal roots. In particular these are counterexamples to the assertions in [Proudfoot–Xu–Young, Conjectures 5.1 and 6.5 (the latter already for the trivial group)](https://arxiv.org/html/1706.05575#S5). This is a stronger placement of the failure than a numerical root plot.

## 2. The three-coefficient separator

For flags of a matroid, write \(f_I\) for the number of strict chains with rank set \(I\), and \(W_i=f_{\{i\}}\). The same local formulas used by z123 can be derived from [palindromicity and the degree bound for KL polynomials](https://arxiv.org/html/1706.05575). For rank three or four the linear coefficient of the KL polynomial is \(W_{r-1}-W_1\), by comparing \(t\) with \(t^{r-1}\) in \(Z\). For rank five, comparing \(t^2\) and \(t^3\) gives
\[
[t^2]P=W_3-W_2+f_{24}-f_{23}-f_{14}+f_{12}.
\]
Indeed the low coefficient has the additional rank-one contraction term \(f_{14}-f_{12}\), while the high coefficient has the rank-two contraction term \(f_{24}-f_{23}\).

For rank six, identify the first three coefficients with the readings at degrees \(5,4,3\). Only the following flag ranks may contribute (all the flags here are those of \(M\)):
\[
\begin{aligned}
z_1&=W_5,\\
z_2&=W_4+f_{35}-f_{34},\\
z_3&=W_3+f_{25}-f_{23}+f_{14}-f_{13}+f_{135}-f_{134}-f_{125}+f_{123}.
\end{aligned}
\]
For a flag of \(M\) of these ranks, the truncation does not alter the corresponding flats of the sum. Split each flag rank between the two summands, monotonically on each side. Count the resulting pair of multichains by the product of the two strict-flag counts after deleting repeated ranks and the bottom/top. This is a literal sum over the data shared at the interface, as in the repo's [marked run-count readout](../chanpak_uniform_and_rooted_readout_equivalences_20261126/README.md); it is not naive multiplication of the polynomials **after** truncating.

For the particularly simple \(U_{3,3}\) case, the rank splits give the following short expressions, where now \(W_i,f_I\) refer to \(G_q\) alone:
\[
\begin{aligned}
z_1={}&W_2+3W_3+3,\\
z_2={}&3f_{13}+3f_{23}-3W_1+4W_2+5W_3+1,\\
z_3={}&-f_{12}+f_{123}+7f_{13}+5f_{23}-7W_1+7W_2+W_3+1.
\end{aligned}
\]
Using Corollary 4.14, their leading terms are \(\frac32q^{11}\), \(\frac{41}2q^{11}\), and \(q^{12}\). In particular the exact reason for the failure is a one-degree **flag spike**: the middle count inherits the \(q^{12}\) chain contribution \(f_{123}\); neither neighbor at indices one and two does. Thus the gap begins with \(-\frac32q^{23}\). Only a single Gabidulin block is needed to generate the spike. The second uniform block changes which rank levels can access it after truncation.

### Exact sign and root certificates

The replay computes the two-variable rational polynomial gap by Gaussian and ordinary binomials. In variables \(N-76,q-181\) all its nonzero coefficients are negative and its constant term is negative. For each \(2\le q\le179\), every coefficient after shifting to \(N-3\) is positive. At \(181\) it evaluates the remaining finite range \(3\le N<76\). For \(N=3\), the gap polynomial in \(q-239\) is coefficientwise negative. The mixed witness has \(195349420523\) ground elements. For comparison, replaying the linked rank-seven family gives a formal integer crossover at \(q=26\), so its first failing prime-power size at this index is 27.

For real roots write \(Z(t)=t^3 F(t+t^{-1})\) with
\[
F(x)=x^3+z_1x^2+(z_2-3)x+z_3-2z_1.
\]
The replay expands the two fixed-\(N\) discriminants at their respective thresholds and finds all coefficients strictly negative. For \(F(-2)\) it finds all coefficients strictly positive. Thus the unique real cubic root lies below \(-2\), giving two negative reciprocal roots of \(Z\); the two nonreal cubic roots give four nonreal roots of \(Z\). This is an exact real-variable certificate related to the [repo's earlier rank-five dual Chow root work](https://github.com/1337isnot1337/solved-conjectures/tree/main/results/caiolo-ferroni-hoster), but for a different invariant and a different mechanism.

## 3. How far the spike survives padding

Direct sum before a truncation is not a product rule for this calculation, but direct sum **without further truncation** does multiply \(Z\), as follows from its definition and KL multiplicativity. For \(M_{q,3}\), adjoin \(k\) coloops after truncation. The effect is \((1+t)^k Z(t)\). For each fixed \(k\) the one-degree spike still defeats log concavity at index two at all sufficiently large \(q\), so the mechanism works in every fixed rank at least six.

Even for \(k/\sqrt q\to u\ge0\), the new index-two gap divided by \(q^{23}\) tends to
\[
\frac98u^2-\frac32 .
\]
Indeed the new first three coefficients are \(z_1+k,\ z_2+kz_1+\binom{k}{2},\ z_3+kz_2+\binom{k}{2}z_1+\binom{k}{3}\); substitute the three leading terms above. Thus the index-two defect survives coloop growth below \(k/\sqrt q=\sqrt{4/3}\), and the leading sign reverses above it (nothing claimed at equality or about other indices). This isolates a quantitative dilution threshold, rather than just another example. It parallels the repo's [exact padding versus target distinction](../solved_conjectures_survey_and_simplex_gap_theorem_20261127/README.md): neutral padding keeps one chosen observation but can change what a neighboring coefficient reads.

Run `python3 replay.py` (SymPy, exact arithmetic) for the symbolic flag equations, threshold signs, and cubic certificates. The new target left open here is rank five or smaller with a different geometry; rank-six minimality is not claimed.
