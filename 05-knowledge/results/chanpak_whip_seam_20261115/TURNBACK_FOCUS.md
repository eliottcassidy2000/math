# A literal turnback: exact gain, an energy subtraction, and matched-height resonance

The following return operation changes actual spine arcs and therefore escapes the ordinal binary observer. It is an analogy with wave focusing, not a mechanical claim about a whip. The randomly sampled [Farey resonance session](../../../07-reflections/cuts-as-farey-geodesics-resonance-and-the-hyperoctahedral-metagraph-opus-20260629.md) suggested checking exact divisibility, not visual closeness of two growing pulses. Here an integer resonance falls out of the count, rather than being a premise.

## Exact slot law for one returning tail vertex

Let R be a nonempty tournament, followed ordinally by a transitive tail t_1,…,t_L. Pick 1 ≤ j ≤ L and a nonempty marked set M ⊆ R. Replace s → t_j by t_j → s precisely for s ∈ M; no other contacts are changed. For m ≥ 1, let (Q_m) count ordered (m)-slot directed path covers (empty slots allowed), and define, writing (|P|) for vertex order,
```
U_m(M) = sum over nonempty R paths P starting in M of Q_m(R-P),
V_m(M) = same sum, with sign (-1)^|P|,
J_m(M) = 2 sum over even-order R paths starting and ending in M of Q_m(R-P).
```
The change from the ordinal baseline (Q_m(R)m^L) is exactly
```
 ΔQ_m=m^{L-j+1}((m+1)^{j-1}U_m(M)+(m-1)^{j-1}V_m(M))-m^L J_m(M).
```
Here (0^0=1). This is a literal simultaneous law at every slot number, not a formal estimate from its first evaluations.

Proof. Use the [odd-cycle packing expansion](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20), with two choices per cycle, then assign each cycle and remaining singleton to one of the (m) slots. A new crossing cycle must contain (t_j), and a packing contains at most one such cycle. Its form is a path (P) in (R) starting in (M), an increasing subsequence (Y) of earlier tail vertices, then (t_j), then the return to the start of (P). When (Y ≠ ∅) all forward arcs are available. For (Y=∅) the result is forbidden precisely when (P) ends in (M). Retain even (|P|+|Y|) so the cycle is odd. The unaffected residual tournament is ordinal and supplies (Q_m(R-P)m^{L-1-|Y|}); the new cycle supplies (2m). Summing with the even indicator gives the two binomial factors, and the forbidden empty subsequences contribute exactly (m^LJ_m). A two-vertex noncycle cannot enter this filter since it would have odd (|P|+|Y|).

For (m=1,j≥2), the formula becomes
```
 H=H(R)+2^{j-1}U_1(M)-J_1(M).
```
So a return past the first tail vertex necessarily raises (H): it creates a triangle (s	o t_i	o t_j → s) for every s ∈ M, (i<j). Equivalently its counted extra odd packings are strictly positive. This particular one-vertex turnback cannot close a window **at the old h**, though both arms can turn back to a new common (h).

## Exact resonance and a physical split of higher readouts

For two arms on a fixed (R), with masks (M,M') and depths (j,j'ge2), equality of (H) forces
```
 2^{j-1}U_1(M)-2^{j'-1}U_1(M')=J_1(M)-J_1(M').
```
Once (2^{min(j,j')-1}>|J_1(M)-J_1(M')|), it is equivalent to **equal return energy** (J_1) and equal doubled pulse (2^{j-1}U_1). In particular a fixed base admits only finitely many depth offsets at which two late arms can resonate, decided by the odd parts and 2-adic valuations of their positive (U_1)'s. It is not enough for two large gains to be approximately equal. For single-contact returns (J_1=0), and the same arithmetic criterion holds at every depth (j,j'ge2).

The higher readouts can still split rapidly through the same matched (H). Use the actual five-vertex core with rows ((16,21,17,7,8)), (H(R)=9,Q_2(R)=84), and masks 9 and 20 (indices from zero). Exact kernels give
```
M     U1 J1  U2  V2  J2
9     25  4 122 -34  20
20    25  4 111 -31  20
```
Place the returning vertex at the tip (L=j≥2). Both new spines have (H=5+25·2^{j-1}), but
```
 Q_2(9)-Q_2(20)=22·3^{j-1}-6.
```
This is a literal family with the first channel matched and a faster channel separated, not an assertion that their complete painted continuation states match or that the width-eight window closes. It illustrates what to measure in a new turnback: enforce the slow count resonance and then control, rather than discard, the faster slot channels. A single returning tip is itself a distinguished-vertex cut over the remaining ordinal spine; no claim of bypassing the known fixed-source scalar capacity estimates is made.

Replay: `python replay_turnback_focus.py` computes every displayed small kernel from directed path endpoints, then independently counts the turned-back tournaments by subset DP with ordered-slot convolution. It tests different depths, masks, and slot counts, not just the displayed pair.
