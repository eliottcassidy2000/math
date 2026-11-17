# An optimal cube of lower-only Hadamard switches

In the [two-tiling chart](../third_cover_binary_tilings_20261117/README.md), a diagonal-normalized Hadamard matrix gives two tournaments under a common directed frame: upper signs for the first and opposite transposed lower signs for the second. Earlier localized switches changed both readouts. Here a single fixed upper tiling admits **2^(n/4−1) different lower tilings** at every dyadic order n≥8. The count refers to literal tilings, not isomorphism classes or distinct path totals. Each independent switch uses the minimum possible number n of entries. It is a Walsh-character packing of the closed-quadruple switches of [Orrick](https://arxiv.org/abs/math/0507515).

## Construction and optimality inside this seam

Let V=F₂^m, n=2^m, m≥3. Start with M_xy=(−1)^(x·y+y·y), with diagonal +1. Use the row subgroup R=span(e₀,e₁+e₂) and the two sectors distinguished by l(x)=x₁+x₂. Place sector 1 first, sector 0 last, alternating Hamming parity odd/even within each sector. Both sectors are balanced because adding e₀ preserves the sector and reverses parity. Hence joint row/column sign changes can make every upper adjacent sign positive and lower adjacent sign negative: before the changes the two signs multiply to (−1)^wt(x+y)=−1.

Take each row coset of R contained in sector 0, and each of the two column fields in sector 1 specified by a fixed value of y₀. Each field is a coset of R^⊥. There are n/8 row cosets and two fields, making n/4 disjoint rectangles, all strictly below the diagonal. Omit the unique rectangle containing the first-subdiagonal cell at the sector boundary. Every other rectangle lies strictly below the subdiagonal, and its sign flip leaves the upper half, the diagonal, and the shared directed frame unchanged.

Choose the flips arbitrarily. Row pairs in the same R coset get identical signs on each field and keep their inner product. For rows from different cosets, each changed field contributes zero to their inner product, because x+x' lies outside R and its character sums to zero on each coset of R^⊥. Rows left entirely unchanged count as choosing no flips. Thus every bit pattern remains Hadamard, with independent disjoint supports of 4·(n/4)=n entries.

For optimality in this half-by-half seam, any two distinct Hadamards differ in at least n entries. Indeed, if k rows change and one changes in d entries, its difference vector is orthogonal to all unchanged old rows and has inner product of magnitude at most 2d with each of the remaining k old rows. Parseval gives 4nd≤4kd², so d≥n/k. Summing gives the bound n. Only n²/4−1 cross-seam entries are available after preserving the frame, so no family of pairwise disjoint nonzero Hadamard switches supported in that seam can have more than n/4−1 bits. This is not a bound on non-disjoint switches in all of the lower triangle.

## Exact readings and where the apparent independence ends

At n=8 take the frame (2,3,4,5,1,0,7,6) and switch rows {0,1,6,7} against columns {2,4} (unpermuted labels). The upper count stays 233. The lower count moves from prime 233 to 249=3·83. In the two-core count H=ad+be−cf this particular switch cannot preserve H with a fixed other profile even allowing either switched port to move. For any old/new port pair of the switched arm, let Δb and Δc denote changes in the cover and start counts. Divide the seam difference by the other arm's weighted cover count B. The result is Δb+16 A/B−Δc C/B with 0≤C≤A≤B. Its three possible extreme values on this triangle are Δb,Δb+16,Δb+16−Δc. The complete 64-pair replay verifies that all three have the same sign and each has magnitude at least 6, giving gap ≥6B. This is not a claim of obstruction for every frame, order or coordinated switch.

At n=16 the canonical three-bit cube has eight distinct lower counts with the same upper tiling, so these eight lower tournaments are even distinct up to converse. In mask order:

    612097097, 612081469, 612051653, 612047197,
    606270393, 606264149, 606257381, 606246749.

The mask-1 and mask-2 mixed difference is 11172. Thus independent matrix bits do not become additive count digits. This is still useful structural freedom: one can alter the lower core without moving the upper at every dyadic order, but there is no theorem that an arbitrary desired seam count can be obtained by choosing the bits. Unlike a bounded patch, the minimum move itself grows proportionally to n.

Run `python replay.py` for direct Gram checks on all masks at 8,16,32, the literal order-eight profiles and inequalities, and the order-sixteen counts by independent subset DP. No numerical computation is used in the all-order orthogonality or optimality proofs.
