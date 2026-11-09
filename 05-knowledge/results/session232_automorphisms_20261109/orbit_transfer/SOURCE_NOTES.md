# Source pins and exact scope

Checked against the linked v1 HTML on 2026-11-09. Each paper-derived paragraph below is under 200 words. The finite results in PROOF.md have elementary proofs independent of these papers. Do not expand the source summaries elsewhere by duplicating additional paper-derived prose.

## Xinyu Shao, Cohomology of integral structure sheaves for open varieties

Source: https://arxiv.org/html/2610.31289v1

Pins: Theorem 1.2; section 2.2, after equation (2.1), weight decomposition and uniform-lag argument; Proposition 3.8 and its proof.

Shao works over a mixed-characteristic perfectoid field containing all roots of unity. Theorem 1.2 identifies the cohomology tower of quasi-compact opens exhausting the trivial-log locus with a constant tower, up to almost isomorphism. The section 2.2 calculation separates the persistent coefficient part from additional weights. For coefficients modulo p^m, lag max(m,2) kills every additional weight's cohomology: sufficiently large weights die already on coefficients; smaller weights are killed by a Koszul annihilator whose valuation the transition scalar dominates. Uniformity over all weights is essential. Proposition 3.8 obtains almost vanishing of the first derived inverse limit under finite-generation assumptions on successive images. Its proof replaces decreasing images by an almost stationary subsystem. These are geometric/module assertions, not finite combinatorial compactness theorems.

## Yuanning Zhang, The stacky approach to de Rham–Witt cohomology

Source: https://arxiv.org/html/2610.28149v1

Pins: Definitions 2.16 and 2.18; Proposition 2.20; Remark 2.27; Proposition 2.30, equation (2.31) and proof; Lemma 4.9.

Zhang's naive r-Nygaard stack uses the relation ut=p^r; the genuine stack is a chain of Nygaard charts with specified Frobenius gluing. Proposition 2.20 gives their contraction. Remark 2.27 compares taking even stacks before and after the relevant iterated fiber product. Proposition 2.30 computes pushforward as a fiber of chart sections mapping to adjacent-overlap sections, with the gluing maps retained; its proof uses adjacent intersections, absence of triple intersections, descent, and affine restriction of scalars. Lemma 4.9 obtains the Verschiebung product identities from trace projection formulas. The distinction is not an assertion that one stack is the ordinary coarse orbit set of the other. Our finite analogy keeps actual stabilizer/transport data rather than treating these geometric constructions as tournament constructions.

## Repository inheritance checks

Before deriving the finite recurrence, searched canon and recent results for pointwise stabilizers, prefix orbits, endpoint/subset quotients, weighted quotient DP, and genuine/naive quotient terminology. Closest inherited mechanisms are the paths' free action (LEM-003), all-order presentation gaps (THM-4812), SCC/OCF sheafification (THM-4892), and Session231's seven-coordinate carries and content filter. No literature-priority claim follows from this bounded repository search.
