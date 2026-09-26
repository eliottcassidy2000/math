---
id: THM-4503
title: "Two-chain Collatz parity posets, exact swap carry, and the obstruction and repair for arithmetic height selection"
status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED
source: crossroads-poset-20260926
depends_on: []
related:
  - 05-knowledge/results/crossroads_poset_20260926_bridge.md
  - 05-knowledge/results/crossroads_poset_20260926_width.md
scripts:
  - 04-computation/experiments/crossroads_poset_20260926_bridge.py
  - 04-computation/experiments/crossroads_poset_20260926_width.py
outputs:
  - 05-knowledge/results/crossroads_poset_20260926_bridge.out
  - 05-knowledge/results/crossroads_poset_20260926_width.out
audit: "Root and automata independently derived the two-chain bridge; geometry audited carrier, scope, and height-law obstruction."
---

# THM-4503 — An exact poset bridge with its arithmetic sidecar

Fix counts a,b of odd and even letters, ell=a+b, and require every nonempty
prefix with i odd and j even letters to satisfy 3^i>2^(i+j). Set
`m_j=min{i:3^i>2^(i+j)}`. For b>0,a>=m_b these words are exactly the linear
extensions of the width-at-most-two poset

    A1<...<Aa, B1<...<Bb, A_(m_j)<B_j (1<=j<=b).

For b=0 there is one all-odd word; for b>0,a<m_b there are no words, not
a poset with no extensions. This is a multiplicative prefix condition,
sufficient but not necessary for actual no descent of a specified source.

For a word w, with indices j starting at zero, its carry and source are

    C_w=sum_(j:w_j=1) 2^j 3^(number of ones after j),
    T_w(n)=(3^a n+C_w)/2^ell,
    r_w=-C_w*3^(-a) mod 2^ell.

Swapping adjacent 10 to 01 at j, with s odd letters after the pair, changes
the carry by 2^j3^s and the source residue by -2^j3^(s-a). In particular the
residue displacement has exact 2-adic valuation j. Rearrangement preserves
the endpoint multiplier but changes the integer being iterated.

For ell=6,a=5 the four extensions 110111,111011,111101,111110 have least
residues 27,39,47,31. Restricting sources to at most 31 selects the first
and last orders, which are not the full extension set of any poset on
these same six labels. The resulting order-statistic law also reverses
XYZ covariance: with Z=7 times the unit order coordinates,

    Cov(Z_A2-Z_A3,Z_B1-Z_A2)=1/4

instead of -5/32 for the full uniform law. Its support is nonconvex.
These are obstructions to a transfer, not counterexamples to a poset theorem.

The repair in the many-period regime is exact. For m distinct selected
residues modulo N, write X=qN+s and let t of them occur in 1,...,s. The
residue law of a uniform selected positive integer at most X has distance

    TV=t(m-t)/(m(qm+t)) <= 1/(4q), q>=1,

from uniform. Any full-law comparison balance b survives at least as b-TV.
For parity words N=2^ell, so 2^ell=o(X) restores uniformity; this estimate
does not reach a fixed source as ell tends to infinity.

Proofs and explicit positive/hostile controls:
[bridge note, sections 4–6](../../05-knowledge/results/crossroads_poset_20260926_bridge.md),
[width audit, sections 4–6](../../05-knowledge/results/crossroads_poset_20260926_width.md).
No result here imports the unverified submitted poset papers or proves Collatz.
