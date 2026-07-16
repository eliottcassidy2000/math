# HYP-7097 — The circulant max-cut theorem for the parallel-class book (general odd n)

**Status:** RESOLVED (death-star-2026-07-16-S28) — **THM-913 upgraded to PROVED for all
odd n**: (L1) ξ(d) = (d−1)(n−1−d)/2 refereed exactly (all odd n ≤ 31); (L2) the three-line
Faulhaber-cancellation identity F(contiguous) = Z(n) (refereed m ≤ 60 + independent arc
bookkeeping); (L3) optimal over ALL 2-page drawings via AAFRS — coloring-optimality not
needed. REFINEMENT: ŵ(k) = (n−1)/2 − n/(4sin²(πk/n)) identified exactly ⟹ coloring
max-cut ⟺ arc-maximality of the cycle Green's form (⟺ arc minimizes Σ_{pairs}d(n−d)) —
exhaustively verified all odd n ≤ 19; the general-n arc lemma is the one named residual
(affects only coloring-universality).

Targets: (1) PROVE ξ(d) = (d−1)(n−1−d)/2 (the class-crossing profile); (2) PROVE the main
identity F(contiguous split) = Z(n) — hand proof found at claim time: with M = 2m−1,
F = Σ_{e=1}^{m−1} e(M−2e)(M−e)/2; the M²S₁ and −3M·S₂ Faulhaber terms cancel exactly
(S₂ = m(m−1)M/6), leaving 2S₃/2 = m²(m−1)²/4 = Z(2m+1) — THREE LINES; (3) the coloring
max-cut: Fourier-reduce (ŵ(k) closed form via the cycle Green's function 1/(4sin²(πk/n)));
prove the compression cases; extended exhaustive verification n ≤ 19–21; name what remains.
Then THM-913 upgrades to PROVED-construction-optimal (via AAFRS lower bound — which needs
NO coloring-optimality: any drawing achieving Z(n) is optimal among ALL 2-page drawings).

-> THM-913, THM-906(II), boxeph-S39 (2Z(n) refutation noted), opus THM-900; death-star-S28.
