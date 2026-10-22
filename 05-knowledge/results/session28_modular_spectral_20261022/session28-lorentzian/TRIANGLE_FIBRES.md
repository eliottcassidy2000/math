# Exact Hamiltonian response along score-invisible triangle circulations

Status: PROVED by path contraction, FINITE-EXACT independent subset-DP controls. A specialization of the classical/previously studied even-path-forest expansion, not a claim of a new general expansion.

Let n≥3k. On k vertex-disjoint triples, give each of the three edges a common cyclic bias t_i∈[-1,1]; choose all remaining edges independently with probability1/2. Thus each triangle's expected circulation is oriented cyclically, and its bias contribution has zero incidence Bx=0. Every vertex's expected outdegree remains(n−1)/2, regardless of all t_i.

Write F_n(t)=E H(T). Then

    F_n(t_1,...,t_k)
      = 2^(-(n−1)) Σ_(J⊆{1,...,k}) 6^|J| (n−2|J|)! ∏_(i∈J)t_i².

Equivalently,

    F_n(t)=2^(-(n−1)) ∫_0^∞ e^(−u) u^(n−2k) ∏_(i=1)^k(u²+6t_i²) du.

The integral is elementary Gamma integration of a polynomial, not a modular or spectral claim. It gives a useful positive representation of the exact finite polynomial.

## Proof

Expand the expected path count as

    2^(-(n−1)) Σ_(Hamiltonian vertex orders P) ∏_(e along P)(1+signed bias_e).

A selected set of biased edges lies in disjoint triangles. A Hamiltonian path can contain at most two edges from any one triangle. If a selected triangle contributes a single edge, reversing that selected path component pairs its sign with its negative, so its total coefficient is0. The only surviving selections use either0 or2 edges from every triangle.

For a selected set J of j triangles, there are3^j choices of their two-edge paths. Their signs in cyclic orientation are positive: each is a coherent directed two-edge path, and reversing a two-edge component preserves the product of its two signs. Each of the j components has2 possible traversal directions. Contracting these length-two paths leaves n−2j blocks, which can be ordered in(n−2j)! ways. The coefficient is consequently3^j·2^j·(n−2j)!=6^j(n−2j)!, as asserted.

## What survives, and what does not

- F_n is strictly increasing in each t_i² on this family when the other biases are fixed.
- The expected score vector, hence the score-Lorentzian point formed from expected scores, is constant. At n25 it is always c=(12^25;70).
- The expected cyclic-triangle count is

      E c3(T)=binom(n,3)/4 + (3/4) Σ_i t_i².

  The source is edge variance/covariance, not a change in expected scores.
- The score lift of an expectation is not the expectation of its quadratic norm. In particular Q(E v)≠E Q(v) in general. Jensen/variance is the missing sidecar.
- These are probabilistic ensembles, not a theorem that every individual score fiber has a Hamiltonian improvement operation, and not a Hamiltonian-maximizer classification.

For one triangle,

    F_n(t)=n!/2^(n−1) · (1+6t²/[n(n−1)]).

At n25 a fully biased triangle increases the expected Hamiltonian count by exactly1% while keeping all expected scores12. For eight disjoint fully biased triangles at n25,

    F_25(1,...,1)/F_25(0,...,0)
      =248148042263449/229042013200000
      ≈1.0834171373038255.

The eight-triangle source has E c3=581, whereas the completely unbiased source has575. The Lorentzian lift of their common expected score vector would instead encode650 if incorrectly treated as an actual deterministic regular tournament. This concrete discrepancy exposes the extension mistake.

## Checks

`python3 check_triangle_fibres.py > check_triangle_fibres.json` compares the contraction formula to an independent exact rational subset dynamic program for six configurations at n3,4,6,7,9, including negative and fractional biases and three simultaneous circles. It separately computes expected cyclic triples and verifies the constant expected scores. The n25 ratio uses the proved closed polynomial, not25-vertex path enumeration.
