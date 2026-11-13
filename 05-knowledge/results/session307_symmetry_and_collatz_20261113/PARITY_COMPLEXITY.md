# Parity factor complexity controls the occupation of a putative divergent orbit

For a nonrepeating positive shortcut Collatz orbit x_i, let P_w(k) be the number of distinct length-k consecutive factors of its infinite parity word w, and N(X)=#{i : x_i≤X}. Then for every X≥1,
```
N(X) ≤ P_w(ceil(log₂(X+1))).
```
Indeed length-k parity vectors distinguish residue classes modulo 2^k (Terras parity bijection, inherited in [THM-4634](../../../01-canon/theorems/THM-4634-collatz-parity-number-dictionary-irrationality-exponents-of-parity-vectors-and-jarnik-strata.md)). Two distinct positive orbit values at most X are distinct residues for this k; hence their next k parities are different. The argument uses actual consecutive factors beginning at visited nodes, not all words admissible under an observer. It also works for every fixed odd b in the shortcut (3x+b)/2 at odd x, as long as the trajectory is positive and nonrepeating.

Consequences:

1. If w has linear (respectively polynomial) factor complexity then N(X)=O(log X) (respectively O((log X)^d) for degree d). In particular a hypothetical orbit coded by a Sturmian word would occupy at most ceil(log₂(X+1))+1 positions below X. This is a necessary sparsity constraint, not proof that such an orbit exists.
2. If w has zero factor entropy then for every ε>0 the orbit has N(X)=O_ε(X^ε). Its Dirichlet sum ∑ x_i^(−s) converges for every s>0, by dyadic shells. More generally if its entropy in bits is η, the abscissa of convergence is at most η (and separately bounded by h from THM-4499). The entropy exists by submultiplicativity P_w(k+l)≤P_w(k)P_w(l).
3. Conversely, if x_i≤C(i+1)^α for all i, a nonrepeating orbit requires η≥1/α. For all large k take all i with C(i+1)^α<2^k to obtain this lower bound on P_w(k), and take k→∞. Thus a hypothetical divergent orbit with zero factor entropy cannot have even a global polynomial upper envelope. This is independent of whether its average parity slope equals log₃2.

These observations convert a size estimate into a necessary symmetry-breaking constraint: at scale 2^k no two distinct orbit nodes share the same continuation word of length k. Contrast the [low-complexity hard-class work](../collatz_procgen_20260922_hard_class.md), where repetitiveness controls Diophantine approximations to real and 2-adic coding series; here distinct forward words directly bound an individual orbit's spatial occupancy. The proof does not infer low entropy for all integer trajectories, or rule out superpolynomial escape.

