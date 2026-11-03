# Feedback-core gap formula and exponential spectra

## Exact theorem
Let T be a tournament on [n] with reference order 1<...<n. Let B be the pairs i<j for which the arc is j→i, and let S be the set of vertices incident with B. Put U=V−S and s=|S|≤2|B|. For s>0, let π=(p1,...,ps) range over permutations of S. There are s+1 gaps, numbered 0,...,s. For u in U define its eligible gaps

Eπ(u)={j in {0,...,s}: (j=0 or p_j<u) and (j=s or u<p_{j+1}}.

Define Dπ={j in {1,...,s−1}: the tournament arc is p_{j+1}→p_j}. Then the number H(T) of directed Hamiltonian paths is

H(T)=Σ_π Σ_{F⊆Dπ} (−1)^|F| ∏_{u∈U} |Eπ(u)\F|.

The empty product is 1. If s=0, H(T)=1. Complexity O(s! 2^(s−1) n s), assuming the input order is given and counting arithmetic operations (integer bit costs extra).

Proof: In any directed Hamilton path restrict to S to get π. Between special vertices, all vertices are ordinary. Every edge touching an ordinary vertex points upwards in the reference order. Hence ordinary vertices in a gap occur in increasing order, and u may be assigned to precisely the eligible gaps above. Each unrestricted assignment produces one sequence. For j in Dπ, the gap j must be nonempty, otherwise its two boundary special vertices are in the wrong orientation. Inclusion-exclusion on these empty-gap events yields the stated product. Every Hamilton path has exactly one special permutation and gap assignment, so no overcounting occurs.

## Consequence for constructive spectrum families
Fix the ordered special set S and its internal tournament, and stipulate that all edges involving newly inserted ordinary vertices point forward in the reference order. Fix nonnegative integers c_0,...,c_s. For each integer t≥1 insert exactly t c_i ordinary vertices into the i-th interval of the order around the special vertices. Then

H(T_t)=Σ_{π,F} (−1)^|F| λ_{π,F}^t,
λ_{π,F}=∏_{i=0}^s |Eπ(u_i)\F|^{c_i},

where u_i is a notional point in the i-th interval, Eπ(u_i) is well-defined for every vertex in that interval, and x^0=1 including 0^0. All bases λ are nonnegative integers. Therefore H(T_t), t≥1, satisfies a constant coefficient integer linear recurrence with characteristic polynomial dividing ∏_{distinct λ}(X−λ). The recurrence is an identity of finite exponential sums, not a claim about all tournaments or all spectrum values. Fixed offsets of ordinary counts give integer coefficients multiplying λ^t.

Example: Take m+2 vertices in natural order and reverse only the arc between the first and last. The forward special ordering contributes 1 for m≥1 (all m middle vertices must occupy the interior gap); the reverse special ordering contributes 2^m (each middle vertex may independently go to either outer gap). Thus H=2^m+1 for m≥1; m=0 is exceptional and has H=1.

## Checks and literature
`python check.py` exhaustively compares this formula with permutation-based direct Hamilton-path counts for all labelled tournaments n≤5 (1,099 total). This is a formula implementation check, not a proof of the complexity bound or novelty. Source context: Bar-Noy and Naor, *Sorting, Minimal Feedback Sets, and Hamilton Paths in Tournaments*, SIAM J. Discrete Math. 3 (1990), 7–20, https://epubs.siam.org/doi/10.1137/0403002; their Hamilton-path/minimal-feedback-set bijection supplies the feedback framework. See also https://doi.org/10.1016/j.dam.2019.09.005 on exact tournament path counts. We make no novelty claim. Search of local canon for exact phrases around single backward edges found no direct matching title, not a comprehensive novelty check.
