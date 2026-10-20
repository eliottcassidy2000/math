# Maximum degree of the tournament line metagraph

Status: **PROVED conditional on the stated interpretation of the symbols**;
the Hamiltonian-path upper bound is **CITED**. No exact all-order extremal
formula or leading asymptotic constant is claimed.

## 1. Definition and inherited counting object

The user's clarification identifies the tournament metagraph, but no unique
definition of D(z) or M_k was recovered. We explicitly choose the degree
interpretation:

    D_n(z) = sum over converse-merged nodes v of z^deg(v),
    M_n = deg D_n = max_v deg(v).

This is the **complement-line multigraph**, not the single-arc-flip graph.
Keep parallel lines; count a loop twice in the degree. Write N=binom(n,2),
m=binom(n-1,2), a_n for ordinary tournament classes and s_n for self-converse
classes. Its vertex count is V_n=(a_n+s_n)/2, and its total degree is 2^m.
These are inherited from [THM-1430, tiling/class
dictionary](../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md)
and [THM-4731, twins and
depth](../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).

For an ordinary class T let h(T) be its number of directed Hamiltonian
paths and g(T)=|Aut(T)|. Its merged degree is h/g when self-converse,
and 2h/g otherwise. Thus

    D_n(z) = sum_(T self-converse) z^(h(T)/g(T))
           + (1/2) sum_(T not self-converse) z^(2h(T)/g(T)).

The sums here run over ordinary isomorphism classes. This is an exact
general coefficient formula, not a solution of the extremal optimization.
The companion [metagraph note](boundary_support_20261020_metagraph.md)
gives the finite coefficient tables and the bridge/wiring results.

## 2. Elementary class-count asymptotics, with a proof

We rederive the estimates needed here rather than importing an unexamined
asymptotic claim from a historical reflection:

    a_n ~ 2^N/n!,       s_n=o(a_n),       V_n ~ 2^(N-1)/n!.       (1)

Burnside's lemma gives a_n=(1/n!) sum_pi Fix(pi).
Only permutations of odd order can be tournament automorphisms: an
automorphism of even order has a nontrivial involution power, which
interchanges two vertices and reverses their arc.

Suppose a nonidentity odd-order permutation moves k vertices (k>=3).
Exactly binom(n-k,2) unordered pairs are fixed individually. Every other
pair orbit has odd length at least three. Therefore its number of pair
orbits is at most

    N - (2/3)(N-binom(n-k,2)) = N-k(2n-k-1)/3.

Each pair orbit supplies at most one orientation bit. The number of such
permutations is at most n^k, so the ratio of the nonidentity Burnside terms
to the identity term is bounded above by

    sum_(k=3)^n n^k 2^(-k(2n-k-1)/3)
    <= sum_(k=3)^n (n 2^(-(n-1)/3))^k = o(1).

This proves the first assertion of (1).

For the second, a labelled self-converse tournament has a permutation pi
that reverses all arcs. Such pi fixes at most one vertex. Individually
fixed unordered pairs must be transpositions of pi, at most floor(n/2);
every other pair orbit has length at least two. There are at most
(N+floor(n/2))/2 pair orbits. Consistency may reduce the count further, so
for any fixed pi there are at most 2^((N+floor(n/2))/2) such tournaments.
Taking a union over n! permutations, and bounding the number of
self-converse classes by the number of their labelled representatives,

    s_n/a_n <= (n!)^2 2^(-N/2+floor(n/2)/2) -> 0.

This deliberately loose bound suffices. It proves the remaining assertions.

## 3. The asymptotic scale of M_n/n!

The maximum degree is at least the average degree. From (1),

    M_n >= 2^m/V_n = (1+o(1)) n!/2^(n-2).                      (2)

For the upper bound use g(T)>=1 and the published Hamiltonian-path bound:

    max_T h(T) <= C n^(3/2) n!/2^(n-1).

**CITED input:** Noga Alon, *The Maximum Number of Hamiltonian Paths in
Tournaments*, Combinatorica 10 (1990), 319–324, Theorem 1.1; checked in the
[author's manuscript](https://www.tau.ac.il/~nogaa/PDFS/hamilton.pdf),
printed page 3, with the proof on printed pages 4–6.
This is a sufficient upper bound, not a claim that it is the sharpest
currently available polynomial factor. No lower bound for unweighted
Hamiltonian-path maxima is automatically transferred through |Aut(T)|.

Since every merged degree is at most twice h(T),

    M_n <= C n^(3/2) n!/2^(n-2).                              (3)

Consequently the following are rigorous:

    M_n/n! = 2^(-n+O(log n)),
    lim_(n->infinity) (M_n/n!)^(1/n) = 1/2,
    lim_(n->infinity) M_n/n! = 0,
    M_n^(1/n)/n -> 1/(2e).                                   (4)

Equations (2)–(3) leave a polynomial-factor gap. In particular they do NOT
prove M_n ~ c n!/2^n for any c, and they do not determine the maximizing
classes at all n. The external maximum h(T) is not M_n: the latter retains
the automorphism divisor and the self-converse/paired distinction.

## 4. Why this bears on boundary support

The average-degree lower bound uses full state mass before quotienting:
each tiling is exactly one line endpoint. The maximum-degree question
requires the sizes of the individual path-presentation fibres, not just
the total or the number of isomorphism classes. This is the same loss that
prevents replacing a cubic pole by another one merely because their
boundary supports have equal cardinality. The support-preserving strip
reduction retains every boundary word; a degree enumerator does not retain
which fibres the complement involution joins.

No reduction from arbitrary tournament wiring to cubic noncolourability,
and no Petersen-minor consequence of (4), is asserted.
