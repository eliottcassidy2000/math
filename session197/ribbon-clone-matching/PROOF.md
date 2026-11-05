# Clone matching ranks and the finite-seed obstruction

Session197. **All-order structural and counting statements.** These results isolate limitations of clone-only inverse constructions; they do not exclude other physical births or prove any new global Hamilton-path gap. No literature-priority claim is made.

Let T be a nontrivial strong tournament, let r be a vertex, and let T_m replace r by a transitive module of m vertices, for m>=1. Every module vertex has the old incidences of r toward the other vertices. Write H and C for the numbers of Hamilton paths and directed Hamilton cycles modulo rotation.

Partition the old vertices into O=N+(r) and I=N-(r). Define the bipartite graph B_r on O and I by an edge oi exactly when o->i. Its matching number is nu. Every edge is precisely a directed triangle r,o,i,r. Strongness implies nu>=1: otherwise every I vertex would dominate every O vertex, preventing any path from O back to r.

## 1. Two insertion facts, with their scope

Hamilton-path count is monotone under arbitrary induced extension. For a path P and a new vertex x, insert x first if x beats the first vertex, last if the last vertex beats x, and otherwise at the first adjacent transition p_i->x->p_(i+1). Such a transition exists. The chosen rule is deterministic; deleting x recovers P, so the map is injective. Iterate for all extra vertices.

Hamilton-cycle count is monotone between STRONG induced tournaments S contained in T. The following stronger form is useful: a cycle of S extends to a cycle of T with its old cyclic order preserved.

If an outside vertex x has both signs toward S, the old cycle has an edge u->v with u->x->v, permitting deterministic insertion. If no outside vertex is mixed, partition the outside vertices into those Q dominated by S and those P dominating S. Strongness of T forces both sides and an arc q->p with q in Q and p in P. Otherwise S could not reach P or be reached from Q. Replace an old cycle edge u->v by u->q->p->v. Each step enlarges a strong induced subtournament by one or two vertices. Continue until T. Every insertion preserves the old cyclic order, and deleting the added vertices recovers the original cycle, giving the injection.

The second statement is not monotonicity under arbitrary weak intermediate extensions. Its initial and final tournaments must be strong. All uses below meet that requirement: T_m is strong, and each rooted-triangle subgraph with all m clones is strong.

## 2. Rooted matching controls exact clone growth

Let S(m,j) denote a Stirling number of the second kind. There are nonnegative integers a_j and b_j, depending only on (T,r), such that for m>=1,

    H(T_m) = sum_(j=1)^(nu+1) a_j j! S(m,j),
    C(T_m) = sum_(j=1)^nu b_j j! S(m,j).             (1)

Moreover a_(nu+1)>0 and b_nu>0. Consequently the exact leading exponential bases are nu+1 for H and nu for C.

To see the path formula, split a Hamilton path into its maximal clone runs. Delete the clones while retaining the ordered old-vertex fragments and the positions of the j nonempty clone runs. This is a finite skeleton independent of m. Every assignment of the m distinct clone vertices to the j ordered nonempty runs works, because the vertices in each run must be increasing and all external incidences are identical. There are j! S(m,j) assignments. The same construction for a cycle, rooted at one fixed old vertex to remove rotation ambiguity, gives the cycle formula with nonnegative integral skeleton counts.

A path with j clone runs has j-1 disjoint internal old fragments, each starting in O and ending in I. Each such fragment contains a consecutive O->I transition, and therefore contains a rooted triangle. The transitions give a matching of size j-1 in B_r. Thus j<=nu+1. A cycle with j clone runs has j such old return fragments, so j<=nu.

Conversely, choose nu disjoint rooted triangles and nu+1 clones. Alternating the clones with the nu two-vertex return paths gives a Hamilton path of that induced subgraph with nu+1 clone runs. Insert every remaining old vertex by Section1. The existing clone runs remain separated, proving a_(nu+1)>0. For cycles, use nu clones and the nu return paths in cyclic order. Their union is strong. The cycle extension in Section1 preserves the cyclic order and therefore the nu separated clone runs, proving b_nu>0.

There is also the explicit lower bound

    C(T_m) >= (nu-1)! nu! S(m,nu).                  (2)

On the selected rooted-triangle subgraph, put the nu return paths in any cyclic order and assign the clones to nu ordered nonempty increasing blocks. Fixing the first return path leaves (nu-1)! orders, and there are nu! S(m,nu) block assignments. These cycles are distinct, and the strong-induced injection transfers them to T_m. For m<nu the right side is zero, so the same bound remains valid.

## 3. Bounded-cycle clone rays are binary

The following are equivalent:

1. The numbers C(T_m) are bounded as m varies.
2. nu=1.
3. Every edge of B_r has a common endpoint.
4. C(T_m)=C(T) for every m>=1.

The matching-one/star equivalence follows directly: two edges with different endpoints would form a matching; if two edges share one endpoint, every further edge must use that same endpoint. Section2 proves the other equivalences.

In this case there are exact integers A>=C(T)>0 and B such that

    H(T_m)=A 2^m+B,
    A=(H(T_2)-H(T))/2,
    B=2H(T)-H(T_2).                                (3)

Indeed (1) has only j=1,2 terms, a_1=H(T), and A=a_2. Every original Hamilton cycle, cut at r, gives one of the two-run path skeletons with empty outside fragments, so A>=C(T).

If nu>=2, select two rooted triangles in (2), obtaining

    C(T_m)>=2^m-2.

In particular C(T_3)>=6. Thus a transitive-clone ray that stays within C<=2, or even within C<=5, must have the binary form (3). It cannot have a genuinely ternary leading term. This does not say its offset B is a unit, nor does it exclude a ternary construction using non-clone births.

## 4. Finite strong seeds and arbitrary choices of cloned vertices

Fix finitely many nontrivial strong seed tournaments S_i of orders k_i. Starting from any seed, allow any finite sequence of replacements of any current vertex by a nonempty transitive module. The chosen vertex may change at every step; there is no persistent-root assumption.

Every descendant is still a transitive blow-up of the ORIGINAL seed: it has positive multiplicities m_1,...,m_(k_i), and replacing a vertex inside one transitive module only enlarges that module, preserving its transitive order. This is an induction on the construction sequence.

Every vertex of a strong seed lies on a directed triangle. Select the full module of multiplicity m_j at that vertex and one vertex from each of the other two triangle modules. The induced tournament is a cyclic triangle with one transitive part of size m_j and two singleton parts. Its Hamilton-path count is exactly

    2^(m_j)+1.                                     (4)

For completeness, the paths with the two singleton vertices in their forward order split the module into an arbitrary before/after pair of increasing subsequences, giving 2^(m_j) choices. In their reverse order the module must lie between them, giving one additional path. These cases exhaust all paths because one singleton has only the other singleton as an out-neighbor within this induced graph.

Induced-path monotonicity from Section1 now gives H>=2^(m_j)+1 for every j. Therefore the number of distinct output values H<=X from the finite seed bank is at most

    sum_i floor(log_2(X-1))^(k_i),                  (5)

for X>=3. This bound counts multiplicity vectors, so graph isomorphisms and equal counts can only reduce it. An optional singleton seed contributes only H=1, since all its descendants are transitive.

In particular this entire grammar has zero density among the positive integers and among the odd integers. It cannot cover a cofinite odd spectrum. A binary inverse proof based only on transitive cloning of finitely many fixed graph seeds must therefore fail somewhere, even if it changes roots freely and observes many unit affine moves. It needs genuinely non-clone operations or an unbounded family of new seed types.

The strong-seed hypothesis is used in (4) for EVERY seed vertex. We do not silently apply the logarithmic multiplicity bound to singleton SCCs of a weak seed, whose multiplicity can grow without increasing H.

## 5. Verification and boundary examples

The standalone verifier checks all labelled strong tournaments through order five and every root, using subset path/cycle DP on literal clone graphs. It extracts the Stirling-basis coefficients and verifies the predicted exact degree, positivity, star cases, and the two-triangle lower bound. A regular seven-vertex example supplies matching rank three. It independently checks (4) through clone size ten, and verifies the original-seed blow-up invariant under explicit replacements at changing vertices.

All-order conclusions come from the proofs above, not extrapolation of these controls. The program does not enumerate a global spectrum or treat gcd one as a realization certificate.
