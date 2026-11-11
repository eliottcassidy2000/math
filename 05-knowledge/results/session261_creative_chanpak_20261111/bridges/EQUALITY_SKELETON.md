# Equality is hereditary and has a two-component chain skeleton

**PROVED necessary structure and exact local obstruction formula.** This continues `EXTENSION4.md`. No sufficiency or complete finite forbidden-pattern classification is claimed.

Write

    S(T,z)=H(T)−H(T−z)−d(G_z)+1 ≥ 0.

## 1. Monotonicity under rooted induced extension

If (T,z) is an induced rooted subtournament of (T',z), then

    S(T',z) ≥ S(T,z).                                    (1)

For each old even nonempty vertex set W, its p(W), ε(W), and κ(W) in equation(6) of `EXTENSION4.md` are unchanged. The complementary Hamilton count can only increase: inherited OCF includes every old packing as a packing of the larger tournament with new vertices unused. All newly indexed W contribute nonnegative fibres. Summing proves (1).

Thus the equality class S=0 is hereditary under rooted induced subgraphs, and each strict local witness remains strict after arbitrary induced extension. This does not say S is monotone under old-edge changes.

## 2. The first possible forbidden rooted patterns have order five

If |T|≤4 then S(T,z)=0. The only nonempty even subsets of U=T−z have size two, and each A→B graph edge contributes 2 both to the newborn and domination expansions. There are no other terms.

For |T|=5, let U have four vertices and let p be its number of spanning A-to-B Hamilton paths. Every size-two term again cancels, so exactly

    S(T,z)=2p−ε(V(U))2^{κ(V(U))}.                        (2)

Consequently:

- if G is connected, S=2(p−1);
- if G is 2K₂, S=2(p−2);
- otherwise S=2p.

The first case applies to connected bipartite graphs of even order four; in particular it includes a star, P₄ and K₂,₂. Equation(2), rather than a broad tournament census, is an exact classifier of all minimal-order strict witnesses. Every strict rooted order-five witness is automatically minimal by the preceding paragraph.

## 3. A wrong nesting direction creates a strict P₄ witness

Suppose a,a'∈A have common neighbour b∈B, a' has a neighbour c∈B not adjacent to a, and a→a'. On W={a,a',b,c}, G[W] is P₄. There are at least two A-to-B Hamilton paths:

    a',c,a,b,

and, according to the orientation between b and c, either

    a,a',c,b       if c→b,
    a,a',b,c       if b→c.

Thus p(W)≥2, although this connected even fibre needs only one. Its contribution is at least2, so S(T,z)≥2. The statement is invariant under simultaneously reversing all tournament edges and interchanging the shore names.

This is a literal five-vertex rooted obstruction, not a comparison of numerical spectra.

## 4. Equality forces G to be a disjoint union of chain graphs

An induced P₅ in G has two same-shore internal vertices with a common neighbour and with an exclusive neighbour each. Whichever direction their tournament edge takes, one choice of four vertices satisfies the wrong-nesting obstruction in §3. Therefore S=0 forbids induced P₅ in G.

For completeness, a connected bipartite graph without induced P₅ is a chain graph. If two vertices a,a' in the same shore had incomparable neighbourhoods, choose exclusive neighbours b,b'. If a,a' have a common neighbour c, the vertices b,a,c,a',b' induce P₅. If they have no common neighbour, a shortest connecting path has even length at least four, and its first five vertices induce P₅. Hence all neighbourhoods in a shore are nested. The same argument applies componentwise.

Therefore every nontrivial connected component of an equality graph is a chain graph. This is a necessary equality statement; the universal inequality itself has no chain restriction.

## 5. Each nontrivial graph component must induce a transitive tournament

Fix one connected chain component C of G, assuming S=0.

For A vertices in C with strictly nested neighbourhoods, the larger neighbourhood must point to the smaller one in U. Otherwise their nonempty common neighbourhood and an exclusive neighbour give §3's obstruction. Dually, for B vertices the smaller neighbourhood must point to the larger one.

Vertices with equal opposite-shore neighbourhood form a twin class. Each such class induces a transitive tournament. Otherwise choose a directed triangle in the class and one common opposite-shore neighbour (nonempty because C has no isolates). The resulting four-vertex induced graph is K₁,₃. Its three same-shore vertices form a directed triangle, so there are three spanning A-to-B paths, or the reversed dual. The local fibre contributes 2·3−2=4, contradicting S=0.

The forced directions between strict neighbourhood classes and transitivity within each class exclude all same-shore directed triangles. A mixed triangle with two A vertices would require an exclusive B neighbour and precisely the forbidden direction between their nested neighbourhoods. The B-dual excludes the other mixed triangles. Thus U[C] has no directed triangle and is transitive.

## 6. At most two nontrivial components; every graph edge is a feedback pair

Three nontrivial graph components contain three edges inducing 3K₂. Their six endpoints have κ=3; the universal fibre contributes at least2·3!−2³=4. Hence an equality graph has at most two nontrivial components, in addition to any isolated vertices.

For any edge ab of G, its size-two fibre is exactly

    2[H(U−{a,b})−1].                                    (3)

Equality therefore forces U−{a,b} transitive for **every** graph edge. In particular any graph edge disjoint from a directed triangle in U yields a strict rooted order-six witness, with slack at least4.

Combining the results, if G has an edge and S=0, then:

- G is the disjoint union of at most two nontrivial chain components and isolated vertices;
- each nontrivial component induces a transitive subtournament of U;
- deleting the endpoints of any graph edge leaves U transitive.

For G empty, S=0 always: U is the ordinal sum U[B]→U[A], and inserting z between these blocks changes no Hamilton count. Cycles inside the shores are then allowed. This exceptional case must not be erased.

## 7. Why this is still only a necessary skeleton

The three bullets are not sufficient. Take U on vertices a,x,b,c with A={a,x}, B={b,c}, and orientations

    x→a,  a→b, a→c,  b→c,  b→x, c→x.

Then G is K₁,₂ on {a,b,c} plus the isolated vertex x. Its nontrivial component is transitive, there is only one such component, and deleting either graph-edge pair leaves two vertices, hence a transitive tournament. Nevertheless

    x,a,b,c

is a spanning A-to-B Hamilton path. G has odd components of orders3 and1, so ε=0 and equation(2) gives S≥2 (in fact S=2). This is another minimal rooted order-five obstruction.

A useful local consequence is that an isolated A vertex cannot point to a nonisolated A vertex of degree at least two in an equality example; two of its B neighbours produce this witness. Dually, a degree-at-least-two B vertex cannot point to an isolated B vertex.

The remaining architecture is now concrete: examine cross-component and isolated-vertex seams within a two-chain skeleton, retaining the edge-feedback condition. A claim that rooted patterns of order at most seven completely characterize equality would still need a proof. The present results identify explicit minimal patterns and bounded witnesses for each stated necessary condition, not that stronger classification.
