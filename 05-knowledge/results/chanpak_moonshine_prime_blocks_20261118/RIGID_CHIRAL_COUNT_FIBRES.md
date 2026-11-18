# Infinitely many rigid nonconverse count-fibre pairs with coprime parents

**PROVED for the original two-source recontact family, j>=1.** The free-action divisibility principle is inherited; the all-order combination with the explicit recontact and deletion signatures is the result here. These are facts about count fibres, not a claim that a particular merged-metagraph edge joins the representatives.

For every j>=1, take the two literal recontact graphs from `UNIT_GAIN_FAMILIES.md`, with

    M=16^j, B=4M, u=(7M-2)/5,
    V=(299M+36-10j)/10,
    WL=B+4BV, WR=u+4BV.

Mark each final added vertex. Let the common path count be N.

**Theorem.** The two graphs are asymmetric, nonisomorphic and not converse-isomorphic. Their distinguished deletion counts dL,dR are coprime and satisfy

    dL-dR=32*16^j,  dL<N/2, dR<N/2.

Both marked weighted starts equal N+5, and N=4 modulo5. Thus the scalar count fibre contains at least two distinct classes not identified by converse with genuinely different smaller physical parents at every level of this infinite construction.

## Why the symmetry group is trivial

The displayed V lies between24M and32M. Its first two binary digits are therefore11 and its bit length is4j+5. The final two additions are both letter1. Every old vertex other than the old root beats both new vertices, the old root already had outdegree3, and every earlier new vertex beats the two final vertices. Hence the last new vertex is the unique vertex of outdegree1.

Rerooting gives T=N+eta with eta=5, since the penultimate letter is1. The R address is even because u is even for j>=1, so its count formula gives N=9+25WR-5K(WR)=4 modulo5. Consequently gcd(N,T)=gcd(N,N+5)=1.

Any automorphism fixing the root acts freely both on Hamiltonian paths and on the ordered spanning path-cover objects counted by T. Its order divides both N and T and is therefore1. The unique-outdegree1 property forces every unrooted automorphism to fix this root. Thus the full automorphism group is trivial. This is stronger than merely finding no automorphism in a finite search.

## Parent counts distinguish the two classes

The reroot theorem gives dL-dR=2^s, with s=4j+5, hence the displayed difference and coprimality. Since each graph has a unique outdegree1 vertex, any isomorphism must match the distinguished vertices and therefore their deletion counts. The unequal counts exclude an isomorphism.

For the sharper bound, in any prefix of L's word put q=t-h. Its initial value is14 and its first1 can lower it by only4, giving q>=10 thereafter. For R the analogous bound is q>=12. Indeed leading1s lower q only at the first step; each later run0^k1 has net nonnegative change, since0 increases q by t and the next1 subtracts at most d=5, while t>=5. At the final step,

    N-2d_parent=q_pre-e_pre,

and e_pre is0 or5. Thus N-2dL>=5 and N-2dR>=7.

## A directed deletion signature rules out converse equivalence

In both graphs the original vertex4 is the unique indegree1 vertex. Its original outdegree was3 and it beats every new vertex. The other original outdegree3 vertex is the old marked root, which loses to at least two added1 vertices. Other original vertices started with lower degree; new vertices cannot reach outdegree n-2 because each beats at most one original vertex.

Deleting original vertex4 from L leaves the four-vertex rooted profile (3,0,12,3). Consequently its deletion count is

    D4_L=3+12WL-3K(WL).

Since WL is even,

    2D4_L-N=WL-K(WL)-3>0.

The last inequality follows from K(WL)=K(V)+2V, K(V)<=V and WL=B+4BV. Thus D4_L>N/2, whereas deleting R's unique outdegree1 vertex gives dR<N/2. A converse-isomorphism would have to match L's unique indegree1 vertex to R's unique outdegree1 vertex and preserve their deletion counts. The opposing inequalities exclude it. The same comparison with dL also excludes a self-converse isomorphism of L. Consequently L, its converse L*, and R are three distinct isomorphism classes with the same count N at the same order. In particular the count-fibre multiplicity is at least3 at every level. No all-order claim that R is non-self-converse is needed.

## What this constrains in the metagraph

Any compression of a count fibre that retains a representative's continuation obligations must distinguish these two classes despite their identical H and marked T. The marked deletion count supplies an exact separating coordinate, and the directed high/low deletion signature also detects the converse mismatch. At these levels, count equality is not caused by a hidden automorphism or by identifying a converse pair.

This does not assert anything about adjacency in the merged metagraph, about every count fibre, or about every asymmetric tournament. In particular it is not a proof of Chan–Pak. It gives an all-order test family that a proposed recursive metagraph encoding must represent faithfully.

Replay: `python3 finite_edit_recontact/rigid_fibre_replay.py`. It checks literal degree certificates and deletion-source profiles, exact transfer and separation inequalities through j=12, and direct isomorphism/automorphism controls for the first two levels using a standalone degree-refined backtracking checker.
