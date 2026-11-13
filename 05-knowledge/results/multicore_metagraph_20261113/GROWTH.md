# Exact growth with a marked final vertex

**PROVED:** the forward transfer and orbit-weight identities below, for every
finite tournament. **FINITE-EXACT:** the implementation through order eight.
This refines the failure of bare-node growth in
[THM-796, three-sorted recursive incidence, section 8](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md).
It does not reconstruct complement-line incidence from a rooted state.

## A presentation and its continuation state

A presentation is a tournament with an ordered Hamiltonian path, up to an
isomorphism preserving that path. Name the path 0,1,...,n−1. Its consecutive
arcs are fixed; its other binomial(n−1,2) arcs are free bits. Thus there are
2^binomial(n−1,2) presentations. This is the inherited tiling dictionary,
not an identification of tilings with unmarked tournament classes.

An **end state** is an isomorphism class (T,p) in which some Hamiltonian path
ends at p. It forgets the preceding path order. To append a new last vertex
x, force p→x and choose all other incident arcs freely. Let

    F_n((T,p),(U,x))

count the 2^(n−1) choices that yield the indicated child end state. This
integer depends only on (T,p): a rooted isomorphism transports the attachment
choices bijectively. Hence end states form an exact forward quotient of the
presentation growth process. Each row of F_n sums to 2^(n−1). No assumption
that presentations are uniformly distributed inside an unmarked node is made.

Write e_T(p) for the number of Hamiltonian paths ending at p and A(T,p) for
the root stabilizer in Aut(T). A path fixes every vertex, so this stabilizer
acts freely on the paths. The number of presentations in the end state is

    mu(T,p) = e_T(p)/|A(T,p)|.

Starting with the singleton of weight one, the exact integer recursion is

    mu_(n+1)(U,x) = sum_(T,p) mu_n(T,p) F_n((T,p),(U,x)).       (G1)

Every child presentation has exactly one parent presentation, obtained by
deleting its last vertex. Thus (G1) counts each once. Its total mass is
2^binomial(n,2). Forgetting the root gives the ordinary class weight
H(T)/|Aut(T)|. Every tournament appears because it has a Hamiltonian path.

This is a sufficient state theorem, not a minimal-state theorem. The full
rooted tournament is retained. An endpoint count vector alone is not asserted
to suffice. Nor does (G1) make backward deletion Markovian on end states:
the penultimate vertex has been forgotten.

## A second description by rooted deletion cards

Fix representatives (T,p),(U,x), with orders n,n+1. Let I be the number of
vertices v of U−x such that v→x and (U−x,v) is isomorphic to (T,p). Then

    F_n((T,p),(U,x)) = I |A(T,p)| / |A(U,x)|.                 (G2)

For each such v there are |A(T,p)| rooted isomorphisms into the deletion
card. Each realized attachment is counted |A(U,x)| times by isomorphisms
of its rooted child with (U,x). This proves (G2), including integrality.
Substitution in (G1) is precisely the Hamiltonian endpoint recurrence
e_U(x)=sum_(v→x) e_(U−x)(v), with automorphism factors retained.

## An all-size counting formula for end states

Let A(z)=sum a_n z^n count ordinary tournament classes, including a_0=1.
Let S(z) count nonempty strong tournament classes (including the singleton),
S_bullet(z) their vertex-rooted classes, B(z) end states, and R(z) all
vertex-rooted tournaments. Strong components form a unique ordered word;
automorphisms preserve each component. A vertex can end a Hamiltonian path
exactly when it belongs to the last component. Necessity follows from the
component order. Sufficiency follows by traversing Hamiltonian paths in the
earlier components and ending a Hamiltonian cycle of the last strong
component at the chosen vertex; singleton components need no cycle.
Consequently, as formal power series,

    A=1/(1−S),       B=A S_bullet,       R=A² S_bullet=A B.   (G3)

There is also a direct Burnside expression. For a partition lambda of n
whose parts are all odd, put m_j=the multiplicity of j,
z_lambda=product j^m_j m_j!, and

    o(lambda) = sum_i (lambda_i−1)/2 + sum_(i<j) gcd(lambda_i,lambda_j).

An even permutation cycle cannot fix a tournament: its halfway power swaps
two vertices of that cycle. For odd cycles, each unordered-edge orbit can
be oriented in two ways, and there are o(lambda) such orbits. Thus

    a_n = sum_(lambda odd) 2^o(lambda)/z_lambda,
    r_n = sum_(lambda odd) m_1 2^o(lambda)/z_lambda,
    b_n = r_n − sum_(i=1)^(n−1) a_i b_(n−i).                 (G4)

The second formula marks a vertex fixed by the permutation. G3 explains
why formal division by A gives nonnegative integers, which is less obvious
from G4 alone. These are formulas for counts of states, not their transitions
or complement matches. The standard a_n are [A000568](https://oeis.org/A000568);
deleting the distinguished root identifies r_(n+1) with binary-coloured
unlabelled n-vertex tournaments, the interpretation listed for
[A093934](https://oeis.org/A093934). No new OEIS term is claimed.

|Order n|a_n|Merged nodes|End states b_n|Strong rooted classes|Tiling presentations|
|---:|---:|---:|---:|---:|---:|
|1|1|1|1|1|1|
|2|1|1|1|0|1|
|3|2|2|2|1|2|
|4|4|3|7|4|8|
|5|12|10|33|24|64|
|6|56|34|236|194|1024|
|7|456|272|2681|2395|32768|
|8|6880|3528|50709|47588|2097152|

At order eight this replaces 2,097,152 path presentations by 50,709 end
states for forward counting. It does not claim to outperform specialized
unmarked tournament generators, or to extend A000568's known census.

Reproduction: `python -B 05-knowledge/results/multicore_metagraph_20261113/rooted_growth.py 8`.
Every transfer attachment is enumerated. Independent endpoint subset
counting checks every state; literal path permutations through five and
all tilings through six check the weights by different routes.
