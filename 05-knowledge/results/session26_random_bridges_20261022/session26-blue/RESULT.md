# The one-apex stratum: exact masses, blue incidences, and forced wiring

22 October 2026. Status: PROVED statements with explicit proofs below.
Independent review is a separate artifact. No claim of priority is made.

Let tau(T) be the smallest number of vertices whose deletion makes a
tournament transitive. The one-apex stratum means tau(T)<=1; it includes
the transitive tournament. All metagraph statements use n>=3 and the
fundamental converse-merged multigraph, with loops counted twice in degree.
A line is blue when both pinned-path presentations are reflection-fixed;
nonfixed lines are paired with their reflections to give black edges.

## 1. A feedback-number constraint on EVERY metagraph edge

Suppose T'=kappa_P(T), where P is a directed Hamiltonian path and kappa
reverses exactly the arcs not on P. Then

    tau(T)+tau(T') >= floor(n/3).                         (1)

Partition the ordered vertices of P into floor(n/3) disjoint consecutive
triples, discarding at most two terminal vertices. On each triple the two
path arcs are retained and the chord is reversed. Exactly one endpoint
tournament therefore has a directed triangle on that triple. A feedback
set for that endpoint must meet the triple. The triples are disjoint, so
summing their required deletions proves (1).

Isomorphism and converse preserve tau. Consequently (1) descends to the
actual merged graph, independently of the edge color. In particular:

- The nodes tau<=d form an independent set whenever n>=6d+3.
- Every neighbor of a tau<=d node has tau>=floor(n/3)-d.
- A loop can occur only at tau>=ceil(floor(n/3)/2).

The bound is sharp for every n. Start with the transitive tournament, whose
Hamiltonian path is unique. In its path complement, a triple is cyclic
exactly when its sorted positions are consecutive. If neither or just one
of the two adjacent gaps has length one, its orientation is transitive;
if both do, it is cyclic. Deleting every third path vertex kills all such
triples. A tournament is transitive iff it has no directed triangle, so
the unique neighbor of the transitive node has tau exactly floor(n/3).

A path-sensitive refinement retains the color of EVERY consecutive
three-window, indicating which endpoint owns its triangle. For each
color c, let nu_c be the maximum number of disjoint windows of that color.
Then separately tau(T)>=nu_0 and tau(T')>=nu_1. These are interval-packing
lower certificates, not an exact characterization: nonconsecutive
triangles have been forgotten. For n=8 with alternating window colors,
nu_0=nu_1=2, sharper than the universal sum bound 2.

## 2. Class counts and root-forgetting conventions

We use the deletion-apex classification proved in
[the companion classification note](../session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md).
For clarity its mechanism is repeated here.

A nontransitive tournament's deletion-transitive vertices are the
intersection of all its cyclic triangles. There are at most three.
Three occur only when the sole nontrivial strongly connected component
is C3. If a strong tournament of order n>=4 has two such apices u->v,
the other vertices form a chain. The two apices have insertion cuts in
that chain; strong connectivity forces opposite extreme cuts. The
tournament is exactly C3(1,1,n-2), with the long block transitive.
All other strong one-apex tournaments have one deletion-transitive apex.

An automorphism fixing such an apex fixes the unique order on the
remaining transitive tournament, hence is the identity. Tournament
automorphism groups have odd order: an involution would swap some pair
and reverse its arc, and any even finite group contains an involution.
Thus every strong one-apex tournament of order at least four is rigid.
The exceptional C3 has automorphism order three.

Every nontransitive one-apex tournament consists of a transitive prefix,
one strong core, and a transitive suffix. Its Hamiltonian paths and
automorphisms are precisely those of that core, with the unique outside
orders attached.

Let A_n count ordinary one-apex classes and C_n their self-converse
classes. For n>=2,

    A_n = 2^(n-1) - binom(n-1,2) - 1.                   (2)

For n=2r>=2 and n=2r+1>=3 respectively,

    C_(2r)=r,       C_(2r+1)=2^r+r-1.                  (3)

The number of selected merged nodes is (A_n+C_n)/2.
These are a subclass of A000568, NOT its full enumeration.

Brief derivation: strong orders n>=4 have 2^(n-3) rooted binary words,
with exactly one double-apex class counted twice, giving 2^(n-3)-1
ordinary classes. Order three contributes C3 once. Sum over the
n-k+1 prefix/suffix placements of each order-k core, and add the
transitive class, to obtain (2). For converse, the exceptional double-
apex class is always self-converse. A unique-apex strong class is
self-converse iff its word equals its reverse-complement. Therefore
strong self-converse counts are 1 at even orders >=4 and
2^((n-3)/2)+1 at odd orders >=5, with C3 separately counted once.
Padding is self-converse exactly when prefix and suffix lengths agree;
summing gives (3).

## 3. Exact selected tiling mass

For an ordinary class T write h(T)=H(T)/|Aut(T)|. This is the number of
fixed-Hamiltonian-path presentations (tilings) of that class: automorphisms
act freely on directed Hamiltonian words. Let

    V_n = sum_(ordinary one-apex classes T) h(T).

Then for every n>=2,

    V_n = 4^(n-2) - 3*2^(n-2) - (n^2-7n+4)/2.         (4)

Here is a direct calculation. A rooted strong one-apex word has length
m=n-1, first bit 1 and last bit 0. Number the base vertices 1,...,m;
b_i=1 means apex->base_i. Its path count is
1+sum_(1<=j<i<=m) b_j(1-b_i)w_(i-j), where w_1=2 and w_d=2^(d-1) for d>=2.
For m>=3 the sum over all 2^(m-2) strong words is

    9*2^(2m-6) + 2^(m-3).

Indeed, the mandatory endpoint pair contributes 2^(2m-4), the endpoint-
interior pairs contribute another 2^(2m-4), and the interior pairs
contribute 2^(m-4)(2^(m-2)-2); the constant terms contribute 2^(m-2).
The formula also holds at m=3, where the interior-pair contribution is zero.

At strong order n>=4 all classes are rigid, but the exceptional double-
apex class has two words and H=2^(n-2)+1. Subtract it once. The resulting
strong mass is

    9*4^(n-4) - 3*2^(n-4) - 1.

C3 contributes h=1, not H=3. Summing over all prefix/suffix placements
and adding the transitive class yields (4).

For example, n=3,4,5,6,7 gives selected mass 2,8,43,209,926.
This counts actual tilings in the stratum, not just rooted apex words.

## 4. Blue incidences from the actual reflection

Let B_n be the sum, over self-converse classes in the stratum, of their
reflection-fixed tiling counts b(T). Then

    B_(2r)   = 2^r+r-2             (r>=1),
    B_(2r+1) = 2^(2r-1)+r-1       (r>=1).              (5)

To prove this, first take a rigid strong class with unique apex. It can
be self-converse only at n=2r+1. Its unique anti-isomorphism fixes the
apex and reverses the base order. A symmetric Hamiltonian path has the
apex at its center. Its first half is an increasing set selecting one
vertex from each reversal pair; the second half is forced.

If the first r bits of its reverse-complement-fixed word form the
binary integer J, the number of such paths is

    b(T)=J+1-b_r = 2*floor(J/2)+1.                      (6)

For an explicit count, a transversal with maximum selected index
i>r has 2^(i-r-1) choices; the all-left transversal has maximum r
and contributes once. The path enters the apex exactly when that
maximum has bit 0. Reverse-complement symmetry gives (6).
Here J ranges from 2^(r-1) to 2^r-1, each yielding one rooted class.

Next consider the double-apex class C3(1,1,t), t=n-2>=2.
Write its arcs as u->v->chain->u. Its H=2^t+1 Hamiltonian words are

    S, u, v, chain\S            (all 2^t subsets S, ordered internally)
    v, chain, u                (one extra word).

Its anti-isomorphism swaps u,v and reverses the chain. The extra word
is symmetric. A word of the first kind is symmetric exactly when S
selects one vertex from each chain-reversal pair. Hence

    b(T)=1 if t is odd, and b(T)=2^(t/2)+1 if t is even. (7)

Thus strong blue totals are 2^(n/2-1)+1 for even n>=4, and
3*2^(n-4)+1 for odd n>=5. C3 contributes one blue tiling.
Equal transitive padding preserves these counts; arbitrary unequal
padding is not self-converse. Add the transitive class and sum to
obtain (5).

Within this stratum the points of symmetry (nodes with ONLY blue
incidences) are exactly the transitive class and, when n is odd>=3,
the centrally padded C3 class. Indeed these have h=b=1. Every rigid
strong unique-apex class has the Hamiltonian path apex,base-order,
which cannot be symmetric because its fixed apex is not central.
The same path in the double-apex class ends in the chain rather than
at the anti-image v of its starting apex, so is also nonsymmetric.
Padding preserves such a witness. This is NOT a classification of
points of symmetry outside tau<=1.

## 5. An all-order boundary-edge census

By (1), the entire one-apex stratum is independent for n>=9.
Consequently every line incident to it leaves the stratum, and exact
fundamental colored boundary-edge counts are

    blue boundary edges  = B_n,
    black boundary edges = (V_n-B_n)/2.                (8)

A reflection-fixed tiling supplies one blue incidence. Each pair of
nonfixed reflected tilings supplies one black incidence. Converse
merging preserves these orbit counts. Since no edge has two selected
endpoints at n>=9, incidences equal boundary edges without division
by another factor of two.

At n=9 there are 123 selected merged nodes, 15,989 selected tilings,
131 blue boundary edges, and 7,929 black boundary edges.
At n=10 there are 240 selected merged nodes, 64,751 selected tilings,
35 blue boundary edges, and 32,358 black boundary edges.
Every adjacent node has tau>=floor(n/3)-1.

For n<9 equations (2)-(5) still hold, but (8) cannot be inferred:
internal edges and loops require their own incidence accounting.

## 6. Generating functions and limits

Let a_n=A_n for n>=0, with a_0=a_1=1. Let v_n=V_n with v_0=v_1=1.
Their ordinary generating functions, understood formally, are

    sum a_n z^n =
      1/(1-z) + [z^3 + z^4/((1-z)(1-2z))]/(1-z)^2,

    sum v_n z^n =
      1/(1-z)
      + [z^3 + z^4(9/(1-4z)-3/(1-2z)-1/(1-z))]/(1-z)^2.

They count classes and tilings respectively, not the same objects.
The stratum has A_n~2^(n-1) classes and V_n~4^(n-2) tilings.
Its exact boundary budget is a constraint on the full metagraph,
not a uniqueness theorem for its global wiring, a bridge classification,
or a shortcut to unrestricted tournament enumeration.

