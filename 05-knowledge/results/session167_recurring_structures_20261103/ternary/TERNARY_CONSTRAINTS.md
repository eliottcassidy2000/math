# A profile-invisible ternary switch: capacity constraints, Ramsey recolorings, and cellular clocks

**Status.** Exact constructions and elementary proofs, with self-contained finite certificates. The seed triangle and its minimal tournament order are independently checkable by exhaustive finite enumeration. The arbitrary-size statements follow from the proofs. This packet does not improve a Ramsey bound, solve Friedman's fence-density problem, or classify all tournament metagraphs. No literature-priority claim is made.

The useful common object is a **finite state space with coordinate changes**. A newly certified three-state tournament switch lets us preserve those changes while making every path-cover count identical. Fence capacities determine one allowed-state set, Ramsey exclusions determine another, and a cellular automaton supplies an evolution on the states. These are exact constructions, not identifications of the original problems.

## 1. Definitions and a three-state switch

Let M_n be the loopless graph of isomorphism classes of n-vertex tournaments, adjacent when some representatives differ by reversing one arc. Let S_n be the induced subgraph of strong tournaments. These are the single-arc-flip class graphs, **not** the repository's blue/black complement-line metagraph.

For a tournament T, let p_j(T) count its partitions into j nonempty directed paths, with a directed order on each path but with the paths themselves unordered. Write

    P_T(z) = sum_j p_j(T) z^j,
    F_T(r) = sum_j p_j(T) (r)_j.

Here (r)_j is the falling factorial, and p_1(T)=H(T). Equivalently F_T(r) counts r labeled path slots, allowing empty slots. The same definition works at the empty tournament, with F_empty(r)=1.

Encode a tournament on {0,...,5} by the lexicographically ordered pairs 01,02,...,45, with bit 1 meaning i→j for i<j. Set

    A_0 = mask 153,    A_1 = mask 177,    A_2 = mask 181.

**Theorem 1 (the ternary seed).** These are three pairwise nonisomorphic strong tournaments. They induce K_3 in S_6 and have the same complete profile

    P(z) = 29z + 131z² + 162z³ + 77z⁴ + 15z⁵ + z⁶,
    F(r) = r⁶ + 12r⁴ + 16r².

Tournament order six is minimal for a triangle all of whose vertices have the same complete profile, even without requiring strong connectivity.

**Finite proof certificate.** `check.py` constructs the three literal tournaments, checks reachability, and computes their profiles by two separate algorithms: Hamilton counts on subsets followed by set-partition DP, and permutation/cut enumeration. For each pair it checks all 720 relabellings and finds minimum quotient Hamming distance exactly one. Here are witnesses; permutations map old labels to new labels:

| Source mask | Target mask | Relabelling of target | Relabelled target | Flip |
|---:|---:|---|---:|---|
| 153 | 177 | (4,2,1,5,0,3) | 8345 | {3,5} |
| 153 | 181 | (0,4,2,1,3,5) | 185 | {1,2} |
| 177 | 181 | (0,1,2,3,4,5) | 181 | {0,3} |

For minimality, the same program independently generates **every** labelled tournament for n=1,...,5 and partitions the masks into full relabelling orbits. The class counts are 1,1,2,4,12. It computes every class profile and every one-flip class edge. There are no equal-profile edges through order four; at order five the only such edge joins masks 8 and 10. In particular there is no equal-profile triangle. The minimality claim depends on this complete finite certificate, not a sampling argument.

### Why these are classes rather than three fixed labelled bitstrings

The labelled tournament flip graph is a cube and cannot contain a triangle. Our triangle necessarily uses relabellings. A literal three-flip lift is

    153 → 8345 → 8337 → 8369,

where the successive classes are A_0,A_1,A_2,A_0. The last mask is obtained from 153 by the odd permutation (4,2,1,5,0,3). Thus the switch has precisely the odd relabelling holonomy discussed in S166. It cannot be advertised as three fixed labelled tournaments that are pairwise one bit apart.

## 2. The ternary constraint compiler

The ordinal sum T⊕U orients every cross arc from T to U. Restrict a path in T⊕U to the two blocks: it is a T-path followed by a U-path, allowing either part empty. Doing this in each labeled slot gives

    F_(T⊕U)(r) = F_T(r) F_U(r).

For a word x=(x_1,...,x_k) in {0,1,2}^k, define

    T_x = A_(x_1) ⊕ ... ⊕ A_(x_k).

**Theorem 2 (exact compiler).** For every k≥1 and every set X⊆{0,1,2}^k, the map x↦[T_x] is an induced embedding of the graph on X in which words are adjacent exactly when they differ in one coordinate. Every image has order 6k and

    F_(T_x)(r) = (r⁶+12r⁴+16r²)^k,
    H(T_x) = 29^k.

Consequently every coefficient of the full path-cover profile is independent of x, even after imposing arbitrary restrictions on X.

**Proof.** The SCCs of T_x are exactly its k strong six-vertex blocks, in their uniquely ordered condensation. An isomorphism therefore preserves each block position, proving injectivity. A reversal within one block gives a member of the displayed family only if it changes that one factor to another seed class. A reversal between two blocks merges all SCCs in the interval between them: the new backward arc, the unchanged forward cross arcs, and internal strong paths make that interval mutually reachable. It therefore leaves the family. Theorem 1 makes every distinct pair of values at a coordinate adjacent. This proves the induced-graph statement, including absence of extra edges. The profile identities follow from the slot product. ∎

This is the S166 Cartesian-product mechanism with a new alphabet: the earlier invisible edge yields binary cubes; the triangle yields **ternary Hamming graphs** K_3^□k. It is not a universality theorem for arbitrary abstract graphs.

### Strong outputs, with an explicit order bound

The same graph embeds inducedly in S_(30k+10), still with one common complete path-cover profile. For completeness, the S166 wrapper uses n=6k, q=2n+5, regular q-vertex guards B,C, and cyclic cross orientations

    B → T_x → C → B.

The guard-B degrees equal 2n+2, core degrees lie in [2n+5,3n+4], and guard-C degrees equal 3n+7. These bands are separated by at least three. An isomorphism preserves the bands. After one arc reversal, a touched guard has a degree different from its required constant, and the gaps still identify the bands, so only core reversals stay in the wrapper family. This proves inducedness.

Profile equality survives replacing a module by a same-profile module: cut every global path at exits and entries of the module, match its local j-path covers bijectively, and reconnect their constituent paths to the same external segments. Uniform external adjacency makes every reconnection legal, and the inverse matching reverses it. The global number of paths and absence of cycles are unchanged. Apply this separately in each seed slot. The wrapped Hamilton count is common but is **not** asserted to be 29^k.

## 3. Capacity constraints: retaining the fence bins

The per-fence accounting notes show why angle-only relaxations lose information: the pieces on each side of a unit fence must sum to one. A finite list of candidate pieces with nonnegative resource requirements gives an exact downward-closed **capacity model**:

    Δ = { I⊆[k] : sum_(i in I) w_(b,i) ≤ 1 for every bin b }.

For literal predetermined fence-side pieces, w_(b,i) is their assigned length in bin b. Several bins and pieces used in more than one resource constraint are allowed. This model retains the selected capacities; it does not certify polygon closure, junction angles, a plane embedding, or a full fence configuration.

More generally let Δ be any nonempty downward-closed family of subsets of [k], and let

    f_Δ(t) = sum_(I in Δ) t^|I|.

Allow a word x exactly when its nonzero support belongs to Δ; values 1 and 2 are two tags for an occupied candidate. No geometric interpretation of the tags is assumed.

**Theorem 3 (capacity reconfiguration).** The allowed words form an induced tournament-class graph as in Theorem 2, with

    number of vertices = f_Δ(2),
    number of edges = 3 f_Δ'(2).

Its graph distance between two states is their Hamming distance. If at least one singleton is allowed, its chromatic number is exactly three.

**Proof.** Each support I has 2^|I| tagged states. From any state one can remove each occupied coordinate, giving 2f_Δ'(2) edges in total, counted once at their larger support. Recoloring 1↔2 within a fixed support gives f_Δ'(2) further edges. This proves the edge formula.

To change x to y, first remove the coordinates occupied only in x, recolor coordinates occupied in both as necessary, and then add coordinates occupied only in y. Every intermediate support is a subset of the initial or final support, so is allowed. Each differing coordinate changes exactly once. No path can use fewer coordinate changes, proving Hamming isometry **within the displayed graph**. The sum of coordinates modulo three is a proper three-coloring. An allowed singleton together with its empty and two tagged states gives a triangle. ∎

The isometry statement does not exclude a shortcut through other tournament classes outside this displayed induced subgraph.

### Why pairwise conflicts are insufficient

Three pieces of length 2/5 assigned to one unit bin fit pairwise but not all together. Their capacity family has all proper subsets of {1,2,3} and excludes the full set. Here

    f(t)=1+3t+3t²,    |V|=19,    |E|=45.

A pairwise conflict graph would incorrectly allow all three and give the full 27-state ternary cube. The extra eight forbidden words are precisely those with all three pieces occupied. This example captures the higher-order bin constraint that a conflict graph alone forgets; it does not settle the geometric fence problem.

### A closed form for a strip model

If consecutive positions cannot both be occupied, Δ is the independence complex of a k-vertex path. Its polynomial is

    f_k(t)=sum_(j=0)^floor((k+1)/2) binom(k-j+1,j)t^j,
    f_k(t)=f_(k-1)(t)+t f_(k-2)(t).

Thus the profile-invisible tournament graph has

    a_k=f_k(2)=(2^(k+2)-(-1)^k)/3

states and e_k=3f_k'(2) edges. In particular

    a_0=1, a_1=3;  a_k=a_(k-1)+2a_(k-2),
    e_0=0, e_1=3;  e_k=e_(k-1)+2e_(k-2)+3a_(k-2).

For k=1,...,8 the state counts are 3,5,11,21,43,85,171,341. `check.py` verifies every state, edge, and pairwise distance in these finite controls, and separately checks the formulas and all distances for all 19 nonempty downward-closed families on a three-element ground set. This is a discrete capacity-strip construction, not a new fence-packing density estimate.

## 4. The (3,3,5) and (3,3,6) Ramsey spaces transfer exactly

Fix n≥2 and put k=binom(n,2), indexing the **labelled** edges of K_n. Let X_(n,s) consist of all three-color words with no color-0 triangle, no color-1 triangle, and no color-2 K_s. Its one-coordinate graph is exactly the usual admissible single-edge-recoloring graph for the (3,3,s) constraint.

**Corollary 4.** For s=5 or 6, this entire admissible recoloring graph embeds inducedly in M_(6 binom(n,2)), and in S_(30 binom(n,2)+10), with the complete path-cover profile constant throughout each displayed family. In the unwrapped version every state has H=29^binom(n,2).

This transfers the actual constraints and recoloring adjacency. It does not produce an admissible coloring where none was known, identify colorings under permutations of K_n's vertices, or improve R(3,3,5) or R(3,3,6).

Two exact enumerators provide nonextremal controls. One enumerates disjoint triangle-free red and blue edge sets; the other independently extends every labelled five-vertex coloring by its last vertex:

| n | Edge slots | Admissible (3,3,5) | Admissible (3,3,6) |
|---:|---:|---:|---:|
| 3 | 3 | 25 | 25 |
| 4 | 6 | 543 | 543 |
| 5 | 10 | 29,196 | 29,197 |
| 6 | 15 | 3,730,190 | 3,731,612 |

A separate inclusion-exclusion control explains the n=6 difference 1422. Among red/blue-triangle-free colorings, a specified green K_5 allows arbitrary colors on the five edges to the remaining vertex, giving 3^5 choices. The intersection of two distinct specified green K_5 events allows three choices for their one uncovered edge; three or more such events force all edges green. Thus their union has

    6·3^5 - binom(6,2)·3 + binom(6,3)-binom(6,4)+binom(6,5)-binom(6,6)
    = 1423

colorings. The (3,3,6) condition excludes just the all-green member, so the difference is 1422.

For current external context, Radziszowski's April 2026 survey gives 45≤R(3,3,5)≤57, and Boza's March 2026 paper proves R(3,3,6)≤91. These are cited background, not results of this packet. [Radziszowski, §6.1.2](https://www.cs.rit.edu/~spr/ElJC/sur.pdf); [Boza, arXiv:2603.10851v2](https://arxiv.org/abs/2603.10851v2).

## 5. An exact reversible cellular clock inside the same profile fiber

On a periodic row of k ternary sites, let S be the cyclic shift and define the radius-one linear cellular automaton over F_3 by

    Φ(x)_i = x_i + x_(i+1) mod 3,    Φ=I+S.

**Theorem 5.** If k=3^a, a≥0, then Φ^k=-I, and Φ has exact order 2k. The configuration with one 1 and all other entries 0 has exact period 2k. The tournament states T_(Φ^t x) all retain the same complete path-cover profile; for the one-1 orbit they give 2k distinct classes, all with H=29^k.

**Proof.** In characteristic three, repeated cubing gives

    (I+S)^k = I+S^k = 2I = -I.

Thus Φ^(2k)=I and Φ is invertible. For k=1 the order is two. For k≥3 put q=k/3. Then

    Φ^(2q) = (I+S^q)^2 = I - S^q + S^(2q).

Applied to the one-1 configuration, this is not the original configuration because the three indicated sites are distinct. Its period therefore does not divide 2k/3; it also does not divide k because the k-th iterate is its negative. Every proper divisor of 2·3^a divides one of those two numbers, so its period is exactly 2k. The global map consequently has that order too. Theorem 2 supplies the tournament encoding and profile constancy. ∎

A synchronous CA step can change several coordinates. It is **not** generally one arc-flip edge. It has a path of exactly d coordinate changes in the encoded Hamming graph, where d is the number of changed sites, by updating those sites one at a time. Those intermediate updates require an external schedule; this construction does not claim an intrinsic autonomous tournament rule. Within a block, the physical representative for a transition may require relabelling, as §1 explains.

The exact observation is stronger than an analogy: arbitrarily long reversible local clocks exist in a prescribed constant-profile family. The path-cover profile records none of their phase. This does not make their states isomorphic, and it does not imply that every tournament dynamic is a cellular automaton.

## 6. What the numbers 7,28,126 contribute

The inherited S158 carrier is an unordered set partition in which a block is either a singleton or has an unordered pair of designated terminals. Its EGF is

    exp(x + x² exp(x)/2),

and its counts on 3,4,5 labels are 7,28,126. The new verifier checks them by the direct distinguished-block recurrence.

For a given tournament, weight a nonsingleton terminal block by the number of its directed Hamilton paths having that terminal pair. The total weight of carriers with j blocks is exactly p_j(T). This is a genuine combinatorial route from the carrier to our interface, but **the weights are essential**: 7,28,126 themselves are neither the seed's Hamilton counts nor a new tournament spectrum law.

The new point is that this entire weighted carrier profile is simultaneously equal on all three seed states. It remains equal under every ordinal composition and every allowed-state restriction above. Preserving the carrier totals therefore does not preserve the state, the Ramsey coloring, the selected capacity pieces, or the CA phase.

## 7. Reproduction, inheritance, and open questions

Run:

    python3 check.py
    python3 -O check.py
    g++ -O3 -std=c++17 check_ramsey.cpp -o check_ramsey
    ./check_ramsey > ramsey_counts.jsonl
    python3 independent_ramsey.py
    python3 -O independent_ramsey.py

The Python verifier uses only the standard library and explicit checks that stay active under optimization. Its `certificate.json` contains every seed witness, the complete small-order minimality summary, the common profile of all nine 12-vertex two-slot products, capacity counts, and CA clock controls through k=81. The main infinite statements have proofs above; finite checks are not substituted for those proofs. The Ramsey counts use exhaustive enumeration only through n=6.

Inherited sources:

- [S166 Cartesian flip graphs](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/cartesian_flip/NOTE.md): SCC Cartesian embedding, strong wrapper, module profile replacement, and relabelling holonomy.
- [S158 weighted carriers](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session158_creative_20261103/carrier/NOTE.md): terminal-pair carriers and the path-cover interface. The carrier definition and weighting are restated here so the argument is self-contained.
- [Per-fence accounting](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/procgen_fence2_20260926_per_fence_accounting.md): literal fence-side capacities and their geometric scope. Our capacity model takes that resource constraint as input and does not assert a geometric realization theorem.

Useful next questions are whether four pairwise adjacent strong tournaments can share a full profile, how much the strong-wrapper order can be reduced while retaining every ternary constraint edge, and which additional endpoint-resolved invariants distinguish the three states without losing practical compositionality. None is settled here.
