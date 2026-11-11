# Square-path ports: a sharp Gram sector, its first negative seam, and a primitive-triple census compiler

**Session257. PROVED IN THIS NOTE** for the elementary statements below;
**FINITE-EXACT** for the saved controls. No literature-priority claim, no
canon promotion, and no Chan–Pak coverage theorem. The square-sum/Pythagorean
path theorem and the sorted-spinor tournament formula are inherited, explicitly
named inputs. New work is the restricted gluing geometry, its sharp realized
probe bank and negative boundary, and the color-census decoder connecting them
to the inherited tournament carrier.

## 0. Inheritance and live concepts

- Closest proved source: `01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md`,
  Corollary A3''': primitive `s²+t²=u²`, `s<t`, supplies a path on `[t²−1]`
  using exactly the three reflection matchings with targets `s²,t²,u²`.
- Closest response mechanism: Session256 `gram/PROOF.md`, whose cyclic
  tournament Hamiltonian-count kernel is positive with a finite profile bank.
- Exact older square-sum ports: `05-knowledge/results/session32_square_monster_20261023/ports/PORT_TRANSFER.md`;
  contracting the branched 3–4–5 blocks loses Hamiltonian feasibility.
- Recovered near miss: `05-knowledge/results/session40_marked_ports_square_sum_20261023/square-sum/README.md`,
  especially §6. Its endpoint/two-path birth polynomial is already established.
  Endpoint signatures alone cannot handle vertex arrival: Q24 has no paths,
  but two-component states create all ten Q25 paths. We do not rebrand that
  polynomial or its connectivity sidecar as new.
- Least-used sidecar here: the NUMBER OF SEAM EDGES and the PAIRING of exposed
  ports, rather than just their set or color counts.

Concept board: primitive reflection chains; named endpoint response vectors;
seam-connectedness pairings; positive tournament kernels; spinor graph compilers.
The anchor is a lawful square-sum/Chan–Pak bridge, the niche is the exact seam
positivity boundary, and the wildcard is decoding a primitive triple from
an abstract edge-color census.

**Two different colorings must not be conflated.** Sections 4–5 use colors
meaning the *three ordered square targets*. These are not the user's recovered
three-color partition of Q25 edges (first part containing 23–13,23–2,25–11,
25–24,15–1,15–10,17–8,16–9). Nor are they the Q12 component colors. The Q12
component sizes 3,4,5 alone remain insufficient for a structural triple map.

## 1. Exactly two matching seam edges give a positive kernel

Fix a finite set B of k external port names. A graph G has at least two
vertices and an injectively port-labelled subset of its vertices; it need not
use every name. Arithmetic vertex labels, when present, are separate data.
For distinct ports a,b let h_G({a,b}) be the number of **unoriented** spanning
Hamiltonian paths with those endpoints, zero if either port is absent.

Take disjoint copies of G and H. For each common port a add one seam edge
between its two copies. Keep all original edges. Let C₂(G,H) count unoriented
Hamiltonian cycles in this glued graph that use **exactly two seam edges**.
Then

    C₂(G,H) = Σ_{ {a,b} ⊂ B } h_G({a,b}) h_H({a,b}).             (1)

Proof: deleting the two seam edges leaves one spanning path in each half.
They have the same two port names as endpoints. Conversely such a pair of
paths has a unique completion using the two prescribed seam edges. This is
an edge-set bijection, so no factor2 from orientations is present. Both halves
have at least two vertices, excluding degenerate parallel-edge cycles. QED.

Thus C₂ is a positive-semidefinite kernel, with

    C₂(G,H)² ≤ C₂(G,G) C₂(H,H).                                (2)

Equality holds iff the two endpoint vectors are linearly dependent (including
zero vectors). For two nonzero nonnegative vectors this means positive
proportionality. Its rank is at most binom(k,2).

This is an ACTUAL cycle-count interpretation, not positivity assigned to an
arbitrary numerical feature vector. But the glued graph need not be another
square-sum graph, and C₂ is not the total cycle count without the seam filter.
Those are mandatory type boundaries.

## 2. Primitive reflection chains realize a sharp finite probe bank

For every k≥2, choose one primitive triple s<t with t²−1≥k. Let P be its
three-target reflection graph on `[t²−1]`. By the inherited THM-4510 result,
P itself is a path, and hence has exactly one unoriented Hamiltonian path.

For each unordered pair {a,b}⊂B, give its two endpoints those external port
names; assign all other k−2 names to distinct interior vertices. Call the
result P_ab. Then

    h_(P_ab) = e_ab,       C₂(G,P_ab)=h_G({a,b}).                (3)

Therefore the rank bound binom(k,2) is sharp **even with every probe realized
by the same primitive square-reflection chain, differently port-labelled**.
The binom(k,2) probe responses recover the entire endpoint signature. Two
sources are indistinguishable by every two-seam gluing context iff these
finite responses agree. The dimension is sharp on realized sources, not
merely on a formal vector space.

Port renaming is part of the construction. We are not claiming those renamed
ports remain the arithmetic labels of the square-sum vertices. A rule insisting
that all ports are their arithmetic integers is a narrower context family and
does not inherit this sharpness proof.

For fixed arithmetic port names 1,...,25, the engine finds:

- Q15, Q16, Q17 have single endpoint pairs {8,9}, {8,16}, {16,17};
  their two-seam Gram matrix is I₃, although all have scalar path count1.
- Q23 has three endpoint pairs {2,18},{9,18},{18,22}.
- Q25 has ten, each with multiplicity1, containing those three.
  Their 2-by-2 Gram block is [[3,3],[3,10]], determinant21.

The occurrence of 21 here is an exact determinant, NOT a map to the
Chan–Pak forbidden count21. Its mechanism is ordinary support overlap.

## 3. Four seams destroy positivity: the minimal pairing obstruction

Take four named ports and three source graphs, each a matching of two edges:

    M₁ = {12,34}, M₂ = {13,24}, M₃ = {14,23}.

Each has zero Hamiltonian-path endpoint signature. Glue as above and require
four seam edges. Equal matchings give two disjoint 4-cycles, so no Hamiltonian
cycle. Different matchings give one 8-cycle. Consequently the actual four-seam
Hamiltonian-cycle matrix is

                  [0 1 1]
    C₄(M_i,M_j) = [1 0 1] = J−I.                             (4)
                  [1 1 0]

Its eigenvalues are 2,−1,−1; vector (1,−1,0) has quadratic value−2.
Already the principal 2-by-2 block has determinant−1. The unrestricted
Hamiltonian-cycle gluing kernel is also indefinite on these examples,
because no other seam count can produce a Hamiltonian cycle here.

These source blocks are not exotic: each is the **same induced square-sum
graph on {1,3,7,9}**, whose only edges are 1–3 and 7–9, with three port maps:

    ports (1,2,3,4) -> arithmetic (1,3,7,9), (1,7,3,9), (1,7,9,3).

Thus the obstruction occurs in the port-decorated square-sum block category
itself. It does not assert that the assembled seam edges are square sums.

Why does it fail? With two seams, each half must be one path, and the two
signatures pair by equality. With four seams, each half has two path pieces;
their endpoint pairings must be DIFFERENT for the union to be connected.
Connectivity contributes an off-diagonal compatibility matrix, not a positive
identity matrix. The older port transfer retains exactly the kind of ownership
information needed to see this distinction. The Session256 tournament Gram
proof consequently does not transfer to unrestricted square-sum gluing by
replacing directed profiles with undirected ones.

One seam is impossible across a cut; two is the first nonzero sector and is
positive; four is the first possible higher even sector and already fails.
That is the precise meaning of “first negative seam,” not a census of every
conceivable gluing operation.

## 4. The primitive triple is recoverable from three edge-color counts

Let P(s,t) be the inherited reflection path on `[t²−1]` for a primitive triple
s²+t²=u² with s<t. Color edges 0,1,2 by their ordered targets s²,t²,u².
Let (a,b,c) count edges of the three colors. Then

    ε := b−a−c is either0 or1,
    s² = 2a+1+ε,                  t² = 2b+2−ε.                (5)

Moreover ε=0 iff s is odd, and ε=1 iff s is even.

Proof: put A=s², B=t², N=B−1. The matching of target A has floor((A−1)/2)
edges; that of target B has floor((B−1)/2). For target A+B, eligible vertices
range from A+1 through B−1. Because A+B is odd, there is no fixed point, so
its edge count is (B−A−1)/2. If A is odd and B even, the counts are

    ((A−1)/2, (B−2)/2, (B−A−1)/2),

whereas if A is even and B odd they are

    ((A−2)/2, (B−1)/2, (B−A−1)/2).

Substitution gives (5). QED.

Examples:

    triple (3,4,5)     -> census (4,7,3),   ε=0;
    triple (5,12,13)   -> census (12,71,59), ε=0;
    triple (8,15,17)   -> census (31,112,80),ε=1.

This is a decoder from an abstract three-colored graph invariant to its
primitive triple, provided the source is in the certified reflection-chain
family. No original arithmetic vertex labels or endpoint values are needed.
The order of the three color names IS needed.

A proposed census can be tested without a path search: compute ε,A,B, require
ε∈{0,1}, A,B,A+B positive squares with smaller roots s<t and gcd(s,t)=1.
These checks are necessary and sufficient for that census to belong to SOME
primitive reflection chain: necessity is above; sufficiency constructs the
inherited chain and its counts. This does not certify the proposed colored
path itself. The engine swaps two colors of P(3,4), preserves its census,
and finds no starting arithmetic label whose reflection itinerary visits
1,...,15 exactly once. Its full witness is in results.json.

Thus color census preserves the primitive parameter; it loses the actual
color itinerary and arithmetic embedding. Reconstructing the canonical
reflection carrier from the recovered triple repairs existence of a source,
not equality with an arbitrary submitted source.

## 5. Compile the decoded triple into actual Chan–Pak source tournaments

Given the triple recovered in §4, let o be its odd leg. Define

    m = sqrt((u+o)/2), n = sqrt((u−o)/2), x=m−n, y=n.

For primitive triples these are positive integers, m>n, gcd(m,n)=1, and
m−n is odd. Form the decorated tournament

    T = C₃[TT_x, TT_y, K₁].

This is a directed cyclic triangle with transitive blocks of sizes x,y,1.
The inherited Session256 `spinor/PROOF.md` §8 gives its actual Hamiltonian
path count

    H(T) = 2+(2^x−1)(2^y−1).                                  (6)

The map

    colored reflection chain -> color census -> primitive triple
      -> Euclidean spinor -> decorated tournament

is fully specified, and both endpoints recover the SAME primitive triple.
It is not a numerical identification of Hamiltonian counts: P(s,t) has one
unoriented Hamiltonian path as a path graph, while T has at least3 directed
Hamiltonian paths. The predicate preserved is primitive-triple identity and
realization by explicit graph carriers; Hamiltonian-count equality is
explicitly destroyed.

For the first three triples above the spinors are respectively (2,1),(3,2),
(4,1), giving tournament H=3,5,9. Independent subset DP checks this and all
primitive spinors with m≤10 in the finite receipt.

The three Berggren spinor branches from Session256 induce lawful operations
on these decorated reflection carriers: decode (m,n), apply the branch,
rebuild its reflection path, and compile the same new spinor to T. This is
an exact commuting construction, but it is **reconstruction**, not a bounded
local graph rewrite or a source-preserving extension. It cannot be advertised
as such. The order of the two legs must be restored after each branch before
constructing the reflection carrier.

A useful arithmetic boundary survives the bridge. Since x is odd,
2^x−1 ≡1 mod3. Thus (6) never equals1 mod3. This inherited sorted-family
obstruction holds for the entire compiled square-reflection family too.
The bridge creates genuine source tournaments, but cannot supply all
Chan–Pak targets. Restoring the nonsorted cut itinerary is the already known
Session256 way to escape this congruence obstruction; the reflection census
alone does not prescribe that itinerary.

## 6. Connections, scope, and next gates

1. **New contextual transfer:** primitive square chains supply a realized,
   sharp finite basis for endpoint-sensitive two-seam continuation. This is
   the undirected analogue of the fixed-observer response bank, with a
   different signature and a rigorously delimited gluing operation.
2. **New failure boundary:** four-seam connectivity is already indefinite
   on induced square blocks. Any broader positive kernel needs a different
   observable, a restricted port-pairing class, or an explicit correction.
3. **New arithmetic bridge:** square-target color census recovers primitive
   spinors and hence actual cyclic source tournaments, with its retained
   predicate and missing itinerary stated explicitly.
4. **OPEN:** find a richer source operation that couples this arithmetic
   recovery to a useful Hamiltonian-count transformation such as3H+2 or3H+4.
   The present compiler does not do that, and the seam Gram theorem supplies
   no such formula.
5. **OPEN:** decide which restrictions on multi-seam boundary pairings make
   their connection kernel positive. Replacing connectedness with equality
   trivially creates a Gram matrix, but would count the wrong objects.

No sporadic-group, supersingular-prime, permutohedron, LRC, or universal
ternary-coverage consequence is inferred. Other lanes can use this precise
port interface, but must preserve its naming, seam count and connectivity
sidecars when doing so.

## 7. Reproduction

Run `python3 -B check.py` in this directory. Standard library only.
All validation uses explicit checks raising RuntimeError, not optimizable
assertions. Normal and `python3 -O -B check.py` runs were verified to produce
byte-identical results.json and standard output.
`results.json` records:

- exact endpoint signatures for Q12,Q15,Q16,Q17,Q18,Q23,Q24,Q25;
- 64 independent direct Hamiltonian-cycle recounts, all ordered pairs of
  simple graphs on three named ports, verifying the two-seam identity;
- all nine four-seam matching counts and their arithmetic square embeddings;
- all31 primitive spinors 2≤m≤12, 1≤n<m, gcd=1, opposite parity: explicit
  reflection-chain construction, edge-color census and spinor recovery;
- independent tournament subset DP for the subrange m≤10;
- the equal-census, unrealizable-itinerary hostile.

These controls verify the new claims, not the whole inherited theorem corpus.
The all-size proofs are above; no finite scan is promoted to an infinite
coverage claim. All files in this lane are new; the dirty math checkout was
read only, and publication is the parent session's responsibility.
