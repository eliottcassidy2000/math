# Why square-sum cycles start at 32, not 31

Session 32 bounded extension, 2026-10-23 (UTC).

Let G_n be the **simple undirected** graph on {1,...,n}, with distinct
vertices adjacent exactly when their sum is a square. A cycle has at least
three vertices. A spanning 2-factor is a spanning degree-two subgraph, not
necessarily connected. Cycles and paths are counted as undirected edge sets,
so reversal and cyclic rotation do not create additional objects.

## Status, inheritance, and scope

**PROVED (elementary):** G31 has exactly two spanning 2-factors, of component
sizes 11+20 and 14+17. Hence it is not Hamiltonian, despite minimum degree
two. A displayed cycle proves G32 Hamiltonian. Together with the elementary
low-degree obstructions below, this proves that the first Hamiltonian cycle
occurs at 32.

**FINITE-EXACT:** the independent certificates additionally show that G32 has
exactly five spanning 2-factors and a unique Hamiltonian cycle. Equivalently,
G31 has exactly one Hamiltonian path with prescribed endpoints 4 and 17.
The complete G32 2-factor types are 3+6+23, 3+29, 5+27, 9+23, and 32, each
occurring once. The uniqueness claim is a bounded exact computation, not a
new all-order construction.

**INHERITED CITED, not reproved here:** the existing square-sum report cites
Gerbicz for Hamiltonian cycles at every n>=32. The local threshold proof
below does not establish persistence for larger orders and does not audit
his base/gluing tables. The current extension needs no further literature
claim.

Closest mechanism: the endpoint-aware forced-chain reduction and the marked
ear correspondence in `../session32-square-sum/README.md`, §§3–5. Canonical
hostile: a degree-two vertex can be an endpoint of a Hamiltonian path, so
unconditional degree-two forcing is valid here for **cycles**, not for
arbitrary paths. Corrected near miss: suppressing chains while forgetting
which endpoint port an optional edge actually uses. Least-used sidecar:
residual degree at each named port, rather than only total degree of the
contracted block.

Concept board: leaf window; forced chains; subdivided K4; port compatibility;
marked-ear arrival. The boundary mechanism is not a change in arrival
degree: both 31 and 32 arrive with degree two. The extra ear at 32 exposes
the specific old endpoint pair {4,17}.

## 1. The leaf obstruction ends at 31

For 18<=n<=30, vertex 18 has the single neighbour 7: among possible sums
19,...,18+n<=48, only 25 gives a distinct positive summand at most n.
At n=31, the second neighbour 31 arrives because 18+31=49.

Earlier orders have an immediate degree-at-most-one obstruction too:

- n=1,2: vertex 1 has degree zero.
- 3<=n<=6: vertex 2 has degree zero.
- n=7: vertex 2 has degree one.
- 8<=n<=16: vertex 8 has the single neighbour 1.
- n=17: vertex 16 has the single neighbour 9.

Every vertex of G31 has degree at least two, as follows either from its
listed edges or the inherited degree formula

    deg_Gn(x) = floor(sqrt(x+n)) - floor(sqrt(x)) - 1_{2x is a square}.

Thus 31 is the first order where the elementary minimum-degree condition
allows a cycle. It is not the first order where connectivity of a spanning
2-factor can actually be achieved.

## 2. Exact forcing at 31

Assume a spanning 2-factor C exists. Start with every square-sum edge allowed
and no forced edges. Two implications are valid:

1. If a vertex has exactly two allowed incident edges, both are forced.
2. If a vertex has two forced incident edges, every other edge there is
   excluded.

These implications preserve the entire set of spanning 2-factors; they make
no connectedness assumption. The full graph has 44 edges. The following
short trace reaches a fixed point:

1. Original degree-two vertices 16,17,18,25,26,27,28,29,30,31 force both
   incident edges. There are 19 distinct forced edges at this stage.
2. Saturated vertex 7 excludes 2–7 and 7–9; saturated 8 excludes 1–8;
   saturated 19 excludes 6–19; saturated 20 excludes 5–20.
3. Vertex 2 now has only neighbours 14,23, so forces both edges. This
   saturates 23 and excludes 13–23. Vertex 13 now has only neighbours 3,12,
   so forces both of those edges too.

The result is 23 forced edges and 38 allowed edges, of which 15 are optional.
The machine certificate records every individual inference and edge.

The forced subgraph has these eight path components, including singletons:

    A = [1]
    B = [3,13,12]
    C = [4]
    D = [5,31,18,7,29,20,16,9,27,22]
    E = [6,30,19,17,8,28,21]
    F = [10,26,23,2,14]
    G = [11,25,24]
    H = [15].

A nontrivial forced path needs **one** optional edge at each of its two
specified endpoint ports. A singleton needs two optional edges. Internal
vertices of the forced paths are saturated. These are necessary and
sufficient local degree conditions; connectivity remains a separate test.

The optional edges are

    1–3, 1–15, 1–24,
    3–6, 3–22,
    4–5, 4–12, 4–21,
    5–11, 6–10, 10–15, 11–14, 12–24, 14–22, 15–21.

This is an exact eight-block kernel with original endpoint ports retained.
Expanding each block recovers a unique selected spanning subgraph. It is
not legitimate to replace the two port obligations of a nontrivial block
by the weaker condition “two incident edges anywhere on the block.”

### A concrete hostile to forgetting ports

The bare contracted graph has the spanning cycle

    A–B–C–D–G–F–E–H–A.

But at D the edges C–D and D–G both use original port 5, leaving port 22
unused. At G the edges D–G and G–F both use port 11, leaving port 24 unused.
The quotient cycle cannot expand into a degree-two spanning subgraph.
Thus ordinary Hamiltonicity of the bare eight-vertex quotient gives a
false positive. Retaining the labelled ports repairs exactly this failure.

## 3. A subdivided K4 proves there are only two completions

Look at the optional-edge graph on the ten exposed path ports and three
singletons. Its degree-three junctions are 1,3,4,15; all other vertices
have optional degree two. It is the following subdivision of K4:

    a: 1–3
    b: 1–15
    c: 1–24–12–4
    d: 3–6–10–15
    e: 3–22–14–11–5–4
    f: 4–21–15.

Every internal vertex of these corridors has residual demand one, so
selected and unselected edges alternate along each corridor. Let a,b,c,d,e,f
be the 0/1 indicators of the first displayed edge of each corridor.
For the odd-length c,d,e corridors, the first and last indicators agree.
For the length-two f corridor, the last indicator is 1−f.

The junction residual demands are two at 1,4,15 but **one at 3**, because
3 already uses the forced edge 3–13. Therefore all possibilities are
exactly the binary solutions of

    a+b+c = 2,             (junction 1)
    a+d+e = 1,             (junction 3)
    c+e+f = 2,             (junction 4)
    b+d+1−f = 2.           (junction 15)

From the last equation, f=b+d−1. Substituting into the third gives
b+c+d+e=3. The sum of the first two equations is
2a+b+c+d+e=3, so a=0. Consequently b=c=1, d+e=1, and f=d.
There are exactly two solutions:

    (a,b,c,d,e,f) = (0,1,1,1,0,1) or (0,1,1,0,1,0).

This is an exhaustive algebraic proof, not an inference from the census.
Expanding gives the following two 2-factors; each displayed line closes
back to its first entry.

For d=1, the component lengths are 11 and 20:

    3,6,30,19,17,8,28,21,4,12,13;
    1,15,10,26,23,2,14,22,27,9,16,20,29,7,18,31,5,11,25,24.

For d=0, the component lengths are 14 and 17:

    3,13,12,4,5,31,18,7,29,20,16,9,27,22;
    1,15,21,28,8,17,19,30,6,10,26,23,2,14,11,25,24.

Both are disconnected. Hence **G31 has no Hamiltonian cycle**.

## 4. The new ear 4–32–17 supplies the missing endpoint pair

Vertex 32 has exactly two neighbours, 4 and 17, with sums 36 and 49.
Every spanning cycle of G32 must use both edges 4–32 and 32–17.
Deleting vertex 32 therefore gives a spanning path of G31 with endpoints
4 and 17. Conversely, adjoining that ear to any such spanning path closes
it into a Hamiltonian cycle. These maps are inverse, giving an exact
bijection of undirected edge sets.

An explicit such path of G31 is

    4,21,28,8,1,15,10,26,23,2,14,22,27,9,16,20,29,7,18,31,
    5,11,25,24,12,13,3,6,30,19,17.

Joining its endpoints through 32 gives the cycle

    1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,
    29,20,16,9,27,22,14,2,23,26,10,15,1.

It visits every integer 1,...,32 once before returning, and every adjacent
sum is a square. This alone is a positive certificate at 32.

For an explicit local repair of the first G31 2-factor (11+20), delete

    1–24, 4–12, 8–17,

and add

    1–8, 12–24, 4–32, 17–32.

The result is exactly the displayed 32-cycle. In particular, the arrival
at 32 frees vertex 17 from its formerly forced edge 17–8, making the
previously excluded connection 1–8 usable. It does not simply insert a
vertex into an already existing G31 Hamiltonian cycle, since there is none.

The small degree solver finds five G31 spanning subgraphs with degree one
at 4,17 and degree two elsewhere. Only one is connected. Under the ear
bijection, these five subgraphs correspond to the five 2-factors of G32;
the disconnected cases have a path plus one or more disjoint cycles.
The single connected case establishes the finite-exact uniqueness claim.

The mechanism continues the earlier marked-ear account at 25, where the
compulsory ear was 11–25–24. Here the relevant old state is specifically
the Hamiltonian-path endpoint pair {4,17}; the total number of old paths
or the arrival degree by itself does not identify that compatibility.

## 5. Reproduction, completeness, and limits

From this directory:

    python3 check_cycle_threshold.py --output certificate.json
    python3 -O check_cycle_threshold.py --output certificate.optimized.json
    cmp certificate.json certificate.optimized.json

The checker uses Python's standard library and exact integers only. Its
universe is all simple square-sum edges on the specified vertex sets.
The inherited filters are just the two proved degree-forcing rules. There
is no premature proper-cycle rejection: disconnected 2-factors are kept
and their component lengths are reported.

- At 31, the producer examines exactly C(15,8)=6,435 optional-edge subsets,
  checks degrees, and recovers the same two solutions as the binary proof.
- At 32, the degree-obligation solver exhausts its finite branching tree
  in nine nodes. It enumerates five 2-factors, one connected.
- For G31 with endpoints 4,17, the same general solver uses nine nodes;
  five degree-feasible subgraphs survive, one connected.
- Positive control: every displayed cycle/path is expanded and checked.
- Hostile controls: both disconnected G31 completions survive; the bare
  quotient has a spurious spanning cycle; degree-at-most-one witnesses
  exclude every earlier order. The checker also verifies the exact
  insertion/deletion bijection and the three-edge/four-edge local surgery.
- No assertion can disappear under optimization: all checks raise explicit
  exceptions. Normal and optimized runs produce byte-identical certificates.

The certificate records graph edges, forcing traces, labelled forced chains,
port-preserving kernel edges, all cycle-cover components, the G32 cycle,
the prescribed-endpoint path, and the checker SHA-256. This is a bounded
threshold result, not an all-n Hamiltonicity theorem. No Git state was
changed by the producer.
