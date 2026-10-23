# Complete square-sum color holdouts33--36 and additive palette transport

Session42 continuation. This closes the concrete pair-only stress-test limitation of the earlier square-color packet. The original feature grammar is unchanged and its frozen files are not modified. All numerical conclusions here are finite; no surviving coloring is called an all-size decoder.

## 1. Independent complete carriers

`enumerate.cpp` independently enumerates Hamiltonian paths of Q_n by direct vertex DFS from every starting vertex. Every outgoing unused square-neighbor is explored. A completed path is retained exactly when its starting label is smaller than its ending label, so each unoriented path appears exactly once. There is no heuristic pruning or reliance on earlier path counts.

| n | all paths up to reversal | DFS states |
|---|---:|---:|
|33|669|13,834,243|
|34|4,041|47,358,366|
|35|17,175|145,789,407|
|36|12,960|272,661,152|

These independently computed counts agree with the inherited Session40 count-only census. Every path is subsequently checked in Python for complete vertex coverage, square adjacency, canonical orientation and duplicates. Exact file digests certify the regenerated carriers. Large generated path files are not included in the publication manifest; the short enumerator reproduces them.

## 2. No retuning: all1549 earlier successes face complete holdouts

The six vertex features and seven edge features, the mod2/mod3 coefficient grammars, and deduplication by global color names on the order80 domain are exactly those of the earlier frozen packet. All1584 distinct original colorings are regenerated independently. Only the1549 that separated every endpoint-marked Q32 path enter the test cohort.

A label consists of the actual endpoint pair followed by the ordered vertex- or edge-color word. This is not a proper coloring requirement. The C++ test uses exact string equality; hash-table collisions do not identify unequal strings. A separate Python implementation uses exact byte translation and set equality and rechecks every candidate's count at every held-out order.

| cohort | Q32 successes | first failure33 | first failure34 | first failure35 | first failure36 | survives33--36 |
|---|---:|---:|---:|---:|---:|---:|
| vertex2 |56|1|33|17|1|4|
| vertex3 |359|1|30|37|10|281|
| edge2 |60|1|31|11|2|15|
| edge3 |1074|1|36|46|28|963|

Total:286 of1549 earlier successes fail genuine complete holdouts;1263 survive all four orders. Only19 of116 two-color candidates survive, versus1244 of1433 three-color candidates. These counts describe this predeclared grammar, not all colorings or a probability distribution on natural mathematical models.

All1549 had separated the two old hostile paths at36 and the two at80. Nevertheless each of the four families contains a coloring that fails already at33. This demonstrates exactly why pair-only tests were insufficient.

`ledger.json` records every candidate's fixed coefficients, distinct-label counts at32,...,36, first failed order and an actual pair of path indices witnessing it. `witnesses.json` expands representative witnesses into full paths for every nonempty cohort/first-failure category.

## 3. First failures and genuine finite survivors

The vertex coloring c(v)=v mod3 passes Q32 but fails Q33: two distinct paths with endpoints18 and26 have the same ordered residue word. The edge coloring c(a,b)=a+b mod3 has the identical collision partition once the starting vertex is supplied, for the structural reason below.

The two-color vertex rule

c(v)=v+floor(sqrt(v))+ZeckendorfDigitCount(v) mod2

also first fails at33, with endpoints6 and29. Its additive edge counterpart has the same collision. Both full path pairs are retained in `witnesses.json`.

The four surviving two-color vertex formulas, in the original feature order

(v, popcount(v), floor(sqrt(v)), ZeckendorfDigitCount(v), smallestZeckendorfIndex(v), sumZeckendorfIndices(v)),

have coefficient vectors

(0,0,0,0,0,1), (0,1,0,0,0,1), (1,1,0,1,0,0), (1,1,1,0,0,0).

In particular, parity of the sum of Zeckendorf indices survives every complete carrier tested here. That is a concrete finite positive result, not a proof that it identifies arbitrary large square-sum paths. The other survivors show that this test does not uniquely select Fibonacci structure.

## 4. Exact additive vertex/edge color correspondence

Let G be any graph, let c:V(G)->A take values in an abelian group, and assign

e({u,v})=c(u)+c(v).

Given the actual starting vertex of a path, its ordered vertex-color word and its ordered edge-color word determine each other:

c(v_(i+1))=e({v_i,v_(i+1)})−c(v_i).

Thus, on any fixed family of oriented paths with actual starting vertices retained, vertex-color and induced edge-color observations have **exactly the same collision partition**. This is an all-size structure-preserving map. It explains some paired experimental results instead of treating them as independent coincidences.

For example, mod2 the square-shell label k of a+b=k² obeys k²=k. The vertex formula v+floor(sqrt(v))+ZeckendorfDigitCount(v) induces exactly the tested edge feature combination k+floor(sqrt(a))+floor(sqrt(b))+ZeckendorfDigitCount(a)+ZeckendorfDigitCount(b).

### Odd cycles distinguish two from three colors

Consider the unsigned incidence linear map c->e over a field. Its kernel satisfies c(v)=−c(u) along each edge.

- Over F2, a connected component has a one-dimensional constant kernel.
- Over a field of odd characteristic, a bipartite component has a one-dimensional alternating kernel; a non-bipartite component has zero kernel, since an odd cycle forces2c(root)=0.

Therefore on any graph with b bipartite connected components and k total connected components, the ranks over F3 and F2 are respectively |V|−b and |V|−k. Isolated vertices count as bipartite components.

For Q15 and every later Q_n, connectedness and the retained five-cycle6--3--1--15--10--6 imply that an F3 vertex coloring is uniquely recoverable from all its additive edge colors, whereas an F2 coloring retains one global flip ambiguity. An actual starting-vertex color removes the latter ambiguity for path-word transport. The statement concerns recovery from the entire edge-color assignment, not recovery of an arbitrary path from a compressed palette.

`check_additive_colors.py` independently verifies unsigned-incidence matrix ranks for Q1,...,Q80 and exact edge transport. It uses Gaussian elimination, not the bipartite formula itself, to calculate ranks.

The image of the additive vertex-color map is generally a proper subset of all edge-colorings. This theorem does not claim every procedural edge coloring has such a vertex lift. It supplies a testable common universe and identifies the information lost when endpoints or starting colors are omitted.

## 5. Reproduction and scopes

Run:

- `python3 verify.py --fresh`
- `python3 -O verify.py --fresh`
- `python3 check_additive_colors.py`
- `python3 -O check_additive_colors.py`

The fresh verifier compiles the supplied C++ enumerator with g++, regenerates all four complete path carriers, validates their digests and paths, reconstructs the original feature grammar, and independently recomputes every ledger count and witness. Paired receipts agree byte for byte. `test.cpp` is the separately implemented original C++ scoring pass; the Python verifier does not invoke it or trust its counts without recomputation.

The supplied training paths are exactly the Q32 subset of the pinned Session40/Session42 inherited carrier. Generated full holdout paths are reproducible local outputs, not missing proof assumptions. A working C++17 compiler is required for fresh enumeration.

Remaining obligation: new orders or structurally selected hostile families must be genuinely independent holdouts. No coefficients were selected or modified after observing33--36, and no extrapolation beyond36 is made.
