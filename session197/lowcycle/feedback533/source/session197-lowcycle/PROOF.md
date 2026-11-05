# Every tournament with 67 Hamiltonian paths

Session 197, 2026-11-05. **Complete all-order classification, using a published order bound and a finite exact certificate.** The full `{7,21}` Hamiltonian-path spectrum remains open. No literature-priority claim is made.

## 1. Main conclusions

Let H(T) count directed Hamiltonian vertex words, and C(T) count directed Hamiltonian cycles modulo cyclic rotation. Every tournament with H(T)=67 has exactly one nontrivial strong component. That component is one of **six rigid tournaments of order seven**. Three have C=3 and three have C=4. Every other strong component is a singleton.

Consequently:

1. No strong tournament with one or two Hamiltonian cycles has H=67, at any order. This removes the two-apex qualification from the earlier low-cycle exclusion.
2. The least possible number of Hamiltonian cycles in a **strong** H=67 tournament is exactly **three**. The earlier word `12301` with four cycles was a positive example, not a minimum-cycle witness.
3. Every H=67 tournament has feedback vertex number exactly two. Thus the three-cycle escape does not require leaving the two-apex model.
4. The number of unlabelled order-n tournaments with H=67 is

       0 for n<7;       6(n-6) for n>=7.

   All these tournaments are rigid. The labelled count is therefore `6(n-6)n!` for n>=7.

The strongness qualification in item 1 is essential. Arbitrarily padding a seven-vertex core by universal sources or sinks preserves H=67 but makes the whole tournament weak and gives C=0.

## 2. Sources and advance

The named [Session 196 low-cycle packet](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/lowcycle/PROOF.md) proves the two-apex C=1-or-2 exclusion by exact finite-state methods and supplies the H=67,C=4 word `12301`. Its [guide](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/00-navigation/S196.md) carefully leaves the general low-cycle assertion open.

The imported quantitative theorem is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3, Theorem 1, [DOI](https://doi.org/10.37236/1141), [official paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). The theorem was independently checked in the official paper. We use only its consequence

       H(T)^3 >= 5^(n-1)                              (1)

for every strong tournament T of order n>=3. This is an external theorem, not a new claim or a conclusion inferred from our census.

The new-in-session result is the exhaustive six-core classification, the sharp C=3 threshold, and the resulting all-order unlabelled stratum. The inaccessible general gluing metadata is not used.

## 3. A rigorous order cutoff

If H(T)=67 and T is strong, (1) gives

       5^(n-1) <= 67^3 = 300763 < 390625 = 5^8.

Hence n<=8. There are no strong tournaments of order two; orders one through eight are nevertheless all included in the checker. No arbitrary search cutoff or stabilization assumption is involved.

The complete unlabelled tournament inventories generated at orders 1 through 8 have sizes

       1, 1, 2, 4, 12, 56, 456, 6880.

Exactly six strong classes in this complete range have H=67, all at order seven. There are twelve H=67 classes at order eight, all weak, arising by adding a single source or sink to these six cores.

## 4. The six strong cores

Outgoing row i has bit j equal to one exactly when i points to j. The following canonical representatives are the complete list; the numbering agrees with `certificate.json`.

1. Rows `[2,24,19,21,97,15,47]`, C=3, a two-apex word `22102`.
2. Rows `[2,48,11,19,37,77,31]`, C=4, a two-apex word `23230`.
3. Rows `[2,48,11,35,13,85,31]`, C=3, a two-apex word `13230`.
4. Rows `[4,9,34,21,71,27,47]`, C=4, a two-apex word `12301`.
5. Rows `[4,33,10,19,39,77,31]`, C=4, a two-apex word `13021`.
6. Rows `[16,17,11,35,44,71,31]`, C=3, a two-apex word `23122`.

Here the two-apex convention is exactly the S196 convention: a->b, increasing base x_1,...,x_5, bit 0 for a->x_i and bit 1 for b->x_i. The rows use a different canonical vertex ordering; the certificate records the actual apex pair and base permutation for every listed word. Therefore the word and rows describe isomorphic graphs, not necessarily the same labelled adjacency matrix.

The converse pairs are 1 with 6 and 4 with 5. Cores 2 and 3 are self-converse. Every core has trivial automorphism group. Each has a transitive deletion of two vertices and no transitive deletion of one vertex, giving feedback vertex number exactly two.

All of these claims are checked on the actual adjacency matrices. In particular the three-cycle witnesses are not inferred merely by fitting a transfer formula.

## 5. Why this classifies arbitrary orders

The strong components of a tournament are totally ordered, with every arc between earlier and later components pointing forward. A Hamiltonian path must traverse the components in that order and visits each component in one consecutive run. Therefore

       H(T) = product over strong components D of H(D).       (2)

Every factor is a positive integer. Since 67 is prime, exactly one component has H=67 and all others have H=1. A tournament has H=1 if and only if it is transitive. For the direction needed here, a nontrivial strong component contains a directed cycle and hence a directed triangle; that triangle has three Hamiltonian paths, which inject into Hamiltonian paths of the full component by successive vertex insertion. Thus a nontrivial strong component cannot have H=1. Every remaining strong component is therefore a singleton.

The unique nontrivial component is one of the six cores already classified. Accordingly every H=67 tournament has the unique form

       Tr_p -> Q_i -> Tr_q,       p,q>=0, 1<=i<=6,             (3)

where Q_i is a listed core. Conversely every graph in (3) has H=67 by (2).

For fixed order n, p+q=n-7 has n-6 nonnegative solutions, and the uniquely ordered SCC decomposition distinguishes their isomorphism classes. This proves `6(n-6)`. An automorphism preserves each ordered strong component; the six cores and all singleton components are rigid, so every output is rigid. Deleting the two certified core vertices makes the entire graph transitive, and no one-vertex deletion can make its induced core transitive, proving the feedback-vertex statement.

## 6. Completeness and independent counting in the certificate

`verify.py` is standalone standard-library Python. It starts with the singleton and generates every order-(n+1) tournament by adjoining a vertex in all 2^n incidence patterns to every retained order-n representative. This is complete inductively: delete any vertex of a target tournament, identify its retained representative, and relabel the deleted vertex's incidence pattern accordingly.

The canonicalization uses isomorphism-invariant color refinement, starting with outdegrees and iteratively recording outgoing-neighbor color multisets. It then explicitly tries every vertex permutation within the resulting color classes and chooses the least full adjacency-row tuple. The refinement only restricts to invariant color classes; all remaining permutations are tried. Each canonical tuple is an actual relabelling of its input, so it cannot merge nonisomorphic tournaments. The invariant construction gives the same tuple to isomorphic inputs. Thus the process is an exact isomorphism enumeration, not a hash heuristic.

Hamiltonian paths are counted by independent vertex-subset endpoint dynamic programming. Hamiltonian cycles are counted by depth-first enumeration with vertex zero fixed first, so each cycle is counted once. Strongness is independently checked by transitive closure. All labelled tournaments through order five, 1,099 graphs in total, cross-check both counters against literal permutation enumeration.

For the six target cores the checker explicitly enumerates all 7! permutations to verify rigidity. It enumerates apex pairs and tests the remaining base order, and separately checks every single-vertex deletion. The complete census and all literal graph witnesses are written to `certificate.json`.

Reproduce:

    python3 verify.py > normal.stdout
    cp certificate.json certificate.normal.json
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp certificate.normal.json certificate.json

All checks use explicit exceptions and remain active under `-O`. The finite census plus Busch's theorem proves the all-order conclusion. The broader assertion that every C=1 tournament has a strong bridge, or that the regular order-five tournament is the only bridge-free C=2 tournament, is not proved by this packet. Finite evidence for those questions and the separate exact rooted gluing calculus are recorded in `ROOTED.md`.
