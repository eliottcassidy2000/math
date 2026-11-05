# Independent review: returning three-port ternary clone clock

**PASS**, no corrections requested, 2026-11-05.

Frozen source path `session197-hybrid-constructions/ternary-clock/`:

- PROOF.md: 6ca5ee094d0f13b6e4a48ba513891594e490cc67910776aeb7107943db2fcebf
- core.py: 3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c
- verify.py: c9dda7477ba23725130d9fc4d9a7385dc0dd408d44550487e51d70402e62ebf5
- Both validation receipts: 1818956248ff697fbf0b092783f56c5902d319b254873d1619ba0448ba020c89

Copied-source normal and optimized runs reproduce this receipt hash exactly. The original files remain unchanged.

## Proof and interface

The seed rows define a valid eight-vertex tournament. Its distinguished vertex beats exactly the three fixed ports 0,1,4 and loses to the other four ordinary vertices. The explicit clone rows have the old clones directed forward toward a newly appended clone. Thus adding a vertex that beats exactly those three retained ports and loses to every other old vertex is literally T_m -> T_(m+1), without relabelling a response-equivalent vertex or changing old arcs. The root module grows, but the exposed physical ports do not move.

The supplied seed strongness paths are valid. Nonempty substitution into a strong quotient preserves strongness, including a transitive clone block: two vertices inside that block can reach each other by going around a quotient cycle when the internal direction is wrong.

The maximal-run skeleton decomposition is bijective. Once clone vertices are assigned to the ordered nonempty runs, their increasing order inside each run is forced. The number of assignments is k!S(m,k), not a binomial coefficient. The O-to-I crossing graph has exactly the stated four edges and matching number two; each internal fragment supplies a vertex-disjoint crossing edge. Hence no fourth run is possible.

The exact coefficients a_1=179,a_2=180,a_3=60 simplify to 60*3^m-1 for every positive m. The zero-occurrence count59 is separately H(T-r) and is not used as a cloning-formula term. The offset identity is exactly 3H+2.

The local two-port exclusion is correctly restricted to this particular one-vertex extension: the newborn has three old out-neighbors and at least five old in-neighbors. Relabelling cannot change these two degrees, but a different graph, multiple births, or another interface is not excluded.

## All-order +4 obstruction

Permuting the k-1 internal nonempty fragments is a free action on r-skeletons. Every fragment remains legal between two r occurrences, and disjoint nonempty label sets make the action free. Thus (k-1)! divides a_k, in particular a_3 is even.

For any eventual pure-clone recurrence H_(m+1)=3H_m+4, positivity of H makes the tail coefficient in C*3^m-2 nonzero. The largest nonzero run coefficient a_K is also the coefficient of the largest exponential K^m, so independence of positive exponential sequences forces K=3 and C=a_3. The resulting count is even, contradicting the odd Hamilton-path count of a tournament. This applies at every seed order and any number of exposed ports, but only to pure transitive-root cloning as stated.

## Independent finite controls

The new `ternary-clock-audit/independent.py` does not import either source module.

- It enumerates literal Hamilton words on every subset of the seven ordinary vertices, collecting initial fragments ending in I, final fragments starting in O, and nonempty internal fragments starting in O and ending in I.
- Disjoint-subset convolution of those fragment arrays independently reconstructs the entire vector [59,179,180,60,0,0,0,0,0]. This is distinct from both the source gap-word enumeration and its visited-vertex/run-count DP.
- A separately built adjacency formula is checked by literal Hamilton-path DFS, without subset DP, for m=1,...,6. It gives179,539,1619,4859,14579,43739.
- Full Floyd reachability confirms strongness, and row-by-row comparison confirms all six physical returning operations.

The source replay additionally exhausts 1,285,200 skeleton candidates and all60 three-run involution examples. The all-order clock and obstruction follow from the displayed decompositions and free action, not from extrapolating these finite controls.
