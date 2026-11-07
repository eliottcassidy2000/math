# Square-sum onsets, the 3/10 counts, and the missing port coordinate

**Status: source-pinned synthesis; elementary all-graph bijection; FINITE-EXACT replay through order 25.** No new global Hamiltonicity, tournament-spectrum completeness, or literature-priority claim. The counts below are inherited and independently recomputed, not newly discovered numerical values. The new contribution is their explicit use as transfer obstructions and their comparison through a common boundary state.

Reference main: `214586af6c8803d8b21959951e8f151e64bff36f`. `SOURCES.json` pins the repository source paths and blob IDs. The shared checkout was left unchanged.

## 1. Three different axes, not three interchangeable pairs

Let `Q_n` be the simple undirected graph on `1,...,n` with edge `xy` exactly when `x != y` and `x+y` is a square. Let `p(Q_n)` count its Hamiltonian paths up to reversal, equivalently as unoriented spanning path edge sets. For a tournament `T`, let `H(T)` count directed Hamiltonian paths and put `q(T)=(H(T)-1)/2`.

The cited source trail resolves the numbers as follows.

- **`7,21`: tournament count values.** They are excluded values of `H`, at every order. The assertion that they are the only excluded positive odd values remains open. Use THM-1370 *h-spectrum-omits-7-21-all-n*, not the other theorem sharing that ID, and the corrected Session210 atlas.
- **`3,10`: tournament index values.** In June's HYP-2198 *single-core-signature-gap-complete* they are missing weighted inversion sums for a binary apex word. This is a restricted arithmetic construction, also missing `6,14,17,...`. In Session210 they are the globally excluded cardinalities of the signed odd-cycle-packing quotient. This is a literal combinatorial interpretation of `q`, not a separately classified number-theory spectrum.
- **An arithmetic observation, not a transport.** The October2 *camion_busch_gaps_polyhedra_collatz* note §5 explicitly labels `3=T_2`, `10=T_4`, and `2q+1=q_field^2+q_field+1` at field sizes 2 and 4 as numerology in its attempted Collatz correspondence. The same section gives the genuine one-apex weighted-inversion formula and its extra gaps. No independent arithmetic family whose complete spectrum is precisely `N_0 minus {3,10}` was identified in this targeted inheritance pass.
- **`15,25`: square-sum vertex orders.** The first nontrivial chain is at 15. Gerbicz's all-orders chain existence for every `n>=25` is CITED by the September22 source, not reproved here. The exact finite onsets are **15,23,25**, not merely the pair 15,25: a chain exists at `15,16,17,23`, and at every order from 25 in the cited theorem. The trivial graph `Q_1` is a separate convention.
- **The old `15,25` coincidence.** The October22 *monster_levels_square_sum_7_21* note §5 already rejects the identification of these square-sum orders with the strong-tournament minimum values `f(6)=15,f(7)=25`. Those are different axes and no map is supplied by the shared numbers. Its endpoint correction and MISTAKE-668 must be retained.

Thus one should not describe the present result as discovering a common exceptional-number law. It identifies exactly which structure can and cannot cross between the problems.

## 2. A sharp obstruction to an affine count transport

The previous square-sum census, independently replayed by three formulations here, gives

    p(Q15)=p(Q16)=p(Q17)=1,
    p(Q23)=3,  p(Q24)=0,  p(Q25)=10.

Consequently no assignment defined on all square-sum graphs can produce a tournament `F(Q_n)` satisfying

    q(F(Q_n)) = p(Q_n),
    equivalently H(F(Q_n)) = 2p(Q_n)+1.

**Proof.** At `n=23` it would give `H=7`, and at `n=25` it would give `H=21`, both impossible. This is a nonexistence statement for this exact scalar-preserving transport, independent of what tournament orders or auxiliary gadgets it uses. It does not rule out maps preserving existence, weighted data, marked-path subcounts, or another statistic. It does not disprove tournament-spectrum completeness.

There is an informative symmetry distinction. Reversal acts freely on the directed traversals of an undirected path whenever the graph has at least two vertices, so the symmetric digraph associated with `Q_n` has `2p(Q_n)` directed Hamiltonian paths. Tournament reversal is a map to the converse tournament, not a second path in the same tournament. Adding a vacuum `+1` to the former count is not justified by either symmetry.

For the general recognizable-gate construction `a -> b -> D -> a`, cycles are in bijection with directed spanning paths of the core `D`. Hence for the symmetric digraph of `Q23` or `Q25` the gate-cycle counts are 6 or 20, not 7 or 21. The resulting strong digraph is not a tournament: core edges come in both directions and core nonedges remain missing. The session's common path-cover packet supplies the general gate/cover proof. Its exact map does not turn arbitrary core relations into tournament relations.

## 3. The correct microscopic state at a degree-two birth

Here is the boundary-state operation shared by square sums and the path-cover viewpoint.

**Proposition (two-port insertion).** Let `G` be a simple graph and let `G+x` add one vertex whose distinct neighbors are exactly `u,v`. Suppose `G` has no Hamiltonian path. Hamiltonian paths of `G+x`, up to reversal, are in bijection with spanning linear forests of `G` having exactly two components, each a path (singletons allowed), such that `u` and `v` lie in different components and are endpoints of their respective components. The map deletes `x`; its inverse adds the edges `ux,xv`.

**Proof.** The new vertex cannot be an endpoint of a spanning path, since deleting an endpoint would give a Hamiltonian path of `G`. Thus both its incident edges are used. Deleting it produces exactly two path components with the stated marked endpoints. Conversely the two connecting edges concatenate any such forest into one spanning path. Both maps are inverse on unoriented edge sets. A singleton has one vertex port, not two artificial orientation choices. No assumption that a graph-degree-two vertex is interior has been used except where the absence of a predecessor path proves it.

This is an elementary special case of the September endpoint-aware reductions and Session32's named-port transfer theorem, not a novelty claim for boundary-state dynamic programming.

For the three square-sum onsets the ports and forest counts are exactly

    n=15: predecessor Q14, ports {1,10}, count 1;
    n=23: predecessor Q22, ports {2,13}, count 3;
    n=25: predecessor Q24, ports {11,24}, count 10.

Each pair is determined by the actual square equations: `15+1=16`, `15+10=25`; `23+2=25`, `23+13=36`; `25+11=36`, `25+24=49`.

The script independently computes the number of these forests without constructing paths in `Q_n`: for each allowed vertex subset `S` containing `u` but not `v`, multiply the number of paths starting at `u` spanning `S` by the number starting at `v` spanning its complement, and sum. Every component is oriented uniquely from its distinguished port, including singleton components. This gives 1,3,10.

This explains the productive micro/macro relationship: new chains are assembled from compatible **two-path** states of a graph that itself has **zero one-path** states. Information invisible to the one-path scalar controls the next scale.

### Scalar recurrence obstruction

All the additions at

    n=13,14,15,19,20,21,22,23,25

have the same coarse input `(p(Q_(n-1)), degree_new)=(0,2)`, but their outputs are

    0,0,1,0,0,0,0,3,10.

Therefore no universal recurrence depending only on the old path count and the new vertex's degree can recover these births. The missing coordinate is the designated-port two-component cover count. This statement does not exclude formulas also using `n`, the port labels, or a fuller graph state. In particular it is not a proof that all finite-state or arithmetic recurrences fail.

## 4. Endpoint inventories and the genuine arithmetic scaling

All three `Q23` chains have endpoint 18. Their other endpoints are `2,9,22`, each once. All ten `Q25` chains have endpoint 18. Their other endpoints are

    2,3,4,8,9,10,11,13,22,23,

each once. `RESULT.json` includes all 14 onset paths and the exact two-component forests after deleting their newly inserted vertex. These named endpoints are the required sidecar; their cardinality alone does not explain the next birth.

The arithmetic mechanism at the first onset is already much stronger than a numeric coincidence. THM-4510 *sum-graphs-as-reflection-orbits-pythagorean-zigzags* interprets each target sum `s` as the reflection `x -> s-x`; two reflections compose to translation by their target difference. The tight triple `9,16,25` yields the first window `15,16,17`, and the primitive Pythagorean-triple theorem produces chains at orders `t^2-1,t^2` from `s^2+t^2=u^2`, with the specified primitive and ordering hypotheses. The next cited examples are `(5,12,13)` at order 143 and `(8,15,17)` at order 224. This is an actual cross-scale arithmetic construction, unlike the coincidences with 7 and 21. A fresh finite replay in `verify_pythagorean_scaling.py` verifies the restricted three-square paths at both orders for all 11 primitive triples with `1<=s<t<=60`, through order 3600. It does not enumerate all paths of the larger square-sum graphs; normal and optimized results agree.

A second inherited macro mechanism is Gerbicz's 25-fold block blow-up: within a block `25a_k+(-1)^(k+1)c`, adjacent sums equal `25(a_k+a_(k+1))`, hence remain squares. The full gluing requires its odd-order, prescribed endpoint, and junction hypotheses. Those cannot be replaced by the scalar `p(Q_n)`. The cited global threshold theorem requires further constructions; this packet does not claim to verify them.

## 5. Connection and loss ledger

- **Square paths -> two-port covers:** exact deletion/insertion bijection; preserves the entire path object; requires the two named ports and component compatibility. Coarse degree or mere number of components loses the answer.
- **Square paths -> gate cycles:** exact after replacing each undirected edge by both directions; preserves directed paths and doubles unoriented counts. Loses tournament status, so neither tournament oddness nor its forbidden `H` values transports.
- **Square paths -> tournament odd-cycle index:** proposed scalar transport is refuted by `Q23,Q25`.
- **Pythagorean triples -> square chains:** inherited exact reflection/rotation construction with hypotheses. Preserves square-edge validity; does not preserve tournament Hamilton counts.
- **Triangular values `3,10` -> forbidden tournament indices:** a numerical identity, not an independent classification. The actual counting carrier is the Session210 odd-cycle quotient.

The promising shared object is therefore the path-cover algebra with named boundary data, not equality of the three small number pairs.

## 6. Reproduction

    python3 -B verify_square_bridges.py
    python3 -O -B verify_square_bridges.py
    python3 -B verify_pythagorean_scaling.py
    python3 -O -B verify_pythagorean_scaling.py

The complete universe is every `Q_n`, `1<=n<=25`; no random sampling or order filter. The three implementations are simple-path DFS with justified leaf-root symmetry reduction, memoized subset endpoint recurrence, and designated-port two-path forest convolution. The first two agree at every order; the third agrees at all nine zero-count/degree-two birth controls. The normal and optimized result bytes agree. An independent auditor separately enumerated all paths without reading this script.
