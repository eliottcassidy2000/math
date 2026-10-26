# Production ordinary counting with a reverse valuation potential

2026-10-25. **IMPLEMENTED, with PROVED recurrence equivalence:** build the
ordinary A000568 structural DAG, compute its reverse minimum-cost potential,
and evaluate the normalized count a(N)/2^B modulo2^K using only K-bit
residues. Safe reduced-cost pruning is included. **FINITE-EXACT:** the
declared universes in `check.py` and both receipts. This closes the
production-integration gap left by the preceding algebraic construction.
It does not reduce the ordinary structural state graph or establish a
runtime or full bit-complexity advantage.

## 1. Inheritance, scope, and the implemented change

The nearest proved mechanism is
[the reverse valuation construction](../../port_residual_integral_20261025/counting/REVERSE_POTENTIAL.md)
and its [independent review](../../port_residual_integral_20261025/counting/REVERSE_REVIEW.md).
That packet deliberately implemented only a generic DAG control. The
[integer gauge](../../port_residual_integral_20261025/counting/README.md)
supplies this production graph's nonnegative edge shifts and odd factors.
The [bounded-cycle precision method](../../port_residual_integral_20261025/precision/README.md)
supplies an independent normalized-residue oracle.

The canonical hostile is a zero-cost initial skip followed by a costly
continuation: subtracting the whole baseline at the first edge gives a
negative exponent. The corrected near miss is using a residue-pruned
absolute graph to construct the normalized computation. The least-used
sidecar is the minimum remaining valuation of every structural state,
computed before any residue is discarded. The live concepts are structural
reachability / exact defect merges / odd denominator units / remaining
valuation / normalized observations / modular edge pruning.

This implementation is the ordinary, zero-excess transfer only. It does
not modify the signed clone-excess recurrence, its semantic support, or
the exact clone-cut interpretation in the preceding packets. Positive
arithmetic weights in the ordinary graph do not provide an objectwise
positive model for clone cancellation.

## 2. The structural DAG is independent of all coefficient values

Fix N. Process odd lengths ell in descending order. A graph vertex is

    (ell,S,D),      u=N-S,

where S is remaining degree and D(d)=d*m-f(d) is retained for the odd
coordinates d<=min(ell,S). Here m is the actual prefix cycle count and
f(d)=sum gcd(c,d) over its past cycles. The retained state does not store
m. The ordinary gauge makes every continuation depend only on this state.
Coordinates discarded because d>S or d>ell never become affordable later:
both remaining degree and stage decrease.

For each k>=0 with ell*k<=S, choose k cycles of length ell. At ell=1
force k=S. Set

    S'=S-ell*k,       u'=u+ell*k,
    D'(d)=D(d)+k*(d-gcd(ell,d)),
    b=max_(retained d>=3) ceil(D(d)/(d-1)),
    b'=the corresponding bound at the destination.

Empty maxima are0. The edge goes to the next stage, even when k=0 or
S=0. These stages remain distinct vertices. At the last stage all paths
reach the one vertex(0,0,()). For N=0 this is already the initial vertex.

The old integer scalar coordinate is the summed weight

    Y=O(u)*2^[S*(m-b)]*w,

where O(u) is the odd part of u!. Each edge has

    Y_destination += T*2^J*Y_source,
    T=O(u')/[O(u)*ell^k*O(k)] >0 odd,
    J=A_ell(k)-v2(k!)+k*(ell*b-D(ell))
       +S'*(k+b-b') >=0,
    A_ell(k)=k*(ell-1)/2+ell*binom(k,2).              (1)

When ell>S only k=0 is possible and D(ell) is not read. The preceding
integer-gauge proof establishes both exact divisibility defining T and
nonnegativity of J. Graph construction checks the nonnegative costs;
the optional exact-integer replay checks the odd-factorial divisions.
The modular production path uses their proved odd-unit identity directly.

`OrdinaryDAG(N)` creates vertices and edges solely from these formulas.
No exact weight, sign, residue, or nonzero-support test participates in
graph generation. All generated vertices are reachable by construction.
Every remaining positive degree can be completed with singletons at the
last stage, so every vertex is terminal-reaching. Every edge decreases
the stage, establishing acyclicity without an infinite-orbit assumption.

The genuine merger at N=26 illustrates the retained information. Before
stage3, histories(23) and(9,9,5) both reach S=3,D=(0,2), although their
past cycle counts differ. The implementation assigns them the same vertex
and preserves their normalized contributions. Omitting the stage from the
key would make additional, unjustified identifications; the checker counts
and guards against those lost distinctions.

## 3. Reverse pass and pathwise normalization

After the entire structural graph is known, compute

    V(terminal)=0,
    V(v)=min_(v->w) [J(v,w)+V(w)],
    B=V(initial).

The reverse pass performs one minimum relaxation per edge and stores one
integer V per vertex. Define

    Z_v=2^[V(v)-B]*Y_v,
    J'(v,w)=J(v,w)+V(w)-V(v)>=0.                    (2)

For any prefix path p to v, appending a minimizing suffix gives a complete
path of cost cost(p)+V(v)>=B. Its contribution to Z_v is therefore an odd
integer times a nonnegative power of two. This proves integrality before
merging. All histories at a vertex use the same V(v), so addition commutes
with the transformation. The ordinary defect quotient preserves every
available continuation, which makes concatenating the prefix and suffix
legitimate even at merged vertices.

The production recurrence is now

    Z_initial=1,
    Z_w += T(v,w)*2^J'(v,w)*Z_v.                    (3)

It runs from its initial state. It never divides an already truncated
absolute Y value by a power of two. This is a gauge on the reachable
path lattice; arbitrary independent integer Y values need not have
integral transformed coordinates.

A complete graph path chooses exactly one odd cycle partition lambda.
The intermediate gauges telescope, leaving O(N) times its Burnside
summand. The product of T factors is the odd integer

    O(N)/product_d[d^k_d*O(k_d)].

Hence its total original edge cost is exactly

    e(lambda)-sum_d v2(k_d!),

the rational summand's binary valuation. The termwise theorem gives
B=floor(N/2) for N>=3, and B=0 at N=0,1,2. Equality is attained by(N) for
odd N and by(1,N-1) for even N>=4. The producer computes V and checks this
identification; it does not replace the reverse pass with the known B.

B need not be v2(a(N)). At N=50, B=25 but the final count has valuation26.
Binary carries among positive minimum-valuation terms can raise the
valuation of their sum.

## 4. K-bit evaluation, odd units, and safe pruning

At termination Z_terminal=O(N)*a(N)/2^B. Compute all residues modulo2^K,
then multiply by the inverse of odd O(N). No B additional residue digits
are needed in this coordinate system. Ordinary arithmetic temporaries can
be wider than K bits; the receipt measures those widths separately from
stored residue precision.

The production modular path constructs O(j) modulo2^K from
O(j)=O(j-1)*oddpart(j). It obtains T modulo2^K by inverting its odd
denominator factors. No exact factorial, rational Burnside summand, or
large exact count is formed in this path.

If J'>=K, the entire edge coefficient vanishes modulo2^K because T is an
integer. The evaluator rejects that edge before computing its modular T
or multiplying a weight. This is safe for every source residue. It leaves
the structural graph and V intact: the edge is omitted only from this
precision's coefficient work.

The evaluator also skips an edge whose already normalized source residue
is zero. This is safe after all predecessor contributions have merged at
that vertex. Layer-by-layer evaluation ensures that ordering. It is not
the invalid earlier operation of constructing the structural graph from
nonzero absolute residues. Disabling both skips gives a full edge-by-edge
modular reference; disabling just source-zero skipping isolates the effect
of reduced-cost pruning.

Residue storage uses two consecutive graph layers. The full structural
graph, its defect tuples, edges, and V remain stored. A smaller residue
precision therefore does not imply a smaller structural state graph or
smaller total memory in every implementation.

## 5. Production API and verification

`reverse.py` exposes:

* `OrdinaryDAG(N)`: generate the structural graph and its reverse potential.
* `graph.residue(K)`: normalized a(N)/2^B modulo2^K and work counters.
* `graph.residue(B+K, normalized=False)`: the original absolute gauge on
  the same graph; dividing its final residue by2^B gives the same target.
* `graph.exact(audit_original=True)`: optional full-integer replay, checking
  Y_v and Z_v at every intermediate vertex, and returning the full a(N).
* `graph.structure()`: state, edge, merger, mask, and storage counters.

The optional exact replay is a checker and usable exact-count route;
`residue` does not call it. Stage and defect data are explicit in every
stored state, and each edge retains its chosen cycle multiplicity for
partition-path audits.

Run from the repository root:

    python 05-knowledge/results/lca_ports_reverse_20261025/counting/check.py --receipt 05-knowledge/results/lca_ports_reverse_20261025/counting/normal.json
    python -O 05-knowledge/results/lca_ports_reverse_20261025/counting/check.py --receipt 05-knowledge/results/lca_ports_reverse_20261025/counting/optimized.json

Both modes use explicit exceptions and produce identical receipts. The
controls cover:

* Full counts and every intermediate Y-to-Z identity for N=0..50, with
 12,589 checked coordinates, compared against the prior integer response
 engine. Additional full-count controls at N=80,100 also compare that engine.
* 408 normalized residue comparisons, N=0..50 and K=1..8, against both
 the exact count and the separately implemented bounded-cycle method.
* 153 comparisons with fully unpruned evaluation and with only edge-cost
 pruning enabled, verifying both the answer and edge-work accounting.
* 3,725 complete partition paths at N=0..34, checking each path's original
 cost, odd factor, and the unmerged rational count independently.
* The actual N=26 merger, a zero-cost initial skip at N=7, the distinction
 between absolute and normalized low bits, and N=50's extra final factor
 of two. Invalid orders and precisions are rejected.

The special-order counts0,1,2 are included. In particular, the two-singleton
partition at N=2 has cost0, not the invalid baseline N/2=1.

## 6. Measured work, with explicit exclusions

These are deterministic operation observations from the receipts, not
timings. Each comparison uses this same DAG and targets six normalized
bits; the absolute route uses B+6 stored residue bits.

| N | Vertices including terminal | Edges / reverse relaxations | Normalized modular updates | Absolute modular updates | Normalized multiply-operand bit sum | Absolute multiply-operand bit sum |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 927 | 2,732 | 402 | 1,116 | 2,899 | 46,518 |
| 80 | 5,120 | 18,109 | 1,245 | 6,061 | 9,094 | 403,618 |
| 100 | 13,445 | 50,044 | 2,382 | 15,331 | 17,648 | 1,274,624 |

At N=100,K=6, reduced costs discard31,103 edges before coefficient work;
16,559 additional edges have zero normalized source residues, leaving
2,382 modular updates. All50,044 structural edges are still scanned, and
the reverse potential already cost50,044 relaxations. The full graph has
not become smaller.

The operand-bit sum adds the bit lengths of the source residue and odd
factor for each performed weight multiplication. It is not multiplication
bit complexity. The receipt also records maximum product, shifted-product,
sum, and stored-residue widths; cache entries; odd inverse calls; odd-
factorial setup products; and the peak number of scalar weight slots in
two layers. The layer-residue bit sum is accumulated across layers, not
peak memory. These counters exclude Python object overhead and do not
combine modular inversion, graph generation, and instrumentation into an
end-to-end cost model. No fastest-engine or runtime claim is made.

## 7. Exact connection and remaining boundary

Source: the existing ordinary structural defect DAG with odd integral
factors and nonnegative costs. Target: the same graph with reduced costs
and a normalized K-bit scalar at each vertex. The map is the state-local
reachable-lattice gauge(2). Preserved: every path contribution after its
known common baseline is removed, exact state mergers, and the final
normalized ordinary count. Deliberately lost: binary digits above the
chosen precision. Needed sidecars: stage, N, exact defect masks, and V.

The production integration is complete. The next boundary is whether a
less expensive state-dependent lower potential can replace materializing
the full reverse graph while preserving nonnegative reduced costs and
the same baseline. No such replacement is claimed here. The signed clone
transfer remains governed by its own already audited integral and modular
recurrences; this ordinary result does not bypass its signed semantics.
