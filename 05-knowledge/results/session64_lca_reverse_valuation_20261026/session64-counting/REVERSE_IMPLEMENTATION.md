# Ordinary A000568: the reverse valuation pass on its actual state graph

2026-10-26. **IMPLEMENTED / FINITE-EXACT:** `reverse_ordinary.py` integrates
the inherited proved reverse potential (`05-knowledge/results/port_residual_integral_20261025/counting/REVERSE_POTENTIAL.md`)
into the unrestricted ordinary tournament recurrence. It builds the actual
structurally reachable `(stage, S, D)` DAG, computes the suffix potential,
and evaluates normalized binary digits directly with K-bit stored residues.
**PROVED:** the optional valuation-budget filter below discards exactly
those graph edges that cannot lie on a complete path of valuation below
`B+K`. This is an arithmetic working-subgraph reduction after a full DAG
pass, not a new quotient of the ordinary structural recurrence.

## Inheritance, scope, and map

Closest mechanism: the integral one-channel defect transfer and the finite-DAG
reverse potential, including its independent `REVERSE_REVIEW.md`. The latter
left integration into the actual ordinary state graph open. This closes that
implementation boundary, not a new proof of the already established gauge.
The concept board is structural continuation / integer valuation / exact
mergers / normalized precision / cost accounting. The canonical hostile is
order 50: the absolute residue modulo 4 is zero, while
`a(50)/2^25 == 2 (mod 4)`. The least-used sidecar is the *minimum prefix cost*
after reweighting; it supplies a safe per-edge arithmetic budget test.

Source: the structurally reachable ordinary odd-cycle prefix DAG, retaining
stage, remaining degree and every affordable defect coordinate. Target: the
same DAG with reweighted integer coefficients; optionally its precision-active
edge list. Preserved: every ordinary continuation contribution and every state
merger, exact integer count, and the requested normalized residue. Lost already
upstream: literal tournament arcs, clone owners, ports and individual
isomorphism identities. Sidecars: target N, suffix potential V and requested
precision K; the optional filter also uses minimum prefix potential P.

**Ordinary means excess e=0.** All unrestricted tournaments, including classes
with clones, remain counted. No signed core-excess transition is silently
replaced by an ordinary one. No positive clone object model or black-core
connectivity result is inferred. Literal owner/port and signed-excess carriers
remain separate. The independent full class census below imposes no prime,
pair-free, connectedness, automorphism, or clone restriction.

## 1. Exact graph and integer edge transfer

For each target N, process all odd lengths ell in decreasing order. A vertex
is `(ell,S,D)`. Include zero-cycle skips and early `S=0` vertices at their
actual stages; all paths end at one `(-1,0,())` vertex. Every structurally
reachable vertex has an all-singleton completion. No coefficient, residue,
or modular nonzero test is consulted while constructing this graph.

At a vertex set `u=N-S`, `b=b(D)`. Choose k actual ell-cycles, with
`S'=S-ell*k`, and force `S'=0` when ell=1. The inherited equations give

    D'(d)=D(d)+k*(d-gcd(ell,d)),   b'=b(D'),
    T=O(u+ell*k)/(O(u)*ell^k*O(k)),
    J=k*(ell-1)/2+ell*binom(k,2)-v2(k!)
      +k*(ell*b-D(ell))+S'*(k+b-b').

Unavailable `D(ell)` is never read when ell>S; then k=0. The exact factor T
is a positive odd integer, and J>=0 by the inherited integral-gauge theorem.
The implementation verifies both statements on every tested edge. It stores
edge parameters, not large exact T values, so the modular production pass
computes T using only modular odd factorials and odd-unit inverses.

## 2. Reverse potential and reachable integrality

Set `V(terminal)=0` and, backwards in topological order,

    V(v)=min_(v->w) (J(v,w)+V(w)),   B=V(initial),
    J'(v,w)=J(v,w)+V(w)-V(v) >= 0.

The existing theorem identifies `B=floor(N/2)` for N>=3 and B=0 otherwise.
The code computes B from the actual graph and checks, rather than inserts,
this identity. Initialize `Z(initial)=1`, then use

    Z(w) += T(v,w) * 2^J'(v,w) * Z(v).

For every reachable partial sum,

    Z(v)=2^(V(v)-B)*Y(v).

Although its displayed exponent may be negative, a prefix path of cost c can
be completed by a suffix of minimum cost V(v). Hence `c+V(v)>=B`, proving
integrality path by path before state mergers. The state-defined potential
makes additions commute. At the terminal,

    Z(terminal)=O(N)*a(N)/2^B.

Thus the modular pass directly computes `a(N)/2^B mod 2^K` using only K-bit
stored residue precision and a final odd-unit inverse. It never divides an already
truncated absolute residue by a power of two. The exact large-integer passes
in the checker are independent validation work, not prerequisites of the
modular API `normalized_modular(build_dag(N),K)`. Multiplication and shift
temporaries in Python may exceed K bits before reduction; the precision claim
concerns the stored residues, not machine-operation widths. No B+K-bit absolute
state or full exact count is stored by this modular path.

## 3. A precision-active arithmetic subgraph

Compute forward minimum reduced costs

    P(initial)=0,
    P(w)=min_(v->w) (P(v)+J'(v,w)).

The minimum reduced suffix cost from every vertex is zero: any original
minimum-cost suffix telescopes to `V(v)-V(v)=0`. Therefore the minimum reduced
complete-path cost among paths traversing edge v->w is exactly

    P(v)+J'(v,w).

An edge with this value >=K lies on no complete path contributing modulo
`2^K`. Conversely, an edge with value <K lies on at least one such path.
So the optional filter keeps exactly the union of complete paths with total
original valuation <B+K. In particular a vertex can be omitted from arithmetic
if P(v)>=K. This uses structural min-cost bounds, never the accidentally zero
current residue. It is exact before sums and carries, not a probabilistic
zero test or a claim that kept contributions cannot cancel modulo `2^K`.

The implementation retains the complete DAG for auditable construction and
potential computation. The active count below measures modular edge updates;
it does not pretend the rejected edges were never generated or visited by
preprocessing. The code also compares every intermediate modular value with
an unfiltered pass and with the exact transformed integer value.

## 4. Verified controls and costs

Run from this directory or pass `--repo /path/to/repository`:

    python reverse_ordinary.py > reverse-normal.json
    python -O reverse_ordinary.py > reverse-optimized.json
    cmp reverse-normal.json reverse-optimized.json

Both receipts are identical. Explicit exception gates remain active under -O.

* Exact ordinary counts agree with the inherited integer engine for every
  N=0..50 and additionally N=60,80,100.
* All 32,875 reachable prefix coordinates satisfy the exact Y/Z gauge.
* All 432 normalized modular comparisons (54 targets, K=1..8) agree at every
  intermediate vertex with exact Z and with the unfiltered modular pass.
* Independent unmerged odd partitions agree for N=0..27, including each
  graph-computed minimum summand valuation B.
* The known A000568 prefix through N=12 is checked.
* Full labelled tournament enumeration and exact S_n orbit marking gives
  class counts `1,1,1,2,4,12,56` for N=0..6. This exhausts the entire labelled
  universes, including all 32,768 tournaments at N=6.
* Distinct-stage equal `(S,D)` keys, absolute-zero/normalized-nonzero states,
  negative naive initial-edge shifts, and invalid input types are exercised.
* Order 50 gives the required normalized residue 2 modulo 4 with K=2 directly.

| N | All input states | All edges | K=2 active states | K=2 updates | K=8 active states | K=8 updates |
|---:|---:|---:|---:|---:|---:|---:|
|30|212|472|41|44|115|143|
|50|926|2732|149|158|364|440|
|60|1720|5499|149|156|583|704|
|80|5119|18109|359|374|1178|1408|
|100|13444|50044|549|568|2180|2576|

Terminal vertices are excluded from the input-state columns. These all-state
counts agree exactly with the inherited one-channel graph, not with the larger
pre-defect f-state graph. The reduced forward arithmetic is a genuine exact
operation-count saving for the stated normalized precision. Construction and
two potential passes still process all states/edges; the reference evaluator
also allocates a full vertex array. No elapsed-time, peak-memory, polynomial
state bound, asymptotic speedup, or end-to-end bit-complexity claim follows.

A cost hostile is important: exact reweighted integers can be *larger* at some
intermediate vertices. At N=100 the sum of snapshot value bit lengths is
5,193,606 for Y but 5,534,143 for Z, even though the final maximum drops from
4,854 to 4,804 bits. Those sums are not peak memory. The gain is the ability
to truncate Z at the requested normalized K-bit precision, not a guarantee
that exact Z arithmetic is always cheaper.

## 5. Boundaries retained

The reverse pass is on a finite decreasing-stage DAG, not an infinite-orbit
termination argument. It is a reachable-lattice gauge, not an integral map
on arbitrary independent vertex assignments. Subtracting B on the first edge
alone fails on zero-cost skips. Dropping structural states because their
absolute residues vanish loses normalized information. Forgetting stage can
identify vertices with different available continuations. Applying the
ordinary formula to signed excess without constructing and auditing its own
weighted graph is outside this implementation.

There is no demonstrated further structural quotient for unrestricted
ordinary A000568 here: the existing degree-gauge quotient is inherited,
not rediscovered. The new delivered result is its reverse-potential
implementation and a provably exact precision-active arithmetic filter.
