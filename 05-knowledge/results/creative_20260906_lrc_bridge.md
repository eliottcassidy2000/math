# Divisor components compile the exact decoder-walk cancellation spectrum

**Status: PROVED ELEMENTARILY + INDEPENDENTLY AUDITED; FINITE-EXACT controls.** This is an exact
compiler for the inherited conditional endpoint-walk criterion. It replaces
an unbounded walk search by at most 512 gcd states per endpoint in an eleven
vertex core, and produces witnesses of at most 19 edges. It does not establish
a new unconditional safe class or prove that a qualifying component always
exists. Actual entry and LRC(14) remain **OPEN**.

The [independent root audit](creative_20260906_root_audit.md) accepts the
proof, inherited consumer scope, and matching replay without repair.

## Inheritance, portfolio, and the retained object

The closest mechanism is the independently audited endpoint-walk identity in
[overnight14_20260906_lrc_endpoint_walk.md](overnight14_20260906_lrc_endpoint_walk.md):
for a positive-integer-labelled walk, its exact cancellation is

    J = gcd(endpoint labels) / gcd(all visited labels).

The parent endpoint-gcd supplier is the audited
[overnight12 signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md),
which uses **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
and **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`.
The connected-gcd descent supplier is
[overnight12_20260906_lrc_decoder_descent.md](overnight12_20260906_lrc_decoder_descent.md).
The current hereditary caps come from the completed, audited
[lrc14_joint_shadow_empty_core_next_sep06.md](lrc14_joint_shadow_empty_core_next_sep06.md),
not its reserved theorem namespace.

The canonical hostile is the actual three-vertex atlas graph `6--2--3`:
its endpoint gcd is three while the collective gcd is one. The corrected
near miss is to erase that factor, or to bound its product-content surrogate
instead of the exact cancellation. The least-used sidecar is the entire
connected component after imposing a divisibility condition. Targeted searches
for gcd/connected-component, divisor-component, prime-barrier, endpoint-gcd,
and walk cancellation on current canon and the overnight reports found no
matching spectrum/compiler statement; this is no external priority claim.
The recent endpoint-walk warning in `01-canon/MISTAKES.md` was read.

Anchor: the actual eleven-core entry obstruction. Niche: the divisor lattice
of a labelled graph. Wildcard: eliminate all walk chronology while retaining
its exact attainable cancellation and maximum distinct support. The concept
board has five objects: endpoint gcd, divisor-induced component, collective
gcd closure, distinct support size, and physical coefficient/phase consumer.
The first four now have an exact correspondence; the fifth retains all its
inherited actual-entry hypotheses.

The sparse-interval transport report
[synthesis_20260905_lrc_sparse_transport.md](synthesis_20260905_lrc_sparse_transport.md)
was also read. Its component closure suggests testing which optimization is
actually present, but there is no claimed map from an interval contact graph
to a decoder graph. No star-forest property or physical mass conclusion is
transferred between these different graph carriers.

## 1. Exact attainable cancellation factors

Let `G` be a finite undirected graph on `n>=2` vertices with positive integer
labels `w(v)`. Distinct vertices may have equal labels. Fix distinct endpoint
vertices `x,y`, and put `D=gcd(w(x),w(y))`. A walk may repeat vertices and
edges. Its support is its set of distinct vertices.

For each positive divisor `h|D`, let `C_h` be the connected component of `x`
in the graph induced by vertices whose labels are divisible by `h`. Both
endpoints satisfy the divisibility condition, but `y` need not belong to
that component. If it does, put `d_h=gcd(w(v):v in C_h)`.

**Theorem 1 (exact spectrum and support maximum).** An integer `J>0` occurs
as the exact cancellation of some `x`-to-`y` walk if and only if

    J|D,    y in C_(D/J),    d_(D/J)=D/J.                 (1)

For every such `J`, the greatest attainable number of distinct visited
vertices is exactly `|C_(D/J)|`. Hence a walk with cancellation `J` and at
least `r` distinct vertices exists exactly when (1) holds and that component
has at least `r` vertices.

**Proof.** If a walk has support `S` and cancellation `J`, the inherited
identity gives `gcd(w(S))=h=D/J`. Its support is connected, contains both
endpoints, and all its labels are divisible by `h`. Therefore `S subset C_h`.
All labels in `C_h` are divisible by `h`, so `h|d_h`. Conversely the gcd of
`C_h` divides the gcd of its subset `S`, which equals `h`. Thus `d_h=h`.
This also bounds the support size by `|C_h|`.

Conversely, when (1) holds, a walk from `x` to `y` can visit every vertex of
the finite connected graph `C_h`; for example traverse a spanning tree,
returning along each branch and finishing along the tree path to `y`.
Its collective gcd is `d_h=h`, so its cancellation is `D/h=J`. This attains
the upper bound on support. QED.

For a budget `B>=1`, the equality `d_h=h` may be omitted if the question is
whether some `J<=B` occurs with at least `r` distinct vertices: it is
equivalent to the existence of `j|D`, `j<=B`, for which `y in C_(D/j)` and
`|C_(D/j)|>=r`. Indeed the full-component walk has actual factor
`D/d_(D/j)`, which divides `j` and is therefore at most `j`. This distinction
between an exact factor and a budget is necessary.

This is not a theorem about exactly `r` distinct vertices. A long induced
path can require more than `r` vertices just to connect the endpoints.

## 2. Finite certificates independent of integer height

**Theorem 2 (bounded witnesses and factorization-free state space).** Every
attainable factor `J` has a witness with at most `2n-3` edges. All factors and
their maximum support sizes can be computed by inspecting at most `2^(n-2)`
gcd states, without factoring any integer or enumerating walks.

**Proof.** For a qualifying component of size `s`, choose a spanning tree and
let its endpoint path have `ell>=1` edges. Traverse each edge off that path
in both directions and each edge on that path once, ending at `y`. The walk
uses `2(s-1)-ell<=2s-3<=2n-3` edges and visits the full component.

Every attainable `h=D/J` is the gcd of a connected vertex subset containing
`x,y`. It therefore belongs to the possibly redundant list

    { gcd(D, w(v):v in T) : T subset V(G) minus {x,y} }.  (2)

The list has at most `2^(n-2)` entries. Start with the state set `{D}` and,
for each non-endpoint vertex `v`, add `gcd(h,w(v))` for each current state
`h`. For each distinct final state inspect its induced component and apply
(1). Theorem 1 proves completeness and soundness. QED.

The factor set is invariant under multiplying all labels by the same
positive integer. Components are relabeling-covariant. Walk chronology,
edge multiplicity, and the larger accumulated numerator/denominator content
are discarded. The graph and the induced-component gcd are the sidecars
that preserve exact cancellation, endpoints, and maximum distinct support.

## 3. The precise LRC consumer

Retain the full actual `11+2` entry hypotheses of the inherited signed-box
result: primitive thirteen-speed physical row, `sum speeds<=Q^2`, actual
THM-3818 decoder graph, and `W=V_dec`. Write the eleven physical core as
`tV`, where `gcd(V)=1` and `K=max V`, and start at `x=tK`.
The signed-box endpoint criterion uses

    Q=91^6,  H=floor(Q/(42*177))=76,388,115.

For a hypothetical strict failure the inherited physical subset-gcd caps
are `M5=11,342,250`, `M6=31,950`, `M7=90`, and `t<=2`. For a fixed other
endpoint `y`, Theorem 1 gives a complete finite test for the existence of
an inherited sufficient walk with at least `r` distinct vertices and

    J <= floor(tH/M_r),    r in {5,6,7}.                (3)

The budget formulation after Theorem 1 may be used, or the exact spectrum
may be filtered. A walk visiting more than `r` vertices still has collective
gcd at most `M_r`, because its gcd divides the gcd of any `r` visited
vertices. Therefore (3) implies

    gcd(K,y/t) = D/t = J gcd(visited physical labels)/t <= H,

and the inherited endpoint theorem supplies a phase with clearance at least
`1/14`. The closure proof uses a strict-failure contradiction; scalar gcd
caps are not asserted for every safe actual entry.

There are ten endpoints, at most 512 gcd states per endpoint, and witnesses
with at most 19 edges. This is an exact finite compiler at a given core,
not a finite enumeration of all core labels. It establishes no new
unconditional safe class beyond the inherited walk criterion. Proving
that every remaining core contains a qualifying component is still OPEN.

The new restriction `t<=2` increases the eligibility range for the signed
coefficient gate, but does not by itself replace `H` by `Q/2`: the phase
comparison still needs `g>QK/(D p)>=42K`. The old `H` was limited by that
comparison, not merely by the coefficient gate `tD/delta<=Q`.

## 4. Hostiles, controls, and reproduction

On the actual atlas graph `6--2--3`, the endpoint pair `6,3` has `D=3`.
At `h=3`, the divisible vertices `6,3` are disconnected. At `h=1`, the
component has gcd one and three vertices. Thus the spectrum is exactly
`{J=3: maximum support 3}`. Erasing connectivity falsely admits `J=1`;
erasing cancellation falsely turns collective gcd one into endpoint gcd one.
The first failed implication is respectively induced divisibility to a
connecting walk, or collective gcd to endpoint gcd. The repaired form is
Theorem 1. The missing coordinate is the component through which the walk
must pass, and its exact collective gcd.

A five-vertex star with central label one and leaf labels `2,3,5,7`, with
endpoints seven and five, has spectrum `{J=1: maximum support 5}`. Every
simple endpoint path visits only three vertices. A six-edge repeated walk
visits all five. This generic labelled-graph control is not asserted to be
an induced actual decoder graph. It shows why replacing arbitrary walks by
simple paths loses valid support certificates.

The [standalone source](../../04-computation/creative_20260906_lrc_bridge.py)
compares three independent descriptions: all connected vertex subsets,
divisor-induced components, and the factorization-free gcd state machine.
It separately reconstructs each spanning-tree walk and calculates its
cancellation from rational edge ratios and prefix denominator lcms.
The finite universe contains every four-vertex graph, all ordered positive
label vectors in `{1,2,3,4}^4`, and all six distinct endpoint pairs, including
disconnected graphs and repeated labels; it adds the path, cycle, star,
and complete five-vertex graphs with every distinct label subset of
`{1,...,8}` and every endpoint pair. No actual-entry filter is silently
inherited by these generic graph tests.

    python3 -B 04-computation/creative_20260906_lrc_bridge.py
    python3 -B -O 04-computation/creative_20260906_lrc_bridge.py

The matching [frozen output](creative_20260906_lrc_bridge.out) records the
case count and always-active gate count. The theorem is proved above;
the finite bank tests its normalization, both hostiles, and common scaling
at a factor of `10^80`. No theorem ID, external priority, universal phase
closure, or mathematical conclusion from a failed sufficient test is claimed.

Normal and optimized runs are byte-identical: **100,544** labelled
Graph/endpoint cases and **559,903** always-active gates. Source SHA256:
`1dc2a4441fd3b343cb174c13dee9312f7aada94cb14ad2de3b19e762f6d31a3a`.
Output SHA256:
`c935c2e12178d578a749a8b2100e2d53739c93a80af485e02f6a73219288fe23`.

**Incoming freshness check.** The concurrently fetched
`overnight15_20260906_lrc_larger_unit.md` extends actual unit-containing
component closure to further split types. Its introduction and scope were
read through `git show origin/main:...`; it does not replace the nonunit
endpoint obligation compiled here. No new dependency on that result is used.
