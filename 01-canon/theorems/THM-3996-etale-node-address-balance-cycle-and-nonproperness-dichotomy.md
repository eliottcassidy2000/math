---
id: THM-3996
title: "Etale node-address balance forces cycles or nonproperness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a dominant etale polynomial plane map and an
  integral target curve with a node outside the nonproperness locus, normalize
  the finite pullback and direct each source-node address from its plus-branch
  component to its minus-branch component. Every vertex has indegree equal to
  outdegree, both equal to its finite-etale covering degree. Hence every edge
  lies on a directed cycle; a complete forest forces the target node into the
  nonproperness/Jelonek locus. Applied to THM-3992, distinct owners of the two
  known companion germs require an additional node address unless the node is
  nonproper. If the two known addresses are the full connected packet, the
  companion germs have one owner and form a two-edge cycle with t=0.
source: root + nodal_companion_completion / THM-3992 completion audit, 2026-08-24
audit: >
  PASS (root + thm3996_hostile_audit, 2026-08-24). The graph proof was
  independently checked through finite-etale base change to normalization and
  a condensation-DAG balance argument. The companion verifies the THM-3992
  normalization, node addresses, cycles of lengths 2..12, the deleted-address
  path hostile, Riemann-Hurwitz degree-one control, and an embedded two-cycle
  model. Normal and optimized runs match.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - MISTAKE-469
script: 04-computation/jc2_etale_node_address_balance_thm3996.py
output: 05-knowledge/results/jc2_etale_node_address_balance_thm3996.out
script_sha256: f3a7560606de94cfe0ee705459f491b27a6b25f7c2d51ac4e35e6cb9734f16c0
output_sha256: 998952528ef0dae5c370f216846328323d12de483235eb8147669a713f3ede05
hash_basis: raw LF bytes
---

# THM-3996 -- etale node-address balance forces cycles or nonproperness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let

```text
F:X=A2 -> Y=A2
```

be a dominant etale polynomial morphism. Let `N subset Y` be an integral
affine curve with an ordinary node `o`, smooth near every other point in a
chosen neighborhood of `o`. Write

```text
nu:N_tilde -> N,              nu^(-1)(o)={o_+,o_-}      (1)
```

for its normalization at the node. Denote by `S_F` the nonproperness locus of
`F` (over `C`, the usual Jelonek set).

The theorem first treats any node `o notin S_F`, then applies the result to
THM-3992's forced nodal cubic.

## 1. The complete oriented node-address graph

Assume `o notin S_F`. Choose a connected Zariski-open neighborhood `V` of `o`
on which `F` is proper. Since an etale morphism is quasi-finite, the restriction

```text
F^(-1)(V) -> V                                              (2)
```

is finite etale. Put

```text
N^o=N intersect V,              N_tilde^o=nu^(-1)(N^o),
Z^o=F^(-1)(N^o),
W=Z^o times_(N^o) N_tilde^o.                              (3)
```

The map `W->N_tilde^o` is finite etale. The normalized curve
`N_tilde^o` is connected, so each connected component `W_i` is a finite
etale cover of a constant positive degree `r_i`.

Define a directed multigraph `Gamma_o(F,N)` as follows:

```text
vertices: connected components W_i;
edges:    source points q in F^(-1)(o);
q: i->j  when the lift of q over o_+ lies in W_i and the
          lift of q over o_- lies in W_j.                 (4)
```

Loops and parallel edges are retained. This is the **complete** address graph:
it uses every source point over `o`, not only intersections with one selected
component.

## 2. Finite-etale degree is the exact balance law

For every vertex `W_i`, the fibre over `o_+` has exactly `r_i` points. Each is
the plus lift of a unique source address `q`, so

```text
outdegree(W_i)=r_i.                                      (5)
```

The same argument over `o_-` gives

```text
indegree(W_i)=r_i.                                       (6)
```

Thus every vertex is balanced.

This local equality has a global graph consequence. Contract every strongly
connected component of a weak component of `Gamma_o` to form its condensation
DAG. Summing `(5)--(6)` over each contracted block shows that its total
incoming and outgoing edge counts agree. A nontrivial finite DAG has a source;
balance would force that source to have no outgoing edge, while weak
connectivity forces an outgoing edge, a contradiction.
Hence no edge joins distinct strongly connected blocks. Every weak component
is strongly connected, and every directed edge lies on a directed cycle.

Consequently a nonempty complete address graph cannot be a forest, where a
loop is a one-cycle and two parallel oppositely directed edges form a
two-cycle. The contrapositive is the exact nonproperness certificate:

```text
Gamma_o(F,N) nonempty, complete, and cycle-free
                                      ==>  o in S_F.      (7)
```

## 3. Application to the forced THM-3992 node

In THM-3992 put

```text
N_a: C^2-A^3+(3a^2/4)A+a^3/4=0,             a!=0,
o=(-a/2,0).                                             (8)
```

Its normalization is the affine line

```text
X |-> (A,C)=(X^2+a, X^3+(3a/2)X),
nu^(-1)(o): X^2=-3a/2.                                 (9)
```

The source line `L={t=0}` maps to `(9)` with degree one. THM-3992 proves that
the two roots in `(9)` are simple source clutch points and that locally

```text
F^(-1)(N_a): t Q(x,t)=0;                               (10)
```

at each clutch, one branch is `L` and the other is a companion germ.

Orient `(1)` so that the degree-one normalization vertex `W_L` supplies the
plus branch at one clutch and the minus branch at the other. The two known
addresses are therefore one edge leaving `W_L` and one edge entering it.

If their companion germs have distinct global owners and `o notin S_F`, the
first companion owner has an unmatched incoming edge and the second an
unmatched outgoing edge. Balance `(5)--(6)` forces at least one further source
address over `o`; a single edge between those owners makes the sharp triangle.
Equivalently,

```text
distinct companion owners
  ==>  an additional node address or o in S_F.          (11)
```

If the two known addresses are the entire connected packet containing `L`
and `o notin S_F`, balance forces their other endpoints to be the same vertex
`W_C`. The packet is

```text
W_L -> W_C -> W_L,                                     (12)
```

a two-edge multicycle. Thus a full two-address audit forces one companion
owner; a two-owner audit without further addresses forces the node itself
into the nonproperness locus.

This is the correction in MISTAKE-469. The roots of `Q(x,0)` enumerate only
clutches lying on `L`; they are not a census of `F^(-1)(o)`. Moreover these
curves lie inside the Keller source. THM-3951's forest theorem concerns
completion-boundary primes and cannot be applied to this interior graph.

## 4. Proper packets are literal directed cycles

There is a stronger classification when the connected pullback packet
`Z_L->N_a` containing `L` is proper. It is then finite etale. Base change to
the normalization `(9)` makes every `W_i` a connected finite etale cover of
`A1`.

Such a cover has degree one in characteristic zero. Indeed, compactify it to
a degree-`d` map of smooth projective curves to `P1`. If `r` points lie over
infinity, all ramification is there and Riemann--Hurwitz gives

```text
2g-2=-2d+sum(e_j-1)=-2d+(d-r)=-d-r.                   (13)
```

Since `g>=0`, `d>=1`, and `r>=1`, equation `(13)` forces `d=r=1` and `g=0`.
Thus every vertex of the connected packet has indegree and outdegree one, so
the packet is one directed cycle.

If its total degree is `d`, it has `d` vertices and `d` node addresses. In
the THM-3992 packet, `d=2` forces the common companion owner in `(12)`.
Distinct companion owners require `d>=3` and exactly `d-2` further node
addresses beyond the two on `L`.

## 5. Sharp models and failure boundary

The cycle conclusion is sharp in every degree. Put

```text
R={h in k[z]: h(1)=h(-1)}
 =k[z^2,z(z^2-1)],
Spec R: v^2=u(u-1)^2.                                  (14)
```

For `d>=2`, with indices modulo `d`, define

```text
B_d={(h_0,...,h_(d-1)) in k[z]^d:
     h_i(1)=h_(i+1)(-1) for every i}.                  (15)
```

The diagonal map `R->B_d` is finite of degree `d`: `B_d` is an `R`-submodule
of the finite normalization product. Away from the node it is the disjoint
union of `d` normalization sheets. At each of its `d` source nodes, the two
completed branches map isomorphically to the two target branches, so
`Spec B_d->Spec R` is finite etale. Its source is connected and its address
graph is the directed `d`-cycle. The case `d=3` has distinct companions at a
chosen vertex and exactly one additional address, proving sharpness of `(11)`.

Delete that third source node from the `d=3` cover. The remaining morphism is
etale and quasi-finite, but it is nonproper over `o`; its address graph is the
two-edge path. This is the minimal hostile showing why properness/completeness
cannot be omitted.

For `d=2` there is an embedded plane-curve control:

```text
Z: t(t-x^2+1)=0,
u=x^2,                    v=x(x^2-1-2t).               (16)
```

Both components map to `(14)`, meet transversely at `(x,t)=(+/-1,0)`, and
their two addresses form the directed two-cycle. This is only a curve model;
it is not an ambient Keller map.

## 6. Exact residual and scope

THM-3996 does not determine which alternative in `(11)` occurs for THM-3992.
It does not prove that the two known addresses exhaust the node fibre, that
the node belongs to `S_F`, that a companion divisor is irreducible, or that a
cycle carries a nonzero class. It proves that these are the exhaustive next
sidecars:

```text
complete node-address census,
component ownership and orientation,
nonproperness of the target node,
cycle class only after the graph is complete.           (17)
```

A cycle is not a Keller obstruction by itself, and a Jelonek value is expected
for a hypothetical nonproper Keller map. The gain is that the former
forest-versus-cycle guess has become the exact dichotomy `(11)`, with one
additional-address minimum and a sharp hostile. **QED.**

## Reproduction

```bash
python3 04-computation/jc2_etale_node_address_balance_thm3996.py
python3 -O 04-computation/jc2_etale_node_address_balance_thm3996.py
sha256sum 04-computation/jc2_etale_node_address_balance_thm3996.py \
  05-knowledge/results/jc2_etale_node_address_balance_thm3996.out
python3 agents/check_docs.py
```
