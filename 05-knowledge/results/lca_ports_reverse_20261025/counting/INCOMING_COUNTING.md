# Incoming ordinary skip contraction: exact integration map

2026-10-25. **RECOVERED CONCURRENT RESULT:** commit `9e79465c2c`,
[reverse implementation](../../session64_lca_reverse_valuation_20261026/session64-counting/REVERSE_IMPLEMENTATION.md)
and [direct skip contraction](../../session64_lca_reverse_valuation_20261026/session64-counting/SKIP_CONTRACTION.md),
with separate [reverse](../../session64_lca_reverse_valuation_20261026/session64-audit-counting/REVERSE_AUDIT.md)
and [contraction](../../session64_lca_reverse_valuation_20261026/session64-audit-counting/CONTRACTION_AUDIT.md)
audits. These supplied implementations are reused, not duplicated here.
The small executable comparison below is **FINITE-EXACT** and passes
normally and under `python -O`.

## Same normalization, two additional operations

The incoming full DAG has the same `(stage,S,D)` states, odd T factors,
original costs J, suffix potential V, and normalized scalar
Z_v=2^[V(v)-B]*Y_v as this packet's `OrdinaryDAG`. Our terminal `(0,0,())`
is its `(-1,0,())`; this terminal convention is the only state-map change.
At N=100 both full graphs have13,445 vertices including the terminal and
50,044 edges. The inherited exact ordinary merger remains unchanged.

Incoming P is the minimum **already-reduced prefix cost**: P(initial)=0
and P(w)=min(P(v)+J'(v,w)), where J'=J+V(w)-V(v). It is not the minimum
sum of original J values. Since every vertex has a zero reduced-cost
suffix, retaining P(v)+J'(v,w)<K keeps exactly the union of complete
paths with original valuation below B+K. This is a structural valuation
filter after full construction. It differs from our local J'>=K pruning
followed by suppression of actually zero source residues; modular carries
can cause the latter to skip further edges. Update counts from those
two evaluators must not be compared as if their pruning rules were equal.

The second addition is a genuine reduction of the constructed graph:

    c(ell,0,())=terminal,
    c(ell,S,D)=(min(ell,largest odd<=S),S,D), S>0.

Removed vertices have only k=0 available. Their affordable D coordinates,
u and b stay unchanged, giving T=1,J=0. Delete the unique forced chains,
or uniquely reinsert them, to obtain a weighted complete-path bijection.
The incoming builder canonicalizes a destination immediately and never
constructs the removed vertices. At N=100 it has10,548 vertices including
the terminal and47,147 edges. This removes administrative skips, not
active defect distinctions, and is implemented only for ordinary e=0.

## Quantities that do and do not descend

V is constant along each forced chain; consequently B and all retained
reduced edge costs agree under c. P and accumulated Y/Z do not descend
pointwise: a canonical destination can receive additional histories.
They agree at surviving full aggregates, not at every removed source.
At N=6, the removed `(1,0,())` has P=1 whereas its terminal image has P=0.
Summing all old values in a contraction fibre would also count a path
once per skipped vertex. The valid operation is path/edge contraction.

Not every k=0 edge is removable. At N=8 the active `(3,3,(0,2))` has
two choices; its k=0 edge has J=3 because b drops from1 to0. Both
choices remain. Active stage is retained throughout, so this is not
indiscriminate stage erasure or a new active-state equivalence.

There is also an evaluator boundary: contracted edges can jump several
stages. Our `residue` method's two-consecutive-layer buffers cannot accept
the contracted graph unchanged. Use the incoming audited full-vertex
accumulator or a correctly scheduled collection of pending stage buckets.
No integration of that alternative storage scheme is claimed here.

## Scope and direct comparison

The published audits already supply unmerged prefix reconstruction,
partition path bijections, unrestricted permutation-action Burnside, and
large exact-count controls. The single block below does not repeat those
large audits. For N=0..20 and N=100, it compares this packet's actual
structural graph with the incoming direct constructor, independently
following forced chains and comparing every surviving labelled edge.
It checks V everywhere, P at surviving aggregates only, and both exact
Y/Z at survivors for the small orders. Normalized residues are compared
at K=1,2,6 for the small orders and K=6 at100. For the final arithmetic
comparison it applies the same incoming P filter to both graphs.

Both replay modes pass 64 normalized comparisons, 14,036 full-vertex
potential checks and 47,860 projected-edge checks. At N=100, K=6 the common
residue is 44; using the same P filter gives 2,318 full-graph updates and
1,009 contracted-graph updates. The 440 strict prefix-loss cases at that
order independently confirm the pointwise-descent boundary.

Run the embedded block from the repository root; it imports the published
implementations but neither of their main audit programs. The same block
works under `python -O`. Graph counts include the common terminal unless
explicitly labelled otherwise. No timing, memory, asymptotic recurrence,
or signed-clone conclusion is inferred.

```python
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import json
import sys

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def load(name,path):
    spec=spec_from_file_location(name,path)
    module=module_from_spec(spec); sys.modules[name]=module
    spec.loader.exec_module(module)
    return module

root=Path('05-knowledge/results')
local=load('local_reverse_integration',root/'lca_ports_reverse_20261025/counting/reverse.py')
incoming=root/'session64_lca_reverse_valuation_20261026/session64-counting'
base=load('reverse_ordinary',incoming/'reverse_ordinary.py')
contract=load('incoming_direct_contraction',incoming/'contracted_ordinary.py')

def vertex(state):
    return (-1,0,()) if state.stage==0 else (state.stage,state.remaining,state.defect)

def reduced_prefix(full):
    P=[None]*len(full.states); P[0]=0
    for i,row in enumerate(full.edges):
        need(P[i] is not None,'structural reachability')
        for edge in row:
            value=P[i]+full.reduced_cost(i,edge)
            P[edge.target]=value if P[edge.target] is None else min(P[edge.target],value)
    return P

orders=list(range(21))+[100]
vertex_checks=edge_checks=coordinate_checks=modular_checks=0
removed_total=strict_prefix_total=0; sample=None
for n in orders:
    full=local.OrdinaryDAG(n); small=contract.build_contracted_dag(n)
    fv=[vertex(s) for s in full.states]; fi={v:i for i,v in enumerate(fv)}
    si={v:i for i,v in enumerate(small.vertices)}; P=reduced_prefix(full)
    mapping=[]; forced=0; strict=0
    for i,state in enumerate(full.states):
        at=i
        while at!=full.terminal and (full.states[at].remaining==0 or
                                     full.states[at].stage>full.states[at].remaining):
            need(len(full.edges[at])==1,'removed vertex has one edge')
            edge=full.edges[at][0]; target=full.states[edge.target]
            need(edge.multiplicity==0 and edge.cost==0,'removed edge is a zero-cost unit')
            need((target.remaining,target.defect)==(full.states[at].remaining,full.states[at].defect),
                 'forced skip retains all response data')
            at=edge.target
        canonical=fv[at]
        need(canonical in si,'independently followed canonical image exists')
        image=si[canonical]; mapping.append(image)
        need(full.potential[i]==small.potential[image],'V descends at every original vertex')
        need(P[i]>=small.prefix_cost[image],'canonical aggregate can improve reduced prefix minimum')
        forced+=at!=i; strict+=P[i]>small.prefix_cost[image]
        vertex_checks+=1
    need({fv[i] for i in range(len(fv)) if fv[i] in si}==set(small.vertices),'exact surviving set')
    need(len(full.states)-len(small.vertices)==forced,'exact contraction count')
    for j,v in enumerate(small.vertices):
        i=fi[v]
        need(P[i]==small.prefix_cost[j],'P agrees at surviving aggregates')
        expected=[(mapping[e.target],full.states[i].stage,e.multiplicity,n-full.states[i].remaining,e.cost)
                  for e in full.edges[i]]
        actual=[(e.target,e.ell,e.k,e.u,e.shift) for e in small.edges[j]]
        need(expected==actual,'each retained weighted edge is exact')
        edge_checks+=len(actual)
    # Reuse the incoming evaluator on the local full graph with the same P filter.
    converted=base.OrdinaryDAG(n,fv,
        [[base.Edge(e.target,full.states[i].stage,e.multiplicity,n-full.states[i].remaining,e.cost)
          for e in row] for i,row in enumerate(full.edges)],full.layers,full.potential,P)
    if n<=20:
        total,normalized,Y,Z,_=base.exact_audit(converted)
        stotal,snormalized,Ys,Zs,_=base.exact_audit(small)
        need((total,normalized)==(stotal,snormalized),'small exact counts')
        for j,v in enumerate(small.vertices):
            need((Ys[j],Zs[j])==(Y[fi[v]],Z[fi[v]]),'surviving exact aggregates')
            coordinate_checks+=1
    for K in ((1,2,6) if n<=20 else (6,)):
        ours,_=full.residue(K)
        incoming_full,fw,_=base.normalized_modular(converted,K)
        contracted,cw,_=base.normalized_modular(small,K)
        need(ours==incoming_full==contracted,'same normalized residue')
        modular_checks+=1
        if n==100:
            sample={'n':n,'K':K,'residue':ours,
                    'full_vertices_including_terminal':len(full.states),
                    'contracted_vertices_including_terminal':len(small.vertices),
                    'full_edges':sum(map(len,full.edges)),
                    'contracted_edges':sum(map(len,small.edges)),
                    'same_P_filter_full_updates':fw['edge_updates'],
                    'same_P_filter_contracted_updates':cw['edge_updates'],
                    'removed_vertices_with_strict_prefix_loss':strict}
    removed_total+=forced; strict_prefix_total+=strict

six=local.OrdinaryDAG(6); P6=reduced_prefix(six)
i=next(i for i,s in enumerate(six.states) if (s.stage,s.remaining,s.defect)==(1,0,()))
need(P6[i]==1 and P6[six.terminal]==0,'pointwise P descent hostile')
eight=local.OrdinaryDAG(8)
i=next(i for i,s in enumerate(eight.states) if (s.stage,s.remaining,s.defect)==(3,3,(0,2)))
zero=next(e for e in eight.edges[i] if e.multiplicity==0)
need(len(eight.edges[i])==2 and zero.cost==3,'active zero-choice must remain')
need(any(small.vertices[e.target][0]<small.vertices[i][0]-2
         for i,row in enumerate(small.edges) for e in row),'contracted edges jump stages')
print(json.dumps({'status':'PASS','small_orders':[0,20],'additional_order':100,
                  'full_vertex_potential_checks':vertex_checks,'projected_edge_checks':edge_checks,
                  'small_surviving_exact_coordinates':coordinate_checks,
                  'normalized_comparisons':modular_checks,'forced_vertices_removed':removed_total,
                  'strict_prefix_loss_cases':strict_prefix_total,'order100_sample':sample},indent=2))
```
