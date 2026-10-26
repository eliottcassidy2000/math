"""Direct ordinary DAG with forced unaffordable/terminal skips contracted.

This is elementary administrative contraction, not a further defect quotient.
The unchanged full engine in reverse_ordinary.py is the independent oracle.
"""
from collections import defaultdict
from math import gcd
import json
import reverse_ordinary as base


need=base.need


def canonical_vertex(ell,s,D):
    if not s:
        need(not D,'terminal defect is empty')
        return (-1,0,())
    affordable=s if s%2 else s-1
    active=min(ell,affordable)
    need(active>=1,'positive remainder has an active stage')
    need(len(D)==(active+1)//2,'forced skip preserves exact affordable defect')
    return (active,s,D)


def build_contracted_dag(n):
    """Build only active vertices, without first building the full DAG."""
    base.valid_n(n)
    if not n: return base.build_dag(0)
    top=n if n%2 else n-1
    initial=(top,n,(0,)*((top+1)//2))
    buckets=defaultdict(set); buckets[top].add(initial)
    raw={}
    for ell in range(top,0,-2):
        for vertex in sorted(buckets[ell]):
            _,s,D=vertex
            need(ell<=s,'no forced skip vertex is constructed')
            b=base.lower_bound(D); selected=D[(ell-1)//2]; outgoing=[]
            for k in range(s//ell+1):
                rest=s-ell*k
                if ell==1 and rest: continue
                Dp=tuple(D[(d-1)//2]+k*(d-gcd(ell,d))
                         for d in range(1,min(ell-2,rest)+1,2))
                bp=base.lower_bound(Dp)
                shift=(k*(ell-1)//2+ell*k*(k-1)//2-base.vfact(k)
                       +k*(ell*b-selected)+rest*(k+b-bp))
                need(shift>=0,'direct integral contracted edge')
                target=canonical_vertex(ell-2,rest,Dp)
                buckets[target[0]].add(target)
                outgoing.append((target,ell,k,n-s,shift))
            raw[vertex]=outgoing
    terminal=(-1,0,()); raw[terminal]=[]
    vertices=sorted(raw,key=lambda v:(-v[0],v[1],v[2]))
    need(vertices[0]==initial and vertices[-1]==terminal,'contracted topological endpoints')
    index={v:i for i,v in enumerate(vertices)}
    edges=[[base.Edge(index[target],ell,k,u,J) for target,ell,k,u,J in raw[v]] for v in vertices]
    V=[None]*len(vertices);V[-1]=0
    for v in range(len(vertices)-2,-1,-1):
        need(all(e.target>v for e in edges[v]),'strict active-stage descent')
        V[v]=min(e.shift+V[e.target] for e in edges[v])
    P=[None]*len(vertices);P[0]=0
    for v,row in enumerate(edges):
        need(P[v] is not None,'direct structural reachability')
        for e in row:
            cost=e.shift+V[e.target]-V[v]
            need(cost>=0,'contracted reduced edge cost')
            candidate=P[v]+cost
            P[e.target]=candidate if P[e.target] is None else min(P[e.target],candidate)
    layers=[[index[v] for v in vertices if v[0]==ell] for ell in sorted(buckets,reverse=True)]
    need(V[0]==(n//2 if n>=3 else 0),'contracted baseline')
    return base.OrdinaryDAG(n,vertices,edges,layers,V,P)


def audit_contraction(full,small):
    """Check actual forced chains and every surviving local transition."""
    fi={v:i for i,v in enumerate(full.vertices)}
    si={v:i for i,v in enumerate(small.vertices)}
    odd=base.odd_factorials(full.n)
    mapping=[]; forced=0; strict_prefix_loss=0
    for v,vertex in enumerate(full.vertices):
        canonical=canonical_vertex(*vertex) if vertex[1] else (-1,0,())
        need(canonical in si,'every full vertex has a surviving image')
        target=si[canonical]; mapping.append(target)
        if vertex != canonical:
            need(len(full.edges[v])==1,'a contracted vertex is deterministic')
            e=full.edges[v][0]
            need(e.k==0 and e.shift==0 and base.edge_factor(e,odd)==1,'forced unit edge')
            need(full.vertices[e.target][1:]==vertex[1:],'forced skip retains S,D')
            need(full.potential[v]==full.potential[e.target],'suffix potential constant along skip')
            forced+=1
        need(full.potential[v]==small.potential[target],'V descends along every forced chain')
        strict_prefix_loss+=full.prefix_cost[v]>small.prefix_cost[target]
    need(len(full.vertices)-len(small.vertices)==forced,'exact skipped vertex count')
    for sv,vertex in enumerate(small.vertices):
        fv=fi[vertex]
        need(full.prefix_cost[fv]==small.prefix_cost[sv],'P agrees at surviving aggregate only')
        projected=[(mapping[e.target],e.ell,e.k,e.u,e.shift) for e in full.edges[fv]]
        actual=[(e.target,e.ell,e.k,e.u,e.shift) for e in small.edges[sv]]
        need(projected==actual,'every direct edge is its full-DAG contraction')
    return mapping,{'forced_vertices_removed':forced,
                    'skipped_prefix_minimum_strictly_above_aggregate':strict_prefix_loss}


def controls():
    samples=[];count_checks=modular_checks=coordinate_checks=0
    prefix_hostile=False
    for n in list(range(51))+[60,80,100]:
        full=base.build_dag(n); small=build_contracted_dag(n)
        mapping,detail=audit_contraction(full,small)
        prefix_hostile|=bool(detail['skipped_prefix_minimum_strictly_above_aggregate'])
        a,z,Y,Z,_=base.exact_audit(full)
        b,w,Ys,Zs,_=base.exact_audit(small)
        need((a,z)==(b,w),'exact unrestricted count and normalized count')
        full_indices={v:i for i,v in enumerate(full.vertices)}
        for i,vertex in enumerate(small.vertices):
            j=full_indices[vertex]
            need((Ys[i],Zs[i])==(Y[j],Z[j]),'surviving exact merged weights')
            coordinate_checks+=1
        count_checks+=1
        for K in range(1,9):
            x,work,_=base.normalized_modular(full,K)
            y,smallwork,_=base.normalized_modular(small,K)
            need(x==y,'contracted normalized residue')
            modular_checks+=1
            if n in (30,50,60,80,100) and K in (2,8):
                samples.append({'n':n,'K':K,'full_input_states':len(full.vertices)-1,
                    'contracted_input_states':len(small.vertices)-1,
                    'full_edges':sum(map(len,full.edges)),
                    'contracted_edges':sum(map(len,small.edges)),
                    'full_active_input_states':full.statistics(K)['active_input_states'],
                    'contracted_active_input_states':small.statistics(K)['active_input_states'],
                    'full_active_updates':work['edge_updates'],
                    'contracted_active_updates':smallwork['edge_updates'],**detail})
    need(prefix_hostile,'skip-source prefix minimum can exceed canonical aggregate minimum')
    # Forced skip contraction is valid even though bare stage erasure makes a loop.
    six=base.build_dag(6); source=six.vertices.index((3,1,(0,)))
    edge=six.edges[source][0]
    need(six.vertices[edge.target]==(1,1,(0,)),'actual stage-erasure loop hostile')
    need(canonical_vertex(*six.vertices[source])==canonical_vertex(*six.vertices[edge.target]),
         'proved skip contraction resolves that specific loop')
    eight=base.build_dag(8); active=eight.vertices.index((3,3,(0,2)))
    zero_cycle=next(e for e in eight.edges[active] if e.k==0)
    need(len(eight.edges[active])==2 and zero_cycle.shift==3,
         'an affordable zero-cycle edge is not a forced zero-cost skip')
    return {'status':'PASS','ordinary_orders':{'consecutive':[0,50],'additional':[60,80,100]},
            'exact_count_comparisons':count_checks,'normalized_modular_comparisons':modular_checks,
            'surviving_prefix_coordinate_checks':coordinate_checks,'normalized_precision_bits':[1,8],
            'samples':samples,'scope':'Direct deterministic zero-cost unit-skip contraction only; stage remains at every active branching vertex. No new defect equivalence or clone-excess model.'}


if __name__=='__main__': print(json.dumps(controls(),indent=2,sort_keys=True))
