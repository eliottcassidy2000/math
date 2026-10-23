"""Finite independent controls for symmetry-transported residual certificates."""
from itertools import permutations, product
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
checks=0
def require(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)

def potential(n, edges, x, costs):
    arcs=[]
    for k,(u,v) in enumerate(edges):
        if x[k]<1: arcs.append((u,n+v,costs[k]))
        if x[k]>0: arcs.append((n+v,u,-costs[k]))
    reach=[[u==v for v in range(2*n)] for u in range(2*n)]
    for u,v,c in arcs: reach[u][v]=True
    for k in range(2*n):
        for u in range(2*n):
            for v in range(2*n): reach[u][v] |= reach[u][k] and reach[k][v]
    h={}
    for root in range(2*n):
        if root in h: continue
        scc={v for v in range(2*n) if reach[root][v] and reach[v][root]}
        h[root]=0
        pending=[root]
        for u in pending:
            for a,v,c in arcs:
                if a!=u or v not in scc: continue
                if v not in h:
                    h[v]=h[u]+c
                    pending.append(v)
                elif h[v]!=h[u]+c:
                    return None
    return h

def table(q,edges):
    return tuple(int(q[u]==v) for u,v in edges)

def main():
    n=3
    edges=list(product(range(n),repeat=2))
    tables=[table(q,edges) for q in permutations(range(n))]
    # Complete independent census of ternary integer costs on K3,3.
    histogram={'constant':0,'nonconstant':0}
    for costs in product((-1,0,1),repeat=9):
        values={sum(a*b for a,b in zip(costs,x)) for x in tables}
        h=potential(n,edges,tables[0],costs)
        require((h is not None)==(len(values)==1),'residual criterion against all six tables')
        histogram['constant' if h is not None else 'nonconstant']+=1
    row_p=(2,-1,3); col_p=(5,0,-2)
    costs=tuple(col_p[v]-row_p[u] for u,v in edges)
    x=tables[0]
    h=potential(n,edges,x,costs)
    require(h is not None,'positive potential control')
    transported=[]
    for shift in range(3):
        phi=lambda u:(u+shift)%3
        edge_map={k:edges.index((phi(u),phi(v))) for k,(u,v) in enumerate(edges)}
        new_cost=[0]*9; new_x=[0]*9
        for k,j in edge_map.items(): new_cost[j]=costs[k]; new_x[j]=x[k]
        new_h={phi(v) if v<n else n+phi(v-n):a for v,a in h.items()}
        for k,(u,v) in enumerate(edges):
            if new_x[k]<1: require(new_cost[k]==new_h[n+v]-new_h[u],'transport forward potential')
            if new_x[k]>0: require(-new_cost[k]==new_h[u]-new_h[n+v],'transport backward potential')
        values={sum(a*b for a,b in zip(new_cost,y)) for y in tables}
        require(values=={sum(a*b for a,b in zip(costs,x))},'common transported baseline')
        transported.append(dict(shift=shift,baseline=next(iter(values))))
    def off(q):return sum(i!=j for i,j in enumerate(q))
    phases=[tuple((i+s)%3 for i in range(3)) for s in range(3)]
    legitimate=[off(q)-off(q) for q in phases]
    union_values=sorted({off(q)-off(r) for q in permutations(range(3)) for r in permutations(range(3))})
    require(legitimate==[0,0,0],'shared odd phase')
    require(-3 in union_values and 3 in union_values,'premature-union hostile')
    # Two unique branches are not enough: projected baselines still differ.
    separate=[]
    for shift in (0,1):
        branch_edges=[(i,(i+shift)%3) for i in range(3)]
        branch_x=(1,1,1)
        branch_costs=tuple(int(u!=v) for u,v in branch_edges)
        require(potential(3,branch_edges,branch_x,branch_costs) is not None,
                'unique branch passes residual test')
        separate.append(sum(a*b for a,b in zip(branch_x,branch_costs)))
    require(separate==[0,3],'literal branch baseline values')
    require(len(set(separate))==2,'baseline hostile')
    output=dict(scope='Unit-capacity bipartite tables; not a native tournament-port census.',
                all_K33_ternary_costs=histogram,transported=transported,
                shared_phase_branches=3,unrestricted_union_tables=36,
                union_projected_values=union_values,distinct_unique_branch_baselines=separate,
                checks=checks)
    text=json.dumps(output,sort_keys=True,indent=2)+'\n'
    target=HERE/'residual.json'
    if '--write' in sys.argv:target.write_text(text,encoding='utf-8')
    else:require(json.loads(target.read_text())==json.loads(text),'frozen receipt')
    print(json.dumps(dict(checks=checks,cost_census=histogram)))

if __name__=='__main__':main()
