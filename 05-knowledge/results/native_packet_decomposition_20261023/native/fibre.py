"""Complete fixed-margin/current fibre, actual28-source pin-complete menu."""
import json,time,sys
from collections import defaultdict,Counter,deque
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent

class DSU:
    def __init__(self,n):self.p=list(range(n));self.z=[1]*n
    def root(self,i):
        while self.p[i]!=i:self.p[i]=self.p[self.p[i]];i=self.p[i]
        return i
    def join(self,i,j):
        i,j=self.root(i),self.root(j)
        if i==j:return False
        if self.z[i]<self.z[j]:i,j=j,i
        self.p[j]=i;self.z[i]+=self.z[j];return True

def main():
    D=json.loads((HERE/'domains.json').read_text());sources=sorted(r['source'] for r in D['sources']);by={r['source']:r for r in D['sources']}
    codes=[(by[s]['code'][0],sum(b<<i for i,b in enumerate(by[s]['code'][1:]))) for s in sources]
    old=json.loads((HERE.parents[1]/'injective_native_blocks_20261023/search/larger16_selected.json').read_text())['circuits'][0]
    initial=sum(1<<sources.index(s) for s in old['negative']);terminal=sum(1<<sources.index(s) for s in old['positive'])
    def charge(mask):
        x=y=0
        for i,c in enumerate(codes):
            if mask>>i&1:x+=c[0];y^=c[1]
        return mask.bit_count(),x%128,y
    wanted=charge(initial)
    exact_X=sum(by[s]['X'] for s in sources if initial>>sources.index(s)&1)
    if charge(terminal)!=wanted:raise RuntimeError('old endpoints have unequal code')
    halves=[]
    for offset in (0,14):
        table=[]
        for mask in range(1<<14):
            x=y=0
            for i in range(14):
                if mask>>i&1:x+=codes[offset+i][0];y^=codes[offset+i][1]
            table.append((mask.bit_count(),x%128,y,mask))
        halves.append(table)
    right=defaultdict(list)
    for k,x,y,mask in halves[1]:right[k,x,y].append(mask)
    states=[]
    for k,x,y,mask in halves[0]:
        for hi in right.get((wanted[0]-k,(wanted[1]-x)%128,wanted[2]^y),[]):states.append(mask|(hi<<14))
    if '--exact' in sys.argv:
        states=[m for m in states if sum(by[s]['X'] for i,s in enumerate(sources) if m>>i&1)==exact_X]
    states.sort();n=len(states);index={m:i for i,m in enumerate(states)}
    if len(index)!=n:raise RuntimeError('MITM duplicate')
    print('fibre',n,'wanted',wanted,flush=True)
    # Prim on the entire complete weighted graph, not a proposed move alphabet.
    distances=[29]*n;parents=[-1]*n;used=[False]*n;distances[index[initial]]=0;tree=[]
    for _ in range(n):
        v=min((i for i in range(n) if not used[i]),key=lambda i:(distances[i],i));used[v]=True
        if parents[v]>=0:tree.append([parents[v],v,distances[v]])
        for j in range(n):
            if not used[j]:
                d=(states[v]^states[j]).bit_count()
                if d<distances[j]:distances[j]=d;parents[j]=v
    maxcost=max(c for _,_,c in tree);adj=defaultdict(list)
    for a,b,c in tree:adj[a].append((b,c));adj[b].append((a,c))
    src,dst=index[initial],index[terminal];todo=deque([src]);prev={src:None}
    while dst not in prev:
        a=todo.popleft()
        for b,c in adj[a]:
            if b not in prev:prev[b]=(a,c);todo.append(b)
    path=[];v=dst
    while prev[v] is not None:
        p,c=prev[v];path.append([states[p],states[v],c]);v=p
    path.reverse();thresholds=[]
    for b in range(2,maxcost+1,2):
        ds=DSU(n)
        for i,j,c in tree:
            if c<=b:ds.join(i,j)
        sizes=Counter(ds.root(i) for i in range(n))
        thresholds.append(dict(budget=b,components=len(sizes),component_sizes=sorted(sizes.values()),
                               endpoints_connected=ds.root(src)==ds.root(dst)))
    endpoint_budget=max(c for _,_,c in path)
    # Shortest event route at the exact minimax budget, over all legal endpoints.
    todo=deque([src]);prev={src:None}
    while dst not in prev:
        i=todo.popleft()
        for j,m in enumerate(states):
            if j not in prev and (states[i]^m).bit_count()<=endpoint_budget:prev[j]=i;todo.append(j)
    route=[];v=dst
    while v is not None:route.append(states[v]);v=prev[v]
    route.reverse()
    outside=((1<<28)-1)^(initial|terminal)
    costs=[sum(u!=v for u,v in zip(by[s]['old_before'],by[s]['old_after'])) for s in sources]
    def work(mask):return sum(c for i,c in enumerate(costs) if mask>>i&1)
    route_debt=[]
    sys.path.insert(0,str(HERE.parents[1]/'injective_native_blocks_20261023/search'))
    import explore as native_data
    target_X=native_data.X[1965]
    target_binary=sum((native_data.Y[j][1965]%2)<<(j-1) for j in range(1,5))
    for m in route:
        oldsum=sum(by[s]['X'] for i,s in enumerate(sources) if (m&~outside)>>i&1)
        extsum=sum(by[s]['X'] for i,s in enumerate(sources) if (m&outside)>>i&1)
        old_count=(m&~outside).bit_count();extra_count=(m&outside).bit_count();delta_old_count=old_count-8
        oldparity=charge((m&~outside))[2]^charge(initial)[2]
        extparity=charge(m&outside)[2]
        route_debt.append(dict(mask=m,old16_V_sum=oldsum,outside_V_sum=extsum,total_V_sum=oldsum+extsum,
                               old16_V_count=old_count,outside_V_count=extra_count,
                               source_charge_X_change_old16=2*(oldsum-exact_X),source_charge_X_change_outside=2*extsum,
                               target_U_degree_delta_old16=-2*delta_old_count,target_U_degree_delta_outside=-2*extra_count,
                               U_current_change_old16=2*(oldsum-exact_X-target_X*delta_old_count),U_current_change_outside=2*(extsum-target_X*extra_count),
                               U_current_change_total=2*(oldsum+extsum-exact_X),
                               V_current_changes='negatives of the U-owned changes',
                               old16_binary_source_charge=oldparity,outside_binary_source_charge=extparity,
                               old16_binary_residual=oldparity^(target_binary if delta_old_count%2 else 0),outside_binary_residual=extparity^(target_binary if extra_count%2 else 0)))
    result=dict(sources=sources,domain='all actual paths retaining best inherited endpoint+indexed arc pins and targetU/V',
                V_count=wanted[0],V_code=list(wanted[1:]),states=states,initial=initial,terminal=terminal,
                complete_graph_MST=tree,global_sharp_slot_budget=maxcost,thresholds=thresholds,
                endpoint_sharp_slot_budget=endpoint_budget,shortest_endpoint_route=route,
                route_supports=[[sources[j] for j in range(28) if (a^b)>>j&1] for a,b in zip(route,route[1:])],
                route_owner_clock_occurrences=[work(a^b) for a,b in zip(route,route[1:])],
                MST_owner_clock_occurrences=[work(states[i]^states[j]) for i,j,c in tree],
                outside_sources=[sources[j] for j in range(28) if outside>>j&1],
                outside_V_sets=[[sources[j] for j in range(28) if (m&outside)>>j&1] for m in route],
                fixed_integer_X_sum=exact_X if '--exact' in sys.argv else None,route_owned_debt_ledger=route_debt)
    (HERE/('fibre_exact_receipt.json' if '--exact' in sys.argv else 'fibre_receipt.json')).write_text(json.dumps(result,indent=2))
    print('global sharp',maxcost,'endpoint sharp',endpoint_budget,'route events',len(route)-1,'thresholds',[(r['budget'],r['components']) for r in thresholds])
    print('route',result['route_supports']);print('outside states',result['outside_V_sets'])
if __name__=='__main__':main()
