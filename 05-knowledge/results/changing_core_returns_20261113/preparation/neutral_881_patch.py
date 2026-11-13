from pathlib import Path
import importlib.util,json
from collections import deque
from heapq import heappop,heappush
from itertools import permutations,combinations
from math import gcd,lcm
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('preparation_independent_counter',HERE/'check.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
data=json.loads((HERE.parent/'selector/receipt.json').read_text(encoding='utf-8'))
out=[]
def exact_wheel(z):
    H,s,q,F=z;delta=F+s-H-q;d=F-q
    if delta%3==0:return dict(period=1,good=[],reason='three divides imbalance')
    primes=[];remaining=abs(delta);p=3
    while remaining%2==0:remaining//=2
    while p*p<=remaining:
        if remaining%p==0:
            primes.append(p)
            while remaining%p==0:remaining//=p
        p+=2
    if remaining>1:primes.append(remaining)
    period=2;obstructions=[]
    for p in primes:
        residues=[];a=1
        while a not in residues:residues.append(a);a=2*a%p
        period=lcm(period,len(residues))
        bad=[i for i,a in enumerate(residues) if (d*a+q)%p==0]
        obstructions.append((p,len(residues),bad))
    e.need(period<2000000,'explicit wheel size')
    good=[k for k in range(period) if (d*pow(2,k,3)+q-delta)%3==0 and all(k%order not in bad for p,order,bad in obstructions)]
    e.need(all((gcd(e.state_return(z,k)[0],e.state_return(z,k)[1])==1 and e.state_return(z,k)[1]%3==0 and e.state_return(z,k)[0]%3!=0)==(k%period in good) for k in range(1,81)),'independent wheel guard')
    return dict(period=period,good=good,obstructions=obstructions)
for source in data['sources']:
    T=tuple(source['rows']);U=tuple(source['prepared_rows']);r=source['root'];flips=source['arc_flips'];m=len(flips);end=(1<<m)-1
    states=[];graphs=[]
    for mask in range(1<<m):
        rows=list(T)
        for bit,(a,b) in enumerate(flips):
            if mask>>bit&1:rows[a]^=1<<b;rows[b]^=1<<a
        rows=tuple(rows);e.need(e.tournament(rows),'partial patch tournament')
        graphs.append(rows);states.append(e.state(rows,r));e.tables.cache_clear()
    e.need(graphs[-1]==U,'full patch endpoint')
    e.need(states[0]==tuple(source['root_profiles'][r]),'independent original root profile')
    e.need(states[-1]==tuple(source['prepared_profile']),'independent prepared root profile')
    neutral=[mask for mask,z in enumerate(states) if z[0]==881]
    def path_bfs(vertices,width):
        allowed=set(vertices);parent={0:None};queue=deque([0])
        while queue:
            a=queue.popleft()
            if a==end:break
            for b in allowed:
                if b not in parent and 0<(a^b).bit_count()<=width:parent[b]=a;queue.append(b)
        if end not in parent:return None
        path=[];a=end
        while a is not None:path.append(a);a=parent[a]
        return list(reversed(path))
    def minimax(monotone):
        costs={0:0};parents={0:None};queue=[(0,0)]
        while queue:
            cost,a=heappop(queue)
            if cost!=costs[a]:continue
            if a==end:break
            for bit in range(m):
                if monotone and a>>bit&1:continue
                b=a^(1<<bit);candidate=max(cost,abs(states[b][0]-881))
                if candidate<costs.get(b,10**30):costs[b]=candidate;parents[b]=a;heappush(queue,(candidate,b))
        path=[];a=end
        while a is not None:path.append(a);a=parents[a]
        return dict(max_absolute_count_excursion=costs[end],mask_path=list(reversed(path)),count_path=[states[a][0] for a in reversed(path)])
    minwidth=next(w for w in range(1,m+1) if path_bfs(neutral,w) is not None)
    e.need(all((states[a][0]%2)==1 for a in range(1<<m)),'odd Hamilton counts')
    arcs=list(combinations(range(len(T)),2));hcache={}
    def flip(rows,selected):
        U=list(rows)
        for a,b in selected:U[a]^=1<<b;U[b]^=1<<a
        return tuple(U)
    def H(rows):
        if rows not in hcache:
            hcache[rows]=e.tables(rows)[0][-1];e.tables.cache_clear()
        return hcache[rows]
    component={T};todo=deque([T]);neutral_edges=[]
    while todo:
        current=todo.popleft()
        for arc in arcs:
            neighbor=flip(current,[arc])
            if H(neighbor)!=881:continue
            a,b=arc
            twins=all(((current[a]>>v)&1)==((current[b]>>v)&1) for v in range(len(T)) if v not in arc)
            e.need(twins and r not in arc,'every neutral edge swaps nonroot adjacent twins')
            neutral_edges.append(dict(source_rows=current,arc=arc,target_rows=neighbor))
            if neighbor not in component:component.add(neighbor);todo.append(neighbor)
    component_states=[]
    for member in sorted(component):
        z=e.state(member,r);e.tables.cache_clear();e.need(z==states[0],'component fixes full root profile')
        component_states.append(dict(rows=member,state=z))
    two_arc=[];wheelcache={}
    for patch in combinations(arcs,2):
        candidate=flip(T,patch)
        if H(candidate)!=881:continue
        z=e.state(candidate,r);e.tables.cache_clear()
        if z not in wheelcache:wheelcache[z]=exact_wheel(z)
        wheel=wheelcache[z]
        row=dict(arcs=patch,rows=candidate,state=z,strong=e.strong(candidate),root_indegree=sum(bool(candidate[v]>>r&1) for v in range(len(T))),ready_defect=z[2]-z[1],wheel=wheel,entry_possible=bool(wheel['good']),source_sorted_outdegrees=sorted(a.bit_count() for a in T),target_sorted_outdegrees=sorted(a.bit_count() for a in candidate))
        if wheel['good']:
            k=next(k for k in range(1,wheel['period']+1) if k%wheel['period'] in wheel['good'])
            V,rr=e.returned(candidate,r,k);actual=e.state(V,rr)
            e.need(actual==e.state_return(z,k),'independent physical two-arc repair return')
            row['literal_return']=dict(k=k,rows=V,root=rr,state=actual)
        two_arc.append(row)
    global_probe=dict(single_arc_component=component_states,neutral_component_edges=neutral_edges,all_source_single_arc_checks=len(arcs),all_source_two_arc_checks=len(arcs)*(len(arcs)-1)//2,two_arc_H881_endpoints=two_arc,two_arc_entry_repairs=[a for a in two_arc if a['entry_possible']],graphs_counted=len(hcache))
    out.append(dict(global_probe=global_probe,parent_order=source['parent_order'],parent_rank=source['parent_rank'],root=r,arcs=flips,original_indegree=sum(bool(T[v]>>r&1) for v in range(len(T))),prepared_indegree=sum(bool(U[v]>>r&1) for v in range(len(U))),all_partial_patch_states=states,neutral_masks=neutral,neutral_single_arc_path=path_bfs(neutral,1),minimal_neutral_simultaneous_width=minwidth,minimal_width_neutral_path=path_bfs(neutral,minwidth),best_single_arc_excursion=minimax(False),best_monotone_single_arc_excursion=minimax(True)))
# Independently check full profiles of all285 rooted labelled controls via permutations and cuts.
controls=0
for n in range(1,5):
    for T in e.labelled(n):
        for r in range(n):
            e.need(e.state(T,r)==e.permstate(T,r),'independent literal full marked profile')
            controls+=1
        e.tables.cache_clear()
result=dict(status='PASS',scope='All32 partial subsets of each prescribed five-arc patch; exhaustive full H881 single-arc components; all two-arc endpoint patches from each source',literal_rooted_controls=controls,sources=out)
(HERE/('neutral_881_patch.json' if __debug__ else 'neutral_881_patch.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k!='sources'}))
for row in out:print(json.dumps({k:v for k,v in row.items() if k not in ('all_partial_patch_states','arcs','global_probe')}))

for row in out:
    g=row['global_probe'];print(json.dumps(dict(parent_order=row['parent_order'],component_size=len(g['single_arc_component']),two_arc_H881_endpoints=len(g['two_arc_H881_endpoints']),two_arc_entry_repairs=len(g['two_arc_entry_repairs']),graphs_counted=g['graphs_counted'])))
