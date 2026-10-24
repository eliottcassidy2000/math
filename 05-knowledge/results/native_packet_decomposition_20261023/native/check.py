"""Independent native-domain, all-combinations fibre and all-pairs graph audit.

No discovery cache is used. Exact-state completeness is checked directly
over every8-subset of28, independently of the producer's14+14 MITM.
"""
import csv,gzip,json,sys
from collections import Counter,defaultdict,deque
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'injective_native_blocks_20261023/search'))
import explore as e
COUNT=0
def need(b,why):
    global COUNT
    COUNT+=1
    if not b:raise RuntimeError(why)
def matrix(canon):
    A=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=(canon>>(27-k))&1;A[j][i]=1-A[i][j]
    return A
def mask(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(e.CHORDS))
class DSU:
    def __init__(self,n):self.p=list(range(n));self.sz=[1]*n
    def root(self,i):
        while self.p[i]!=i:self.p[i]=self.p[self.p[i]];i=self.p[i]
        return i
    def add(self,a,b):
        a,b=self.root(a),self.root(b)
        if a==b:return
        if self.sz[a]<self.sz[b]:a,b=b,a
        self.p[b]=a;self.sz[a]+=self.sz[b]
def components(states,budget):
    ds=DSU(len(states))
    for i,a in enumerate(states):
        for j in range(i):
            if (a^states[j]).bit_count()<=budget:ds.add(i,j)
    return sorted(Counter(ds.root(i) for i in range(len(states))).values())
def shortest(states,a,b,budget=6):
    prev={a:None};todo=deque([a])
    while todo and b not in prev:
        x=todo.popleft()
        for y in states:
            if y not in prev and (x^y).bit_count()<=budget:prev[y]=x;todo.append(y)
    if b not in prev:return None
    out=[]
    while b is not None:out.append(b);b=prev[b]
    return out[::-1]
def main():
    D=json.loads((HERE/'domains.json').read_text());rows=D['sources'];sources=sorted(z['source'] for z in rows);by={z['source']:z for z in rows};U,V=D['targets']
    old=json.loads((RESULTS/'injective_native_blocks_20261023/search/larger16_selected.json').read_text())['circuits'][0]
    edgefile=e.NATIVE/'merged_edges.tsv.gz';edges=[tuple(int(r[k]) for k in ('source','target','blue','black')) for r in csv.DictReader(gzip.decompress(edgefile.read_bytes()).decode().splitlines(),delimiter='\t')]
    weights={(i,j):(b,k) for i,j,b,k in edges}
    need(len(rows)==28 and len(set(sources))==28,'complete inherited28source menu')
    allports=set();resources=[];source_details=[];hostiles=[];raw_pin_total=0;off_target=[]
    for z in rows:
        s=z['source'];A=matrix(int(e.NODES[s]['canon_a']));op=[[A[j][i] for j in range(8)] for i in range(8)]
        actual=[]
        # Direct permutation enumeration, independent of the source-path DFS.
        aut=0
        for p in permutations(range(8)):
            if all(A[i][j]==A[p[i]][p[j]] for i,j in combinations(range(8),2)):aut+=1
            if not all(A[a][b] for a,b in zip(p,p[1:])):continue
            m=mask(A,p);need(e.MERGED[e.OWNER[m]]==s,'actual full source owner')
            t=e.MERGED[e.OWNER[e.FULL^m]];actual.append(dict(path=list(p),target=t,mask=m,phase=(e.X[s]-e.X[t])%8))
        need(aut==1 and len(actual)==z['H'],'complete actual source paths and trivial kernel')
        need(z['X']==e.X[s] and z['code']==[e.X[s]%128]+[e.Y[j][s]%2 for j in range(1,5)],'full original marked charge coordinates')
        p0,p1=z['old_before'],z['old_after'];first,last,k,a,b=z['pins']
        fixed=[[i,x] for i,(x,y) in enumerate(zip(p0,p1)) if x==y]
        need(fixed==z['fixed_clocks'],'all inherited shared clocks retained')
        unpinned=[r for r in actual if r['target'] in (U,V)]
        pinned=[r for r in unpinned if r['path'][0]==first and r['path'][-1]==last and r['path'][k:k+2]==[a,b]]
        fully=[r for r in unpinned if all(r['path'][i]==v for i,v in fixed)]
        raw_full=[r for r in actual if all(r['path'][i]==v for i,v in fixed)]
        raw_pin_total+=len(raw_full)
        off_target.extend(dict(source=s,**r) for r in raw_full if r['target'] not in (U,V))
        need(raw_full==z['raw_all_shared_clocks_domain'],'raw pin domain before target-margin exclusion')
        need(unpinned==z['unpinned_target_domain'],'complete target domain before pins')
        need(pinned==z['exact_arc_target_domain'] and fully==z['exact_all_shared_clocks_target_domain'],'exact complete pinned domains')
        need(len(pinned)==2 and pinned==fully and {r['target'] for r in pinned}=={U,V},'every actual pinned domain exactly binary withboth targets')
        for r in unpinned:
            if r not in pinned:
                p=r['path'];bad=[]
                if p[0]!=first:bad.append('first owner')
                if p[-1]!=last:bad.append('last owner')
                if p[k:k+2]!=[a,b]:bad.append('indexed arc')
                need(bool(bad),'unused extra path violates actual inherited pins')
                hostiles.append(dict(source=s,path=p,target=r['target'],failed_pins=bad))
        caps=[]
        for t in (U,V):
            b0,k0=weights[min(s,t),max(s,t)];count=sum(r['target']==t for r in actual)
            need(b0==0 and k0==2*count,'all full source/target capacities are literal paired black')
            caps.append(dict(target=t,blue=b0,black=k0,unpinned_ordinary_paths=count,pinned_ordinary_paths=1))
        for r in pinned:
            p=r['path'];m=r['mask'];cm=mask(op,list(reversed(p)))
            need(all(op[a][b] for a,b in zip(p[::-1],p[-2::-1])),'actual reversed converse path')
            need(e.MERGED[e.OWNER[cm]]==s and e.MERGED[e.OWNER[e.FULL^cm]]==r['target'],'actual converse owners')
            K=[[0]*8 for _ in range(8)];pp={frozenset(t) for t in zip(p,p[1:])}
            for i,j in combinations(range(8),2):K[i][j]=A[i][j] if frozenset((i,j)) in pp else A[j][i];K[j][i]=1-K[i][j]
            need(mask(K,p)==e.FULL^m,'literal inverse target normalized port')
            words={m,e.FULL^m,cm,e.FULL^cm}
            need(len(words)==4 and allports.isdisjoint(words),'normalized F/converse resources genuinely disjoint')
            allports.update(words);resources.append(dict(source=s,target=r['target'],key=min(words),closed_words=sorted(words)))
        source_details.append(dict(source=s,H=len(actual),aut=aut,capacities=caps))
    need(len(hostiles)==19 and len(resources)==56 and len(allports)==224,'exact native spare-path and resource counts')
    need(raw_pin_total==57 and len(off_target)==1 and off_target[0]['source']==1923 and off_target[0]['target']==1559,'unique raw pin alternative excluded by target1559 zero occupancy')
    for slot in old['slots']:
        need({tuple(slot[k]['path']) for k in ('before','after')}=={tuple(r['path']) for r in by[slot['source']]['exact_arc_target_domain']},'old16literal paths and commonclock records preserved, not only source labels')
    need(e.X[U]==e.X[V]==146 and all(e.Y[j][U]%2==e.Y[j][V]%2 for j in range(1,5)),'targets equal in complete paired dyadic code')
    need(all((b+k)%2==0 for i,j,b,k in edges if U in (i,j)),'4eU is authentic universal per-slot phase hostile')
    initial=sum(1<<sources.index(s) for s in old['negative']);terminal=sum(1<<sources.index(s) for s in old['positive'])
    X=[e.X[s] for s in sources];Y=[sum((e.Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in sources]
    code_states=[];integer_states=[];subset_count=0
    for ids in combinations(range(28),8):
        subset_count+=1;x=sum(X[i] for i in ids)
        if x%128!=120:continue
        y=0
        for i in ids:y^=Y[i]
        if y:continue
        m=sum(1<<i for i in ids);code_states.append(m)
        if x==1016:integer_states.append(m)
    code_states.sort();integer_states.sort()
    need(subset_count==3108105 and len(code_states)==1569 and len(integer_states)==626,'independent complete8subset enumeration')
    graph_reports=[]
    for name,states in [('fibre_receipt.json',code_states),('fibre_exact_receipt.json',integer_states)]:
        f=json.loads((HERE/name).read_text());need(f['states']==states and f['initial']==initial and f['terminal']==terminal,'independent enumeration equals MITMreceipt')
        n=len(states);ds={b:DSU(n) for b in (2,4,6)};adj=[[] for _ in states];index={m:i for i,m in enumerate(states)}
        for i,m in enumerate(states):
            for j in range(i):
                cost=(m^states[j]).bit_count()
                for b in ds:
                    if cost<=b:ds[b].add(i,j)
                if cost<=6:adj[i].append(j);adj[j].append(i)
        thresholds=[]
        for b,z in ds.items():
            sizes=sorted(Counter(z.root(i) for i in range(n)).values());rec=next(r for r in f['thresholds'] if r['budget']==b)
            need(sizes==rec['component_sizes'],'all-pairs component certificate agrees with MST thresholds')
            thresholds.append((b,len(sizes)))
        need(ds[4].root(index[initial])!=ds[4].root(index[terminal]),'literal budget4component separation of old endpoints')
        need(len({ds[6].root(i) for i in range(n)})==1,'all feasible states connected atbudget6')
        td=DSU(n);need(len(f['complete_graph_MST'])==n-1,'saved spanning tree length')
        for i,j,c in f['complete_graph_MST']:
            need(c==(states[i]^states[j]).bit_count() and c<=6,'every saved native tree edge legal budget6');td.add(i,j)
        need(len({td.root(i) for i in range(n)})==1,'explicit spanning tree reaches every feasible state')
        shortest_path=shortest(states,initial,terminal)
        need(len(shortest_path)==len(f['shortest_endpoint_route'])==5,'exact shortest4event route atsharpbudget6')
        route=f['shortest_endpoint_route'];need(all(m in index for m in route),'all route arrival states authenticated')
        costs=[sum(a!=b for a,b in zip(by[s]['old_before'],by[s]['old_after'])) for s in sources]
        works=[sum(c for i,c in enumerate(costs) if (a^b)>>i&1) for a,b in zip(route,route[1:])]
        need(works==f['route_owner_clock_occurrences'],'actual route owner-clock work retained')
        outside=((1<<28)-1)^(initial|terminal)
        for m,entry in zip(route,f['route_owned_debt_ledger']):
            oldsum=sum(X[i] for i in range(28) if(m&~outside)>>i&1);extrasum=sum(X[i] for i in range(28) if(m&outside)>>i&1)
            need(entry['old16_V_sum']==oldsum and entry['outside_V_sum']==extrasum and entry['U_current_change_total']==2*(oldsum+extrasum-1016),'actual target-owned integer debt at every arrival')
            dc=(m&~outside).bit_count()-8;nc=(m&outside).bit_count()
            need(entry['old16_V_count']==dc+8 and entry['outside_V_count']==nc and entry['target_U_degree_delta_old16']==-2*dc and entry['target_U_degree_delta_outside']==-2*nc,'partial target-degree diagonal is retained')
            need(entry['U_current_change_old16']==2*(oldsum-1016-146*dc) and entry['U_current_change_outside']==2*(extrasum-146*nc),'split physical currents include source offdiagonals AND target diagonal')
            if name.startswith('fibre_exact'):need(entry['U_current_change_total']==0,'integer-current route pays at every arrival')
        graph_reports.append(dict(name=name,states=n,thresholds=thresholds,shortest_events=4))
    c=json.loads((HERE/'compensator_receipt.json').read_text());outside=((1<<28)-1)^(initial|terminal)
    one=[m for m in integer_states if (m&outside).bit_count()<=1]
    need(one==c['one_active_full28_states'] and len(one)==40,'complete one-activecompensator family')
    need([ctrl['source'] for ctrl in c['single_extra_source_controls']]==[sources[i] for i in range(28) if outside>>i&1],'singleton-helper controls cover every one of the12outside sources exactly once')
    for ctrl in c['single_extra_source_controls']:
        keep=1<<sources.index(ctrl['source']);ss=[m for m in integer_states if not m&outside&~keep]
        need(ss==ctrl['states'] and components(ss,6)==ctrl['components_at6'],'complete singleton-helper hostile')
        need(shortest(ss,initial,terminal) is None and not ctrl['endpoints_connected'],'no single supplied helper repairsbudget6')
    successful=[]
    for pair in combinations([sources[i] for i in range(28) if outside>>i&1],2):
        keep=sum(1<<sources.index(s) for s in pair);ss=[m for m in integer_states if not m&outside&~keep]
        if shortest(ss,initial,terminal) is not None:successful.append(list(pair))
    need(successful==[r['sources'] for r in c['eligible_two_source_helper_pairs']] and len(successful)==11,'complete66pair helper boundary')
    selected=c['selected_two_source_helper'];need(selected['sources']==[1532,1923] and len(selected['states'])==16,'selected18source complete fibre')
    helpermask=sum(1<<sources.index(s) for s in [1532,1923])
    need(selected['states']==[m for m in integer_states if not m&outside&~helpermask],'selected18source states are exactly complete restricted integer fibre')
    need(components(selected['states'],6)==[16] and len(components(selected['states'],4))>1,'selected18source wholefibre sharp6')
    one18=[m for m in selected['states'] if (m&outside).bit_count()<=1]
    need(one18==c['selected_one_active_states'] and len(one18)==11 and components(one18,6)==c['selected_one_active_components_at6']==[11],'complete11state oneactive18source family connected')
    path=c['selected_route'];need(path[0]==initial and path[-1]==terminal and len(path)==6,'saved five-event oneactive route')
    for m,entry in zip(path,c['selected_route_debt_ledger']):
        need(m in selected['states'] and (m&outside).bit_count()<=1,'all actual18source arrivals allowed withoneactive helper')
        oldsum=sum(X[i] for i in range(28) if(m&~outside)>>i&1);extrasum=sum(X[i] for i in range(28) if(m&outside)>>i&1)
        need(entry['old16_V_sum']==oldsum and entry['outside_V_sum']==extrasum and oldsum+extrasum==1016,'explicit two-source debts cancel at sameactualtarget')
        dc=(m&~outside).bit_count()-8;nc=(m&outside).bit_count()
        need(entry['old16_V_count']==dc+8 and entry['outside_V_count']==nc and entry['source_charge_X_change_old16']==2*(oldsum-1016) and entry['source_charge_X_change_outside']==2*extrasum,'source-charge ledger distinct from current ledger')
        need(entry['target_U_degree_delta_old16']==-2*dc and entry['target_U_degree_delta_outside']==-2*nc and entry['U_debt_old16']==2*(oldsum-1016-146*dc) and entry['U_debt_outside']==2*(extrasum-146*nc),'exact physical target rows include changed diagonal')
        yp=yq=yb=0
        for i,y in enumerate(Y):
            if (m&~outside)>>i&1:yp^=y
            if (m&outside)>>i&1:yq^=y
            if initial>>i&1:yb^=y
        target_y=sum((e.Y[j][U]%2)<<(j-1) for j in range(1,5))
        need(entry['old16_binary_source_charge']==(yp^yb) and entry['outside_binary_source_charge']==yq,'binary source charges retained separately')
        need(entry['old16_binary_residual']==((yp^yb)^(target_y if dc%2 else 0))==entry['outside_binary_residual']==(yq^(target_y if nc%2 else 0)),'all four physical binary target-row residuals include degree terms')
    for (a,b),ev in zip(zip(path,path[1:]),c['selected_route_events']):
        changed=[sources[i] for i in range(28) if (a^b)>>i&1]
        need(ev['before']==a and ev['after']==b and changed==ev['changed_sources'] and len(changed)==6,'actual marked event support')
        need(ev['ordinary_owner_clock_occurrences']==sum(sum(u!=v for u,v in zip(by[s]['old_before'],by[s]['old_after'])) for s in changed),'actual marked event clock cost')
    import four_move_invariant
    four_move_invariant.main()
    need(json.loads((HERE/'four_move_receipt.json').read_text())['endpoint_values']==[-4,-8],'all four-trade structural invariant separates these actual endpoints')
    inputs=[e.NATIVE/'merged_nodes.tsv',e.NATIVE/'owner_map.u16.gz',edgefile,RESULTS/'metagraph_order8_electrical_20261023/receipt.json',RESULTS/'session45_injective_multisource_moves_20261023/source-circuits/all_menus.json',RESULTS/'injective_native_blocks_20261023/search/larger16_selected.json']
    result=dict(status='FINITE-EXACT native complete-pin decomposition; all-dyadic code theorem inherited explicitly',checks=COUNT,
                source_hashes={p.relative_to(RESULTS).as_posix():sha256(p.read_bytes() if p.suffix=='.gz' else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in inputs},
                sources=source_details,normalized_resources=resources,full_fixed_background_endpoint_reservoir=sorted(allports),
                excluded_spare_paths=hostiles,raw_common_clock_paths=raw_pin_total,zero_target_margin_exclusions=off_target,complete8subset_universe=subset_count,graph_results=graph_reports,
                helper_minimum=2,selected_helpers=[1532,1923],selected18source_states=16,selected_route_events=5,selected_route_clock_work=sum(ev['ordinary_owner_clock_occurrences'] for ev in c['selected_route_events']),
                scope='Complete specified fixed-margin/current fibres of independent28or18source slots; fixed actual pins and normalizedresources; aggregate originaldyadic currents, not universal per-slot phases; nofullatlasrewiring')
    result=json.loads(json.dumps(result));out=HERE/'check_receipt.json'
    if '--write' in sys.argv:out.write_text(json.dumps(result,indent=2))
    else:need(json.loads(out.read_text())==result,'saved complete native audit')
    print('PASS',COUNT,'checks;11828actualpaths;3108105eight-subsets;1569/626states;globalbudget6;2helpers')
if __name__=='__main__':main()
