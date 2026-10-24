#!/usr/bin/env python3
"""Independent exact receiver audit. Standard library; no producer-code imports."""
from pathlib import Path
from collections import defaultdict, Counter
from hashlib import sha256
from itertools import combinations
import argparse, json
P=Path(__file__).resolve().parent
CHECKS=0

def require(value,message):
    global CHECKS
    CHECKS+=1
    if not value: raise RuntimeError(message)

def read(name): return json.loads((P/'native/inputs'/name).read_text())
def digest(value): return sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()
def qadd(a,b): return (16*((a//16+b//16)%128)) | ((a^b)&15)
def qsum(xs):
    q=0
    for x in xs: q=qadd(q,x)
    return q

def inputs():
    pins=json.loads((P/'SOURCEPINS.json').read_text())
    actual={p.relative_to(P).as_posix() for directory in ['native','sources'] for p in (P/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    require(actual==set(pins),'complete original source closure')
    for name,row in pins.items():
        b=(P/name).read_bytes()
        require(len(b)==row['bytes'] and sha256(b).hexdigest()==row['sha256'],'source pin '+name)
    domain=read('native_packet_decomposition_20261023/native/domains.json')['sources']
    electric=read('metagraph_order8_electrical_20261023/receipt.json')['order8']
    X=electric['higher_layers'][-1]['lifts'][0]
    Y=electric['radical_lifts_mod4']
    labels={s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5)) for s in range(len(X))}
    records={r['source']: sorted(r['exact_all_shared_clocks_target_domain'],key=lambda a:a['target']) for r in domain}
    costs={s:sum(a!=b for a,b in zip(r[0]['path'],r[1]['path'])) for s,r in records.items()}
    cheap={s for s in records if costs[s]==2}
    base=sorted(cheap-{176,1425,2071,1457}); triple=[1457,1841,3015]
    require(len(records)==28 and len(cheap)==20 and len(base)==16,'rederive native receiver inventory')
    require(sum(costs[s] for s in triple)==8,'native triple clock cost8')
    require(qsum(labels[s] for s in triple)==1263,'native triple charge1263')
    for s in base+triple:
        require([r['target'] for r in records[s]]==[1965,2234],'literal binary native target classes')
    primitives=[([s],1,costs[s],labels[s]) for s in base]+[(triple,3,8,qsum(labels[s] for s in triple))]
    F=[X]+Y; mods=[256]+[4]*len(Y)
    for s in base+triple:
        require(all(2*(f[s]-f[1965])%m==2*(f[s]-f[2234])%m for f,m in zip(F,mods)),'fixed source-row paired potential equality')
    return base,triple,primitives,F,mods,labels,costs,records

def cost(diff):
    t=diff>>16
    base=(diff&65535).bit_count()
    return (base+3*t,2*base+8*t)

def cost_chain(states):
    values=sorted({cost(a^b) for a,b in combinations(states,2)})
    require(all(a[0]<=b[0] and a[1]<=b[1] for a,b in zip(values,values[1:])),'simultaneous native edge costs form a coordinate chain')
    return values

def prim_mst(matrix):
    """Independent dense Prim minimax spanning-tree threshold on cost pairs."""
    n=len(matrix); alive=[True]*n; best=[(99,999)]*n; best[0]=(0,0); answer=(0,0)
    for _ in range(n):
        u=min((i for i in range(n) if alive[i]),key=lambda i:best[i]); alive[u]=False
        answer=max(answer,best[u])
        for v in range(n):
            if alive[v] and matrix[u][v]<best[v]: best[v]=matrix[u][v]
    return answer

def bottleneck_floyd(matrix):
    """All-pair ordinary connectivity, independently of the Prim threshold."""
    a=[r[:] for r in matrix];n=len(a)
    for k in range(n):
        for i in range(n):
            ik=a[i][k]
            for j in range(n):
                candidate=max(ik,a[k][j])
                if candidate<a[i][j]: a[i][j]=candidate
    return a

def monotone_dag(states,matrix):
    """All roots; each edge increases the exact set of switched coordinates."""
    n=len(states); ans=[]
    for source,x in enumerate(states):
        diffs=[x^y for y in states]
        order=sorted(range(n),key=lambda i:diffs[i].bit_count())
        best=[(99,999)]*n; best[source]=(0,0)
        for v in order:
            if v==source: continue
            candidates=[]
            for u in order:
                if u!=v and diffs[u]&diffs[v]==diffs[u]:
                    candidates.append(max(best[u],matrix[u][v]))
            best[v]=min(candidates)
        ans.append(best)
    require(all(ans[i][j]==ans[j][i] for i in range(n) for j in range(n)),'strict path reversal symmetry')
    return ans

def interval_recursion(states,matrix):
    """Independent endpoint recursion; rejects proper intermediates outside interval."""
    n=len(states); answer={}; indecomposable=[]
    pairs=sorted(combinations(range(n),2),key=lambda ij:(states[ij[0]]^states[ij[1]]).bit_count())
    for i,j in pairs:
        x,y=states[i],states[j]; diff=x^y
        mids=[k for k,z in enumerate(states) if k not in (i,j) and (x^z)&diff==x^z]
        if not mids:
            answer[i,j]=matrix[i][j];indecomposable.append((i,j))
        else:
            answer[i,j]=min(max(answer[tuple(sorted((i,k)))],answer[tuple(sorted((k,j)))]) for k in mids)
    return answer,indecomposable

def reachable_component(states,start,budget):
    reached={start}; todo=[start]; previous={start:None}
    while todo:
        x=todo.pop(0)
        for y in states:
            slots,clocks=cost(x^y)
            if y not in reached and slots<=budget[0] and clocks<=budget[1]:
                reached.add(y);previous[y]=x;todo.append(y)
    return reached,previous

def certify_ordinary_witness(states,start,finish,budget):
    reached,previous=reachable_component(states,start,budget)
    require(finish in reached,'literal ordinary joint-budget witness path exists')
    path=[finish]
    while path[-1]!=start:path.append(previous[path[-1]])
    path.reverse()
    slotcut,_=reachable_component(states,start,(budget[0]-1,999))
    clockcut,_=reachable_component(states,start,(99,budget[1]-1))
    require(finish not in slotcut and finish not in clockcut,'separate below-minimum cuts certify both exact ordinary coordinates')
    return {'path':path,'events':[cost(x^y) for x,y in zip(path,path[1:])],'below_slot_cut':sorted(slotcut),'below_clock_cut':sorted(clockcut)}

def arbitrary_domain_checks():
    """All 65,535 nonempty binary acceptance sets on four coordinates.
    Cost = separate nonuniform slot/clock weights; feasibility is genuinely joint.
    """
    totals={'domains':0,'budget_domains':0,'total_domain_ordered_pairs':0}
    weights=(1,3,2,5)
    costs=[(d.bit_count(),sum(w for i,w in enumerate(weights) if d>>i&1)) for d in range(16)]
    budgets=[(0,0),(1,3),(2,4),(2,6),(3,7),(4,11)]
    # Exact bitset transitive closure in each source-relative subset poset.
    for subset in range(1,1<<16):
        states=[x for x in range(16) if subset>>x&1];n=len(states); primitive=[]
        for i,x in enumerate(states):
            for j in range(i+1,n):
                y=states[j]
                if not any(z not in (x,y) and (x^z)&(x^y)==x^z for z in states): primitive.append(costs[x^y])
        for budget in budgets:
            theorem=all(s<=budget[0] and c<=budget[1] for s,c in primitive)
            actual=True
            for x in states:
                order=sorted(states,key=lambda y:(x^y).bit_count());reachable={x}
                for y in order:
                    if y==x:continue
                    if any((x^z)&(x^y)==x^z and costs[z^y][0]<=budget[0] and costs[z^y][1]<=budget[1] for z in reachable):reachable.add(y)
                if len(reachable)!=n:actual=False;break
            require(theorem==actual,'arbitrary four-bit domain joint-budget theorem')
            totals['budget_domains']+=1
        totals['domains']+=1; totals['total_domain_ordered_pairs']+=n*(n-1)
    return totals

def main(write=False):
    base,triple,components,F,mods,labels,prices,records=inputs()
    fibres=defaultdict(list); grades=[0]*(1<<17); charges=[0]*(1<<17)
    potentials=[tuple(0 for f in F)]*(1<<17)
    cweights=[tuple(sum(f[s] for s in sources)%m for f,m in zip(F,mods)) for sources,_,_,_ in components]
    totalpot=tuple(sum(f[s] for s in base+triple)%m for f,m in zip(F,mods))
    row_by_fibre={}
    for mask in range(1<<17):
        if mask:
            bit=mask&-mask; j=bit.bit_length()-1; prior=mask^bit
            grades[mask]=grades[prior]+components[j][1]
            charges[mask]=qadd(charges[prior],components[j][3])
            potentials[mask]=tuple((a+b)%m for a,b,m in zip(potentials[prior],cweights[j],mods))
        key=grades[mask],charges[mask]
        rows=tuple(((2*((19-grades[mask])*f[1965]-(tot-v)))%m,(2*(grades[mask]*f[2234]-v))%m) for f,m,tot,v in zip(F,mods,totalpot,potentials[mask]))
        if key in row_by_fibre:require(row_by_fibre[key]==rows,'literal receiver target rows match every accepted grade/Q fibre')
        else:row_by_fibre[key]=rows
        fibres[key].append(mask)
    sizes=Counter(map(len,fibres.values()));pair_count=sum(n*(n-1)//2 for n in map(len,fibres.values()))
    require((len(fibres),max(sizes),pair_count)==(23444,22,474896),'complete fibre census')
    pergrade={k:{'grade':k,'fibres':0,'assignments':0,'pairs':0,'max_fibre':0,'ordinary':(0,0),'strict':(0,0),'primitive_pair_occurrences':0,'ordinary_optimum_histogram':Counter(),'strict_optimum_histogram':Counter(),'gap_pairs':0,'ordinary_lower_witness':None,'strict_lower_witness':None} for k in range(20)}
    trades=set();primitive_hist=Counter();pair_hasher=sha256(); witness=None;largest_primitive=None;global_edges=set()
    fibre_summary=['grade\tcharge\tstates\tpairs\tordinary_slots\tordinary_clocks\tstrict_slots\tstrict_clocks\tindecomposable_pairs']
    for (grade,q),states in sorted(fibres.items()):
        n=len(states);r=pergrade[grade];r['fibres']+=1;r['assignments']+=n;r['pairs']+=n*(n-1)//2;r['max_fibre']=max(r['max_fibre'],n)
        matrix=[[cost(a^b) for b in states] for a in states]
        global_edges.update(matrix[i][j] for i,j in combinations(range(n),2))
        ordinary=bottleneck_floyd(matrix); mst=prim_mst(matrix)
        strict=monotone_dag(states,matrix);recursive,indecomp=interval_recursion(states,matrix)
        require(mst==max((ordinary[i][j] for i,j in combinations(range(n),2)),default=(0,0)),'ordinary all-pair Floyd agrees with independently chosen Prim MST')
        primitive_max=max((matrix[i][j] for i,j in indecomp),default=(0,0))
        require(primitive_max==max((strict[i][j] for i,j in combinations(range(n),2)),default=(0,0)),'strict uniform threshold equals interval-indecomposable maximum')
        fibre_summary.append('\t'.join(map(str,[grade,q,n,n*(n-1)//2,*mst,*primitive_max,len(indecomp)])))
        for i,j in combinations(range(n),2):
            require(strict[i][j]==recursive[i,j],'every exact strict pair optimum agrees DAG versus interval recursion')
            require(ordinary[i][j]<=strict[i][j]<=matrix[i][j],'ordinary versus strict versus direct containment')
            r['ordinary_optimum_histogram'][ordinary[i][j]]+=1;r['strict_optimum_histogram'][strict[i][j]]+=1
            if ordinary[i][j]>r['ordinary']:
                r['ordinary_lower_witness']={'charge':q,'states':states,'endpoints':[states[i],states[j]],'cost':ordinary[i][j],**certify_ordinary_witness(states,states[i],states[j],ordinary[i][j])}
            r['ordinary']=max(r['ordinary'],ordinary[i][j]);r['strict']=max(r['strict'],strict[i][j])
            if ordinary[i][j]<strict[i][j]:
                r['gap_pairs']+=1
                if grade==9 and (witness is None or strict[i][j]>witness['strict']):witness={'grade':grade,'charge':q,'states':states,'start':states[i],'finish':states[j],'ordinary':ordinary[i][j],'strict':strict[i][j],'direct':matrix[i][j]}
            pair_hasher.update(json.dumps([grade,q,states[i],states[j],ordinary[i][j],strict[i][j]],separators=(',',':')).encode()+b'\n')
        r['primitive_pair_occurrences']+=len(indecomp)
        for i,j in indecomp:
            x,y=states[i],states[j]; trade=tuple(sorted((x&~y,y&~x)));trades.add(trade);primitive_hist[matrix[i][j]]+=1
            if r['strict_lower_witness'] is None or matrix[i][j]>r['strict_lower_witness']['cost']:
                r['strict_lower_witness']={'charge':q,'states':states,'endpoints':[x,y],'cost':matrix[i][j],'accepted_interval':[x,y]}
            if largest_primitive is None or matrix[i][j]>largest_primitive['cost']:largest_primitive={'grade':grade,'charge':q,'states':states,'endpoints':[x,y],'cost':matrix[i][j]}
    global_edges=sorted(global_edges)
    require(all(a[0]<=b[0] and a[1]<=b[1] for a,b in zip(global_edges,global_edges[1:])),'all 474896 actual pair costs jointly form a single coordinatewise chain')
    require(all(s%2==0 and c in (2*s,2*s+2) for s,c in global_edges),'native even slots and exact clock offset')
    unique_hist=Counter(cost(x^y) for x,y in trades)
    require(len(trades)==2543 and unique_hist[(16,34)]==196,'independent normalized primitive trade count')
    expected=[(0,0),(0,0),(4,8),(6,12),(8,18),(10,22),(12,26),(14,30),(16,34),(14,30),(14,30),(16,34),(14,30),(12,26),(10,22),(8,18),(6,12),(4,8),(0,0),(0,0)]
    for k,r in pergrade.items():
        require(r['ordinary']==expected[k],'per-grade exact ordinary threshold '+str(k))
        require(r['strict']==((16,34) if k in (9,10) else expected[k]),'per-grade exact strict threshold '+str(k))
        for key in ('ordinary_optimum_histogram','strict_optimum_histogram'):
            r[key]=[{'slots':s,'clocks':c,'pairs':n} for (s,c),n in sorted(r[key].items())]
    sharp=fibres[8,53]
    require(len(sharp)==2 and cost(sharp[0]^sharp[1])==(16,34),'sharp frozen grade8 charge53 exactly two-state fibre')
    def actual_v(mask):return sorted(s for j,(sources,_,_,_) in enumerate(components) if mask>>j&1 for s in sources)
    def enrich(w):
        w['endpoint_V_sources']=[actual_v(w[k]) for k in ('start','finish')]
        w['all_fibre_V_sources']=[actual_v(x) for x in w['states']]
        return w
    witness=enrich(witness)
    witness.update(certify_ordinary_witness(witness['states'],witness['start'],witness['finish'],witness['ordinary']))
    witness['path_V_sources']=[actual_v(x) for x in witness['path']]
    changes=Counter(s for x,y in zip(witness['path'],witness['path'][1:]) for s in actual_v(x^y))
    witness['reused_sources']={s:n for s,n in changes.items() if n>1}
    require(bool(witness['reused_sources']),'literal smaller ordinary route reacquires resources')
    witness['accepted_endpoint_interval']=[z for z in witness['states'] if (witness['start']^z)&(witness['start']^witness['finish'])==witness['start']^z]
    require(set(witness['accepted_endpoint_interval'])=={witness['start'],witness['finish']},'strict gap witness is interval-indecomposable')
    summary='\n'.join(fibre_summary)+'\n'
    if write:(P/'FIBRE_SUMMARY.tsv').write_text(summary)
    require((P/'FIBRE_SUMMARY.tsv').read_text()==summary,'exact23444-row frozen fibre summary')
    result={'fibre_summary_sha256':sha256(summary.encode()).hexdigest(),'status':'PASS','scope':'full exterior frozen; exact literal binary receiver records; no all-size charged supply claim','census':{'assignments':1<<17,'fibres':len(fibres),'max_fibre':max(sizes),'unordered_endpoint_pairs':pair_count,'size_histogram':dict(sorted(sizes.items()))},'labels':{s:labels[s] for s in base+triple},'native_clock_costs':{s:prices[s] for s in base+triple},'base_sources':base,'triple_sources':triple,'joint_edge_cost_chain':global_edges,'per_grade':list(pergrade.values()),'all_pair_optimum_sha256':pair_hasher.hexdigest(),'unique_primitive_trades':len(trades),'unique_primitive_trade_digest':digest(sorted(trades)),'unique_primitive_cost_histogram':[{'slots':s,'clocks':c,'trades':n} for (s,c),n in sorted(unique_hist.items())],'primitive_occurrence_cost_histogram':[{'slots':s,'clocks':c,'pairs':n} for (s,c),n in sorted(primitive_hist.items())],'sharp_grade8_charge53':{'masks':sharp,'V_sources':[actual_v(x) for x in sharp],'joint_threshold':[16,34]},'grade9_strict_gap_witness':witness,'all_binary_four_coordinate_domains':arbitrary_domain_checks(),'checks':CHECKS}
    producer=json.loads((P/'sources/threshold-producer-certificate.json.txt').read_text())
    require([[x,y,*cost(x^y)] for x,y in sorted(trades)]==producer['all_primitive_trades'],'exact independent full primitive list matches frozen producer')
    require(result['all_pair_optimum_sha256']==producer['complete_pair_census']['pair_record_sha256'],'exact independent complete pair digest matches frozen producer')
    for own,other in zip(result['per_grade'],producer['grades']):
        require(list(own['ordinary'])==other['ordinary_threshold'] and list(own['strict'])==other['strict_threshold'],'independent grade thresholds match frozen producer')
        require(own['primitive_pair_occurrences']==other['primitive_pairs'],'independent primitive occurrence census matches frozen producer')
        for ownkey,otherkey in [('ordinary_optimum_histogram','ordinary_pair_optimum_histogram'),('strict_optimum_histogram','strict_pair_optimum_histogram')]:
            require([[v['slots'],v['clocks'],v['pairs']] for v in own[ownkey]]==other[otherkey],'exact independent pair histogram matches frozen producer')
    result['producer_comparison']='PASS: complete pair digest, all2543 primitive trades, all grade pair histograms and thresholds'
    result['checks']=CHECKS
    return json.loads(json.dumps(result))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=main(args.write);target=P/'AUDIT_RESULT.json'
    if args.write:target.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    else:require(result==json.loads(target.read_text()),'frozen complete independent audit result')
    print(json.dumps({'status':'PASS','checks':CHECKS,'census':result['census'],'joint_uniform_ordinary':[16,34],'joint_uniform_strict':[16,34],'pair_optimum_digest':result['all_pair_optimum_sha256'],'all_binary_four_coordinate_domains':result['all_binary_four_coordinate_domains']},sort_keys=True))
