"""Literal persistent-vertex reset, complete nearby replacement tests, neutral grid."""
from collections import Counter, deque
import itertools
import json
import check as a
import two_birth as b
from frontier import strong

def flip(rows,i,j):
    out=list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)

def main():
    report={}
    traces=[]
    for k in range(1,6):
        source=a.ready_source((2,3,1),k)
        keep=[v for v in range(len(source)) if v != k+4]
        first=a.induced(source,keep)
        out=7*2**(k+1)
        second=a.attach(first,out)
        third=a.attach(second,out)
        graphs=(first,second,third)
        counts=[a.dp(g)[0][-1] for g in (source,*graphs)]
        expect=[34*2**k+5,7*2**k+5,15*2**k+9,31*2**k+17]
        a.need(counts==expect,'exact reset count family')
        for j,R in enumerate(graphs):
            a.need(a.induced(R,range(len(first)))==first,'persistent old arc ledger')
            apex=k+3
            spine=[k+4]+list(range(k+1))+list(range(k+5,k+5+j))+[k+1,k+2]
            a.need(len(spine)==len(R)-1 and apex not in spine,'literal spine labels')
            a.need(all(R[spine[i]] >> spine[t] & 1 for i in range(len(spine)) for t in range(i+1,len(spine))),'spine transitive')
            word=''.join(str((R[apex] >> v)&1) for v in spine)
            a.need(word=='1'+'0'*k+'1'+'0'*(2+j),'exact inner word')
        record={'k':k,'source_rows':source,'kept_source_labels':keep,'birth_outmask':out,
                'intermediate_rows':graphs,'count_trace':counts}
        if k==2:
            alternate=(254,508,504,240,480,448,384,256,9)
            mapping=(1,2,3,6,7,8,0,4,5)
            a.need(third==a.induced(alternate,mapping),'exact history alternate source isomorphism')
            continuation1=a.attach(third,16)
            continuation2=a.attach(continuation1,512)
            a.need(a.dp(continuation1)[0][-1]==215 and a.dp(continuation2)[0][-1]==425,'actual rebase continuation')
            record.update({'current_to_history_labels':mapping,'continuation_rows':[continuation1,continuation2],
                           'full_trace':counts+[215,425]})
        traces.append(record)
    report['literal_resets']=traces
    replacements=[]
    for k in (1,2):
        T=a.ready_source((2,3,1),k)
        W=a.induced(T,[v for v in range(len(T)) if v != k+4])
        target=34*2**k+5
        hits=[]
        first_pruned=0
        for out1 in range(1 << len(W)):
            U=a.attach(W,out1)
            if a.dp(U)[0][-1]>target:
                first_pruned+=1
                continue
            for out2,h in enumerate(b.all_birth_counts(U)):
                if h != target:
                    continue
                V=a.attach(U,out2)
                a.need(a.dp(V)[0][-1]==target,'replacement hit independent DP')
                modules=[]
                for module in range(k):
                    co=b.module_runs(V,module)
                    good=co[2]==34*2**(k-1) and not any(co[3:])
                    modules.append({'module':module,'coefficients':co,'required_uniform_response':good})
                hits.append({'outmasks':[out1,out2],'rows':V,'strong':strong(V),'module_responses':modules})
        a.need(len(hits)==(34 if k==1 else 22),'replacement full hit census')
        good=[hit for hit in hits if any(m['required_uniform_response'] for m in hit['module_responses'])]
        a.need(not any(x['strong'] for x in good),'no strong uniform reset in tested templates')
        if k==1:
            a.need(len(good)==4,'only padded original uniform returns')
            for hit in good:
                a.need(any(a.iso(a.induced(hit['rows'],[v for v in range(len(hit['rows'])) if v!=r]),T)
                           for r in range(len(hit['rows']))),'uniform hit contains exact old source class')
        replacements.append({'k':k,'parent_rows':W,'target':target,'ordered_universe':2**(2*len(W)+1),
                             'first_cut_pruned_above_target':first_pruned,'hits':hits,
                             'uniform_hits_count':len(good),'strong_uniform_hits_count':0})
    report['replacement_probes']=replacements
    # Broader fixed-complexity bank: delete any nonprefix old vertex and
    # replace it by two vertices, preserving the entire H(k) response.
    census=[]
    uniform_strong=[]
    representative=None
    for tail in ((2,3,1),(2,2,1)):
        T=a.ready_source(tail,1)
        for drop in range(1,len(T)):
            W=a.induced(T,[v for v in range(len(T)) if v!=drop])
            hits=uniform=strong_uniform=0
            for out1 in range(1 << len(W)):
                U=a.attach(W,out1)
                if a.dp(U)[0][-1]>73:
                    continue
                for out2,value in enumerate(b.all_birth_counts(U)):
                    if value!=73:
                        continue
                    hits+=1
                    V=a.attach(U,out2)
                    co=b.module_runs(V)
                    if co[2]!=34 or any(co[3:]):
                        continue
                    uniform+=1
                    if not strong(V):
                        continue
                    strong_uniform+=1
                    a.need(a.dp(V)[0][-1]==73,'broader uniform replacement direct count')
                    if representative is None:
                        representative=V
                    mapping=a.iso(V,representative)
                    a.need(mapping is not None and mapping[0]==0,'one pointed isomorphism type in broader uniform bank')
                    reverse=tuple(sum((i!=j and not(V[i] >> j & 1)) << j for j in range(len(V))) for i in range(len(V)))
                    a.need(all(r['gap'] not in (1,2) for Q in (V,reverse) for r in a.root_profiles(Q)),'new uniform source is not forward or converse ready at base')
                    uniform_strong.append({'source_tail':tail,'deleted_vertex':drop,'outmasks':[out1,out2],
                                           'rows':V,'module_coefficients':co,'pointed_iso_to_representative':mapping})
            census.append({'source_tail':tail,'deleted_vertex':drop,'universe':8192,'count_hits':hits,
                           'uniform_hits':uniform,'strong_uniform_hits':strong_uniform})
    a.need(len(uniform_strong)==8,'broader complete strong uniform bank')
    source=b.inflate(representative,0,2)
    a.need(a.dp(source)[0][-1]==141 and strong(source),'new representative at k2')
    layer=[{'rows':source,'h':141,'path':[]}]
    layers=[]
    reached=set()
    while layer:
        children=[]
        examined=0
        for parent in layer:
            counts=b.all_birth_counts(parent['rows'])
            examined+=len(counts)
            for out,value in enumerate(counts):
                if value>425:
                    continue
                child=a.attach(parent['rows'],out)
                if not strong(child):
                    continue
                a.need(a.dp(child)[0][-1]==value,'broader source frontier direct count')
                reached.add(value)
                children.append({'rows':child,'h':value,'path':parent['path']+[out]})
        layers.append({'children_count':len(children),'cuts_examined':examined,
                       'histogram':dict(sorted(Counter(x['h'] for x in children).items())),
                       'children':children})
        layer=children
    a.need([x['children_count'] for x in layers]==[27,11,0] and 425 not in reached,'new uniform bank all-order425 obstruction')
    report['broader_uniform_replacement_bank']={'census':census,'strong_uniform_hits':uniform_strong,
        'representative_rows':representative,'representative_k2_rows':source,'frontier_bound':425,
        'frontier_layers':layers,'reached_counts':sorted(reached)}
    A=a.ready_source((2,3,1),2)
    B=a.ready_source((2,2,1),2)
    queue=deque([A])
    seen={A:0}
    states=[A]
    edges=[]
    cuts=0
    while queue:
        R=queue.popleft()
        for i,j in itertools.combinations(range(len(R)),2):
            cuts+=1
            U=flip(R,i,j)
            if a.dp(U)[0][-1] != 141:
                continue
            if U not in seen:
                seen[U]=len(states)
                states.append(U)
                queue.append(U)
            if seen[R]<seen[U]:
                edges.append({'from':seen[R],'to':seen[U],'edge':[i,j]})
    a.need((len(states),len(edges),cuts)==(8,10,224),'neutral component closure')
    classified=[]
    for i,R in enumerate(states):
        pa,pb=a.iso(R,A),a.iso(R,B)
        a.need(pa is not None or pb is not None,'all neutral states in two obstructed classes')
        classified.append({'index':i,'rows':R,'class':'A' if pa is not None else 'B','isomorphism':pa if pa is not None else pb})
    # Coordinates on K2 x P4 in this deterministic enumeration.
    grid={0:(0,0),1:(1,0),2:(0,1),3:(1,1),4:(0,2),5:(1,2),6:(0,3),7:(1,3)}
    a.need(all(sum(abs(x-y) for x,y in zip(grid[e['from']],grid[e['to']]))==1 for e in edges),'literal grid adjacency')
    report['neutral_edge_component']={'states':classified,'edges':edges,'all_edge_flips_tested':cuts,'grid_coordinates':grid}
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__=='__main__':
    main()
