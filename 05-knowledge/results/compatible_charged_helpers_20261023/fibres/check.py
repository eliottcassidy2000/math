"""Independent small-fibre BFS and literal native-resource/current verification."""
import csv,json,sys
from collections import defaultdict,deque,Counter
from itertools import combinations
from hashlib import sha256
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'injective_native_blocks_20261023/search'))
import explore as e
COUNT=0
def need(b,why):
    global COUNT
    COUNT+=1
    if not b:raise RuntimeError(why)
def domains():return sorted(json.loads((RESULTS/'native_packet_decomposition_20261023/native/domains.json').read_text())['sources'],key=lambda r:r['source'])
def enumerate_fibres(rows,k):
    out=defaultdict(list)
    def rec(start,left,m,x,y):
        if not left:out[x%128,y].append(m);return
        for i in range(start,len(rows)-left+1):
            r=rows[i];yy=sum(v<<j for j,v in enumerate(r['code'][1:]))
            rec(i+1,left-1,m|(1<<i),x+r['code'][0],y^yy)
    rec(0,k,0,0,0);return out
def bfs_components(masks,k,budget=6):
    sets={m:{i for i in range(28)if m>>i&1}for m in masks};left=set(masks);out=[]
    while left:
        a=min(left);left.remove(a);todo=deque([a]);cc=[a]
        while todo:
            m=todo.popleft()
            found=[n for n in left if len(sets[m]&sets[n])>=k-budget//2]
            for n in found:left.remove(n);todo.append(n);cc.append(n)
        out.append(sorted(cc))
    return sorted(out,key=lambda c:(len(c),c))
def adj(s):
    A=[[0]*8 for _ in range(8)];c=int(e.NODES[s]['canon_a'])
    for i,(u,v) in enumerate(combinations(range(8),2)):A[u][v]=(c>>(27-i))&1;A[v][u]=1-A[u][v]
    return A
def word(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j)in enumerate(e.CHORDS))
def main():
    rows=domains();src=[r['source']for r in rows];need(len(rows)==28,'actual full source supply')
    for r in rows:need(r['X']==e.X[r['source']] and r['code']==[e.X[r['source']]%128]+[e.Y[j][r['source']]%2 for j in range(1,5)],'original full paired code')
    groups={};components={};report=[]
    for k in (4,5):
        g=enumerate_fibres(rows,k);groups[k]=g;cc={q:bfs_components(mm,k)for q,mm in g.items()};components[k]=cc
        old=json.loads((HERE/f'census_k{k}.json').read_text())
        bad=[dict(charge=list(q),states=sorted(g[q]),components=cc[q])for q in sorted(g)if len(cc[q])>1]
        need(bad==old['disconnected'],'complete independent DFS/BFS hostile census')
        need(sum(map(len,g.values()))==old['states'] and len(g)==2048,'independent full state and fibre counts')
        need(dict(sorted(Counter(map(len,cc.values())).items()))=={int(k):v for k,v in old['component_counts'].items()},'all component distributions')
        need(all(len(bfs_components(z['states'],k,8))==1 for z in bad),'every residual small-margin fibre connected at sharp uniformbudget8')
        report.append(dict(margin=k,states=old['states'],component_counts=old['component_counts']))
    r=json.loads((HERE/'stabilization_receipt.json').read_text());labels={}
    for q,cc in components[5].items():
        for i,c in enumerate(cc):
            for m in c:labels[m]=(q,i)
    actual_pairs=[]
    for q,cc in sorted(components[4].items()):
        for i,c in enumerate(cc):
            for d in cc[i+1:]:
                for a in c:
                    for b in d:
                        hs=[src[j]for j in range(28)if not(a|b)>>j&1 and labels[a|1<<j]==labels[b|1<<j]]
                        need(len(hs)>=18,'occupied-context repair always exists')
                        actual_pairs.append(dict(before=a,after=b,eligible_sources=hs))
    need(actual_pairs==r['all_disconnected_margin4_state_pairs'] and len(actual_pairs)==1300,'all actual cross-component statepairs, not representatives only')
    bank=r['minimum_choice_bank'];need(len(bank)==4 and all(set(bank)&set(p['eligible_sources'])for p in actual_pairs),'four-source choice bank covers every hostile pair')
    recs=r['all_three_choice_bank_hostiles'];need(len(recs)==3276 and [list(c)for c in combinations(src,3)]==[z['bank']for z in recs],'all triple choicebanks present')
    for z in recs:need(not set(z['bank'])&set(actual_pairs[z['pair_index']]['eligible_sources']),'literal obstruction to each triple choicebank')
    h=r['minimal_hostile'];a,b=h['states'];core=a|b;need(h['charge']==[15,3] and len(groups[4][15,3])==2 and(a^b).bit_count()==8,'smallest two-state margin4 hostile')
    need(sum(row['X']for i,row in enumerate(rows)if a>>i&1)==sum(row['X']for i,row in enumerate(rows)if b>>i&1)==527,'hostile even retains integer current')
    outside=[i for i in range(28)if not core>>i&1];need(len(outside)==20,'all potential single helpers')
    need([z['occupied_source']for z in r['all_twenty_single_source_hostiles']]==[src[i]for i in outside],'singleton helper hostile coverage')
    two=[];exact_two=[]
    for j in outside:
        c=rows[j]['code'];q=((15+c[0])%128,3^sum(v<<i for i,v in enumerate(c[1:])))
        ss=[m for m in groups[5][q]if not m&~(core|1<<j)]
        need(len(ss)==2 and len(bfs_components(ss,5))==2,'one extra actual source never repairs the8source core')
        for t in outside:
            if t==j:continue
            ss=[m for m in groups[5][q]if not m&~(core|1<<j|1<<t)]
            if any(a|1<<j in c and b|1<<j in c for c in bfs_components(ss,5)):two.append((src[j],src[t]))
            ss=[m for m in ss if sum(row['X']for i,row in enumerate(rows)if m>>i&1)==527+rows[j]['X']]
            if any(a|1<<j in c and b|1<<j in c for c in bfs_components(ss,5)):exact_two.append((src[j],src[t]))
    need(two==[(z['occupied_source'],z['unoccupied_source'])for z in r['successful_two_source_contexts']] and len(two)==18,'complete380 oriented helperpair tests')
    need(exact_two==[(z['occupied_source'],z['unoccupied_source'])for z in r['successful_exact_integer_contexts']] and len(exact_two)==12,'complete exact-current helperpair controls')
    chosen=r['selected_ten_source_exact_context'];keep=core|1<<src.index(1480)|1<<src.index(3000);path=chosen['route'];U,V=1965,2234
    need(e.X[U]==e.X[V]==146,'original chosen per-slot phase')
    need(len(chosen['states'])==3 and len(path)==3 and path[0]==a|1<<src.index(1480) and path[-1]==b|1<<src.index(1480),'actual stabilized endpoints')
    j=src.index(1480);c=rows[j]['code'];q=((15+c[0])%128,3^sum(v<<i for i,v in enumerate(c[1:])))
    need(sorted(chosen['states'])==sorted(m for m in groups[5][q]if not m&~keep),'entire ten-source modular fibre equals these three exact-integer states')
    need(len(bfs_components(chosen['states'],5,4))==3,'selected stabilized fibre has sharpbudget6')
    reservoir=set();literal=[];lookup={}
    for i,row in enumerate(rows):
        if not keep>>i&1:continue
        s=row['source'];A=adj(s);op=[[A[v][u]for v in range(8)]for u in range(8)]
        for item in row['exact_all_shared_clocks_target_domain']:
            p=item['path'];m=word(A,p);cm=word(op,p[::-1]);target=item['target'];closed={m,cm,e.FULL^m,e.FULL^cm}
            need(all(A[x][y]for x,y in zip(p,p[1:])),'actual path adjacency')
            need(all(p[j]==v for j,v in row['fixed_clocks']),'all actual owner-clock pins')
            need(m==item['mask'] and e.MERGED[e.OWNER[m]]==s and e.MERGED[e.OWNER[e.FULL^m]]==target,'literal normalized source and inverse target')
            need(e.MERGED[e.OWNER[cm]]==s and e.MERGED[e.OWNER[e.FULL^cm]]==target,'actual converse normalized inverse ports')
            need(len(closed)==4 and not reservoir&closed,'complete F/converse resource reservoir disjoint')
            reservoir|=closed;lookup[s,target]=p;literal.append(dict(source=s,target=target,path=p,closed_words=sorted(closed)))
    need(len(literal)==20 and len(reservoir)==80,'two actual resources per10sources')
    for m in path:
        need(m in chosen['states'] and m.bit_count()==5 and not m&~keep,'all native arrival states inside exact10source supply')
        need(sum(row['X']for i,row in enumerate(rows)if m>>i&1)==599,'integer source sum constant')
        code=[sum(row['code'][0]for i,row in enumerate(rows)if m>>i&1)%128]+[sum(row['code'][j]for i,row in enumerate(rows)if m>>i&1)%2 for j in range(1,5)]
        if m==path[0]:basecode=code
        need(code==basecode,'complete aggregate dyadic code at every arrival')
    for (x,y),ev in zip(zip(path,path[1:]),chosen['events']):
        changed=[src[i]for i in range(28)if(x^y)>>i&1]
        need(len(changed)==6 and changed==ev['changed_sources'],'literal six-slot event')
        work=sum(sum(a!=b for a,b in zip(lookup[s,U],lookup[s,V]))for s in changed)
        need(work==ev['ordinary_owner_clock_work'],'actual owner-clock work')
    inputs=[RESULTS/'native_packet_decomposition_20261023/native/domains.json',e.NATIVE/'owner_map.u16.gz',e.NATIVE/'merged_nodes.tsv',RESULTS/'metagraph_order8_electrical_20261023/receipt.json']
    out=dict(status='FINITE-EXACT independent small-fibre and native stabilization audit',checks=COUNT,small_census=report,
        sources=src,stabilized_native_resources=literal,full_closed_resource_reservoir=sorted(reservoir),
        source_hashes={str(p.relative_to(RESULTS)).replace('\\','/'):sha256(p.read_bytes()if p.suffix=='.gz'else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()for p in inputs})
    file=HERE/'check_receipt.json'
    if '--write' in sys.argv:file.write_text(json.dumps(out,indent=2))
    else:need(json.loads(file.read_text())==out,'saved audit receipt')
    print('PASS',COUNT,'checks; complete margins4/5;1300stabilizedpairs;3state10source native route')
if __name__=='__main__':main()
