"""Exact native relative-block helper and constrained54-bank repair."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import csv,gzip,json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
SUPPLY=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
FULL=(1<<21)-1;U,V,S=1965,2234,2125;CHECKS=0
PAIRS=list(combinations(range(8),2));CHORDS=[(i,j)for i in range(8)for j in range(i+2,8)]
REV=[CHORDS.index((7-j,7-i))for i,j in CHORDS]
def need(v,m):
    global CHECKS
    CHECKS+=1
    if not v:raise RuntimeError(m)
def matrix(c):
    A=[[0]*8 for _ in range(8)]
    for k,(i,j)in enumerate(PAIRS):A[i][j]=(c>>(27-k))&1;A[j][i]=1-A[i][j]
    return A
def mask(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j)in enumerate(CHORDS))
def reflect(m):return sum(((m>>k)&1)<<REV[k]for k in range(21))
def tile(m):return{m,reflect(m),FULL^m,FULL^reflect(m)}
def guard(p):
    a={v:i for i,v in enumerate(p)}
    return a[6]==a[5]+1 and a[2]==a[6]+1 and a[4]==a[3]+1
def encode(a):return(a[0]<<4)|sum(v<<j for j,v in enumerate(a[1:]))
def add(a,b):return((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def sub(a,b):return((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def charge(labels,ids):
    q=0
    for i in ids:q=add(q,labels[i])
    return q
def main():
    supply=json.loads(SUPPLY.read_text());nodes=list(csv.DictReader((RESULTS/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
    merged={int(z[k]):int(z['merged_id'])for z in nodes for k in('class_a','class_b')}
    owner=array('H');owner.frombytes(gzip.decompress((RESULTS/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder!='little':owner.byteswap()
    electrical=json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    X=electrical['higher_layers'][-1]['lifts'][0];Y=electrical['radical_lifts_mod4']
    label=lambda s:[X[s]%128]+[Y[j][s]%2 for j in range(1,5)]
    A=matrix(int(nodes[S]['canon_a']));op=[list(r)for r in zip(*A)];actual=[];auts=[]
    for p in permutations(range(8)):
        if all(A[i][j]==A[p[i]][p[j]]for i,j in PAIRS):auts.append(p)
        if all(A[a][b]for a,b in zip(p,p[1:])):
            m=mask(A,p);need(merged[owner[m]]==S,'entire actual canonical path source')
            actual.append((p,m,merged[owner[FULL^m]]))
    need(len(actual)==401,'complete401 actual Hamiltonian paths')
    multiplicities=Counter(m for p,m,t in actual)
    need(set(multiplicities.values())=={len(auts)},'all source-word kernel fibres have automorphism size')
    domain=sorted([r for r in actual if r[2]in(U,V)],key=lambda r:r[2])
    need(len(domain)==2 and [r[2]for r in domain]==[U,V],'full source U/V domain exactly two')
    P,Q=[r[0]for r in domain]
    need(list(P)==[5,6,2,0,7,3,4,1] and list(Q)==[7,1,0,3,4,5,6,2],'literal two actual records')
    common=[[i,a]for i,(a,b)in enumerate(zip(P,Q))if a==b]
    need(not common,'no individual indexed owner clock retained')
    perm=[Q.index(a)for a in P];cycle=[];i=0
    while i not in cycle:cycle.append(i);i=perm[i]
    need(i==0 and len(cycle)==8,'relative owner-clock permutation is one eight-cycle')
    common_sets=[m for m in range(1,255)if {P[i]for i in range(8)if m>>i&1}=={Q[i]for i in range(8)if m>>i&1}]
    need(not common_sets,'no nonempty proper unordered clock-set pin exists')
    shared_arcs=set(zip(P,P[1:]))&set(zip(Q,Q[1:]));need(shared_arcs=={(5,6),(6,2),(3,4)},'exact common consecutive owner arcs')
    raw=[r for r in actual if guard(r[0])]
    need([r for r in domain if guard(r[0])]==domain,'moving-block guard retains both targets')
    need([r for r in raw if r[2]in(U,V)]==[r for r in actual if r[2]in(U,V)],'complete guarded target domain')
    old=json.loads((RESULTS/supply['old_domain']).read_text())['sources'];reserved=set()
    for z in old:
        for r in z['exact_all_shared_clocks_target_domain']:
            words=tile(r['mask']);need(reserved.isdisjoint(words),'old actual normalized menus disjoint');reserved|=words
    for z in supply['helpers']:
        for name in('before','after'):
            r=z['selected'][name];words=tile(r['mask'])
            need(reserved.isdisjoint(words),'old helper normalized menus disjoint');reserved|=words
    need(len(reserved)==648,'all old28 plus53 helper endpoint words retained')
    weights={}
    with gzip.open(RESULTS/'metagraph_order8_20261023/merged_edges.tsv.gz','rt')as f:
        for z in csv.DictReader(f,delimiter='\t'):
            a,b=int(z['source']),int(z['target'])
            if S in(a,b):weights[a,b]=[int(z['blue']),int(z['black'])]
    records=[]
    for p,m,t in domain:
        at={v:i for i,v in enumerate(p)}
        K=[[0 if i==j else A[i][j]if abs(at[i]-at[j])==1 else A[j][i]for j in range(8)]for i in range(8)]
        cm=mask(op,p[::-1]);words=tile(m)
        need(cm==reflect(m) and mask(K,p)==FULL^m,'literal converse and K complement words')
        need(merged[owner[cm]]==S and merged[owner[FULL^cm]]==t,'converse owners')
        need(guard(p) and all(op[a][b]for a,b in zip(p[::-1],p[-2::-1])),'ordinary block and actual converse Hamiltonian path')
        cp=p[::-1];ca={v:i for i,v in enumerate(cp)}
        need(ca[6]==ca[2]+1 and ca[5]==ca[6]+1 and ca[3]==ca[4]+1,'converse reversed-block guard')
        need(len(words)==4 and reserved.isdisjoint(words),'new source/inverse/converse resources injective')
        reserved|=words
        units={min(tile(mm))for _,mm,tt in actual if tt==t};w=weights[min(S,t),max(S,t)]
        need(w==[0,2*len(units)] and len(units)==1,'exact pure-black unit capacity')
        lifts=[list(pp)for pp,mm,tt in actual if mm in words]
        need(lifts==[list(p)],'complete resource kernel lift is unique')
        need((X[S]-X[t])%8==0,'actual chosen X mod8 phase zero')
        records.append(dict(path=p,mask=m,target=t,phase=0,resource=min(words),closed_words=sorted(words),
                            kernel_lifts=lifts,capacity=dict(blue=w[0],black=w[1],selected_weight=2)))
    need(len(reserved)==656,'old inventory plus8 new closed words, without resource reuse')
    bank=next(b for b in supply['banks']if len(b['source_ids'])==53);ids=bank['source_ids']
    labels={s:encode(label(s))for s in ids+[S]};a,b,c=368,776,3139
    before=set();after=set();repair_witnesses=[]
    for ss in combinations(ids,4):
        if (a in ss)!=(b in ss) and (a in ss)!=(c in ss):before.add(charge(labels,ss))
    need(set(range(2048))-before=={1509},'exact sole missing charge for two XOR constraints')
    for ss in combinations(ids+[S],4):
        if (a in ss)!=(b in ss) and (a in ss)!=(c in ss):
            q=charge(labels,ss);after.add(q)
            if q==1509:repair_witnesses.append(ss)
    need(after==set(range(2048)),'new54 relative-block supply repairs every constrained charge')
    need((776,1596,3139,2125)in repair_witnesses,'specified repaired constrained charge witness')
    residual_triple=sub(1509,labels[a]);residual_pair=sub(sub(1509,labels[b]),labels[c])
    triples=[ss for ss in combinations(ids,3)if charge(labels,ss)==residual_triple]
    need(triples==[(776,1182,2638),(1391,1948,3139),(1400,3111,3139),(2439,2728,3139)],'exact obstructed triple branch')
    fullids=[z['source']for z in supply['all_new_source_labels']];fulllabels={s:encode(label(s))for s in fullids}
    pairs=[ss for ss in combinations(fullids,2)if charge(fulllabels,ss)==residual_pair]
    need(pairs==[(1596,2125)] and residual_pair==239 and residual_triple==305,'unique full72 pair repair')
    out=dict(status='FINITE-EXACT native relative-block repair of coupled pair-predicate obstruction',
             source=S,canon=int(nodes[S]['canon_a']),H=len(actual),aut=len(auts),label=label(S),X=X[S],
             full_UV_count=len(domain),fixed_clocks=common,relative_clock_cycle=cycle,nontrivial_common_clock_sets=common_sets,
             retained_blocks=[[5,6,2],[3,4]],converse_blocks=[[2,6,5],[4,3]],actual_clock_cost=8,
             raw_block_guard_count=len(raw),raw_block_guard_target_counts=dict(Counter(t for p,m,t in raw)),
             records=records,all_reserved_closed_words=len(reserved),
             xor_tree=dict(center=a,leaves=[b,c],occupancy=4,old_count=53,missing=[1509],new_count=54,
                           repaired_charge_count=len(after),repair_witnesses=repair_witnesses,
                           residual_triple=residual_triple,all_triples=triples,residual_pair=residual_pair,full72_pairs=pairs),
             scope='Actual full native path menus, kernels, inverse/converse capacities and target-margin filter; XOR tree is a record predicate, not an asserted shared-owner host embedding.',checks=CHECKS,
             supply_hash=sha256(SUPPLY.read_bytes().replace(b'\r\n',b'\n')).hexdigest())
    out=json.loads(json.dumps(out));path=HERE/'native_2125_receipt.json'
    if '--write'in sys.argv:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',newline='\n')
    else:need(out==json.loads(path.read_text()),'frozen native2125 receipt')
    print(f'PASS {CHECKS}; raw guard{len(raw)}, full U/V2, constrained54 covers{len(after)}')
if __name__=='__main__':main()
