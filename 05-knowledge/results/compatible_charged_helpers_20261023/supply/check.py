#!/usr/bin/env python3
"""Independent native supply and fixed-cardinality charge-section certificate."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import argparse,csv,gzip,json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
CHECKS=0;U,V=1965,2234;FULL=(1<<21)-1
CHORDS=[(i,j) for i in range(8) for j in range(i+2,8)]
PAIRS=list(combinations(range(8),2));REV=[CHORDS.index((7-j,7-i)) for i,j in CHORDS]
def need(v,m):
    global CHECKS
    CHECKS+=1
    if not v:raise RuntimeError(m)
def matrix(c):
    A=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(PAIRS):A[i][j]=(c>>(27-k))&1;A[j][i]=1-A[i][j]
    return A
def mask(A,p):return sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(CHORDS))
def reflect(m):return sum(((m>>k)&1)<<REV[k] for k in range(21))
def tile(m):return {m,reflect(m),FULL^m,FULL^reflect(m)}
def partner(A,p):
    at={v:i for i,v in enumerate(p)}
    return [[0 if i==j else A[i][j] if abs(at[i]-at[j])==1 else A[j][i] for j in range(8)] for i in range(8)]
def charge(labels,indices):
    x=0;y=0
    for i in indices:
        a=labels[i];x+=a[0]
        y^=sum(v<<j for j,v in enumerate(a[1:]))
    return ((x%128)<<4)|y
def group_sum(a,b):return ((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    cert=json.loads((HERE/'selected_supply.json').read_text())
    nodes=list(csv.DictReader((RESULTS/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
    merged={int(z[k]):int(z['merged_id']) for z in nodes for k in ('class_a','class_b')}
    owner=array('H');owner.frombytes(gzip.decompress((RESULTS/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder!='little':owner.byteswap()
    D=json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4']
    label=lambda s:[X[s]%128]+[Y[j][s]%2 for j in range(1,5)]
    edges=[];weights={};neighbors={U:set(),V:set()}
    with gzip.open(RESULTS/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as f:
        for z in csv.DictReader(f,delimiter='\t'):
            a,b,blue,black=(int(z[k]) for k in ('source','target','blue','black'))
            edges.append((a,b,blue+black));weights[a,b]=(blue,black)
            if a in neighbors:neighbors[a].add(b)
            if b in neighbors:neighbors[b].add(a)
    # Direct Laplacian authentication of the actual finite potentials.
    actions=[[0]*len(nodes) for _ in range(6)]
    for a,b,w in edges:
        for k,f in enumerate([X]+Y[:5]):
            d=w*(f[a]-f[b]);actions[k][a]+=d;actions[k][b]-=d
    need(all(v%256==0 for v in actions[0]),'original primitive X harmonic modulo256')
    need(all(all(v%4==0 for v in z) for z in actions[1:]),'all five original radical potentials harmonic modulo4')
    need(label(570)==label(U)==label(V) and X[570]==X[U]==X[V]==146,'actual three-target paired fibre and chosen X phase')
    need([i for i in range(len(nodes)) if label(i)==label(U)]==[570,U,V],'full equal-paired-target fibre')
    # Harmonic4eU is an explicit hostile to universal ordinary per-slot phases.
    haction=[0]*len(nodes)
    for a,b,w in edges:
        d=w*(4*int(a==U)-4*int(b==U));haction[a]+=d;haction[b]-=d
    need(all(v%8==0 for v in haction),'4eU is genuine original harmonic mod8')
    need((4-0)%8==4,'ordinary U/V target phases really differ')
    old=json.loads((RESULTS/cert['old_domain']).read_text())['sources']
    oldids={z['source'] for z in old};allwords=set();oldtiles=set()
    for z in old:
        A=matrix(int(nodes[z['source']]['canon_a']))
        fixed=[[i,a] for i,(a,b) in enumerate(zip(z['old_before'],z['old_after'])) if a==b]
        need(fixed==z['fixed_clocks'],'old exact indexed pins retained')
        for r in z['exact_all_shared_clocks_target_domain']:
            pp=r['path'];m=mask(A,pp)
            need(all(A[a][b] for a,b in zip(pp,pp[1:])) and m==r['mask'],'old actual full source record')
            need(all(pp[i]==a for i,a in fixed),'old shared owner values retained')
            words=tile(m);need(len(words)==4 and allwords.isdisjoint(words),'old normalized resource inventory')
            allwords|=words;oldtiles.add(min(words))
    need(len(oldids)==28 and len(oldtiles)==56 and len(allwords)==224,'unchanged old28 resources and inverse/converse words')
    actualnew=sorted((neighbors[U]&neighbors[V])-oldids)
    need(actualnew==[z['source'] for z in cert['all_new_source_labels']],'complete72 new source alphabet from actual atlas')
    for z in cert['all_new_source_labels']:need(z['label']==label(z['source']),'all72 exact source charges')
    details=[];by={z['source']:z for z in cert['helpers']};helpertiles=set()
    for s,z in sorted(by.items()):
        need(s in actualnew and s not in oldids and s not in (U,V),'new helper physical source row distinct')
        need(z['canon']==int(nodes[s]['canon_a']) and z['label']==label(s) and z['X']==X[s],'actual helper matrix and full paired charge')
        A=matrix(z['canon']);op=[list(x) for x in zip(*A)];actual=[];auts=[]
        for pp in permutations(range(8)):
            if all(A[i][j]==A[pp[i]][pp[j]] for i,j in PAIRS):auts.append(pp)
            if all(A[a][b] for a,b in zip(pp,pp[1:])):
                m=mask(A,pp);need(merged[owner[m]]==s,'actual source owner from entire canonical enumeration')
                actual.append((pp,m,merged[owner[FULL^m]]))
        need(len(actual)==z['H'],'independent8! complete source path census')
        multiplicity=Counter(m for pp,m,t in actual)
        need(set(multiplicity.values())=={len(auts)},'full kernel multiplicity equals source automorphism order')
        selected=z['selected'];a=selected['before'];b=selected['after']
        need(a['target']==U and b['target']==V,'binary bit0 is literal U and bit1 is literal V')
        fixed=[[i,x] for i,(x,y) in enumerate(zip(a['path'],b['path'])) if x==y]
        need(fixed==selected['common_clocks'] and fixed,'explicit nonempty independently chosen helper pins')
        raw=[r for r in actual if all(r[0][i]==x for i,x in fixed)]
        domain=[r for r in raw if r[2] in (U,V)]
        need({(tuple(r['path']),r['mask'],r['target']) for r in (a,b)}==set(domain) and len(domain)==2,
             'complete two-target pin domain has exactly the two supplied full records')
        need(len(raw)==selected['complete_pinned_count'],'raw domain includes all excluded-target paths')
        need(selected['cost']==8-len(fixed),'actual owner-clock replacement cost')
        capacities=[];lifts=[]
        for r in (a,b):
            pp=r['path'];m=mask(A,pp);t=r['target'];cm=mask(op,pp[::-1]);K=partner(A,pp)
            need(cm==reflect(m) and mask(K,pp)==FULL^m,'literal K and converse normalized endpoints')
            need(merged[owner[m]]==s and merged[owner[cm]]==s and merged[owner[FULL^m]]==t and merged[owner[FULL^cm]]==t,'all four actual source/target owners')
            words={m,cm,FULL^m,FULL^cm}
            need(len(words)==4 and m!=cm,'selected tile is black, with literal paired conductance TWO')
            need(allwords.isdisjoint(words),'no source inverse converse normalized resource recycling')
            allwords|=words;helpertiles.add(min(words))
            need(sorted(words)==r['orbit'] and min(words)==r['resource'],'recorded exact normalized unit resource')
            blue,black=weights[min(s,t),max(s,t)]
            need(blue==0 and black>=2 and black%2==0,'original edge has pure black integer multiplicity')
            exacttiles={min(tile(q)) for _,q,tt in actual if tt==t}
            need(black==2*len(exacttiles),'edge multiplicity equals distinct normalized tiles times TWO')
            tilelifts=[list(p) for p,q,tt in actual if q in words]
            pinlifts=[p for p in tilelifts if all(p[i]==x for i,x in fixed)]
            need(pinlifts==[pp],'every full resource lift retained; exact pins select unique witness')
            need((X[s]-X[t])%8==r['phase'] and r['phase']==(X[s]-146)%8,'same chosen original phase, not all ordinary phases')
            capacities.append(dict(target=t,blue=blue,black=black,normalized_units=len(exacttiles),selected_weight=2))
            lifts.append(dict(target=t,normalized_resource=min(words),full_canonical_resource_lifts=tilelifts,
                              selected_pin_lift=pp,source_word_kernel_size=multiplicity[m]))
        details.append(dict(source=s,H=len(actual),aut=len(auts),self_converse=nodes[s]['class_a']==nodes[s]['class_b'],
                            raw_pinned_count=len(raw),excluded_target_counts=dict(Counter(t for _,_,t in raw if t not in (U,V))),
                            fixed_clocks=fixed,cost=selected['cost'],capacities=capacities,kernel_lifts=lifts))
    need(len(by)==53 and len(helpertiles)==106 and len(allwords)==648,'old224 plus424 genuinely new closed endpoint words')
    banks=[]
    for bank in cert['banks']:
        ids=bank['source_ids'];k=bank['V_occupancy'];h=len(ids);labels=[label(s) for s in ids]
        need(h==bank['helper_count'] and len(set(ids))==h and set(ids)<=set(by),'bank actual sources retained')
        image=set();hist=Counter()
        for indices in combinations(range(h),k):
            q=charge(labels,indices);image.add(q);hist[q]+=1
        need(image==set(range(2048)),'complete fixed-k subset charge image covers all2048')
        section=bank['charge_section'];need(set(map(int,section))==set(range(2048)),'complete section keys')
        for text,m in section.items():
            need(m>=0 and m<(1<<h) and m.bit_count()==k,'each section keeps exact helper cardinality')
            need(charge(labels,[i for i in range(h) if m>>i&1])==int(text),'independent section witness charge recomputation')
        need(bank['atomic_support_bound']==2+2*min(k,h-k),'conditional paired-swap plus bank support bound')
        banks.append(dict(helpers=h,V_occupancy=k,charge_values=len(image),subset_count=sum(hist.values()),
                          min_witness_multiplicity=min(hist.values()),max_witness_multiplicity=max(hist.values()),
                          atomic_support_bound=bank['atomic_support_bound'],normalized_tile_reservation=2*h,
                          closed_word_reservation=8*h))
    mainbank=cert['banks'][cert['main_bank_index']];ids=mainbank['source_ids'];labels=[label(s) for s in ids]
    deletion=[]
    for i,s in enumerate(ids):
        image={charge(labels,j) for j in combinations([j for j in range(len(ids)) if j!=i],3)}
        need(len(image)<2048,'main40 bank is inclusion-minimal under individual deletion')
        deletion.append(dict(source=s,remaining_coverage=len(image),first_missing=min(set(range(2048))-image)))
    all_labels=[z['label'] for z in cert['all_new_source_labels']]
    pairs={charge(all_labels,p) for p in combinations(range(72),2)}
    repeated=pairs|{charge(all_labels,(i,i)) for i in range(72)}
    missing=sorted(set(range(2048))-repeated)
    need(len(pairs)==1468 and len(repeated)==1478 and missing[0]==2,'universal two-slot charge section impossible even granting duplicate source labels')
    need(missing==cert['two_sum_obstruction']['missing_even_with_repeated_source_labels'],'complete two-sum obstruction set')
    files=['metagraph_order8_20261023/merged_nodes.tsv','metagraph_order8_20261023/owner_map.u16.gz',
           'metagraph_order8_20261023/merged_edges.tsv.gz','metagraph_order8_electrical_20261023/receipt.json',
           cert['old_domain'],'injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md']
    def digest(p):
        b=p.read_bytes()
        if p.suffix in ('.json','.md','.tsv','.py'):b=b.replace(b'\r\n',b'\n')
        return sha256(b).hexdigest()
    out=dict(status='FINITE-EXACT authenticated native supply and surjective charge sections',
        selected_helpers=len(by),old_slots=28,old_closed_words=224,all_reserved_closed_words=len(allwords),
        helper_details=details,banks=banks,main_bank_deletion_controls=deletion,
        two_sum=dict(new_source_alphabet=72,distinct_pairs=len(pairs),with_repeated_labels=len(repeated),first_missing=2),
        guard_scope='all original dyadic AGGREGATE currents by cited complete paired-code theorem; chosen X mod8 phases; NOT all ordinary per-slot phases',
        source_sha256={f:digest(RESULTS/f) for f in files},certificate_sha256=digest(HERE/'selected_supply.json'),
        checks=CHECKS)
    out=json.loads(json.dumps(out))
    output=HERE/'receipt.json'
    if args.write:output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    else:need(json.loads(output.read_text())==out,'frozen native supply receipt matches')
    print(f'PASS {CHECKS:,}; '+('wrote receipt' if args.write else 'receipt matches'))
if __name__=='__main__':main()
