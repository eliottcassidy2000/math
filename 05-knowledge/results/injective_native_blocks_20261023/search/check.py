"""Independent literal replay of selected paths, shared coordinates and current circuits.

This checker uses no discovery cache. All paths of the seven actual source
tournaments are enumerated again. --write regenerates the small audit receipt.
"""
import csv,gzip,json,sys
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,permutations,product
import explore as e
CHECKS=0
def need(b,why):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(why)

def adjacency(canon):
    A=[[0]*8 for _ in range(8)]
    for bit,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
    return A

def mask_on(A,p):return sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(e.CHORDS))

def main():
    doc=json.loads((e.HERE/'selected_blocks.json').read_text());slots=doc['slots'];targets=doc['targets']
    edges=[tuple(int(r[k]) for k in ('source','target','blue','black')) for r in csv.DictReader(gzip.decompress((e.NATIVE/'merged_edges.tsv.gz').read_bytes()).decode().splitlines(),delimiter='\t')]
    weights={(i,j):(b,k) for i,j,b,k in edges}
    def apply(x):
        y=[0]*3528
        for i,j,b,k in edges:
            v=(b+k)*(x[i]-x[j]);y[i]+=v;y[j]-=v
        return y
    need(all(v%256==0 for v in apply(e.X)),'original full Laplacian certifies frozen X mod256')
    for y in e.Y:need(all(v%4==0 for v in apply(y)),'original modulo4 harmonic lift')
    need(all(v%8==0 for v in apply(e.Z)),'original modulo8 primitive lift')
    sparse=[0]*3528;sparse[targets[0]]=4
    need(all(v%8==0 for v in apply(sparse)),'actual sparse target separator')
    need(sparse[targets[0]]!=sparse[targets[1]],'individual universal phase changes')
    images=json.loads((e.HERE.parent/'electrical/image_witnesses.json').read_text())
    need(images['input_edges_sha256']==sha256((e.NATIVE/'merged_edges.tsv.gz').read_bytes()).hexdigest(),'image witnesses pin exact compressed original edge bytes')
    need(images['modulus']==256 and len(images['solutions'])==3,'all three exact image certificates present')
    expected_charges=[{1425:-1,1477:1,1921:-1,2837:1},{1965:1,2234:-1},{1394:-1,1477:1,1841:1,2185:-1}]
    for sol,charge in zip(images['solutions'],expected_charges):
        need(dict(sol['charge'])==charge and sol['coefficient']==2,'exact paired source or target image charge')
        z=sol['potential_mod256'];need(len(z)==3528 and z[3527]==0 and all(type(a) is int and 0<=a<256 for a in z),'full grounded image potential')
        az=apply(z);rhs=[2*charge.get(i,0) for i in range(3528)]
        need(all((a-b)%256==0 for a,b in zip(az,rhs)),'literal full native Lz=paired charge mod256')
        need([[i,(a-b)//256] for i,(a,b) in enumerate(zip(az,rhs)) if a!=b]==sol['integer_correction_div256_nonzero'],'full integral image correction retained')
    za,_,zb=[s['potential_mod256'] for s in images['solutions']]
    zc=[(a-b)%256 for a,b in zip(za,zb)];az=apply(zc)
    wc={s:expected_charges[0].get(s,0)-expected_charges[2].get(s,0) for s in set(expected_charges[0])|set(expected_charges[2])}
    need(all((a-2*wc.get(i,0))%256==0 for i,a in enumerate(az)),'derived six-source full image certificate')
    sourceinfo=[];sourcepaths={}
    for source in sorted({s['source'] for s in slots}):
        rows=e.paths(source);sourcepaths[source]=rows;A=adjacency(int(e.NODES[source]['canon_a']))
        aut=sum(all(A[i][j]==A[p[i]][p[j]] for i,j in combinations(range(8),2)) for p in permutations(range(8)))
        need(aut==1,'actual selected source has trivial automorphism group')
        caps=[]
        for t in targets:
            paths=[r['path'] for r in rows if r['target']==t]
            b,k=weights[tuple(sorted((source,t)))]
            need(b==0 and k==2*len(paths),'all candidate capacities are literal pure-black paired paths')
            caps.append(dict(target=t,blue=b,black=k,paths=paths))
        sourceinfo.append(dict(source=source,canon=int(e.NODES[source]['canon_a']),H=len(rows),aut=aut,x=e.X[source],sigma=e.SIG[source],capacities=caps))
    for slot in slots:
        source=slot['source'];A=adjacency(slot['canon']);need(slot['canon']==int(e.NODES[source]['canon_a']),'canonical source preserved')
        for key in ('before','after'):
            row=slot[key];p=row['path'];need(sorted(p)==list(range(8)),'actual owner permutation')
            need(all(A[a][b] for a,b in zip(p,p[1:])),'actual Hamiltonian path')
            mask=mask_on(A,p);need(mask==row['mask'],'exact actual clock forward mask')
            need(e.MERGED[e.OWNER[mask]]==source,'actual canonical source owner')
            # Build K_P(A) literally, preserving only the seven selected path arcs.
            pathpairs={frozenset(z) for z in zip(p,p[1:])}
            K=[[0]*8 for _ in range(8)]
            for i,j in combinations(range(8),2):
                K[i][j]=A[i][j] if frozenset((i,j)) in pathpairs else A[j][i]
                K[j][i]=1-K[i][j]
            km=mask_on(K,p);need(km==(e.FULL^mask),'literal K_P target has complementary free word')
            need(e.MERGED[e.OWNER[km]]==row['target'],'actual K_P target owner')
            KK=[[0]*8 for _ in range(8)]
            for i,j in combinations(range(8),2):
                KK[i][j]=K[i][j] if frozenset((i,j)) in pathpairs else K[j][i]
                KK[j][i]=1-KK[i][j]
            need(KK==A,'literal K_P involution restores source on same owner coordinates')
            op=[[A[j][i] for j in range(8)] for i in range(8)];rev=list(reversed(p))
            need(all(op[a][b] for a,b in zip(rev,rev[1:])),'literal converse path is legal')
            cm=mask_on(op,rev)
            need(e.MERGED[e.OWNER[cm]]==source and e.MERGED[e.OWNER[e.FULL^cm]]==row['target'],'converse partners preserve source and target')
            need(row['phase']==(e.X[source]-e.X[row['target']])%8,'original chosen per-slot phase')
            need(any(r==row for r in sourcepaths[source]),'path occurs in complete independent source enumeration')
        p,q=slot['before']['path'],slot['after']['path'];pins=slot['pins']
        need((p[0],p[-1],slot['before']['phase'])==tuple(pins[:3]) and (q[0],q[-1],slot['after']['phase'])==tuple(pins[:3]),'endpoints and chosen phase fixed')
        k,a,b=pins[3:];need(p[k:k+2]==q[k:k+2]==[a,b],'literal indexed shared arc fixed')
        need(slot['all_unchanged_clocks']==[[i,a] for i,(a,b) in enumerate(zip(p,q)) if a==b],'all actual unchanged clock coordinates retained')
        need(sum(a!=b for a,b in zip(p,q))==slot['cost'],'exact owner-clock occurrence cost')
    def charge(ids):
        w=Counter()
        for i in ids:w[slots[i]['source']]+=slots[i]['sign']
        return {s:c for s,c in sorted(w.items()) if c}
    def current(w,x):
        u,v=targets;out={s:2*c*(x[u]-x[v]) for s,c in w.items()}
        z=2*sum(c*x[s] for s,c in w.items());out[u]=z;out[v]=-z
        return out
    def passes(w):
        return sum(w.values())==0 and all(v%256==0 for v in current(w,e.X).values()) and all(v%8==0 for v in current(w,e.Z).values()) and all(all(v%8==0 for v in current(w,[2*a for a in y]).values()) for y in e.Y)
    blocks=[]
    for name,ids in [('first',list(range(4))),('second_reversed',list(range(4,8))),('six',[1,2,3,5,6,7])]:
        w=charge(ids);debt=current(w,e.X)
        need(all(v==0 for v in debt.values()),'complete original X integer debt vanishes')
        need(passes(w),'authenticated paired character and X tests')
        admissible=[]
        for bits in product(range(2),repeat=len(ids)):
            picked=[i for i,b in zip(ids,bits) if b]
            if passes(charge(picked)):admissible.append(bits)
        need(admissible==[tuple([0]*len(ids)),tuple([1]*len(ids))],'supplied circuit has no proper accepted subblock')
        blocks.append(dict(name=name,slots=ids,charge=w,integer_X_debt=debt,accepted_binary_states=admissible,
                           sharp_slot_budget=len(ids),owner_clock_occurrences=sum(slots[i]['cost'] for i in ids)))
    admissible=[];margins=Counter(s['before']['target'] for s in slots)
    for bits in product(range(2),repeat=8):
        state=[s['after'] if b else s['before'] for s,b in zip(slots,bits)]
        inj=len({(s['source'],tuple(r['path'])) for s,r in zip(slots,state)})==8
        shared=all(state[a]['path'][i]==state[b]['path'][j] for (a,i),(b,j) in doc['shared_owner_equalities'])
        if not inj or not shared or Counter(r['target'] for r in state)!=margins:continue
        if passes(charge([i for i,b in enumerate(bits) if b])):admissible.append(bits)
    expected=[(0,0,0,0,0,0,0,0),(0,1,1,1,0,1,1,1),(1,0,0,0,1,0,0,0),(1,1,1,1,1,1,1,1)]
    need(admissible==expected,'complete eight-slot accepted state space')
    normalized_states=[]
    for bits in ([0]*8,[1]*8):
        state=[s['after'] if b else s['before'] for s,b in zip(slots,bits)]
        need(len({(s['class_a'],tuple(r['path'])) for s,r in zip(slots,state)})==8,'ordinary before/after actual selected inventory injective')
        ports={(s['class_a'],tuple(r['path'])) for s,r in zip(slots,state)}|{(s['class_b'],tuple(reversed(r['path']))) for s,r in zip(slots,state)}
        need(len(ports)==16,'converse-closed before/after inventories injective')
        sourcewords=[]
        for s,r in zip(slots,state):
            A=adjacency(s['canon']);op=[[A[j][i] for j in range(8)] for i in range(8)]
            sourcewords.extend([r['mask'],mask_on(op,list(reversed(r['path'])))])
        targetwords=[e.FULL^m for m in sourcewords]
        need(len(set(sourcewords))==16,'normalized ordinary/converse source ports injective')
        need(len(set(targetwords))==16,'actual inverse target ports injective without arbitrary target isomorphisms')
        need(set(sourcewords).isdisjoint(targetwords),'all source and inverse target endpoint ports disjoint')
        need(all(e.MERGED[e.OWNER[m]] in targets for m in targetwords),'every inverse target port has correct actual owner')
        normalized_states.append(dict(source_words=sourcewords,inverse_target_words=targetwords))
    reservoir=sorted({m for state in normalized_states for words in state.values() for m in words})
    need(len(reservoir)==56,'full before/after reservoir has56 distinct endpoint ports')
    # Mathematically, the selected response changes by DeltaL. It is not a
    # replacement of the complete all-path inventory, which already has both.
    row={'status':'FINITE-EXACT selected native circuits; direct full image certificates prove universal mod256 current preservation; all-dyadic extension uses electrical exponent theorem',
         'source_file_hashes':{p.relative_to(e.RESULTS).as_posix():sha256(p.read_bytes() if p.suffix=='.gz' else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [e.NATIVE/'merged_nodes.tsv',e.NATIVE/'owner_map.u16.gz',e.NATIVE/'merged_edges.tsv.gz',e.RESULTS/'metagraph_order8_electrical_20261023/receipt.json']},
         'hash_policy':'raw gzip bytes; text bytes normalized CRLF to LF',
         'sources':sourceinfo,'targets':targets,'target_X':[e.X[t] for t in targets],'target_sigma':[e.SIG[t] for t in targets],
         'circuits':blocks,'eight_slot_accepted_states':admissible,'eight_slot_sharp_connected_move_budget':6,
         'normalized_endpoint_inventories':normalized_states,'fixed_background_omit_endpoint_words':reservoir,
         'per_slot_universal_phase_hostile':{'character':'4e_1965','difference':4,'harmonic_mod8':True},
         'checks':CHECKS,'scope':'selected injective converse-closed inventory with an explicitly unoccupied reservoir; unchanged full native atlas; aggregate all-bank current preservation differs from per-slot universal phase preservation'}
    row=json.loads(json.dumps(row));file=e.HERE/'receipt.json'
    if '--write' in sys.argv:file.write_text(json.dumps(row,indent=2))
    else:need(json.loads(file.read_text())==row,'saved exact receipt replays')
    print('PASS',CHECKS,'checks;',sum(r['H'] for r in sourceinfo),'actual source paths; circuit slot budgets4,4,6; eight-menu budget6')

if __name__=='__main__':main()
