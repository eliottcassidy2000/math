"""Literal validation and exhaustive binary-menu test for8/10/12-source circuits."""
import csv,gzip,json,sys
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations,product
import explore as e
from check import adjacency,mask_on
COUNT=0
def need(b,why):
    global COUNT
    COUNT+=1
    if not b:raise RuntimeError(why)
def main():
    stem='larger16' if '--sixteen' in sys.argv else 'larger'
    doc=json.loads((e.HERE/(stem+'_selected.json')).read_text());targets=doc['targets']
    weights={tuple(int(r[k]) for k in ('source','target')):tuple(int(r[k]) for k in ('blue','black')) for r in csv.DictReader(gzip.decompress((e.NATIVE/'merged_edges.tsv.gz').read_bytes()).decode().splitlines(),delimiter='\t')}
    def image_check(sol,charge):
        need(dict(sol['charge'])==charge and sol['coefficient']==2,'literal larger source/target image charge')
        z=sol['potential_mod256'];need(len(z)==3528 and z[3527]==0 and all(type(a)is int and 0<=a<256 for a in z),'grounded full larger image vector')
        out=[0]*3528
        for (i,j),(b,k) in weights.items():
            d=(b+k)*(z[i]-z[j]);out[i]+=d;out[j]-=d
        need(all((v-2*charge.get(i,0))%256==0 for i,v in enumerate(out)),'full native Lz verifies larger charge mod256')
        need([[i,(v-2*charge.get(i,0))//256] for i,v in enumerate(out) if v!=2*charge.get(i,0)]==sol['integer_correction_div256_nonzero'],'full larger integer correction')
    images=json.loads((e.HERE.parent/'electrical/larger_image.json').read_text())
    need(images['input_edges_sha256']==sha256((e.NATIVE/'merged_edges.tsv.gz').read_bytes()).hexdigest(),'larger images pin actual original edge bytes')
    image_check(json.loads((e.HERE.parent/'electrical/image_witnesses.json').read_text())['solutions'][1],{1965:1,2234:-1})
    rowsby={};sourceinfo={};reports=[]
    for source in sorted({s['source'] for c in doc['circuits'] for s in c['slots']}):
        rows=e.paths(source);rowsby[source]=rows;A=adjacency(int(e.NODES[source]['canon_a']))
        aut=sum(all(A[i][j]==A[p[i]][p[j]] for i,j in combinations(range(8),2)) for p in permutations(range(8)))
        need(aut==1,'trivial actual source automorphism group')
        capacities=[]
        for t in targets:
            paths=[r['path'] for r in rows if r['target']==t];b,k=weights[tuple(sorted((source,t)))]
            need(b==0 and k==2*len(paths),'pure-black literal source and target capacities')
            capacities.append(dict(target=t,blue=b,black=k,ordinary_paths=paths))
        sourceinfo[source]=dict(canon=int(e.NODES[source]['canon_a']),H=len(rows),aut=aut,X=e.X[source],capacities=capacities)
    for circuit in doc['circuits']:
        slots=circuit['slots'];n=len(slots);need(n==circuit['size'] and len({s['source'] for s in slots})==n,'distinct actual sources')
        normalized=[[],[]]
        for slot in slots:
            source=slot['source'];A=adjacency(slot['canon']);op=[[A[j][i] for j in range(8)] for i in range(8)]
            need(slot['canon']==int(e.NODES[source]['canon_a']),'source adjacency pinned')
            for side,name in enumerate(('before','after')):
                row=slot[name];p=row['path'];need(sorted(p)==list(range(8)),'actual full owner permutation')
                need(all(A[a][b] for a,b in zip(p,p[1:])),'actual Hamiltonian path')
                mask=mask_on(A,p);need(mask==row['mask'] and e.MERGED[e.OWNER[mask]]==source,'actual source word and owner')
                pathpairs={frozenset(z) for z in zip(p,p[1:])};K=[[0]*8 for _ in range(8)]
                for i,j in combinations(range(8),2):
                    K[i][j]=A[i][j] if frozenset((i,j)) in pathpairs else A[j][i];K[j][i]=1-K[i][j]
                km=mask_on(K,p);need(km==e.FULL^mask and e.MERGED[e.OWNER[km]]==row['target'],'literal K_P inverse target port')
                rev=list(reversed(p));need(all(op[a][b] for a,b in zip(rev,rev[1:])),'actual converse path legal')
                cm=mask_on(op,rev);need(e.MERGED[e.OWNER[cm]]==source and e.MERGED[e.OWNER[e.FULL^cm]]==row['target'],'actual converse source and target owners')
                need(row in rowsby[source],'complete source enumeration contains literal path')
                need(row['phase']==(e.X[source]-e.X[row['target']])%8,'chosen original electrical phase')
                normalized[side].extend([mask,cm])
            p,q=slot['before']['path'],slot['after']['path'];a,b,phase,k,u,v=slot['pins']
            need(p[0]==q[0]==a and p[-1]==q[-1]==b,'endpoint pins')
            need(slot['before']['phase']==slot['after']['phase']==phase,'per-slot chosen phase retained')
            need(p[k:k+2]==q[k:k+2]==[u,v],'actual shared indexed arc retained')
            need(slot['all_unchanged_clocks']==[[i,u] for i,(u,v) in enumerate(zip(p,q)) if u==v],'every declared shared clock coordinate retained')
            need(sum(a!=b for a,b in zip(p,q))==slot['cost'],'actual owner-clock cost')
        allports=set()
        for sourcewords in normalized:
            targetwords=[e.FULL^m for m in sourcewords]
            need(len(set(sourcewords))==2*n and len(set(targetwords))==2*n,'source and inverse-target ports both injective')
            need(set(sourcewords).isdisjoint(targetwords),'source and inverse target endpoint resources disjoint')
            allports.update(sourcewords);allports.update(targetwords)
        need(len(allports)==8*n,'whole before-after endpoint reservoir is distinct')
        def state_code(bits):
            a=[slot for bit,slot in zip(bits,slots) if bit];margin=sum(s['sign'] for s in a)
            x=sum(s['sign']*e.X[s['source']] for s in a)
            y=[sum(s['sign']*e.Y[j][s['source']] for s in a)%2 for j in range(1,5)]
            return margin,x,y
        allstates=[];level_counts={m:0 for m in (4,8,16,32,64,128,256)}
        for bits in product((0,1),repeat=n):
            margin,x,y=state_code(bits)
            if margin or any(y):continue
            for m in level_counts:
                if x%(m//2)==0:level_counts[m]+=1
            if x%128==0:allstates.append(bits)
        need(allstates==[(0,)*n,(1,)*n],'exhaustive supplied binary menu has exactly two full-dyadic-current states')
        if n==16:need(list(level_counts.values())==[378,194,86,38,18,4,2],'complete lower-dyadic hierarchy on sixteen slots')
        # Any failure of the effective source code is witnessed by an actual
        # dyadic harmonic; completeness and top exponent are audited separately.
        x=sum(s['sign']*e.X[s['source']] for s in slots)
        need(x==0 and e.X[targets[0]]==e.X[targets[1]],'full integer X current debt exactly zero')
        need(sum(s['sign'] for s in slots)==0,'all actual target degrees retained')
        if n in (12,16):image_check(images['solutions'][0 if n==12 else 1],{s['source']:s['sign'] for s in slots})
        reports.append(dict(size=n,positive=circuit['positive'],negative=circuit['negative'],X_sum=circuit['exact_X_sum'],
                            accepted_binary_states=allstates,states_tested=1<<n,sharp_slot_budget=n,owner_clock_occurrences=sum(s['cost'] for s in slots),
                            accepted_states_by_dyadic_modulus=level_counts,
                            normalized_source_inventories=normalized,fixed_background_omit_endpoint_words=sorted(allports),
                            source_charge=[[s['source'],s['sign']] for s in slots]))
    out={'status':'FINITE-EXACT literal larger source circuits; universal dyadic-current sufficiency via electrical complete-code theorem',
         'source_table':sourceinfo,'targets':targets,'circuits':reports,'checks':COUNT,
         'scope':'sharp budgets only for the supplied binary path menus; actual pin schemas retained; aggregate currents, not universal per-slot phases; selected inventories, not complete atlas rewiring'}
    out=json.loads(json.dumps(out));file=e.HERE/(stem+'_receipt.json')
    if '--write' in sys.argv:file.write_text(json.dumps(out,indent=2))
    else:need(json.loads(file.read_text())==out,'saved exact receipt')
    print('PASS',COUNT,'checks;',sum(r['H'] for r in sourceinfo.values()),'actual source paths;',sum(1<<r['size'] for r in reports),'binary states; sharp budgets',[r['size'] for r in reports])
if __name__=='__main__':main()
