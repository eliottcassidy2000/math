import json, pathlib, collections, math, itertools, hashlib, re, os
BASE=pathlib.Path(__file__).resolve().parent
choices=[pathlib.Path(os.environ['MATH_REPO'])] if 'MATH_REPO' in os.environ else [BASE.parent/'math',BASE,*BASE.parents]
REPO=next(p for p in choices if (p/'01-canon').is_dir() and (p/'05-knowledge/results').is_dir())
R=REPO/'05-knowledge/results'
O=pathlib.Path(__file__).parent
sources=[R/'tournament_tiling_metagraph_address_codex_S4.json',R/'tournament_tiling_metagraph_address_n7_codex_S4.json']
atlases=json.loads(sources[0].read_text())['sizes']+[json.loads(sources[1].read_text())]
def graph(n,tile_order,mask):
    a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for b,(v,u) in enumerate(tile_order):
        if mask>>b&1:a[u-1]^=1<<(v-1);a[v-1]^=1<<(u-1)
    return a
def paths(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not s>>v&1:continue
            t=a[v]&~s
            while t:
                b=t&-t;w=b.bit_length()-1;t-=b;dp[s|b][w]+=dp[s][v]
    return sum(dp[-1])
def strong(a):
    n=len(a);reach=[x|(1<<i) for i,x in enumerate(a)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(x==(1<<n)-1 for x in reach)
def brute_aut_sc(a):
    n=len(a);degrees=[x.bit_count() for x in a];auts=antis=0
    for p in itertools.permutations(range(n)):
        if all(degrees[v]==degrees[p[v]] for v in range(n)):
            auts+=all(bool(a[u]>>v&1)==bool(a[p[u]]>>p[v]&1) for u in range(n) for v in range(u+1,n))
        if all(degrees[v]+degrees[p[v]]==n-1 for v in range(n)):
            antis+=all(bool(a[u]>>v&1)!=bool(a[p[u]]>>p[v]&1) for u in range(n) for v in range(u+1,n))
    return auts,antis
rows=[];summaries={};atoms=collections.defaultdict(collections.Counter);centers=collections.defaultdict(collections.Counter)
for atlas in atlases:
    n=atlas['n'];counts=collections.defaultdict(collections.Counter); labelsum=0;fibres=0
    for node in atlas['nodes']:
        a=graph(n,atlas['tile_order'],node['tiling_masks'][0]);H=paths(a);sc=node['self_converse'];mult=1 if sc else 2;h=node['tiling_count']//mult
        assert node['tiling_count']%mult==0 and H%h==0
        aut=H//h;isstrong=strong(a)
        if n<=6:
            aut_b,anti_b=brute_aut_sc(a)
            assert aut==aut_b and sc==bool(anti_b)
        counts[H]['ordinary']+=mult;counts[H]['self_converse']+=sc;counts[H]['merged']+=1;counts[H]['labelled']+=mult*math.factorial(n)//aut
        labelsum+=mult*math.factorial(n)//aut;fibres+=mult*h
        if isstrong:
            atoms[H][n]+=mult
            if sc:centers[H][n]+=1
        rows.append(dict(n=n,node=node['id'],H=H,h=h,aut=aut,self_converse=sc,ordinary_multiplicity=mult,strong=isstrong,adjacency_outmasks=a))
    assert labelsum==1<<(n*(n-1)//2)
    assert fibres==1<<((n-1)*(n-2)//2)
    assert sum(x['ordinary'] for x in counts.values())==atlas['classes']
    summaries[n]={H:dict(v) for H,v in sorted(counts.items())}
# Independently frozen labelled n=8 histogram, not converted to ordinary counts.
p=R/'h_spectrum_n8_exhaustive_monad.out';sources.append(p)
label8={int(h):int(k) for h,k in re.findall(r'^\s*H=\s*(\d+):\s*(\d+)\s*$',p.read_text(),re.M)}
assert sum(label8.values())==1<<28
result=dict(schema=1,source_pins={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},counts_by_order=summaries,strong_atom_inventory_through_order7={h:dict(sorted(c.items())) for h,c in sorted(atoms.items())},self_converse_strong_inventory_through_order7={h:dict(sorted(c.items())) for h,c in sorted(centers.items())},labelled_n8=label8,records=rows)
n7_reps=[]
for row in rows:
 if row['n']!=7:continue
 a=row['adjacency_outmasks'];n7_reps.append(' '.join(map(str,a)))
 if not row['self_converse']:n7_reps.append(' '.join(str(x^127^(1<<i)) for i,x in enumerate(a)))
assert len(n7_reps)==456
(O/'n7_representatives.txt').write_text('\n'.join(n7_reps)+'\n')
(O/'census.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'class_counts':{n:sum(c['ordinary'] for c in d.values()) for n,d in summaries.items()},'atoms_H_below45':{h:dict(c) for h,c in sorted(atoms.items()) if h<45},'centers_H_below45':{h:dict(c) for h,c in sorted(centers.items()) if h<45}},indent=2))
