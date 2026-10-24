#!/usr/bin/env python3
"""Exact balanced order-ten blue carrier and induced native-C multigraph.
Only the frozen almost-regular atlas is input. No full black carrier is visited.
"""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import product,permutations
from math import factorial,prod
import argparse,csv,hashlib,importlib.util,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
CHECKS=Counter()
N=10
PAIRS=tuple((i,j) for i in range(N) for j in range(i+2,N))
FULL=(1<<len(PAIRS))-1

def need(ok,message):
    CHECKS[message]+=1
    if not ok:raise RuntimeError(message)
def read(name):return json.loads((HERE/name).read_text())
def sha(name):return hashlib.sha256((HERE/name).read_bytes()).hexdigest()
def reverse(T):return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))
def normalize(T,order):return tuple(sum(((T[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
def invert(order):
    p=[0]*len(order)
    for i,v in enumerate(order):p[v]=i
    return tuple(p)
def stable(T):
    colors=tuple(r.bit_count() for r in T)
    while True:
        masks={c:sum(1<<v for v,d in enumerate(colors) if c==d) for c in sorted(set(colors))}
        signatures=tuple((colors[v],)+tuple((T[v]&m).bit_count() for m in masks.values()) for v in range(len(T)))
        names={s:i for i,s in enumerate(sorted(set(signatures)))}
        nxt=tuple(names[s] for s in signatures)
        if len(set(nxt))==len(set(colors)):
            return tuple(tuple(v for v in range(len(T)) if nxt[v]==c) for c in range(max(nxt)+1))
        colors=nxt

def canonical(T):
    cells=stable(T);cycles=[]
    for cell in cells:
        k=len(cell);need(k in (1,3,5),'actual stable cells size one three five')
        mask=sum(1<<v for v in cell)
        need(all((T[v]&mask).bit_count()==(k-1)//2 for v in cell),'stable cells internally regular')
        if k==1:cycles.append(cell);continue
        succ={}
        for v in cell:
            candidates=[w for w in cell if (T[v]>>w)&1 and (k==3 or (T[v]&T[w]&mask).bit_count()==1)]
            need(len(candidates)==1,'intrinsic successor unique');succ[v]=candidates[0]
        cycle=[];v=min(cell)
        for _ in cell:cycle.append(v);v=succ[v]
        need(v==cycle[0] and set(cycle)==set(cell),'intrinsic successor one cycle')
        need(all(((T[cycle[i]]>>cycle[j])&1)==int(0<(j-i)%k<=(k-1)//2) for i in range(k) for j in range(k)),'intrinsic cyclic arcs exact')
        cycles.append(tuple(cycle))
    rotations=prod(map(len,cycles));need(rotations<=25,'full-matrix rotation bound')
    best=None;orders=[]
    for shifts in product(*(range(len(c)) for c in cycles)):
        order=tuple(v for c,k in zip(cycles,shifts) for v in c[k:]+c[:k]);value=normalize(T,order)
        if best is None or value<best:best=value;orders=[]
        if value==best:orders.append(order)
    return best,tuple(orders),cells

def fullgroup_error(T):
    key,orders,cells=canonical(T);need(key==T,'canonical class fixed')
    group=tuple(sorted(invert(order) for order in orders));sizes={v:len(c) for c in cells for v in c}
    unseen=set(range(N));orbits=[]
    while unseen:
        v=min(unseen);orb={g[v] for g in group};need(orb<=unseen,'complete group orbit partition');unseen-=orb;orbits.append(tuple(sorted(orb)))
    error=sum(sizes[o[0]]-len(o) for o in orbits)
    return group,error,tuple(orbits),cells

def decode(word):
    rows=[1<<(i+1) if i<N-1 else 0 for i in range(N)]
    for k,(i,j) in enumerate(PAIRS):
        a,b=(i,j) if word>>k&1 else (j,i);rows[a]|=1<<b
    return tuple(rows)
def encode(T,path):return sum(((T[path[i]]>>path[j])&1)<<k for k,(i,j) in enumerate(PAIRS))
def isstrong(T):
    def reach(A):
        seen=1
        while True:
            nxt=seen
            for i in range(N):
                if seen>>i&1:nxt|=A[i]
            if nxt==seen:return seen==(1<<N)-1
            seen=nxt
    return reach(T) and reach(reverse(T))

def frozen_inputs():
    pins=read('SOURCE_PINS.json')['files']
    for pin in pins:need(sha(pin['path'])==pin['sha256'],'source SHA256')
    for tag in sorted({p['path'].split('/')[1].split('-MANIFEST')[0] for p in pins if p['source_artifact']=='MANIFEST.json'}):
        source=read('sources/'+tag+'-MANIFEST.json');byname={r['path']:r['sha256'] for r in source['files']}
        for p in pins:
            if p['path'].startswith('sources/'+tag+'-') and p['source_artifact']!='MANIFEST.json':need(p['sha256']==byname[p['source_artifact']],'input payload linked to frozen manifest')
    need(sha('sources/regular-eleven-MANIFEST.json')=='0462baf160610115c8b2436c2c57cce1964d0439ce49e86b9b43ce00f2bda3f2','exact inherited atlas freeze')
    need(sha('sources/balanced-metagraph-MANIFEST.json')=='7431d29c032d92aa6da1807c436bb839ed9c84909872de31d7f64b09d4963d43','exact inherited scalar freeze')

def atlas_and_selfconverse():
    tokens=iter((HERE/'sources/regular-eleven-eleven-O2-parents.txt').read_text().split());count=int(next(tokens));atlas=[]
    for _ in range(count):
        T=tuple(int(next(tokens)) for _ in range(N));order=int(next(tokens));group=tuple(tuple(int(next(tokens)) for _ in range(N)) for _ in range(order));atlas.append((T,group))
    need(next(tokens,None) is None,'atlas consumed exactly')
    index={T:i for i,(T,g) in enumerate(atlas)};need(len(index)==len(atlas),'distinct input class rows')
    sc={};converse=[];patterns=Counter();labelmass=0
    for i,(T,g) in enumerate(atlas):
        need(all(0<=r<(1<<N) and not (r>>v&1) for v,r in enumerate(T)),'loopless input rows')
        need(all(((T[a]>>b)&1)+((T[b]>>a)&1)==1 for a in range(N) for b in range(a)),'input tournament arcs')
        need(Counter(r.bit_count() for r in T)=={4:5,5:5},'input balanced scores')
        key,orders,cells=canonical(T);need(key==T,'input already exact canonical class')
        exactgroup=tuple(sorted(invert(p) for p in orders));need(exactgroup==tuple(sorted(g)),'input full automorphism group independently recovered')
        labelmass+=factorial(N)//len(g)
        revkey,revorders,_=canonical(reverse(T));need(revkey in index,'converse remains in atlas');mate=index[revkey];converse.append(mate)
        if mate==i:
            antis=tuple(invert(p) for p in revorders)
            need(all(all(((T[a]>>b)&1)!=((T[q[a]]>>q[b])&1) for a in range(N) for b in range(a)) for q in antis),'anti maps reverse every arc')
            involutions=tuple(q for q in antis if all(q[q[v]]==v for v in range(N)))
            need(involutions and all(all(q[v]!=v for v in range(N)) for q in involutions),'SC has fixed-free anti-involution')
            need(isstrong(T),'all SC balanced source classes strong')
            exactg,error,orbits,cells=fullgroup_error(T)
            sc[i]={'rows':T,'group':exactg,'anti_involutions':involutions,'error':error,'vertex_orbits':orbits,'stable_cells':cells}
            patterns[tuple(sorted(map(len,cells)))]+=1
    need(all(converse[converse[i]]==i for i in range(len(atlas))),'independent converse involution')
    receipt=read('sources/regular-eleven-eleven-O2.json')
    need(len(atlas)==receipt['almost_regular_parent_classes'] and len(sc)==receipt['self_converse_parents'],'independent atlas and SC counts match frozen receipt')
    need(labelmass==receipt['labelled_almost_regular_mass'],'atlas full labelled mass')
    states=list(csv.DictReader((HERE/'sources/regular-eleven-eleven-O2-states.tsv').open(),delimiter='\t'))
    comps=list(csv.DictReader((HERE/'sources/regular-eleven-eleven-O2-components.tsv').open(),delimiter='\t'))
    need(len(states)==len(atlas),'one frozen root state per atlas class')
    need(all(int(s['state'])==i and int(s['converse_state'])==converse[i] for i,s in enumerate(states)),'independent converse matches all frozen states')
    rootmap={i:int(states[i]['component']) for i in sc};selfcomponents={int(c['component']) for c in comps if c['component']==c['converse_component']}
    need(len(set(rootmap.values()))==len(sc) and set(rootmap.values())==selfcomponents,'unique SC root class bijection onto SC regular eleven components')
    for i,cid in rootmap.items():
        members=tuple(map(int,comps[cid]['states'].split(',')))
        need(i in members and set(members)&set(sc)=={i},'each SC regular component has precisely one SC balanced root orbit')
        sc[i]['regular_eleven_component']=cid
    return atlas,sc,patterns,labelmass

def gray_blue_carrier():
    pos={pair:k for k,pair in enumerate(PAIRS)};unused=set(range(len(PAIRS)));orbits=[]
    while unused:
        k=min(unused);i,j=PAIRS[k];orb=tuple(sorted({k,pos[N-1-j,N-1-i]}));unused-=set(orb);orbits.append(orb)
    masks=[sum(1<<k for k in orb) for orb in orbits]
    delta=[]
    for orb in orbits:
        d=[0]*N
        for k in orb:i,j=PAIRS[k];d[i]+=1;d[j]-=1
        delta.append(tuple((v,x) for v,x in enumerate(d) if x))
    deg=[r.bit_count() for r in decode(0)];word=0;gray=0;answer=set()
    for step in range(1<<len(orbits)):
        if step:
            bit=(step&-step).bit_length()-1;sign=1 if not (gray>>bit&1) else -1;gray^=1<<bit;word^=masks[bit]
            for v,x in delta[bit]:deg[v]+=sign*x
        if all(d in (4,5) for d in deg):answer.add(word)
    return answer,len(orbits)

def graph_statistics(vertices,adjacency):
    neighbors={v:set() for v in vertices}
    for a,b in adjacency:
        if a!=b:neighbors[a].add(b);neighbors[b].add(a)
    def components(omit=None,omitpair=None):
        unseen=set(vertices)-({omit} if omit is not None else set());parts=[]
        while unseen:
            start=min(unseen);seen={start};todo=[start]
            while todo:
                v=todo.pop()
                for w in neighbors[v]-seen:
                    if w==omit or (omitpair is not None and tuple(sorted((v,w)))==omitpair):continue
                    seen.add(w);todo.append(w)
            unseen-=seen;parts.append(sorted(seen))
        return parts
    parts=components();baseline=len(parts);clock=0;discovery={};low={};bridges=[];cuts=set()
    def visit(v,parent=None):
        nonlocal clock
        clock+=1;discovery[v]=low[v]=clock;children=0
        for w in sorted(neighbors[v]):
            if w==parent:continue
            if w not in discovery:
                children+=1;visit(w,v);low[v]=min(low[v],low[w])
                if low[w]>discovery[v] and adjacency[tuple(sorted((v,w)))]==1:bridges.append(tuple(sorted((v,w))))
                if parent is not None and low[w]>=discovery[v]:cuts.add(v)
            else:low[v]=min(low[v],discovery[w])
        if parent is None and children>1:cuts.add(v)
    for v in sorted(vertices):
        if v not in discovery:visit(v)
    # Independent exhaustive deletion controls; pair deletion is sufficient
    # for bridge testing precisely when the pair has multiplicity one.
    brute_cuts={v for v in vertices if len(components(omit=v))>baseline}
    need(cuts==brute_cuts,'Tarjan articulations match exhaustive vertex deletion')
    brute_bridges={pair for pair,m in adjacency.items() if pair[0]!=pair[1] and m==1 and len(components(omitpair=pair))>baseline}
    need(set(bridges)==brute_bridges,'Tarjan multigraph bridges match exhaustive single-edge deletion')
    return {'connected_component_count':baseline,'component_sizes':sorted(map(len,parts),reverse=True),'components':parts,'bridge_pairs':sorted(bridges),'articulation_nodes':sorted(cuts),'simple_neighbor_degree_histogram':sorted(Counter(len(neighbors[v]) for v in vertices).items()),'cycle_rank':sum(adjacency.values())-len(vertices)+baseline}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'normal.json'));ap.add_argument('--cost-only',action='store_true');args=ap.parse_args()
    frozen_inputs();atlas,sc,patterns,labelmass=atlas_and_selfconverse()
    anti_count=sum(len(r['anti_involutions']) for r in sc.values());trials=anti_count*2**(N//2)*factorial(N//2)
    cost={'atlas_classes':len(atlas),'SC_balanced_classes':len(sc),'anti_involutions':anti_count,'halfpath_candidate_ceiling':trials,'independent_blue_cube_candidates':1<<((N-1)**2//4),'full_black_carrier_enumerated':False}
    need(trials<=10000000,'halfpath cost ceiling ten million')
    if args.cost_only:
        Path(args.output).write_text(json.dumps(cost,sort_keys=True,indent=2)+'\n');print(json.dumps(cost,sort_keys=True));return
    owner={};paths={};node_counts={};path_counts={}
    for idx,rec in sorted(sc.items()):
        T=rec['rows'];multiplicity=Counter();witnesses={}
        for q in rec['anti_involutions']:
            pairs=tuple((v,q[v]) for v in range(N) if v<q[v])
            for chosen in product(*pairs):
                for half in permutations(chosen):
                    path=half+tuple(q[v] for v in half[::-1])
                    if not all(T[a]>>b&1 for a,b in zip(path,path[1:])):continue
                    word=encode(T,path);multiplicity[word]+=1;witnesses.setdefault(word,path)
        need(multiplicity and set(multiplicity.values())=={len(rec['group'])},'blue paths quotient freely by complete automorphism group')
        for word in multiplicity:
            need(word not in owner,'ordinary source classes have disjoint blue words');owner[word]=idx;paths[word]=witnesses[word]
            need(decode(word)==normalize(T,paths[word]),'every blue word has actual source path witness')
        node_counts[idx]=len(multiplicity);path_counts[idx]=sum(multiplicity.values())
    cube,blue_dimension=gray_blue_carrier();need(cube==set(owner),'independent full blue reflection cube equals anti-involution carrier exactly')
    print(json.dumps({'stage':'carrier_complete',**cost,'blue_ports':len(owner)},sort_keys=True),flush=True)
    internal=Counter();boundary=Counter();adjacency=Counter();endpoint=Counter();boundary_keys={};boundary_classes={};boundary_word_target={}
    for word,idx in sorted(owner.items()):
        T=decode(word);target=FULL^word;U=decode(target);ds=[r.bit_count() for r in T];dt=[r.bit_count() for r in U]
        need(target!=word and (FULL^target)==word,'C free involution')
        flipped=list(T)
        for i,j in PAIRS:flipped[i]^=1<<j;flipped[j]^=1<<i
        need(tuple(flipped)==U,'C equals literal nonpath arc reversal')
        need(all(((T[i]>>j)&1)==((T[N-1-j]>>(N-1-i))&1) for i in range(N) for j in range(N)),'literal word is blue')
        need(all(dt[v]==N-1-ds[v]+int(v==0)-int(v==N-1) for v in range(N)),'endpoint degree transport exact')
        isinternal=target in owner;need(isinternal==(ds[0]==5 and ds[-1]==4),'balanced partner iff HL')
        endpoint[ds[0],ds[-1]]+=1
        if isinternal:
            other=owner[target];internal[idx]+=1
            need(U==normalize(sc[other]['rows'],paths[target]),'internal C partner has exact target path witness')
            if word<target:adjacency[tuple(sorted((idx,other)))]+=1
        else:
            boundary[idx]+=1
            need((ds[0],ds[-1])==(4,5),'all blue boundary is LH')
            need(sorted(dt)==[3]+[4]*4+[5]*4+[6],'every boundary target has exact forced score profile')
            key,orders,cells=canonical(U);boundary_word_target[word]=(key,orders[0])
            if key not in boundary_classes:
                group,error,orbits,cells=fullgroup_error(key)
                need(error==0,'every boundary target directly zero error')
                need(canonical(reverse(key))[0]==key,'every boundary target class self converse')
                boundary_classes[key]={'canonical_rows':key,'automorphism_order':len(group),'stable_cells':cells,'error':error,'vertex_orbits':orbits}
    edgecount=sum(adjacency.values());need(sum(internal.values())==2*edgecount,'each internal C line exactly two ports including loops')
    need(len(owner)==sum(internal.values())+sum(boundary.values()),'source port conservation')
    degree=Counter()
    for (a,b),m in adjacency.items():degree[a]+=m;degree[b]+=m
    need(all(degree[i]==internal[i] and internal[i]+boundary[i]==node_counts[i] for i in sc),'every node degree and port conservation')
    # Independent scalar polynomial control, not a generation constraint.
    spec=importlib.util.spec_from_file_location('frozen_blue_engine',HERE/'sources/balanced-metagraph-engine.py');engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    scalar,peak=engine.blue_half_polynomial(N//2)
    need((len(owner),sum(internal.values()),sum(boundary.values()))==(scalar['balanced_blue_ports'],scalar['internal_blue_ports'],scalar['boundary_blue_ports']),'all three generated port totals match independent scalar polynomial')
    inherited=next(r for r in read('sources/balanced-metagraph-normal.json')['counts'] if r['n']==N)
    need(all(scalar[k]==inherited[k] for k in scalar),'scalar rerun matches frozen n10 receipt')
    # Embed the former positive-error induced graph by exact matrix keys.
    old=read('sources/defect-wiring-normal.json');oldrecs=read('sources/converse-errors-normal.json')['representatives'];oldmap={};excluded=[]
    index={T:i for i,(T,g) in enumerate(atlas)}
    for r in oldrecs:
        if not r['self_converse']:continue
        key=canonical(tuple(r['rows']))[0]
        if key in index:
            idx=index[key];need(idx in sc and sc[idx]['error']==40,'old balanced defect node exact class match');oldmap[r['index']]=idx
        else:excluded.append(r['index']);need(Counter(x.bit_count() for x in r['rows'])=={3:5,6:5},'excluded old defect nodes are nonbalanced degree-one family')
    positive={i for i,r in sc.items() if r['error']>0};need(set(oldmap.values())==positive,'all balanced positive-error nodes exactly old embedded controls')
    control_adj=Counter({tuple(sorted((oldmap[a],oldmap[b]))):m for a,b,m in old['adjacency_multiplicities']})
    need(control_adj==Counter({pair:m for pair,m in adjacency.items() if set(pair)<=positive}),'entire induced positive subgraph matches frozen edge multiplicities')
    for e in old['internal_edges']:
        need(owner[e['word']]==oldmap[e['source_index']] and owner[e['complement_word']]==oldmap[e['target_index']],'every old internal line exact word and class match')
    # Global relabel and converse controls on all output classes.
    for T in [r['rows'] for r in sc.values()]+list(boundary_classes):
        for order in [tuple(reversed(range(N))),tuple((7*v+3)%N for v in range(N))]:need(canonical(normalize(T,order))[0]==T,'output normal form invariant under nontrivial relabeling')
    keys=sorted(boundary_classes);keyids={key:i for i,key in enumerate(keys)};badj=Counter();targetcounts=Counter()
    ports=[]
    for word,idx in sorted(owner.items()):
        target=FULL^word
        if target in owner:kind='internal';targetid=owner[target];order=paths[target]
        else:
            kind='boundary';key,order=boundary_word_target[word];targetid=keyids[key];badj[idx,targetid]+=1;targetcounts[targetid]+=1
        ports.append([word,idx,list(paths[word]),target,kind,targetid,list(order)])
    for i,key in enumerate(keys):boundary_classes[key]['boundary_class_id']=i;boundary_classes[key]['edges_from_balanced_layer']=targetcounts[i]
    graph=graph_statistics(sc,adjacency)
    nodes=[{'atlas_index':i,'canonical_rows':r['rows'],'automorphism_order':len(r['group']),'anti_involutions':r['anti_involutions'],'stable_cells':r['stable_cells'],'vertex_orbits':r['vertex_orbits'],'refinement_error':r['error'],'regular_eleven_component':r['regular_eleven_component'],'blue_hamiltonian_paths':path_counts[i],'blue_presentation_words':node_counts[i],'internal_degree':internal[i],'boundary_edges':boundary[i]} for i,r in sorted(sc.items())]
    result={'status':'PASS','scope':'complete induced native-C blue multigraph on all self-converse almost-regular order-ten classes, plus exact boundary target classes; no black graph','cost_gate':cost,'blue_cube_dimension':blue_dimension,'atlas_labelled_mass':labelmass,'SC_balanced_nodes':len(sc),'blue_ports':len(owner),'internal_blue_ports':sum(internal.values()),'internal_blue_edges':edgecount,'internal_blue_loops':sum(m for (a,b),m in adjacency.items() if a==b),'internal_blue_nonloops':sum(m for (a,b),m in adjacency.items() if a!=b),'boundary_blue_edges':sum(boundary.values()),'boundary_classes':len(keys),'boundary_score_profile':[3]+[4]*4+[5]*4+[6],'scalar_polynomial':scalar,'scalar_polynomial_peak_states':peak,'endpoint_histogram':[[list(k),v] for k,v in sorted(endpoint.items())],'SC_stable_pattern_histogram':[[list(k),v] for k,v in sorted(patterns.items())],'SC_automorphism_histogram':sorted(Counter(len(r['group']) for r in sc.values()).items()),'SC_error_histogram':sorted(Counter(r['error'] for r in sc.values()).items()),'induced_degree_histogram':sorted(Counter(internal[i] for i in sc).items()),'blue_presentation_histogram':sorted(Counter(node_counts.values()).items()),'boundary_degree_histogram':sorted(Counter(boundary[i] for i in sc).items()),'adjacency_multiplicity_histogram':sorted(Counter(adjacency.values()).items()),'boundary_target_automorphism_histogram':sorted(Counter(r['automorphism_order'] for r in boundary_classes.values()).items()),'boundary_target_stable_pattern_histogram':[[list(k),v] for k,v in sorted(Counter(tuple(sorted(map(len,r['stable_cells']))) for r in boundary_classes.values()).items())],'known_positive_control':{'old_to_atlas_index':sorted(oldmap.items()),'excluded_nonbalanced_old_indices':sorted(excluded),'embedded_nodes':len(oldmap),'induced_edges':sum(control_adj.values())},**graph,'adjacency_multiplicities':[[a,b,m] for (a,b),m in sorted(adjacency.items())],'boundary_adjacency':[[a,b,m] for (a,b),m in sorted(badj.items())],'nodes':nodes,'boundary_target_classes':[boundary_classes[k] for k in keys],'port_record_schema':['word','source_atlas_index','source_path','C_word','target_kind','target_class_index','target_path_or_canonical_order'],'ports':ports,'checks':dict(CHECKS)}
    Path(args.output).write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    omit={'nodes','boundary_target_classes','ports','checks','adjacency_multiplicities','boundary_adjacency','components'}
    print(json.dumps({k:v for k,v in result.items() if k not in omit},sort_keys=True,indent=2))
if __name__=='__main__':main()
