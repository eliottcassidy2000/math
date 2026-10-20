#!/usr/bin/env python3
"""Independent finite audit. No imports from producer scripts; no producer writes.
Universe: all 2^21 labelled orientations, built by Gray-code sweep.
Run: python3 independent_audit.py [producer directory]
"""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter, defaultdict, deque
from math import gcd
import json, sys

HERE=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'session16-sequence-knots'
def check(test,msg):
    if not test: raise RuntimeError(msg)

edges=list(combinations(range(7),2)); ei={e:k for k,e in enumerate(edges)}
ALL=(1<<21)-1
steps={1:(1,0),2:(-1,1),3:(0,1),4:(0,-1),5:(1,-1),6:(-1,0)}
def edgebit(i,j): return 1<<ei[tuple(sorted((i,j)))]
def follows(t,i,j): return bool(t&edgebit(i,j))==(i<j)
def word_masks(word,closed=False):
    mask=req=0
    pairs=list(zip(word,word[1:]))
    if closed:pairs.append((word[-1],word[0]))
    for i,j in pairs:
        bit=edgebit(i,j);mask|=bit
        if i<j:req|=bit
    return mask,req
def accepts(t,mr):
    mask,req=mr
    return t&mask==req

def lattice_displacement(pairs):
    x=y=0
    for i,j in pairs:
        dx,dy=steps[(j-i)%7];x+=dx;y+=dy
    return x,y

def homology_closed(pairs):
    x,y=lattice_displacement(pairs)
    check((x+3*y)%7==0,'closed lattice displacement')
    return ((x+3*y)//7,y)

def flux(t):
    return homology_closed((i,j) if t>>k&1 else (j,i) for k,(i,j) in enumerate(edges))

def source_mask(rows): return sum((rows[i]>>j&1)<<k for k,(i,j) in enumerate(edges))

# Independent exhaustive Gray-code universe; all 2^21 orientations are visited.
deg=list(range(7));prev=0;states=[]
for k in range(1<<21):
    t=k^(k>>1)
    if k:
        bit=t^prev;j=bit.bit_length()-1;u,v=edges[j]
        change=1 if t&bit else -1
        deg[u]+=change;deg[v]-=change
    if deg==[3]*7:states.append(t)
    prev=t
states.sort();S=set(states)
C=json.loads((SOURCE/'face_flip_certificate.json').read_text())
controls=json.loads((SOURCE/'euler_torus_controls.json').read_text())
check(S=={source_mask(r) for r in C['orientations']},'independent exhaustive universe agrees with certificate')
check(len(S)==2640,'2640 regular orientations')

# Build oriented triangles from actual positive-area lattice triangles, not stored face order.
ccw=[(i,(i+1)%7,(i+3)%7) for i in range(7)]
ccw += [((i+1)%7,(i+4)%7,(i+3)%7) for i in range(7)]
check({tuple(sorted(f)) for f in ccw}=={tuple(f) for f in C['faces']},'face sets agree')
fm=[word_masks(f,True) for f in ccw]
inc=defaultdict(list)
for f,tri in enumerate(ccw):
    check(homology_closed(list(zip(tri,tri[1:]+tri[:1])))==(0,0),'each face is contractible in lattice')
    for i,j in zip(tri,tri[1:]+tri[:1]):
        inc[ei[tuple(sorted((i,j)))]].append((f,1 if i<j else -1))
check(len(inc)==21 and all(len(v)==2 and v[0][1]==-v[1][1] for v in inc.values()),'two distinct oppositely oriented faces per edge')
for v in range(7):
    link=defaultdict(list)
    for f in ccw:
        if v in f:
            p,q=[w for w in f if w!=v];link[p].append(q);link[q].append(p)
    check(len(link)==6 and all(len(x)==2 for x in link.values()),'vertex link degree two')
    seen={next(iter(link))};todo=list(seen)
    while todo:
        for w in link[todo.pop()]:
            if w not in seen:seen.add(w);todo.append(w)
    check(len(seen)==6,'vertex link one six-cycle')
check(7-21+14==0,'connected oriented closed surface Euler characteristic zero')
# Dual adjacency retained for exact integer face potentials.
dual=defaultdict(list)
for k,((f,s),(g,_)) in inc.items():
    dual[f].append((g,k,s));dual[g].append((f,k,-s))

def potential(t,u):
    # q=(J(t)-J(u))/2 is the Zhang-Qian oriented difference t-u.
    z={0:0};todo=[0]
    while todo:
        f=todo.pop()
        for g,k,s in dual[f]:
            value=z[f]-(((t>>k)&1)-((u>>k)&1))*s
            if g in z:check(z[g]==value,'integral potential consistency')
            else:z[g]=value;todo.append(g)
    check(len(z)==14,'dual connected')
    mn=min(z.values())
    return sum(v-mn for v in z.values())

paths=[(w,word_masks(w)) for w in permutations(range(7))]
cycles=[(w,word_masks(w,True),homology_closed(list(zip(w,w[1:]+w[:1])))) for p in permutations(range(1,7)) for w in [(0,)+p]]
H={t:sum(accepts(t,mr) for w,mr in paths) for t in states}
check(Counter(H.values())=={175:720,171:1680,189:240},'raw H populations independent all-permutation check')
sectors=defaultdict(list)
for t in states:
    f=flux(t);check(all(x%2==0 for x in f),'checkerboard even flux')
    sectors[f].append(t)
    for reverse in (False,True):
        reached={0};todo=[0]
        while todo:
            i=todo.pop()
            for j in range(7):
                if j!=i and follows(t,j,i) if reverse else j!=i and follows(t,i,j):
                    if j not in reached:reached.add(j);todo.append(j)
        check(len(reached)==7,'every regular orientation strongly connected')
check(Counter(map(len,sectors.values()))=={660:1,273:6,42:6,14:6,1:6},'flux size multiset')

# Both un-oriented move graph and actual ccw-only directed graph are rebuilt.
adj={};di={};computed_sector_rows=[]
for t in states:
    adj[t]=[];di[t]=[]
    for mask,req in fm:
        if t&mask in (req,mask^req):
            u=t^mask;check(u in S and flux(u)==flux(t),'all legal face flips preserve regularity and flux');adj[t].append(u)
        if t&mask==req:di[t].append(t^mask)
for f,ss in sorted(sectors.items()):
    row={'flux':list(f),'size':len(ss),'H_population':dict(sorted(Counter(H[t] for t in ss).items())),
         'edge_count':sum(len(adj[t]) for t in ss)//2,'degree_histogram':dict(sorted(Counter(len(adj[t]) for t in ss).items()))}
    old=next(r for r in controls['flux_sectors'] if r['flux']==list(f))
    for key in ('size','edge_count'):check(row[key]==old[key],'sector '+key)
    for key in ('H_population','degree_histogram'):check({str(k):v for k,v in row[key].items()}==old[key],'sector '+key)
    computed_sector_rows.append(row)

# Exhaustive source-target check of the cited directed distance formula in this finite universe.
checked_pairs=0;largest_distance=0;sector_directed_diameters={}
for f,ss in sorted(sectors.items()):
    diameter=0
    for t in ss:
        distance={t:0};todo=deque([t])
        while todo:
            v=todo.popleft()
            for u in di[v]:
                if u not in distance:distance[u]=distance[v]+1;todo.append(u)
        check(set(distance)==set(ss),'ccw flip SCC equals flux sector')
        for u in ss:
            check(distance[u]==potential(t,u),'Zhang-Qian directed distance equals normalized integral face potential')
            checked_pairs+=1
        diameter=max(diameter,max(distance.values()))
    sector_directed_diameters[str(f)]=diameter;largest_distance=max(largest_distance,diameter)

model_rows={};class_orbits=[]
for name,data in sorted(controls['models'].items()):
    t=source_mask(data['out_masks'])
    hist=Counter(h for w,mr,h in cycles if accepts(t,mr))
    check(hist=={tuple(r['slope']):r['count'] for r in data['cycle_slope_counts']},'model marked cycle polynomial')
    check(H[t]==data['H'] and flux(t)==tuple(data['Euler_flux']),'model counts and flux')
    check(all(gcd(*h)==1 for h in hist),'all model Hamilton cycles nonzero primitive')
    hp=[w for w,mr in paths if accepts(t,mr)]
    R=sum(not follows(t,w[-1],w[0]) for w in hp)
    check(H[t]==7*sum(hist.values())+R==data['H'],'H=nC+R')
    # All relabelled copies prove independent three-class coverage, rather than assuming it.
    orbit=set()
    directed=[(i,j) if t>>k&1 else (j,i) for k,(i,j) in enumerate(edges)]
    for perm in permutations(range(7)):
        relabel=0
        for i,j in directed:
            a,b=perm[i],perm[j]
            if a<b:relabel|=edgebit(a,b)
        orbit.add(relabel)
    class_orbits.append(orbit)
    # All legal same-endpoint paired complements for each model.
    pairs_checked=0;flux_changes=Counter()
    for p in hp:
        pm,pr=word_masks(p);mid=t^(ALL^pm)
        check(all(follows(mid,i,j) for i,j in zip(p,p[1:])),'first operation preserves P')
        middle=[v for v in range(7) if v not in (p[0],p[-1])]
        pdisp=lattice_displacement(zip(p,p[1:]))
        for tail in permutations(middle):
            q=(p[0],)+tail+(p[-1],);qm,qr=word_masks(q)
            if mid&qm!=qr:continue
            u=mid^(ALL^qm);check(u in S,'same-endpoint second complement is regular')
            # J(u)-J(t)=2(Q-P), checked edgewise, not only after homology.
            for k,(i,j) in enumerate(edges):
                pv=(1 if pr>>k&1 else -1) if pm>>k&1 else 0
                qv=(1 if qr>>k&1 else -1) if qm>>k&1 else 0
                check(((u>>k)&1)-((t>>k)&1)==qv-pv,'paired path edge-chain identity')
            qdisp=lattice_displacement(zip(q,q[1:]));dx,dy=qdisp[0]-pdisp[0],qdisp[1]-pdisp[1]
            check((dx+3*dy)%7==0,'Q-P closed integer cycle')
            dh=((dx+3*dy)//7,dy);ft,fu=flux(t),flux(u)
            check(fu==(ft[0]+2*dh[0],ft[1]+2*dh[1]),'paired path flux law')
            pairs_checked+=1;flux_changes[dh]+=1
    model_rows[name]={'H':H[t],'C':sum(hist.values()),'R':R,'flux':flux(t),'slope_counts':[{'slope':h,'count':c} for h,c in sorted(hist.items())],
        'labelled_isomorphism_orbit_size':len(orbit),'paired_path_presentations_checked':pairs_checked,'paired_return_half_flux_changes':[{'slope':h,'count':c} for h,c in sorted(flux_changes.items())]}
check(set.union(*class_orbits)==S and sum(map(len,class_orbits))==len(S),'three independently relabelled model orbits partition complete universe')

# Torus rotation and marked knot hostile, including all slopes, not only three axis cycles.
paley=source_mask(controls['models']['paley']['out_masks'])
def M(h):return (2*h[0]-h[1],7*h[0]-3*h[1])
for w,mr,h in cycles:
    if accepts(paley,mr):
        w2=tuple(2*x%7 for x in w);mr2=word_masks(w2,True)
        check(accepts(paley,mr2),'multiply by two preserves Paley direction')
        h2=homology_closed(zip(w2,w2[1:]+w2[:1]));check(h2==M(h) and M(M(M(h)))==h,'torus automorphism homology action')
check([M((1,0)),M((2,7)),M((-3,-7))]==[(2,7),(-3,-7),(1,0)],'unknot and two distinct torus knot slopes in one automorphism orbit')

# Natural matching graph independently generated as edge subsets of path components.
match=[]
for a in range(1<<11):
    selected=[k for k in range(11) if a>>k&1]
    if all(not (k+1 in selected and k!=4) for k in selected):match.append(a)
ms=set(match);mh=Counter(sum((a^(1<<k)) in ms for k in range(11)) for a in match)
check(len(ms)==273 and sum(k*v for k,v in mh.items())//2==914 and max(mh)==11,'independent matching graph 273/914/11')
check(sum(len(adj[t]) for t in sectors[(0,2)])//2==770,'torus graph 273/770')

# Distinguishable arcs, loops, parallel arcs, and a balance-failure hostile.
def direct_line_test(name,arcs):
    m=len(arcs)
    H=sum(all(arcs[w[i]][1]==arcs[w[i+1]][0] for i in range(m-1)) for w in permutations(range(m)))
    E=sum(all(arcs[w[i]][1]==arcs[w[(i+1)%m]][0] for i in range(m)) for tail in permutations(range(1,m)) for w in [(0,)+tail])
    return {'name':name,'m':m,'H_line':H,'E_cyclic':E}
line_tests=[direct_line_test('three distinguishable loops',[(0,0)]*3),direct_line_test('two parallel arcs each direction',[(0,1),(0,1),(1,0),(1,0)]),direct_line_test('balanced bouquet of two cycles',[(0,1),(1,0),(0,2),(2,0)]),direct_line_test('unbalanced directed path hostile',[(0,1),(1,2)])]
for row in line_tests[:-1]:check(row['H_line']==row['m']*row['E_cyclic'],'balanced general line-digraph identity')
check(line_tests[-1]['H_line']==1 and line_tests[-1]['E_cyclic']==0,'balance hypothesis cannot be dropped')

result={'status':'PASS','producer_imports':False,'orientation_masks_exhaustively_examined':1<<21,'regular_orientations':len(S),
        'embedding':{'vertices':7,'edges':21,'faces':14,'closed_orientable_genus':1,'face_incidence':'every edge in two distinct faces with opposite signs','vertex_links':'seven six-cycles'},
        'flux_sectors':len(sectors),'H_populations':dict(Counter(H.values())),'sector_size_multiplicities':dict(Counter(map(len,sectors.values()))),
        'ccw_directed_distance_pairs_checked':checked_pairs,'maximum_ccw_directed_distance':largest_distance,
        'sector_ccw_directed_diameters':sector_directed_diameters,'models':model_rows,'line_graph_boundary_checks':line_tests,
        'matching_hostile':{'vertices':len(ms),'edges':914,'maximum_degree':max(mh),'torus_273_vertices_edges':770},
        'finite_scope':'Specified triangular torus K7, three model path-return checks, and named small line-digraph boundary examples only.'}
(HERE/'independent_results.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('models','sector_ccw_directed_diameters')},sort_keys=True,indent=2))
