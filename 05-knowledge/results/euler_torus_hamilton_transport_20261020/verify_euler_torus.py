#!/usr/bin/env python3
"""Exact, stdlib-only controls for the marked toroidal K7 Euler/Hamilton packet.
Run from any directory. All checks remain active under python -O.
"""
from collections import Counter, defaultdict, deque
from itertools import combinations, permutations, product
from math import factorial, gcd
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
N = 7
STEPS = {1:(1,0), 2:(-1,1), 3:(0,1), 4:(0,-1), 5:(1,-1), 6:(-1,0)}
CAROUSEL = (14,28,56,112,97,67,7)
PALEY = (22,44,88,49,98,69,11)
MIDDLE = (26,28,41,112,100,67,7)
FACES = sorted({tuple(sorted((i,(i+1)%7,(i+3)%7))) for i in range(7)} |
               {tuple(sorted(((i+1)%7,(i+3)%7,(i+4)%7))) for i in range(7)})

def require(test, message):
    if not test: raise RuntimeError(message)

def arcs(out):
    return [(i,j) for i in range(len(out)) for j in range(len(out)) if out[i]>>j&1]

def regular_orientations(n=7):
    """Degree-pruned row generation, not a 2^21 mask sweep or class import."""
    require(n%2==1, 'odd n required')
    r=(n-1)//2; out=[0]*n; known_degree=[0]*n
    def visit(i):
        if i==n:
            yield tuple(out); return
        future=list(range(i+1,n)); need=r-known_degree[i]
        if not 0<=need<=len(future): return
        for wins in combinations(future,need):
            wins=set(wins); old=out[i]
            for j in future:
                if j in wins: out[i]|=1<<j
                else: out[j]|=1<<i; known_degree[j]+=1
            unknown=len(future)-1
            if all(known_degree[j]<=r<=known_degree[j]+unknown for j in future):
                yield from visit(i+1)
            for j in future:
                if j not in wins: out[j]^=1<<i; known_degree[j]-=1
            out[i]=old
    yield from visit(0)

def hp_dp(out):
    n=len(out); dp=[{} for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v,count in dp[mask].items():
            candidates=out[v]&~mask
            while candidates:
                bit=candidates&-candidates; candidates-=bit; w=bit.bit_length()-1
                d=dp[mask|bit]; d[w]=d.get(w,0)+count
    return sum(dp[-1].values())

def edge_homology(i,j):
    a,b=STEPS[(j-i)%7]
    # Lift vertex i to (i,0).  Endpoint discrepancy is (i+a-j,b).
    require((i+a-j+3*b)%7==0,'edge lift in kernel')
    return ((i+a-j+3*b)//7,b)

def flux(out):
    labels=[edge_homology(i,j) for i,j in arcs(out)]
    return tuple(sum(z[k] for z in labels) for k in (0,1))

def cycle_hist_permutations(out):
    hist=Counter(); examples={}; cycles=[]
    for rest in permutations(range(1,7)):
        c=(0,)+rest
        edges=list(zip(c,c[1:]+c[:1]))
        if not all(out[i]>>j&1 for i,j in edges): continue
        h=tuple(sum(edge_homology(i,j)[k] for i,j in edges) for k in (0,1))
        require(gcd(*h)==1,'simple cycle homology must be primitive here')
        hist[h]+=1; examples.setdefault(h,c); cycles.append(c)
    return hist,examples,cycles

def cycle_hist_dp(out):
    """Independent DP with actual lattice-step displacement, not edge labels."""
    dp=[Counter() for _ in range(1<<7)]; dp[1][(0,0,0)]=1
    for mask in range(1<<7):
        for (last,a,b),count in dp[mask].items():
            for nxt in range(1,7):
                if mask>>nxt&1 or not out[last]>>nxt&1: continue
                da,db=STEPS[(nxt-last)%7]
                dp[mask|(1<<nxt)][(nxt,a+da,b+db)]+=count
    hist=Counter()
    for (last,a,b),count in dp[-1].items():
        if not out[last]&1: continue
        da,db=STEPS[-last%7]; a+=da; b+=db
        require((a+3*b)%7==0,'closed displacement in lattice')
        hist[((a+3*b)//7,b)]+=count
    return hist

def det_bareiss(matrix):
    a=[row[:] for row in matrix]; n=len(a); sign=1; previous=1
    for k in range(n-1):
        if not a[k][k]:
            j=next(j for j in range(k+1,n) if a[j][k])
            a[k],a[j]=a[j],a[k]; sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                require(numerator%previous==0,'Bareiss exact division')
                a[i][j]=numerator//previous
        for i in range(k+1,n): a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]

def arborescence_determinant(out,root=0):
    vertices=[v for v in range(len(out)) if v!=root]
    return det_bareiss([[out[i].bit_count() if i==j else -int(out[i]>>j&1)
                          for j in vertices] for i in vertices])

def arborescence_direct(out,root=0):
    vertices=[v for v in range(len(out)) if v!=root]
    choices=[[w for w in range(len(out)) if out[v]>>w&1] for v in vertices]
    count=0
    for targets in product(*choices):
        parent=dict(zip(vertices,targets)); good=True
        for v in vertices:
            seen=set()
            while v!=root and v not in seen: seen.add(v); v=parent[v]
            if v!=root: good=False; break
        count+=good
    return count

def euler_transition_count(out):
    """Independent enumeration of local incoming/outgoing matchings (6^7)."""
    edges=arcs(out); incoming=[[] for _ in out]; outgoing=[[] for _ in out]
    for k,(i,j) in enumerate(edges): outgoing[i].append(k); incoming[j].append(k)
    options=[list(permutations(row)) for row in outgoing]
    count=0; m=len(edges)
    for local in product(*options):
        successor=[0]*m
        for v in range(len(out)):
            for e,f in zip(incoming[v],local[v]): successor[e]=f
        e=successor[0]; length=1
        while e!=0: length+=1; e=successor[e]
        count+=length==m
    return count

def line_graph(out):
    edges=arcs(out)
    return tuple(sum(1<<f for f,(x,y) in enumerate(edges) if j==x) for i,j in edges)

def flip(out,tri):
    changed=list(out)
    for i,j in combinations(tri,2): changed[i]^=1<<j; changed[j]^=1<<i
    return tuple(changed)

def directed_triangle(out,tri):
    return all(sum(out[i]>>j&1 for j in tri)==1 for i in tri)

def list_record(counter):
    return [{'slope':list(k),'count':v} for k,v in sorted(counter.items())]

def main():
    output={'scope':'All labelled regular orientations of the specified triangular-torus K7; no larger-order claim.',
            'lattice_basis':[[7,0],[-3,1]],'faces':[list(f) for f in FACES], 'models':{}}
    for name,out in [('carousel',CAROUSEL),('paley',PALEY),('middle',MIDDLE)]:
        hist,examples,cycles=cycle_hist_permutations(out)
        require(hist==cycle_hist_dp(out),name+' two cycle counters')
        hp=[p for p in permutations(range(7)) if all(out[p[i]]>>p[i+1]&1 for i in range(6))]
        H=hp_dp(out); require(H==len(hp),name+' two HP counters')
        nonclosing=[p for p in hp if not out[p[-1]]>>p[0]&1]
        require(H==7*sum(hist.values())+len(nonclosing),'closing decomposition')
        tau=arborescence_determinant(out)
        require(tau==arborescence_direct(out),'two arborescence counters')
        E=euler_transition_count(out); require(E==tau*128,'independent BEST check')
        # Gauge changes of vertex lifts preserve every closed-cycle winding.
        shifts=[(v*v-4*v,(-1)**v*v) for v in range(7)]
        for c in cycles:
            edges=list(zip(c,c[1:]+c[:1]))
            before=tuple(sum(edge_homology(i,j)[k] for i,j in edges) for k in (0,1))
            after=tuple(sum(edge_homology(i,j)[k]+shifts[i][k]-shifts[j][k] for i,j in edges) for k in (0,1))
            require(before==after,'vertex-lift gauge invariance')
        output['models'][name]={'out_masks':out,'H':H,'directed_Hamilton_cycles':len(cycles),
            'nonclosing_Hamilton_paths':len(nonclosing),'first_nonclosing_path':nonclosing[0],
            'arborescences':tau,'Euler_circuits_mod_rotation':E,'Euler_length':21,
            'Hamilton_paths_of_line_digraph':21*E,'Euler_flux':flux(out),
            'cycle_slope_counts':list_record(hist),
            'cycle_examples':[{'slope':list(k),'vertices':list(v)} for k,v in sorted(examples.items())]}
    require(flip(CAROUSEL,(0,2,4))==MIDDLE,'specified winding triangle')
    triangle=[(0,2),(2,4),(4,0)]
    h=tuple(sum(edge_homology(i,j)[k] for i,j in triangle) for k in (0,1))
    require(h==(1,3) and flux(MIDDLE)==(4,8),'triangle flux update')
    output['winding_triangle']={'vertices':[0,2,4],'homology':h,'before':[6,14],'after':[4,8]}
    # Directly audit Aut(Paley), then its action on nonclosing paths.
    autos=[p for p in permutations(range(7)) if all(bool(PALEY[i]>>j&1)==bool(PALEY[p[i]]>>p[j]&1) for i in range(7) for j in range(7))]
    affine={tuple((a*x+b)%7 for x in range(7)) for a in (1,2,4) for b in range(7)}
    require(set(autos)==affine,'Paley exact automorphism group')
    witness=tuple(output['models']['paley']['first_nonclosing_path'])
    orbit={tuple(p[v] for v in witness) for p in autos}
    require(len(orbit)==21 and all(not PALEY[p[-1]]>>p[0]&1 for p in orbit),'one nonclosing orbit')
    output['paley_automorphism_order']=21
    # Small independent direct line-graph Hamilton counts.
    output['small_line_graph_checks']=[]
    for n in (3,5):
        out=tuple(sum(1<<((i+d)%n) for d in range(1,(n+1)//2)) for i in range(n))
        E=euler_transition_count(out); Hline=hp_dp(line_graph(out)); m=n*(n-1)//2
        require(Hline==m*E,'direct line digraph identity')
        output['small_line_graph_checks'].append({'n':n,'Euler_circuits':E,'H_line':Hline})
    all_out=sorted(regular_orientations()); require(len(all_out)==2640 and len(set(all_out))==2640,'independent regular census')
    by_flux=defaultdict(list); Hvalues={}
    for out in all_out:
        require(all(row.bit_count()==3 for row in out),'regular degree')
        f=flux(out); require(f[0]%2==0 and f[1]%2==0,'checkerboard even flux')
        by_flux[f].append(out); Hvalues[out]=hp_dp(out)
    require(set(Hvalues.values())=={171,175,189},'raw H support')
    output['all_regular_H_counts']=dict(sorted(Counter(Hvalues.values()).items()))
    index={out:i for i,out in enumerate(all_out)}; trees=[]; sector_rows=[]
    for f,states in sorted(by_flux.items()):
        remaining=set(states); component_sizes=[]; degrees=Counter(); sector_edges=0
        while remaining:
            start=min(remaining); remaining.remove(start); queue=deque([start]); size=0
            trees.append([index[start],-1,-1])
            while queue:
                out=queue.popleft(); size+=1; degree=0
                for face_id,tri in enumerate(FACES):
                    if not directed_triangle(out,tri): continue
                    neighbour=flip(out,tri); degree+=1
                    require(flux(neighbour)==f,'face flip flux conservation')
                    require(neighbour in index,'face flip stays regular universe')
                    if neighbour in remaining:
                        remaining.remove(neighbour); queue.append(neighbour)
                        trees.append([index[neighbour],index[out],face_id])
                degrees[degree]+=1; sector_edges+=degree
            component_sizes.append(size)
        require(len(component_sizes)==1,'one component per flux in K7 finite universe')
        sector_rows.append({'flux':f,'size':len(states),'H_population':dict(sorted(Counter(Hvalues[t] for t in states).items())),
            'component_sizes':component_sizes,'edge_count':sector_edges//2,'degree_histogram':dict(sorted(degrees.items()))})
    output['flux_sectors']=sector_rows
    require(len(sector_rows)==25,'25 flux sectors')
    size_counts=Counter(row['size'] for row in sector_rows)
    require(size_counts==Counter({660:1,273:6,42:6,14:6,1:6}),'sector size multiset')
    output['sector_size_multiplicities']=dict(sorted(size_counts.items()))
    zero_examples={}
    for out in by_flux[(0,0)]: zero_examples.setdefault(Hvalues[out],out)
    output['zero_flux_H_counterexamples']=dict(sorted(zero_examples.items()))
    # Cheap hostile against identifying a 273-state sector with the natural matching-flip graph.
    matching_states={(a,b) for a in range(32) if not a&(a<<1) for b in range(64) if not b&(b<<1)}
    md=Counter(sum((a^(1<<i),b) in matching_states for i in range(5))+sum((a,b^(1<<i)) in matching_states for i in range(6)) for a,b in matching_states)
    matching_edges=sum(d*c for d,c in md.items())//2
    torus_row=next(row for row in sector_rows if tuple(row['flux'])==(0,2))
    output['273_matching_hostile']={'matching_states':len(matching_states),'matching_edges':matching_edges,
        'matching_degree_histogram':dict(sorted(md.items())),'torus_flux':[0,2],
        'torus_edges':torus_row['edge_count'],'torus_degree_histogram':torus_row['degree_histogram']}
    require(matching_edges!=torus_row['edge_count'],'natural 273-state graphs not isomorphic')
    # A spanning forest is an independently replayable connectivity certificate.
    certificate={'orientations':[list(t) for t in all_out],'forest_child_parent_face':trees,'faces':[list(f) for f in FACES]}
    (HERE/'face_flip_certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    (HERE/'euler_torus_controls.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'regular_orientations':len(all_out),'flux_sectors':len(sector_rows),
       'sector_sizes':output['sector_size_multiplicities'],'models':{k:{q:v[q] for q in ('H','directed_Hamilton_cycles','nonclosing_Hamilton_paths','arborescences','Euler_circuits_mod_rotation','Euler_flux')} for k,v in output['models'].items()},
       '273_hostile':output['273_matching_hostile']},indent=2,sort_keys=True))

if __name__=='__main__': main()
