#!/usr/bin/env python3
"""Session130 independent n5 colored-multigraph audit.
No producer imports or producer data reads. Build complete labelled tournament
orbits first, then locate the fixed-path cube in that independent universe.
All gates use explicit exceptions, including negative mutation controls.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd
from functools import reduce
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def rejected(fn):
    try:
        fn()
    except RuntimeError:
        return True
    raise RuntimeError('Hostile control was not rejected')


class Universe:
    def __init__(self,n):
        self.n=n
        self.pairs=list(combinations(range(n),2))
        self.weights={e:1 << (len(self.pairs)-1-k) for k,e in enumerate(self.pairs)}
        self.full=(1 << len(self.pairs))-1
        self.perms=list(permutations(range(n)))
        self.classof={}
        self.orbits={}
        for code in range(self.full+1):
            if code not in self.classof:
                orbit={self.rename(code,p) for p in self.perms}
                canon=min(orbit)
                need(canon==code,'Orbit traversal skipped a smaller representative')
                self.orbits[canon]=orbit
                for c in orbit:
                    need(c not in self.classof,'Different isomorphism orbits overlap')
                    self.classof[c]=canon
        self.ordinary=sorted(self.orbits)
        self.ordinary_ids={c:i for i,c in enumerate(self.ordinary)}
        self.merged=sorted({min(c,self.classof[c^self.full]) for c in self.ordinary})
        self.node_ids={c:i for i,c in enumerate(self.merged)}
        self.free=[(n-x,n-y) for y in range(1,n-1) for x in range(n,y+1,-1)]
        self.bitfull=(1<<len(self.free))-1
        self.graphs={m:self.fixed_graph(m) for m in range(self.bitfull+1)}
        self.node={m:self.owner(g) for m,g in self.graphs.items()}
        self.fibers=defaultdict(list)
        for m,u in self.node.items(): self.fibers[u].append(m)
        need(set(self.classof[g] for g in self.graphs.values())==set(self.ordinary),
             'Fixed-path atlas missed a labelled tournament class')

    def edge(self,g,i,j):
        need(i!=j,'Loop used as tournament arc')
        return bool(g & self.weights[(i,j)]) if i<j else not bool(g & self.weights[(j,i)])

    def rename(self,g,p):
        # Pull back the full labelled relation through p; no fixed-path producer used.
        return sum(w for (i,j),w in self.weights.items() if self.edge(g,p[i],p[j]))

    def fixed_graph(self,m):
        return self.full ^ sum(self.weights[e] for k,e in enumerate(self.free) if m>>k&1)

    def owner(self,g):
        return self.node_ids[min(self.classof[g],self.classof[g^self.full])]

    def reflect(self,m):
        r=self.rename(self.graphs[m]^self.full,tuple(reversed(range(self.n))))
        matches=[t for t,g in self.graphs.items() if g==r]
        need(len(matches)==1,'Converse plus vertex reversal left fixed-path atlas')
        return matches[0]

    def grid_reflect(self,m):
        f={e:k for k,e in enumerate(self.free)}
        return sum(1<<f[(self.n-1-j,self.n-1-i)] for k,(i,j) in enumerate(self.free) if m>>k&1)

    def induced(self,g,vertices,lower):
        return sum(w for (i,j),w in lower.weights.items() if self.edge(g,vertices[i],vertices[j]))

    def census(self,lower=None):
        rows=[]
        for u,canon in enumerate(self.merged):
            cls=sorted({self.classof[self.graphs[m]] for m in self.fibers[u]})
            paths=[p for p in self.perms if all(self.edge(canon,p[i],p[i+1]) for i in range(self.n-1))]
            aut=factorial(self.n)//len(self.orbits[canon])
            syms=[p for p in paths if all(self.edge(canon,p[i],p[j]) != self.edge(canon,p[-1-i],p[-1-j]) for i,j in self.pairs)]
            blue=[m for m in self.fibers[u] if self.reflect(m)==m]
            need(len(self.fibers[u])*aut==len(cls)*len(paths),'Hamiltonian/aut fiber law')
            need(len(blue)*aut==len(syms),'Blue anti-palindromic-path/aut law')
            row={'node':u,'canon':canon,'classes':[self.ordinary_ids[c] for c in cls],
                 'masks':self.fibers[u],'H':len(paths),'Aut':aut,'symmetric_H':len(syms),
                 'degree':len(self.fibers[u]),'blue_degree':len(blue),'black_degree':len(self.fibers[u])-len(blue),
                 'type':'pure_blue' if len(blue)==len(self.fibers[u]) else 'pure_black' if not blue else 'mixed'}
            if lower:
                def profile(blue_filter,verts):
                    counts=Counter(lower.owner(self.induced(self.graphs[m],verts,lower)) for m in self.fibers[u]
                       if blue_filter is None or ((m in blue)==blue_filter))
                    return [counts[v] for v in range(len(lower.merged))]
                raw=profile(None,list(range(self.n-1)))
                need(raw==profile(None,list(range(1,self.n))),'Endpoint deletion mismatch')
                bd=profile(True,list(range(self.n-1)))
                kd=profile(False,list(range(self.n-1)))
                need(all(a==b+c for a,b,c in zip(raw,bd,kd)),'Colored deletion does not sum')
                d=reduce(gcd,raw)
                parent_pairs=Counter()
                for m in self.fibers[u]:
                    r=self.reflect(m)
                    if m<r:
                        endpoints=tuple(sorted(lower.owner(self.induced(self.graphs[x],list(range(self.n-1)),lower)) for x in (m,r)))
                        parent_pairs[endpoints]+=1
                pair_degrees=Counter()
                for endpoints,multiplicity in parent_pairs.items():
                    for x in endpoints:pair_degrees[x]+=multiplicity
                need([pair_degrees[k] for k in range(len(lower.merged))]==kd,'Black parent-pair incidence mismatch')
                row.update(deletion=raw,primitive=[v//d for v in raw],blue_deletion=bd,black_deletion=kd,
                   black_parent_pairs=[{'parents':list(e),'multiplicity':v} for e,v in sorted(parent_pairs.items())])
            rows.append(row)
        lines=[]
        for m in range((self.bitfull+1)//2):
            c=m^self.bitfull; r=self.reflect(m)
            need(r==self.grid_reflect(m),'Two independent reflection implementations disagree')
            need(self.reflect(c)==r^self.bitfull,'C and R fail to commute')
            need(self.node[r]==self.node[m],'R changed merged node')
            need(r!=c,'Unexpected R=C antipalindromic mask')
            lines.append({'masks':[m,c],'nodes':sorted([self.node[m],self.node[c]]),
                          'color':'blue' if r==m else 'black','reflected_line':min(r,r^self.bitfull)})
        edges=Counter((tuple(l['nodes']),l['color']) for l in lines)
        counts=Counter(l['color'] for l in lines)
        degree=Counter(); blue_degree=Counter(); black_degree=Counter()
        for l in lines:
            for u in l['nodes']:
                degree[u]+=1
                (blue_degree if l['color']=='blue' else black_degree)[u]+=1
        for row in rows:
            u=row['node']
            need(degree[u]==row['degree'] and blue_degree[u]==row['blue_degree'] and black_degree[u]==row['black_degree'],
                 'Loop-aware incidence degree mismatch')
        need(all(k%2==0 for (uv,c),k in edges.items() if c=='black'),'Black edge multiplicity not even')
        self_converse=sum(self.classof[c^self.full]==c for c in self.ordinary)
        need(2*len(rows)==len(self.ordinary)+self_converse,'Converse Burnside mismatch')
        # Connected components in the literal undirected multigraph.
        unseen=set(range(len(rows))); comps=[]
        while unseen:
            todo=[min(unseen)]; comp=set()
            while todo:
                a=todo.pop()
                if a in comp: continue
                comp.add(a)
                for l in lines:
                    if a in l['nodes']: todo.extend(v for v in l['nodes'] if v not in comp)
            unseen-=comp; comps.append(sorted(comp))
        return {'n':self.n,'labelled':self.full+1,'ordinary':len(self.ordinary),'self_converse':self_converse,
                'merged':len(rows),'mask_count':self.bitfull+1,'nodes':rows,'lines':lines,
                'edges':[{'nodes':list(uv),'color':c,'multiplicity':k} for (uv,c),k in sorted(edges.items())],
                'line_counts':dict(counts),'components':comps,'cycle_rank':len(lines)-len(rows)+len(comps)}


def component_count(v,edges):
    parent=list(range(v))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    for u,w in edges:
        a,b=find(u),find(w)
        if a!=b:parent[a]=b
    return len({find(i) for i in range(v)})


def fundamental_cycles(v,edges):
    """Integral Z-basis via a spanning forest, with explicit oriented coefficients."""
    parent=list(range(v));tree=[];chords=[];adj=defaultdict(list)
    def find(a):
        while parent[a]!=a:a=parent[a]
        return a
    for k,(u,w) in enumerate(edges):
        a,b=find(u),find(w)
        if a==b:chords.append(k)
        else:
            parent[a]=b;tree.append(k)
            adj[u].append((w,k,1));adj[w].append((u,k,-1))
    out=[]
    for k in chords:
        u,w=edges[k];vec=[0]*len(edges);vec[k]=1
        todo=[(w,[])];seen=set();path=None
        while todo:
            a,p=todo.pop()
            if a==u:path=p;break
            if a in seen:continue
            seen.add(a)
            todo.extend((b,p+[(e,s)]) for b,e,s in adj[a] if b not in seen)
        need(path is not None,'Spanning forest misses chord endpoints')
        for e,s in path:vec[e]+=s
        boundary=[0]*v
        for (a,b),z in zip(edges,vec):boundary[a]-=z;boundary[b]+=z
        need(not any(boundary),'Fundamental chain is not a cycle')
        out.append(vec)
    return out,chords


def rank_mod2(vectors):
    pivots={}
    for vec in vectors:
        word=sum((z&1)<<k for k,z in enumerate(vec))
        while word:
            k=word.bit_length()-1
            if k in pivots:word^=pivots[k]
            else:pivots[k]=word;break
    return len(pivots)


def det_integer(a):
    if not a:return 1
    a=[r[:] for r in a];sign=1;den=1;n=len(a)
    for k in range(n-1):
        if a[k][k]==0:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=a[i][j]*pivot-a[i][k]*a[k][j]
                need(val%den==0,'Bareiss division not exact')
                a[i][j]=val//den
            a[i][k]=0
        den=pivot
    return sign*a[-1][-1]


def lattice_audit(census):
    lines=census['lines'];v=census['merged'];seen=set();orbits=[]
    for k,l in enumerate(lines):
        if k in seen:continue
        r=l['reflected_line'];orb=sorted({k,r});seen.update(orb)
        need(all(lines[j]['nodes']==l['nodes'] for j in orb),'Reflection changed edge endpoints')
        need(len(orb)==(1 if l['color']=='blue' else 2),'Color/orbit-size mismatch')
        orbits.append(orb)
    qedges=[lines[o[0]]['nodes'] for o in orbits]
    full_edges=[l['nodes'] for l in lines]
    qcycles,qchords=fundamental_cycles(v,qedges)
    gcycles,gchords=fundamental_cycles(v,full_edges)
    blacks=[k for k,o in enumerate(orbits) if len(o)==2]
    blues=[k for k,o in enumerate(orbits) if len(o)==1]
    # Lift each Q-cycle by its selected orbit representatives, and append
    # each representative-minus-twin. This is an integral basis, not merely Q-basis.
    split=[]
    for cycle in qcycles:
        row=[0]*len(lines)
        for z,o in zip(cycle,orbits):row[o[0]]=z
        split.append(row)
    for k in blacks:
        row=[0]*len(lines);a,b=orbits[k];row[a]=1;row[b]=-1;split.append(row)
    need(len(split)==len(gcycles),'Cycle splitting rank mismatch')
    transition=[[row[k] for row in split] for k in gchords]
    determinant=det_integer(transition)
    need(abs(determinant)==1,'Selected-orbit cycle splitting not unimodular')
    # If x is a quotient cycle coordinate and y a twin-difference coordinate,
    # R(x,y)=(x,-A*x-y). Hence integral even lifts exist iff A*x is even.
    a=[[cycle[k] for cycle in qcycles] for k in blacks]
    image_rank=rank_mod2(a)
    cb=component_count(v,[qedges[k] for k in blues]);cq=component_count(v,qedges)
    beta_blue=len(blues)-v+cb
    need(image_rank==len(qcycles)-beta_blue,'Parity obstruction rank disagrees with blue-cycle kernel')
    need(image_rank==len(blacks)+cq-cb,'Components/index identity failed')
    expected_power={'3':0,'4':0,'5':10}[str(census['n'])]
    need(image_rank==expected_power,'Small-order integral index boundary mismatch')
    return {'literal_cycle_rank':len(gcycles),'quotient_cycle_rank':len(qcycles),
            'blue_cycle_rank':beta_blue,'blue_components':cb,'quotient_components':cq,
            'black_pairs':len(blacks),'even_rank':len(qcycles),'odd_rank':len(blacks),
            'even_plus_odd_index':1<<image_rank,'parity_rank':image_rank,
            'split_basis_determinant':determinant,'black_restriction_matrix':a,
            'quotient_orbits':orbits,'quotient_cycle_basis':qcycles,
            'orientation':'R-equivariant; loop twins chosen compatibly, not by independent mask order'}


def main():
    universes={n:Universe(n) for n in (3,4,5)}
    out={str(n):universes[n].census(universes.get(n-1)) for n in (3,4,5)}
    need([(out[str(n)]['ordinary'],out[str(n)]['self_converse'],out[str(n)]['merged']) for n in (3,4,5)]==[(2,2,2),(4,2,3),(12,8,10)],'Census boundary mismatch')
    c=out['5']; n5=universes[5]
    need(len({(r['H'],r['blue_degree']) for r in c['nodes']})==10,'n5 (H,blue) classifier collision')
    n6=Universe(6);c6=n6.census(n5)
    need((c6['labelled'],c6['ordinary'],c6['self_converse'],c6['merged'])==(32768,56,12,34),'Independent full-labelled n6 census')
    hostiles=[]
    for masks,expected in [((1,77),(17,5)),((50,43),(29,7))]:
        rows=[c6['nodes'][n6.node[m]] for m in masks]
        need(rows[0]['node']!=rows[1]['node'],'n6 hostile nodes collapsed')
        need(all((r['H'],r['blue_degree'])==expected and r['Aut']==1 for r in rows),'n6 (H,blue) hostile changed')
        need(rows[0]['deletion']!=rows[1]['deletion'],'n6 incidence sidecar failed to separate')
        hostiles.append({'masks':list(masks),'nodes':rows})
    out['n6_classifier_boundary']={'labelled':32768,'ordinary':56,'self_converse':12,'merged':34,'hostiles':hostiles}
    conflict_hostiles=[]
    for m in (1,26):
        g=n5.graphs[m];cycles=[]
        for root in range(5):
            todo=[(root,)]
            while todo:
                path=todo.pop()
                if len(path)>=3 and len(path)%2 and n5.edge(g,path[-1],root):cycles.append(path)
                todo.extend(path+(v,) for v in range(root+1,5) if v not in path and n5.edge(g,path[-1],v))
        def norm(c):
            k=c.index(min(c));return c[k:]+c[:k]
        signatures=sorted((len(c),norm(tuple(4-v for v in reversed(c)))==c) for c in cycles)
        need(signatures==[(3,False),(3,False),(3,True),(5,True)],'Equivariant conflict hostile changed')
        need(all(set(a)&set(b) for k,a in enumerate(cycles) for b in cycles[k+1:]),'Hostile conflict graph is not K4')
        conflict_hostiles.append({'mask':m,'node':n5.node[m],'cycles':cycles,'length_fixed_signature':signatures,
                                 'common_intersection':sorted(set.intersection(*(set(c) for c in cycles))),
                                 'blue_degree':c['nodes'][n5.node[m]]['blue_degree']})
    need([h['blue_degree'] for h in conflict_hostiles]==[1,3],'Equivariant conflict hostile blue distinction')
    out['equivariant_conflict_hostiles']=conflict_hostiles
    need(c['line_counts']=={'blue':8,'black':24},'Blue/black line counts mismatch')
    need([r['degree'] for r in c['nodes']]==[1,2,10,1,9,9,11,5,13,3],'Degree vector mismatch')
    need([r['blue_degree'] for r in c['nodes']]==[1,0,0,1,1,3,3,1,3,3],'Blue degree mismatch')
    need(len({tuple(r['primitive']) for r in c['nodes']})==10,'Primitive coordinates not distinct')
    need([r['node'] for r in c['nodes'] if r['type']=='pure_blue']==[0,3,9],'Pure-blue node boundary mismatch')
    # Ratio-gain construction and exact rational decoder, with selected sidecar.
    b=Fraction(55)*Fraction(1,7)*Fraction(1,17)*Fraction(17,91)
    cc=Fraction(19,51)*Fraction(77,19)*Fraction(77,29)*Fraction(29,33)
    d=Fraction(78,85)*Fraction(23,38)*Fraction(95,23)*Fraction(15,2)
    need(b*d/cc==Fraction(2025,4802),'PRIMEGAME triangle holonomy changed')
    enc=[b**r['primitive'][0]*cc**r['primitive'][1]*d**r['primitive'][2] for r in c['nodes']]
    gains=[]
    for l in c['lines']:
        u,v=l['nodes']; gains.append(enc[v]/enc[u])
    # All edge gains from a vertex potential have trivial cycle holonomy; this
    # exact edgewise identity is stronger than a finite cycle-basis test.
    need(all(g*enc[u]==enc[v] for l,g in zip(c['lines'],gains) for u,v in [l['nodes']]),'Vertex potential gain identity failed')
    need(all(g==1 for l,g in zip(c['lines'],gains) if l['nodes'][0]==l['nodes'][1]),'Potential gain on loop nontrivial')
    def test_swapped_color():
        bad=[dict(l) for l in c['lines']];bad[0]['color']='black'
        need(all((l['color']=='blue')==(n5.reflect(l['masks'][0])==l['masks'][0]) for l in bad),'Mutated color')
    def test_removed_parallel_line():
        bad=c['lines'][:-1]
        need(sum(len(l['nodes']) for l in bad)==64,'Dropped line')
    def test_loop_single_incidence():
        bad=Counter()
        for l in c['lines']:
            for u in set(l['nodes']):bad[u]+=1
        need([bad[u] for u in range(10)]==[r['degree'] for r in c['nodes']],'Loops counted once')
    def test_missing_aut():
        need(all(r['degree']==len(r['classes'])*r['H'] for r in c['nodes']),'Forgot automorphism denominator')
    def test_missing_converse_sheet():
        need(all(r['degree']*r['Aut']==r['H'] for r in c['nodes']),'Forgot nonself converse sheet')
    out['integral_cycle_audit']={str(n):lattice_audit(out[str(n)]) for n in (3,4,5)}
    out['gain_audit']={'block_products':list(map(str,[b,cc,d])),'primegame_triangle_holonomy':str(b*d/cc),
                       'potential_gain_has_trivial_holonomy':True,'potential_gains':[str(g) for g in gains]}
    out['negative_controls']={name:rejected(fn) for name,fn in [
       ('swapped_color',test_swapped_color),('removed_literal_parallel_line',test_removed_parallel_line),
       ('loop_counted_once',test_loop_single_incidence),('omitted_Aut',test_missing_aut),
       ('omitted_converse_sheet',test_missing_converse_sheet)]}
    out['method']={'producer_imports':False,'producer_data_reads':False,'universe':'all labelled tournaments n=3,4,5,6; full n6 retained internally with two exported classifier hostiles',
                   'class_method':'explicit permutation orbits before fixed-path atlas',
                   'optimized_safe':True,'source_baseline':'5277b4505b6426f2381a01eb9186f0779566ce4e'}
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=='__main__':
    import sys
    if '--test-failing-gate' in sys.argv:need(False,'Intentional optimized-safe negative control')
    main()
