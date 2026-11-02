#!/usr/bin/env python3
"""Exact finite controls, not a re-run of AMP's extremal enumeration. Stdlib only."""
from fractions import Fraction as F
from itertools import combinations, product
import json


def need(p, message):
    if not p:
        raise ValueError(message)

# Q(sqrt(3)): a+b sqrt(3), exact rational coefficient pairs.
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def neg(a): return (-a[0],-a[1])
def sub(a,b): return add(a,neg(b))
def mul(a,b): return (a[0]*b[0]+3*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(c,a): return (c*a[0],c*a[1])
Z=(F(0),F(0)); ONE=(F(1),F(0))
TRI=((Z,Z),(ONE,Z),((F(1,2),F(0)),(F(0),F(1,2))))
ROT=((F(1),F(0)),(F(3,5),F(4,5)),(F(5,13),F(12,13)),(F(8,17),F(15,17)))
def rotate(p,r):
    c,s=r
    need(c*c+s*s==1,'rotation not unit')
    return (sub(scale(c,p[0]),scale(s,p[1])),add(scale(s,p[0]),scale(c,p[1])))
def padd(p,q): return add(p[0],q[0]),add(p[1],q[1])
def d2(p,q):
    x,y=sub(p[0],q[0]),sub(p[1],q[1])
    return add(mul(x,x),mul(y,y))
def unit_edges(pts):
    need(len(set(pts))==len(pts),'collision')
    return {(i,j) for i,j in combinations(range(len(pts)),2) if d2(pts[i],pts[j])==ONE}
def hamming(d):
    labels=list(product(range(3),repeat=d)); pts=[]
    for label in labels:
        p=(Z,Z)
        for i,a in enumerate(label): p=padd(p,rotate(TRI[a],ROT[i]))
        pts.append(p)
    edges=unit_edges(pts)
    expected={(i,j) for i,j in combinations(range(len(labels)),2)
              if sum(a!=b for a,b in zip(labels[i],labels[j]))==1}
    need(edges==expected,'product is not faithful')
    need(len(edges)==d*3**d,'product edge formula')
    need(all(sum(labels[i])%3!=sum(labels[j])%3 for i,j in edges),'three coloring')
    # Each coordinate cyclically advances 0->1->2->0 and supplies exactly one
    # outgoing edge per vertex; all d layers are disjoint and exhaustive.
    index={a:i for i,a in enumerate(labels)}; layers=[]
    for k in range(d):
        layer=set()
        for i,a in enumerate(labels):
            b=list(a); b[k]=(b[k]+1)%3; j=index[tuple(b)]
            layer.add(tuple(sorted((i,j))))
        need(len(layer)==len(labels),'coordinate layer not one cycle per fiber')
        layers.append(layer)
    need(set.union(*layers)==edges and sum(map(len,layers))==len(edges),'layer partition')
    return {'d':d,'n':len(pts),'m':len(edges),'chi':3,'functional_layers_exact':d,
            'degree':2*d,'cycle_rank':len(edges)-len(pts)+1}

# Corrected Eisenstein witness: w^2=w-1, Q(a,b)=a^2+ab+b^2;
# actual Euclidean points are (a+b*w)/sqrt(7).
def emul(x,y):
    a,b=x;c,d=y
    return a*c-b*d,a*d+b*c+b*d
def eadd(x,y): return x[0]+y[0],x[1]+y[1]
def norm(x): a,b=x;return a*a+a*b+b*b
def eisenstein21():
    tri=[(0,0),(1,0),(0,1)]
    wheel=[(0,0),(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    pts=[eadd(emul((2,1),p),emul((3,-1),q)) for p in tri for q in wheel]
    need(len(set(pts))==21,'Eisenstein collisions')
    edges={(i,j) for i,j in combinations(range(21),2)
           if norm((pts[i][0]-pts[j][0],pts[i][1]-pts[j][1]))==7}
    need(len(edges)==57,'Eisenstein21 edge count')
    # This is K3 square W7; the triangular norm-1 cap does not apply at norm7.
    need(min(norm((p[0]-q[0],p[1]-q[1])) for p,q in combinations(pts,2))<7,
         'the penny-graph obstruction was not exposed')
    return {'n':21,'m':57,'unit_squared_norm_before_scaling':7,'coordinates':pts,
            'minimum_squared_norm_before_scaling':min(norm((p[0]-q[0],p[1]-q[1])) for p,q in combinations(pts,2))}

# Independent finite formulations: component deficit and Hall/matching capacity.
def tree_components(n,edges):
    adj=[set() for _ in range(n)]
    for a,b in edges: adj[a].add(b);adj[b].add(a)
    unseen=set(range(n)); t=0
    while unseen:
        stack=[unseen.pop()];vs=[]
        while stack:
            v=stack.pop();vs.append(v)
            for w in adj[v]:
                if w in unseen: unseen.remove(w);stack.append(w)
        if sum(len(adj[v]) for v in vs)//2==len(vs)-1: t+=1
    return t

def matching_size(n,edges,k):
    edges=list(edges);owner={}
    def visit(i,seen):
        for v in edges[i]:
            for color in range(k):
                slot=(v,color)
                if slot in seen: continue
                seen.add(slot)
                if slot not in owner or visit(owner[slot],seen):
                    owner[slot]=i;return True
        return False
    return sum(visit(i,set()) for i in range(len(edges)))
def density_ok(n,edges,k):
    return all(sum((mask>>a&1) and (mask>>b&1) for a,b in edges)<=k*mask.bit_count()
               for mask in range(1<<n))
def graph_controls():
    count=0;subdivisions=0
    for n in range(6):
        possible=list(combinations(range(n),2))
        for code in range(1<<len(possible)):
            edges=[e for i,e in enumerate(possible) if code>>i&1]
            t=tree_components(n,edges); defect=len(edges)-n+t
            need(defect==len(edges)-matching_size(n,edges,1),'functional deficit formula')
            for k in range(1,4):
                need((matching_size(n,edges,k)==len(edges))==density_ok(n,edges,k),'Hall layer equivalence')
            if edges:
                a,b=edges[0];subedges=edges[1:]+[(a,n),(n,b)]
                need(len(subedges)-(n+1)+tree_components(n+1,subedges)==defect,'subdivision defect')
                subdivisions+=1
            count+=1
    return {'all_labeled_simple_graphs_n0_through5':count,'layer_capacities_checked':[1,2,3],
            'single_edge_subdivisions_checked':subdivisions}

EXACT=[0,0,1,3,5,7,9,12,14,18,20,23,27,30,33,37,41,43,46,50,54,57]
BOUNDS={22:(60,61),23:(64,66),24:(68,72),25:(72,78),26:(76,84),27:(81,90),28:(85,96)}
def arithmetic():
    need(next(n for n,x in enumerate(EXACT) if x>n)==4,'first two-layer threshold')
    need(next(n for n,x in enumerate(EXACT) if x>2*n)==11,'first three-layer threshold')
    need(all(x<=3*n for n,x in enumerate(EXACT)),'small crossover exclusion')
    need(all(BOUNDS[n][1]<=3*n for n in [22,23,24]),'floor25')
    need(BOUNDS[28][0]>3*28,'ceiling28')
    rows=[]
    for n in [15,20,24,25,26]:
        lo,hi=(EXACT[n],EXACT[n]) if n<len(EXACT) else BOUNDS[n]
        rows.append({'n':n,'AMP2025_lower':lo,'AMP2025_upper':hi,
                     'functional_edges_at_most':n,'functional_edges_missing_at_least':lo-n,
                     'cycle_rank_at_least':lo-n+1,'maximum_degree_at_least':(2*lo+n-1)//n,
                     'honeycomb_girth6_edge_upper':(3*n-6)//2})
    return {'rows':rows,'first_layer_obstructions':{'one_to_two':4,'two_to_three':11,
             'three_to_four_lower':25,'three_to_four_upper':28}}

def expect_rejection(action,label):
    try: action()
    except ValueError: return label
    raise ValueError('negative control unexpectedly accepted: '+label)

def main():
    controls=[]
    controls.append(expect_rejection(lambda:unit_edges([TRI[0],TRI[0]]),'duplicate points rejected'))
    controls.append(expect_rejection(lambda:rotate(TRI[1],(F(3,5),F(3,5))),'nonunit rotation rejected'))
    controls.append(expect_rejection(lambda:need(len(unit_edges(TRI))==2,'wrong edge count'),'wrong edge count rejected'))
    k4=list(combinations(range(4),2))
    controls.append(expect_rejection(lambda:need(density_ok(4,k4,1),'K4 not functional'),'K4 false one-layer claim rejected'))
    controls.append(expect_rejection(lambda:need(matching_size(4,k4,1)==6,'capacity false'),'K4 false capacity rejected'))
    report={'status':'PASS','arithmetic':arithmetic(),'product_witnesses':[hamming(d) for d in range(1,5)],
            'eisenstein21':eisenstein21(),'finite_graph_controls':graph_controls(),'hostile_controls':controls,
            'scope':'Exact finite controls only. AMP global extremality is cited, not recomputed. No Collatz convergence or plane chromatic improvement.'}
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__': main()
