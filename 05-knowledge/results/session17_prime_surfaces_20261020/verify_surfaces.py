#!/usr/bin/env python3
"""Exact finite controls for the prime-field surface/normalization theorem.
No packages, numerical tolerances, or Python assert statements are used.
F9 is F3[t]/(t^2+1), encoded a+3b for a+b*t.
"""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
import json

CHECKS = 0

def check(condition, description):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(description)

class Field:
    def __init__(self, p, degree=1):
        if degree != 1 and (p, degree) != (3, 2):
            raise ValueError('Only prime fields and the explicit F9 model are used')
        self.p, self.degree, self.q = p, degree, p**degree
        self.elements = range(self.q)
        self.minusone = p-1
    def add(self, x, y):
        if self.degree == 1:
            return (x+y) % self.p
        return ((x%3+y%3)%3) + 3*((x//3+y//3)%3)
    def neg(self, x):
        if self.degree == 1:
            return -x % self.p
        return (-x%3) + 3*((-(x//3))%3)
    def sub(self, x, y):
        return self.add(x, self.neg(y))
    def mul(self, x, y):
        if self.degree == 1:
            return (x*y) % self.p
        a,b,c,d = x%3,x//3,y%3,y//3
        return ((a*c-b*d)%3) + 3*((a*d+b*c)%3)
    def vadd(self, v, w):
        return tuple(self.add(a,b) for a,b in zip(v,w))
    def vneg(self, v):
        return tuple(self.neg(a) for a in v)
    def scale(self, a, v):
        return tuple(self.mul(a,x) for x in v)
    def det(self, v, w):
        return self.sub(self.mul(v[0],w[1]), self.mul(v[1],w[0]))
    def rep(self, v):
        return min(v, self.vneg(v))

def components(vertices, neighbors):
    todo = set(vertices)
    answer = []
    while todo:
        seed = min(todo)
        todo.remove(seed)
        seen, stack = {seed}, [seed]
        while stack:
            x = stack.pop()
            for y in neighbors(x):
                if y in todo:
                    todo.remove(y)
                    seen.add(y)
                    stack.append(y)
        answer.append(frozenset(seen))
    return answer

def field_sanity(F):
    for a in F.elements:
        check(F.add(a,F.neg(a)) == 0, ('additive inverse',F.q,a))
        check(F.mul(a,1) == a, ('multiplicative identity',F.q,a))
        if a:
            check(sum(F.mul(a,b)==1 for b in F.elements)==1,
                  ('unique multiplicative inverse',F.q,a))
    for a,b,c in product(F.elements,repeat=3):
        check(F.mul(F.mul(a,b),c)==F.mul(a,F.mul(b,c)),
              ('associative multiplication',F.q,a,b,c))
        check(F.mul(a,F.add(b,c))==F.add(F.mul(a,b),F.mul(a,c)),
              ('distributivity',F.q,a,b,c))
    if F.q==9:
        check(F.mul(3,3)==2, 'F9 polynomial t^2=-1')
        check(all(F.add(F.mul(x,x),1)!=0 for x in range(3)),
              't^2+1 irreducible over F3')

def verify(F):
    p,q = F.p,F.q
    field_sanity(F)
    vertices = sorted({F.rep(v) for v in product(F.elements,repeat=2) if v!=(0,0)})
    edges = {tuple(e) for e in combinations(vertices,2)
             if F.det(*e) in (1,F.minusone)}
    adj = {v:set() for v in vertices}
    for a,b in edges:
        adj[a].add(b); adj[b].add(a)
    faces = {t for t in combinations(vertices,3)
             if all(tuple(e) in edges for e in combinations(t,2))}
    edgefaces = {e:[] for e in edges}
    for t in faces:
        for e in combinations(t,2):
            edgefaces[e].append(t)
    check(len(vertices)==(q*q-1)//2, ('vertex count',q))
    check(len(edges)==q*(q*q-1)//4, ('edge count',q))
    check(len(faces)==q*(q*q-1)//6, ('face count',q))
    check(all(len(adj[v])==q for v in vertices), ('regular degree',q))
    check(len(components(vertices,lambda v:adj[v]))==1, ('connected graph',q))

    edgeplanes = {}
    for e in edges:
        v,w = e
        if F.det(v,w)!=1:
            w=F.vneg(w)
        expected = {F.rep(F.vadd(v,w)), F.rep(F.vadd(v,F.vneg(w)))}
        actual = {next(x for x in t if x not in e) for t in edgefaces[e]}
        check(actual==expected and len(actual)==2, ('two explicit third vertices',q,e))
        W = frozenset(F.rep(F.vadd(F.scale(a,v),F.scale(b,w)))
                      for a,b in product(range(p),repeat=2) if (a,b)!=(0,0))
        check(len(W)==(p*p-1)//2, ('prime-plane size',q,e))
        edgeplanes[e]=W

    # Select orientation by a+b+c=0 and det(a,b)=1, independently of edge planes.
    orientation = {}
    edge_signs = defaultdict(list)
    for t in sorted(faces):
        oriented_representatives=[]
        for ordered in permutations(t):
            for signs in product((1,F.minusone),repeat=3):
                a,b,c=[F.scale(s,v) for s,v in zip(signs,ordered)]
                if F.vadd(F.vadd(a,b),c)==(0,0) and F.det(a,b)==1:
                    oriented_representatives.append((ordered,(a,b,c)))
        check(len(oriented_representatives)==6, ('orientation representative count',q,t))
        orders={x[0] for x in oriented_representatives}
        selected=min(orders)
        cyclic={selected, selected[1:]+selected[:1], selected[2:]+selected[:2]}
        check(orders==cyclic, ('unique cyclic face orientation',q,t))
        orientation[t]=selected
        a,b,c=selected
        for u,v in ((a,b),(b,c),(c,a)):
            e=tuple(sorted((u,v)))
            edge_signs[e].append(1 if u<v else -1)
        check(len({edgeplanes[e] for e in combinations(t,2)})==1,
              ('face stays in a prime plane',q,t))
    for e in edges:
        check(sorted(edge_signs[e])==[-1,1], ('opposite face orientations',q,e))

    # The neighborhood calculation is checked at every vertex, not just e1.
    link_sizes={}
    for v in vertices:
        link_vertices=adj[v]
        check(all(len(adj[u]&link_vertices)==2 for u in link_vertices),
              ('link degree two',q,v))
        comps=components(link_vertices,lambda u:adj[u]&link_vertices)
        check(sorted(map(len,comps))==[p]*(q//p), ('link cycles',q,v))
        link_sizes[v]=tuple(sorted(map(len,comps)))
    N=F.rep((1,0))
    U={x:F.rep((x,1)) for x in F.elements}
    check(set(U.values())==adj[N], ('explicit north neighbors',q))
    for x,y in combinations(F.elements,2):
        check((U[y] in adj[U[x]])==(F.sub(x,y) in (1,F.minusone)),
              ('north link difference predicate',q,x,y))

    # Independent facet-component computation, then compare it with edge spans.
    face_neighbors={t:set() for t in faces}
    for e,ts in edgefaces.items():
        a,b=ts
        face_neighbors[a].add(b); face_neighbors[b].add(a)
    facets=components(faces,lambda t:face_neighbors[t])
    plane_groups=defaultdict(set)
    for t in faces:
        plane_groups[edgeplanes[tuple(t[:2])]].add(t)
    check({frozenset(ts) for ts in plane_groups.values()}==set(facets),
          ('facet components equal prime-plane surfaces',q))
    count=q*(q*q-1)//(p*(p*p-1))
    check(len(facets)==len(plane_groups)==count, ('normalization component count',q))
    surface_counts=[]
    for W,ts in plane_groups.items():
        es={e for t in ts for e in combinations(t,2)}
        vs={v for t in ts for v in t}
        check(vs==set(W), ('component vertices equal prime plane',q))
        check((len(vs),len(es),len(ts))==((p*p-1)//2,p*(p*p-1)//4,p*(p*p-1)//6),
              ('component surface f-vector',q))
        for v in vs:
            neighbors={u for e in es if v in e for u in e if u!=v}
            comps=components(neighbors,lambda u: {w for w in neighbors if tuple(sorted((u,w))) in es})
            check(len(comps)==1 and len(comps[0])==p, ('component connected link',q))
        surface_counts.append((len(vs),len(es),len(ts)))
    multiplicities=Counter(v for W in plane_groups for v in W)
    check(set(multiplicities.values())=={q//p}, ('components through each vertex',q))

    # Matrix orbit of Fp^2 independently checks the coset denominator.
    matrixplanes=Counter()
    matrix_count=0
    for a,b,c,d in product(F.elements,repeat=4):
        v,w=(a,c),(b,d)
        if F.det(v,w)!=1:
            continue
        matrix_count+=1
        W=frozenset(F.rep(F.vadd(F.scale(x,v),F.scale(y,w)))
                    for x,y in product(range(p),repeat=2) if (x,y)!=(0,0))
        matrixplanes[W]+=1
    check(matrix_count==q*(q*q-1), ('SL2 order',q))
    check(set(matrixplanes)==set(plane_groups), ('SL2 orbit equals plane index',q))
    check(set(matrixplanes.values())=={p*(p*p-1)}, ('SL2Fp stabilizer size',q))

    genus=1+(p*p-1)*(p-6)//24
    check(2-2*genus==(p*p-1)*(6-p)//12, ('surface genus arithmetic',q))
    if q==3:
        check(len(edges)==len(vertices)*(len(vertices)-1)//2, 'p3 graph K4')
        check(len(faces)==4, 'p3 tetrahedron boundary, not filled tetrahedron')
    if q==5:
        # Distance distribution of the icosahedron at every vertex.
        for v in vertices:
            seen={v}; frontier={v}; distribution=[1]
            while frontier:
                frontier={u for w in frontier for u in adj[w]}-seen
                if frontier:
                    distribution.append(len(frontier)); seen.update(frontier)
            check(distribution==[1,5,5,1], ('p5 icosahedral distance profile',v))
    if q==9:
        check(all(s==(4,6,4) for s in surface_counts) and len(surface_counts)==30,
              'F9 normalization is thirty tetrahedron boundaries')
        check(all(len(s)==3 for s in link_sizes.values()), 'F9 is not a surface')

    return {'field_order':q,'characteristic':p,'extension_degree':F.degree,
            'field_model':('F3[t]/(t^2+1); a+3b encodes a+b*t' if q==9 else 'integers modulo p'),
            'vertices':len(vertices),'edges':len(edges),'triangles':len(faces),
            'graph_connected':True,'link_cycle_lengths':[p]*(q//p),
            'normalization_components':count,'component_genus':genus,
            'component_f_vector':surface_counts[0],
            'components_per_original_vertex':q//p,
            'euler_characteristic':len(vertices)-len(edges)+len(faces),
            'SL2_order':matrix_count,'SL2Fp_stabilizer_order':p*(p*p-1)}

def hostile_characteristic_two():
    F=Field(2)
    vs=sorted({F.rep(v) for v in product(range(2),repeat=2) if v!=(0,0)})
    es={e for e in combinations(vs,2) if F.det(*e)==1}
    fs={t for t in combinations(vs,3) if all(e in es for e in combinations(t,2))}
    check((len(vs),len(es),len(fs))==(3,3,1), 'characteristic two is a single filled triangle')
    check(all(sum(set(e).issubset(t) for t in fs)==1 for e in es), 'characteristic two edges have one face')
    return {'field_order':2,'f_vector':[3,3,1],'closed_surface':False,
            'reason':'The sign quotient degenerates and every edge has only one incident triangle.'}

if __name__=='__main__':
    runs=[verify(Field(p)) for p in (3,5,7)] + [verify(Field(3,2))]
    hostile=hostile_characteristic_two()
    print(json.dumps({'status':'PASS','checks':CHECKS,'runs':runs,
                      'hostile_control':hostile,
                      'scope':'Finite-exact controls only; uniform proof is in FINITE_FIELD_SURFACES.md.'},
                     indent=2,sort_keys=True))
