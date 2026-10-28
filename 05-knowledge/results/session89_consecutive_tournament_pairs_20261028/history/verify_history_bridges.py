#!/usr/bin/env python3
"""Independent exact checks for the history lane; no producer imports."""
from collections import Counter
from itertools import permutations
import json

checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)

def positions(m):
    a = [(j-2,i) for i in range(m+2) for j in range(i+2,m+2)]
    b = [(i,j-1) for i in range(m+1) for j in range(i+2,m+1)]
    old_b = [(i,j-2) if i < j-2 else None
             for i in range(m+1) for j in range(i+2,m+1)]
    return a,b,old_b

def pack(word, positions, m):
    out = 0
    for k, rc in enumerate(positions):
        if rc is not None and word >> k & 1:
            r,c = rc
            out |= 1 << (r*m+c)
    return out

def unpack(word, positions, m):
    out=0
    for k,(r,c) in enumerate(positions):
        out |= ((word >> (r*m+c)) & 1) << k
    return out

rows=[]
for m in range(1,13):
    a,b,old_b=positions(m)
    ps=a+b
    require(len(ps)==m*m, ('cell count',m))
    require(set(ps)=={(r,c) for r in range(m) for c in range(m)}, ('bijection',m))
    require(len(set(ps))==len(ps), ('injective',m))
    old=a+old_b
    old_image={rc for rc in old if rc is not None}
    holes={(i,m-1) for i in range(m-1)}
    require(set(ps)-old_image==holes, ('old holes',m))
    require(old.count(None)==m-1, ('lost count',m))
    require(len(old_image)==m*m-m+1, ('old image rank',m))
    for k in range(m*m):
        require(unpack(pack(1<<k,ps,m),ps,m)==1<<k, ('basis inverse',m,k))
    if m<=4:
        hist=Counter()
        for word in range(1<<(m*m)):
            out=pack(word,ps,m)
            require(unpack(out,ps,m)==word, ('roundtrip',m,word))
            hist[pack(word,old,m)]+=1
            if m<=3:
                for k in range(m*m):
                    require((out ^ pack(word^(1<<k),ps,m)).bit_count()==1, ('cube edge',m,word,k))
        require(len(hist)==1<<(m*m-m+1), ('image size',m))
        require(set(hist.values())=={1<<(m-1)}, ('fibres',m))
        rows.append({'side':m,'inputs':1<<(m*m),'old_image_size':len(hist),'old_fibre_size':1<<(m-1)})

# Explicit matrix-to-colored-graph reduction: every relabeling preserves exactly the arcs.
def colored_graph(n,word):
    vertices={(kind,i) for kind in ('u','r','c') for i in range(n)}
    edges=set()
    def edge(a,b): edges.add(frozenset((a,b)))
    for i in range(n):
        edge(('u',i),('r',i));edge(('u',i),('c',i))
    for i in range(n):
        for j in range(n):
            if word>>(i*n+j)&1:
                w=('w',i,j);vertices.add(w)
                edge(('r',i),w);edge(('c',j),w)
    return vertices,edges

graph_checks=0
for n in range(1,4):
    for word in range(1<<(n*n)):
        vertices,edges=colored_graph(n,word)
        decoded=sum(1<<(v[1]*n+v[2]) for v in vertices if v[0]=='w')
        require(decoded==word,('colored decode',n,word))
        N=len(vertices);L=N+2
        deg={v:sum(v in e for e in edges)+1 for v in vertices}
        require(all(2<=d<=N+1 for d in deg.values()),('uncolor original degrees',n,word))
        hubs=[(k+1)*L+sum(v[0]==kind for v in vertices) for k,kind in enumerate(('u','r','c','w'))]
        require(min(hubs)>N+1 and len(set(hubs))==4,('uncolor hub separation',n,word))
        for p in permutations(range(n)):
            def image(v): return (v[0],)+tuple(p[i] for i in v[1:])
            moved=sum(((word>>(i*n+j))&1)<<(p[i]*n+p[j]) for i in range(n) for j in range(n))
            target_v,target_e=colored_graph(n,moved)
            require({image(v) for v in vertices}==target_v,('colored vertex equivariance',n,word,p))
            require({frozenset(image(v) for v in e) for e in edges}==target_e,('colored edge equivariance',n,word,p))
            graph_checks+=1

# Two binary matrices similar over Q (indeed any field) but not permutation similar.
A=((1,0),(0,0));B=((1,1),(0,0));S=((1,-1),(0,1))
def mul(X,Y):
    return tuple(tuple(sum(X[i][k]*Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
require(mul(S,A)==mul(B,S),'GL intertwiner')
require(S[0][0]*S[1][1]-S[0][1]*S[1][0]==1,'GL invertible')
require(all(tuple(tuple(A[p[i]][p[j]] for j in range(2)) for i in range(2))!=B for p in permutations(range(2))),'no permutation intertwiner')

# Ordinary 2x2 bit-matrix conjugacy has 10 orbits.
def perm_matrix_code(w,p,n):
    return sum(((w>>(p[i]*n+p[j]))&1) << (i*n+j) for i in range(n) for j in range(n))
require(len({min(perm_matrix_code(w,p,2) for p in permutations(range(2))) for w in range(16)})==10,'matrix n2 orbits')

# Field and substitution hostiles.
require(1*1-1*1==0,'ones determinant')
require(1*1+1*1==2,'ones permanent')
for t in range(-10,11):
    require(t*t-t**2==0,'dyadic curve polynomial collision')
require(0-1**2!=0,'curve polynomial is not globally zero')
# Binary-exponent substitution is injective only with the multilinear support promise.
for d in range(1,13):
    exponents={sum((1<<i) for i in range(d) if mask>>i&1) for mask in range(1<<d)}
    require(exponents==set(range(1<<d)),('multilinear Kronecker',d))

print(json.dumps({'status':'PASS','always_active_checks':checks,'square_exhaustive':rows,
 'coordinate_basis_sides':'1..12','Hamming_all_cube_edges_sides':'1..3',
 'other_checks':['matrix-to-graph equivariance and decoding through n=3','GL-versus-permutation 2x2 hostile','ordinary matrix n2 orbit count 10','determinant/permanent distinction','dyadic substitution collision and scoped exponent injectivity'],
 'limitations':['No global isomorphism algorithm','No unrestricted PIT derandomization','No unlabelled deck reconstruction theorem','No revalidation of all historical claims']},indent=2))
