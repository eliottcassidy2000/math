from itertools import combinations
from math import isqrt

def check(c, message):
    if not c: raise RuntimeError(message)

def square(n): return n>=0 and isqrt(n)**2==n

def edges(V): return {tuple(sorted((x,y))) for x,y in combinations(sorted(V),2) if square(x+y)}

def triangles_graph(N):
    E=edges(range(1,N+1)); adj={x:set() for x in range(1,N+1)}
    for x,y in E: adj[x].add(y);adj[y].add(x)
    return {(x,y,z) for x in adj for y in adj[x] if y>x for z in adj[x]&adj[y] if z>y}

def triangles_roots(N):
    out=set()
    for a,b,c in combinations(range(1,isqrt(2*N-1)+1),3):
        aa,bb,cc=a*a,b*b,c*c
        if aa+bb>cc and (aa+bb+cc)%2==0:
            v=((aa+bb-cc)//2,(aa+cc-bb)//2,(bb+cc-aa)//2)
            if v[-1]<=N: out.add(v)
    return out

def components(V,E):
    V=set(V); out=[]
    while V:
        stack=[min(V)]; C=set()
        while stack:
            v=stack.pop()
            if v in C: continue
            C.add(v)
            stack.extend(y if x==v else x for x,y in E if x==v or y==v)
        V-=C;out.append(C)
    return out

def common_by_factor(x,z,N):
    d=z-x; out=set()
    for u in range(1,isqrt(d)+1):
        if d%u: continue
        v=d//u
        if (u+v)%2 or u>=v: continue
        a=(v-u)//2; y=a*a-x
        if 1<=y<=N and y not in (x,z): out.add(y)
    return out

A={10,6,3,1,8,17};B={12,4,5,11};C={2,7,9,18,16};V=A|B|C
E=edges(V)
expected={(1,3),(1,8),(2,7),(3,6),(4,5),(4,12),(5,11),(6,10),(7,9),(7,18),(8,17),(9,16)}
check(E==expected,'exact user edge set')
check({frozenset(c) for c in components(V,E)}=={frozenset(A),frozenset(B),frozenset(C)},'components')
check(len(E)==12,'edge count')
# All subgraphs of this forest, including isolated vertices, which have maximum degree two are path covers.
valid=[]
for mask in range(1<<len(E)):
    F={e for i,e in enumerate(sorted(E)) if mask>>i&1}
    deg={v:0 for v in V}
    for x,y in F: deg[x]+=1;deg[y]+=1
    if max(deg.values())<=2: valid.append(F)
max_edges=max(map(len,valid)); winners=[F for F in valid if len(F)==max_edges]
check(max_edges==11 and len(winners)==3,'minimum path covers')
check({next(iter(E-F)) for F in winners}=={(2,7),(7,9),(7,18)},'cut edge list')
for N in (30,50,100,200):
    actual=triangles_graph(N); formula=triangles_roots(N)
    check(actual==formula, f'triangle bijection N={N}')
check(triangles_graph(29)==set(),'first triangle lower bound')
check(triangles_graph(30)=={(6,19,30)},'first triangle unique')
for t in range(1,100):
    x,y,z=2*(t*t-1),2*t*t+4*t+3,2*t*t+8*t+6
    check(x>=0 and x<y<z and [x+y,x+z,y+z]==[(2*t+1)**2,(2*t+2)**2,(2*t+3)**2],f'family {t}')
for N in (12,18,30,100,250):
    E=edges(range(1,N+1))
    check(all((x%4==y%4) for x,y in E if x%2==y%2==0),'even residue split')
    check(all(x%4!=y%4 for x,y in E if x%2==y%2==1),'odd bipartition')
    check(all(((x%8,y%8) if x%2 else (y%8,x%8)) in {(1,0),(3,6),(5,4),(7,2)} for x,y in E if (x+y)%2),'mixed parity residues')
    check({(x//4,y//4) for x,y in E if x%4==y%4==0}==edges(range(1,N//4+1)),'scale four')
    check(all(square((x-2)//4+(y-2)//4+1) for x,y in E if x%4==y%4==2),'shifted 2 mod4')
check((2,6) not in edges(range(1,7)) and (1,3) in edges(range(1,7)),'reject factor-two copy')
for x in range(1,61):
    for z in range(x+1,61):
        actual={y for y in range(1,201) if y not in (x,z) and square(x+y) and square(z+y)}
        check(actual==common_by_factor(x,z,200),f'difference factors {x},{z}')
check(not common_by_factor(2,4,1000000),'2 and 4 have no common neighbor')
check(all(square(a+b) for a,b in zip((2,23,77),(23,77,4))),'distance-three witness')
for m in range(1,100):
    k=isqrt(m)+1; x=2*k*k-m;y=2*k*k+m
    check(x>0 and x!=y and (y-x)%m==0 and square(x+y),'periodic-color obstruction')
print('PASS: exact induced graph; minimum path covers; triangle parameterization to 200; infinite-family instances; four-adic structure; common-neighbor factorization to 60; no periodic residue coloring examples')
