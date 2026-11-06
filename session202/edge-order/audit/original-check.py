import itertools,json,hashlib
from pathlib import Path
BASE=Path('/workspace/scratch/07d5532f2fce/session202-threshold-structure')
def tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for bits in itertools.product((0,1),repeat=len(pairs)):
        yield {(a,b) if d else (b,a) for (a,b),d in zip(pairs,bits)}
def paths(vertices,edges):
    # Literal permutations; independent of the author's subset DP.
    return [p for p in itertools.permutations(vertices) if all((a,b) in edges for a,b in zip(p,p[1:]))]
def stat(n,E):
    V=tuple(range(n)); h=len(paths(V,E)); c=0
    for mask in range(1<<n):
        A=[i for i in V if mask>>i&1];B=[i for i in V if not mask>>i&1]
        c+=len(paths(A,E))*len(paths(B,E))
    return h,c
def physical(p,q,P,Q):
    # Encode edge orientation by vertex types rather than use author's build.
    V=['x','y','b']+[f'p{i}' for i in range(p)]+[f'q{i}' for i in range(q)]
    allowed={('x','y'),('y','b'),('b','x'),('x','p'),('x','q'),('p','y'),('q','y'),('p','b'),('b','q'),('p','q')}
    E=set()
    for a,b in itertools.combinations(V,2):
        ta,tb=a[0],b[0]
        if ta==tb:
            source=P if ta=='p' else Q
            forward=(int(a[1:]),int(b[1:])) in source
        else: forward=(ta,tb) in allowed
        E.add((a,b) if forward else (b,a))
    return V,E
def check(p,q,P,Q,unique=False):
    V,E=physical(p,q,P,Q); ps=paths(V,E)
    hp,cp=stat(p,P);hq,cq=stat(q,Q)
    f=sum(z.index('x')<z.index('y') for z in ps);r=len(ps)-f
    expected=(cp*hq+hp*cq+hp*hq,cp*cq)
    if (f,r)!=expected:raise RuntimeError((p,q,f,r,expected))
    cyc=sum(z[0]=='x' and (z[-1],z[0]) in E for z in ps)
    if unique and cyc!=1:raise RuntimeError('cycle')
    return dict(p=p,q=q,f=f,r=r,H=len(ps),cycles=cyc)
mods=[(n,E) for n in range(1,4) for E in tournaments(n)]
records=[]
for p,P in mods:
    for q,Q in mods:records.append(check(p,q,P,Q))
trans=[]
for p in range(1,4):
    for q in range(1,4):
        P=set(itertools.combinations(range(p),2));Q=set(itertools.combinations(range(q),2))
        trans.append(check(p,q,P,Q,True))
files=['EDGE_ORDER.md','verify_edge_order.py','EDGE_NORMAL.json','EDGE_OPTIMIZED.json']
print(json.dumps(dict(status='PASS',method='Independent edge predicates and literal permutations, not subset DP',arbitrary_pairs=len(records),transitive=trans,hashes={f:hashlib.sha256((BASE/f).read_bytes()).hexdigest() for f in files}),indent=2))
