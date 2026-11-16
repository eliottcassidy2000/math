from itertools import product, permutations, combinations
from fractions import Fraction as F

def arithmetic():
    pts=[(a,b,c) for a,b,c in product(range(3),range(2),range(2)) if (a,b,c) not in [(0,0,0),(2,1,1)]]
    f=lambda a,b,c:(b-c,2*abs(a-1)-1 if b!=c else 0)
    g=lambda a,b,c:(b-c,a-1 if b!=c else 0)
    for a,b,c in pts:
        x,y=f(a,b,c)
        assert f(a,c,b)==f(2-a,1-b,1-c)==(-x,y)
        x,y=g(a,b,c)
        assert g(a,c,b)==(-x,y) and g(2-a,1-b,1-c)==(-x,-y)
    assert len({f(*p) for p in pts})==5
    assert len({g(*p) for p in pts})==7

def phase():
    cycle=[2,8,20,26,11,17]
    delay={8,11}
    for i,x in enumerate(cycle):
        assert (19-x)%27==cycle[5-i]
        assert (x in delay)==((19-x)%27 in delay)
        assert (2**(2+2*(i%2))*x)%27==cycle[(i+1)%6]
    assert sorted((19-x)%27 for x in [5,14,23])==[5,14,23]

def tn():
    n=6
    edges=list(combinations(range(n),2))
    tiles=[x for x in edges if x[1]-x[0]>1]
    def bits(u):
        return {e for i,e in enumerate(tiles) if u>>i&1}
    def mask(D):
        return sum(1<<i for i,e in enumerate(tiles) if e in D)
    def mirror(u):
        return mask({(n-1-b,n-1-a) for a,b in bits(u)})
    def compute(u):
        B=bits(u)
        def beat(a,b):
            return ((a,b) not in B) if a<b else ((b,a) in B)
        paths=sum(all(beat(a,b) for a,b in zip(p,p[1:])) for p in permutations(range(n)))
        code=min(sum(int(beat(p[b],p[a]))<<i for i,(a,b) in enumerate(edges)) for p in permutations(range(n)))
        return code,paths
    seq=[68,70,324,102,356]
    arb=(1<<len(tiles))-1
    ids=[compute(u) for u in seq]
    assert [c for c,h in ids]==[145,152,281,156,282]
    complements=[compute(u^arb) for u in seq]
    hub=bits(955)
    def beats(a,b):
        return ((a,b) not in hub) if a<b else ((b,a) in hub)
    for a,b in edges:
        assert beats(a,b)==(((b-a)%3==1) if (a-b)%3 else a>b)
    all8=[68^z for z in (0,2,256,32,2^256,2^32,256^32,2^256^32)]
    assert len({compute(u)[0] for u in all8})==8
    assert {compute(u^arb)[0] for u in all8}=={408}
    assert {c for c,h in complements}=={408}
    assert [mirror(u) for u in seq]==[68,324,70,356,102]
    print("literal metagraph leaves (code,H)",ids,"hub",complements[0])

def shortcut_children(x,cut=True):
    return [2*x]+([(2*x-1)//3] if x%3==2 and (not cut or x!=2) else [])

def skeleton(x):
    pairs=[]
    for y in shortcut_children(x):
        if y%3==0:continue
        e=1
        while y%9 not in (2,8):
            y=2*y;e+=1
        pairs.append((y,e))
    return pairs

def embed(root,base,steps,children,m):
    # Store by binary address, with intrinsic depth and interval.
    A={}
    def visit(x,key,r,S,v):
        A[key]=(r+S/2,S,v)
        if len(key)==steps:return
        for i,(z,e) in enumerate(children(x)):
            T=S/base**e
            visit(z,key+(i,),r if not i else r+S-T,T,v+e)
    visit(root,(),F(0),F(1),0)
    for (a,(x,y,s)),(b,(z,w,t)) in combinations(A.items(),2):
        l=0
        while l<min(len(a),len(b)) and a[l]==b[l]:l+=1
        D=s+t-2*A[a[:l]][2]
        C=(y*y+w*w+(x-z)**2)/(2*y*w)
        def ge(p):return p<=1 or 2*C>=p+1/p
        assert ge(m*m*F(base**D))
        assert not ge(m*m*F(base**(D+1)))
        assert m*m*base**D <= 2*C < m*m*base**(D+1)
        if l<min(len(a),len(b)):assert 2*C < F(base**D)+F(1,base**D)
        else:
            assert 2*C>F(base**D)+F(1,base**D)
            p=F(5,4)*base**D
            assert 2*C<=p+1/p
    print('rational distance decoder',base,len(A),len(A)*(len(A)-1)//2,'pairs')

arithmetic();phase();tn()
embed(8,3,7,skeleton,F(53,81))
embed(1,5,10,lambda x:[(y,1) for y in shortcut_children(x)],F(3,5))
