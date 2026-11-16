from functools import cache
from itertools import combinations, product

def profile(a):
    n=len(a); F=(1<<n)-1
    for i,j in combinations(range(n),2):
        assert ((a[i]>>j)&1)+((a[j]>>i)&1)==1
    def bit(s):
        while s:
            b=s&-s;s^=b;yield b.bit_length()-1
    @cache
    def p(s,i):
        R=s^(1<<i)
        return sum(p(R,j) for j in bit(R&a[i])) if R else 1
    @cache
    def h(s):
        return sum(p(s,j) for j in bit(s)) if s else 1
    def rooted(q):
        B=F^q
        b=[x|((1<<n) if B>>i&1 else 0) for i,x in enumerate(a)]+[q]
        E=O=0
        for s in range(F+1):
            w=(sum(p(s,i) for i in bit(s&q)) if s else 1)*h(F^s)
            if s.bit_count()%2: E+=w
            else: O+=w
        H=profile(b)[0] if n<=11 else None
        return H,E,O,h(F)
    return h(F),rooted

def phase(n,u):
    a=[0]*n
    def put(i,j,c=1):
        if c:a[i]|=1<<j
        else:a[j]|=1<<i
    for i in (0,3):
        for z in range(3):put(i+z,i+(z+1)%3)
    f=iter(u)
    for d in range(3):
        c=next(f)
        for z in range(3):put(z,3+(z+d)%3,c)
    for i in (0,3):
        for j in range(6,n):
            c=next(f)
            for z in range(3):put(i+z,j,c)
    for i,j in combinations(range(6,n),2):put(i,j,next(f))
    return a

def check():
    examples=[]
    cases=0
    for n in range(6,9):
        for u in product((0,1),repeat=3+2*(n-6)+(n-6)*(n-7)//2):
            a=phase(n,u);_,r=profile(a)
            for k in range(1<<(n-6)):
                vals=[r((k<<6)|1|(1<<(3+d))) for d in range(3)]
                cases+=1
                assert all(v[1:]==vals[0][1:] for v in vals)
                hs=sorted(v[0] for v in vals)
                if hs==list(range(hs[0],hs[0]+5,2)):examples.append((n,u,k,vals,a))
    print("rotational phase controls",cases,"consecutive triples",len(examples))
    for row in examples[:5]:print(row)

if __name__=='__main__':
    A=[74,148,289,22,37,11,190,317,123]
    cyc=lambda i: 3*(i//3)+(i+1)%3
    for i,j in combinations(range(9),2):
        assert (A[i]>>j)&1==(A[cyc(i)]>>cyc(j))&1
    _,r=profile(A)
    expected=[[6875,6815,6923],[7287,7175,7201],[7289,7311,7285]]
    for i,j in product(range(3),repeat=2):
        q=(1<<6)|(1<<i)|(56^(1<<(3+j)))
        assert r(q)==(expected[i][j],5586,6357,1521)
    for q,D in [(92,0),(114,2),(116,4)]:
        T=[x|((1<<9) if not q>>i&1 else 0) for i,x in enumerate(A)]+[q]
        for m in range(1,4):
            k=len(T)
            T=[x|((1<<k) if i!=9 else 0) for i,x in enumerate(T)]+[1<<9]
            assert profile(T)[0]==10422*2**m-2387+D
    for n,num_leaves in [(15,2),(16,2),(17,2),(18,3)]:
        d={i:[j for j in range(1,n+1) if i!=j and i+j in (4,9,16,25)] for i in range(1,n+1)}
        ends={i for i,li in d.items() if len(li)==1}
        assert len(ends)==num_leaves
        if n==15:assert ends=={8,9}
    print("nine exact phase profiles, three padding heights, square endpoints: PASS")
    check()
