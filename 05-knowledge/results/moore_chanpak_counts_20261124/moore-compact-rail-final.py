from itertools import combinations,permutations,product
from math import comb
from collections import Counter
from functools import cache

def rows(A):
    return tuple(sum(int(v)<<i for i,v in enumerate(a)) for a in A)

def read(G,O,U):
    m=len(O);k=len(U);n=2*m+k
    r=[0]*n
    def put(i,j,v=True):
        if v:r[i]|=1<<j
        else:r[j]|=1<<i
    W=list(reversed(range(m+k,n)))
    for L in (list(range(m)),list(range(m,m+k)),W[::-1]):
        for i,j in combinations(L,2):put(i,j)
    for i,v in enumerate(O):
        for j,u in enumerate(O):put(i,W[j],v==u or G[v][u])
        for j,u in enumerate(U):
            put(m+j,W[i],G[u][v])
            put(i,m+k-1-j,G[u][v])
    if not all(r[i]>>(i+1)&1 for i in range(n-1)):return
    return tuple(r)

@cache
def heights(r):
    n=len(r)
    def run(cycle):
        F=[[0]*n for _ in range(1<<n)]
        for i in range(1 if cycle else n):F[1<<i][i]=1
        for s,q in enumerate(F):
            for i,x in enumerate(q):
                if not x:continue
                avail=r[i]&~s
                while avail:
                    b=avail&-avail;avail^=b
                    F[s|b][b.bit_length()-1]+=x
        return sum(x for i,x in enumerate(F[-1]) if not cycle or r[i]&1)
    return run(False),run(True)

def blue(r):
    n=len(r);m=n//2
    def go(p,S):
        if len(p)==m:
            if n%2:
                if p and not (r[p[-1]]>>m&1):return
            elif not (r[p[-1]]>>(n-1-p[-1])&1):return
            yield p
            return
        for j in all_but_middle:
            q=(1<<j)|(1<<(n-1-j))
            if not(S&q) and (not p or r[p[-1]]>>j&1):yield from go(p+(j,),S|q)
    all_but_middle=list(range(n))
    if n%2:all_but_middle.remove(m)
    c=Counter(bool(r[n-1-p[0]]>>p[0]&1) for p in go((),0))
    return c[False],c[True]

@cache
def sig(r,v):
    n=len(r);C=[0]*(n+1)
    for p in permutations(i for i in range(n) if i!=v):
        forced=opt=0
        if r[v]>>p[0]&1:opt+=1
        if r[p[-1]]>>v&1:opt+=1
        for a,b in zip(p,p[1:]):
            arc=r[a]>>b&1
            good=(r[a]>>v&1) and (r[v]>>b&1)
            if good:
                if arc:opt+=1
                else:forced+=1
            elif not arc:break
        else:
            for z in range(opt+1):
                if z+forced:C[z+forced]+=comb(opt,z)
    return tuple(C)

def flip(r,i,j):
    q=list(r);q[i]^=1<<j;q[j]^=1<<i;return tuple(q)

def canon(r):
    n=len(r)
    D=[tuple(i for i,x in enumerate(r) if x.bit_count()==d)
       for d in sorted({x.bit_count() for x in r})]
    return min(sum((r[i]>>j&1)<<k for k,(i,j) in enumerate(combinations(sum(p,()),2)))
               for p in product(*(tuple(permutations(x)) for x in D)))

def sc(r):
    s={0}
    for _ in r:s|={j for i in s for j in range(len(r)) if r[i]>>j&1}
    t={0}
    for _ in r:t|={j for i in t for j in range(len(r)) if r[j]>>i&1}
    return len(s)==len(t)==len(r)


def aut(r):
    n=len(r)
    D=[tuple(i for i,x in enumerate(r) if x.bit_count()==d)
       for d in sorted({x.bit_count() for x in r})]
    q=sum(D,())
    return sum(all((r[i]>>j&1)==(r[v[i]]>>v[j]&1)
                   for i,j in combinations(range(n),2))
               for p in product(*(tuple(permutations(x)) for x in D))
               for v in [dict(zip(q,sum(p,())))])

def inflate(r,v,k):
    n=len(r)
    cells=[list(range(sum(k if t==v else 1 for t in range(i)),sum(k if t==v else 1 for t in range(i+1)))) for i in range(n)]
    Q=[0]*sum(map(len,cells))
    for i,cell in enumerate(cells):
        s=sum(1<<x for j,c in enumerate(cells) if r[i]>>j&1 for x in c)
        for t,u in enumerate(cell):Q[u]=s|sum(1<<x for x in cell[t+1:])
    return tuple(Q)

def main():
    G=[[(i-j)%5 in (1,4) for j in range(5)] for i in range(5)]
    C={}
    total=0
    for p in range(5):
        U=[p]+[q for q in range(5) if G[p][q]]
        O=sorted(set(range(5))-set(U))
        for a in permutations(O):
            for b in permutations(U):
                r=read(G,a,b)
                if r is not None:
                    total+=1
                    C.setdefault(canon(r),(a,b,r))
    assert total==20 and len(C)==2
    R=(106,100,25,82,35,76,20);Z=(114,100,89,19,34,76,24)
    assert {canon(R),canon(Z)}==set(C)
    E=list(combinations(range(7),2))
    signatures=[(109,59),(133,126,32),(133,186,68)]
    expected={R:[4,12,6],Z:[0,8,4]}
    for i,r in enumerate((R,Z)):
        assert sc(r) and aut(r)==1 and heights(r)==((147,12),(135,12))[i]
        assert blue(r)==(5,4)
        def rails(t):
            for p in range(7):
                a=sig(t,p)
                if not any(a[2:]):continue
                for f in E:
                    if p in f:continue
                    u=flip(t,*f)
                    b=sig(u,p)
                    if sc(u) and abs(b[1]-a[1])==2 and a[2:]==b[2:]:
                        yield tuple([min(a[1],b[1])]+[x for x in a[2:] if x])
        assert not list(rails(r))
        found=Counter(x for e in E if sc(t:=flip(r,*e)) for x in rails(t))
        assert found==dict((s,v) for s,v in zip(signatures,expected[r]) if v)
        print('compressed pentagon',heights(r),blue(r),found,flush=True)

    f=lambda k:(59*2**k-9,32*3**k+30*2**k-23,68*3**k-18*2**k-35)
    for idx,(p,e,g) in enumerate([(0,(0,4),(3,4)),(0,(1,2),(2,4)),(6,(1,2),(2,4))]):
        t=flip(R,*e);u=flip(t,*g)
        assert sig(t,p)[1:4]==signatures[idx]+(0,)*(3-len(signatures[idx]))
        for k in (1,2,3):
            assert heights(inflate(t,p,k))[0]==f(k)[idx]
            assert heights(inflate(u,p,k))[0]==f(k)[idx]+2
    M=15015
    assert pow(2,60,M)==1 and pow(3,60,M//3)==1
    lower=[{f(k)[i]%M for k in range(1,61)} for i in range(3)]
    L=set.union(*lower);U={(x+2)%M for x in L}
    assert list(map(len,lower))==[60,48,44] and len(L)==151 and not L&U
    print('phase bank',len(L),len(L|U),flush=True)

    G=[list(map(lambda w:not(set(v)&set(w)),combinations(range(5),2)))
       for v in combinations(range(5),2)]
    O=[2,3,4,5,6,1];U=[0,7,8,9]
    r=read(G,O,U)
    assert r==(43198,28988,57464,21680,36128,7232,923,19254,37421,9279,63943,62286,58069,49656,33649,746)
    assert heights(r)==(834157149,27494558) and blue(r)==(19295,26180) and aut(r)==1
    print('compressed Petersen',heights(r),blue(r),flush=True)

if __name__=='__main__':main()
