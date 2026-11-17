from itertools import product, combinations, permutations
from functools import lru_cache

def rows_ord(parts):
    n=sum(map(len,parts))
    a=[];i=0
    for p in parts:
        a.extend([x<<i | ((1<<n)-(1<<(i+len(p)))) for x in p])
        i+=len(p)
    return tuple(a)

def arrays(D):
    n=len(D);N=1<<n;T=N-1
    cnt=[[0]*(n*n) for _ in range(N)]
    H=[1]+[0]*(N-1)
    for i in range(n):cnt[1<<i][i*n+i]=1
    for s in range(1,N):
        a=cnt[s]
        H[s]=sum(a)
        for j in range(n):
            out=D[j]&~s
            if not out:continue
            for i in range(n):
                h=a[i*n+j]
                if not h:continue
                o=out
                while o:
                    b=o&-o;o-=b
                    cnt[s|b][i*n+b.bit_length()-1]+=h
    K=[0]*(n*n);V=[0]*n
    for s in range(1,N):
        a=cnt[s];h=H[T^s]
        if s.bit_count()&1:
            for i in range(n):V[i]+=h*sum(a[i*n:i*n+n])
        else:
            for i,c in enumerate(a):K[i]+=c*h
    return H[T],K,V

def literal(D,U):
    n=len(D)
    a=[x+((not (U>>i&1))<<n) for i,x in enumerate(D)]+[U]
    N=1<<len(a);H=[]
    for only_root in (False,True):
        d=[[0]*len(a) for _ in range(N)]
        for i in ([n] if only_root else range(len(a))):d[1<<i][i]=1
        for s in range(1,N):
            for i,h in enumerate(d[s]):
                if not h: continue
                o=a[i]&~s
                while o:
                    b=o&-o;o-=b;d[s|b][b.bit_length()-1]+=h
        H.append([1]+[sum(p) for p in d[1:]])
    h,p=H
    v=[sum(x*h[(N-1)^s] for s,x in enumerate(p) if s and s.bit_count()%2==k) for k in (0,1)]
    return (h[-1],h[(1<<n)-1],*v)

def rays(g, C, cap):
    results={}
    tests=0
    for l in range(cap+1):
        for r in range(cap-l+1):
            D=rows_ord([(0,)]*l+[C]+[(0,)]*r)
            G,K,V=arrays(D)
            assert G==g,(D,G)
            n=len(D);N=1<<n;rows=[sum(K[i*n:i*n+n]) for i in range(n)]
            for U in range(N):
                tests+=1
                ids=[i for i in range(n) if U>>i&1]
                k=sum(K[i*n+j] for i in ids for j in ids)
                j=sum(rows[i] for i in ids)-k
                E=sum(V[i] for i in ids)
                O=g+j+k
                if not (k==g and O-E==g):continue
                h=g+2*j
                odd=h+g;t=0
                while odd%2==0:odd//=2;t+=1
                if odd not in results or results[odd][0]>t:results[odd]=t,h,D,U
    print('g',g,'tested',tests,'rays',sorted((a,t,h) for a,(t,h,_,_) in results.items()),flush=True)
    for a,(t,h,D,U) in sorted(results.items()):
        assert literal(D,U)==(h,g,(h+g)//2,(h+3*g)//2)
        print(a,t,h,D,U,flush=True)
    return results

if __name__=='__main__':
    c=(2,4,9,3)
    known={tuple(sum(((c[i]>>j)&1)<<b for b,j in enumerate(p)) for i in p) for p in permutations(range(4))}
    for bits in range(1<<6):
        T=[0]*4
        for b,(i,j) in enumerate(combinations(range(4),2)):
            if not bits>>b&1:i,j=j,i
            T[i]|=1<<j
        if arrays(T)[0]==5:assert tuple(T) in known
    a=rays(3,(2,4,1),5)
    assert {i:v[0] for i,v in a.items()}=={9:1,15:1,27:1}
    b=rays(5,(2,4,9,3),6)
    assert {i:v[0] for i,v in b.items()}=={3:5,5:5,9:2,15:1,19:2,21:1,23:1,25:1,71:1}
    print('g5 first 7p on each ray (bounded diagnostic)')
    def prime(q):return q>=2 and all(q%d for d in range(2,int(q**.5)+1))
    for a,(k,h,D,U) in sorted(b.items()):
        for j in range(k,min(k+12,22)):
            y=(a<<j)-5
            if y%7==0 and prime(y//7):print(a,j,y,D,U);break
