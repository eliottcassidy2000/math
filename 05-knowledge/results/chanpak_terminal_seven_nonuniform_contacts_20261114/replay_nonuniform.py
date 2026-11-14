from itertools import combinations, product
from functools import cache
from math import isqrt

def trans(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def literal(blocks,bits):
    n=sum(map(len,blocks))
    rows=list(trans(n+1))
    t=0
    for B,b in zip(blocks,bits):
        for i,r in enumerate(B):
            rows[t+i]=(r<<t)|sum(1<<j for j in range(t+len(B),n+1))
            if b:
                rows[t+i]^=1<<n
                rows[n]|=1<<(t+i)
        t+=len(B)
    return tuple(rows)

@cache
def FH(rows,m):
    n=len(rows)
    if not m:return int(not n)
    @cache
    def end(S,v):
        t=S^(1<<v)
        return 1 if not t else sum(end(t,u) for u in range(n) if t>>u&1 and rows[u]>>v&1)
    @cache
    def H(S):
        return 1 if not S else sum(end(S,u) for u in range(n) if S>>u&1)
    @cache
    def color(S,k):
        if k==1:return H(S)
        v=S;res=0
        while True:
            res+=H(v)*color(S^v,k-1)
            if not v:return res
            v=(v-1)&S
    return color((1<<n)-1,m)

def prod(a):
    out=1
    for x in a:out*=x
    return out

def formula(B,bits,m):
    h=[FH(T,m) for T in B]
    p=[FH(T,m+1) for T in B]
    r=[FH(T,m-1) for T in B]
    out=prod(h)
    for i,j in combinations(range(len(B)),2):
        if bits[i] and not bits[j]:
            out+=prod(h[:i]+h[j+1:])*((p[i]-h[i])*(p[j]-h[j])*prod(p[i+1:j])+
                                     (h[i]-r[i])*(h[j]-r[j])*prod(r[i+1:j]))
    return m*out

def transfer(B,bits):
    P=K=1; E=a=0
    for X,b in zip(B,bits):
        h,ell=FH(X,1),FH(X,2)
        d=ell-h
        P,E,K,a=h*P,ell*E+b*d*P,h*K+(1-b)*(d*E+a*h*P),b
    return P,E,K,a

def seam(U,V,b,c):
    a,eu,ku,x=transfer(U,b)
    d,_,kv,_=transfer(V,c)
    VR=[tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T))) for T in V[::-1]]
    _,gv,_,_=transfer(VR,tuple(1-u for u in c[::-1]))
    return ku*d+a*kv-a*d+eu*gv+x*(1-c[0])*a*d

def inverse(N,h,l,q,L):
    t=2*N-3*h*q
    if t<=0:return []
    out=[]
    for a in range(1,isqrt(t)+1):
        if t%a:continue
        for f,g in {(a,t//a),(t//a,a)}:
            x,y=f-h,g+q
            if x<=0 or x%(2*l) or y%(2*L):continue
            x//=2*l;y//=2*L
            if (x&(x-1)) or (y&(y-1)):continue
            out.append((x.bit_length()-1,y.bit_length()-1))
    return sorted(out)

def main():
    tests=0
    bank=[trans(1),trans(2),(2,4,1),(6,12,8,1)]
    for k in range(1,5):
        for B in product(bank,repeat=k):
            if sum(map(len,B))>7:continue
            for bits in product(range(2),repeat=k):
                rows=literal(B,bits)
                for m in range(1,4):assert FH(rows,m)==formula(B,bits,m)
                assert transfer(B,bits)[2]==FH(rows,1)
                for j in range(1,k):
                    assert seam(B[:j],B[j:],bits[:j],bits[j:])==FH(rows,1)
                tests+=1
    assert (FH(bank[-1],1),FH(bank[-1],2))==(5,32)
    def val(r,s):return (1<<(r+s+6))+5*(1<<s)-(1<<(r+5))+5
    for r,s in [(1,4),(3,6),(5,4)]:
        N=val(r,s)
        assert (r,s) in inverse(N,5,32,1,1)
        if N%7==0:N//=7
        assert all(N%p for p in range(2,isqrt(N)+1))
    assert val(1,4)==2069 and val(3,6)==7*4691
    for r in range(6):
        for s in range(6):
            assert (val(r,s)%7==0 and val(r,s)%3!=0)==((r,s) in ((3,0),(5,4)))
    print(tests,'literal colored nonuniform apex tests, all seams; factorization, inverse and prime witnesses verified')
if __name__=='__main__':
    main()
