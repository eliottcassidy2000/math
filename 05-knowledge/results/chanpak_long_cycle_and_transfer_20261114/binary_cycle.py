from fractions import Fraction as F
from functools import cache
from math import factorial

def create(K,odd=True,spine_bits=None):
    assert K>=1
    plus=spine_bits is not None and K==1<<spine_bits
    d=[int(c) for c in bin(K-int(plus))[2:]]
    if spine_bits is not None:
        assert spine_bits>=1 and len(d)<=spine_bits
        d=[0]*(spine_bits-len(d))+d
    elif odd and not len(d)&1:d.insert(0,0)
    n=len(d)
    a,b,u,v=range(n,n+4)
    T=[sum(1<<j for j in range(i+1,n))|(1<<a)|(1<<b) for i in range(n)]
    T+=[1<<v,(1<<a)|(1<<u)|(1<<v),((1<<n)-1)|(1<<a)|(1<<v),0]
    for i,s in enumerate(d):
        if s:T[v]|=1<<i
        else:T[i]|=1<<v
    if plus:
        T[b]^=1<<v
        T[v]|=1<<b
    return tuple(T),u

def literal(T):
    n=len(T)
    @cache
    def end(m,u):
        other=m^(1<<u)
        if not other:return 1
        return sum(end(other,j) for j in range(n) if other>>j&1 and T[j]>>u&1)
    return [1]+[sum(end(m,u) for u in range(n) if m>>u&1) for m in range(1,1<<n)]

def cycles(T,r):
    n=len(T)
    @cache
    def walk(m,last):
        if not m:return (T[last]>>r)&1
        m2=m&T[last]
        ans=0
        while m2:
            p=m2&-m2
            m2^=p
            ans+=walk(m^p,p.bit_length()-1)
        return ans
    return walk(((1<<n)-1)^(1<<r),r)

def colors(T,upto):
    counts=literal(T)
    f=[1]+[0]*(len(counts)-1)
    out=[0]
    for _ in range(upto):
        g=[0]*len(f)
        for p in range(len(f)):
            s=p
            while True:
                g[p]+=counts[s]*f[p^s]
                if not s:break
                s=(s-1)&p
        f=g
        out.append(f[-1])
    return out

def round_read(vals):
    M=len(vals)-1
    return sum(F(2*(-1)**(i-1)*factorial(M)**2,
                 i*factorial(M-i)*factorial(M+i))*vals[i] for i in range(1,M+1))

def wrap_transfer(row):
    data={}
    for m,c in row.items():
        if not m:continue
        for j,t in ((m,c*m*m),(m+1,c*m),(m-1,-c*m)):
            if j:data[j]=data.get(j,0)+t
    return data

if __name__=='__main__':
    for K in range(1,81):
        T,r=create(K)
        assert len(T)%2==1 and cycles(T,r)==K
        _,r2=create(K,False)
        assert cycles(create(K,False)[0],r2)==K
    for n in range(1,6):
        P,u=create(1<<n,spine_bits=n)
        assert len(P)==n+4 and cycles(P,u)==1<<n
    for K in (1,2,3,5,9,13,26):
        T,r=create(K)
        L=(len(T)+1)//2
        Fm=colors(T,L)
        assert round_read(Fm)==2*K
        row={1:1}
        W=[]
        for j in range(L):
            assert row[j+1]==factorial(j)
            W.append(sum(c*Fm[m] for m,c in row.items()))
            row=wrap_transfer(row)
        # Pad only twice for independent literal verification at small K.
        if K<=3:
            Q=T
            for j in range(min(3,L)):
                assert literal(Q)[-1]==W[j]
                n=len(Q)
                Q=(sum(1<<q for q in range(1,n+1)),)+tuple((v<<1)|(1<<(n+1)) for v in Q)+(1,)
    print('binary cycle targets 1..80; exact signed color and physical wrap readouts')
