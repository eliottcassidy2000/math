from functools import cache
from math import isqrt

K=[(1,1,1),(1,2,2),(3,3,4)]
def matrix(q):
    while len(K)<=q:
        K.append(tuple(sum(tr) for tr in zip(*K[-3:])))
    return K[q]

def state(gaps):
    a,b=1,1<<gaps[0]
    q=0
    for r in gaps[1:]:
        if not r:q+=1;continue
        e,f,g=matrix(q)
        a,b=e*a+f*b,(f*a+g*b)<<r
        q=0
    assert q==0
    return a,b

def packing(gaps):
    P=[1]
    cum=[0]
    for r in gaps:cum.append(cum[-1]+r)
    for j,x in enumerate(cum):
        P.append(P[-1]+sum((1<<(x-cum[i]))*P[i] if x>cum[i]
                            else 2*P[i]*int((j-i)%2==0) for i in range(j)))
    return P[-1]

def physical_rows(gaps):
    anchors=[0]
    for g in gaps:anchors.append(anchors[-1]+g+1)
    n=anchors[-1]+1
    back=set(zip(anchors,anchors[1:]))
    return tuple(sum(1<<j for j in range(n)
                     if i!=j and ((i<j)!=((min(i,j),max(i,j)) in back)))
                 for i in range(n))

def literal(gaps):
    rows=physical_rows(gaps)
    n=len(rows)
    @cache
    def ending(S,v):
        if S==1<<v:return 1
        T=S^(1<<v)
        return sum(ending(T,u) for u in range(n) if T>>u&1 and rows[u]>>v&1)
    return sum(ending((1<<n)-1,v) for v in range(n))

def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def inverse(N,A,B):
    if N%7:return None
    t=N//7-A
    if t%B:return None
    u=t//B
    return u.bit_length()-1 if u>=2 and u&(u-1)==0 else None

def qword(r,q,i,j): return (r,)+(0,)*q+(1,)+(0,)*i+(1,)+(0,)*j+(1,)+(0,)*7+(2,)
PREFS=((1,), (2,), (1,0,1), (2,0,0,5), (1,0,0,3,0,0,1))
COEFF={10:((3805977,4534839),(3805975,4534839)),
       14:((43557897,51898713),(43557899,51898713))}
if __name__=='__main__':
    for prefix in PREFS:
        a,b=state(prefix)
        for q,coef in COEFF.items():
            outputs=[]
            for (i,j),(A,B) in zip(((2,1),(1,2)),coef):
                gaps=prefix+qword(1,q,i,j)[1:]
                outputs.append(packing(gaps))
                assert outputs[-1]==7*(A*a+B*b)
            assert outputs[0]-outputs[1]==(14 if q==10 else -14)*a
    for q,r,p in ((10,3,40084687),(14,8,13329628427)):
        word=qword(r,q,1,2)
        assert prime(p) and packing(word)==7*p and inverse(7*p,*COEFF[q][1])==r
        print(f'q={q}, r={r}, order {len(physical_rows(word))}: H=7*{p}={7*p}')
    for r in range(1,8):
        for s in range(1,8):
            w1=(1,r,0,1,0,s)
            w2=(r,0,1,1,0,s)
            N=(6*(1<<r)+5)*(6*(1<<s)+5)
            assert (packing(w1),packing(w2))==(N+2,N)
            if r+s<=5: assert (literal(w1),literal(w2))==(N+2,N)
    print('arbitrary physical prefixes; elastic cofactor and exact bivariate literal checks passed')
