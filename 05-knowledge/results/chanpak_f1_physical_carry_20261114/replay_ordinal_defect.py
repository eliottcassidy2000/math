"""Literal support checks for the nonnegative ordinal-defect identity."""
from itertools import product, combinations

def paths(rows):
    n=len(rows); full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)]
    for i in range(n): end[1<<i][i]=1
    for S in range(1,full+1):
        for j in range(n):
            if S>>j&1:
                end[S][j]+=sum(end[S^(1<<j)][i] for i,r in enumerate(rows) if r>>j&1)
    return [1]+list(map(sum,end[1:])),end

def decorate(rows,U):
    n=len(rows)
    return tuple(r | ((not(U>>i&1))<<n) for i,r in enumerate(rows))+(U,)

def profile(rows,U):
    n=len(rows); full=(1<<n)-1
    H,end=paths(rows)
    _,start=paths(tuple(full ^ (1<<i) ^ r for i,r in enumerate(rows)))
    h=H[-1]; P=M=E=0
    for A in range(1,full+1):
        left=sum(end[A][i] for i in range(n) if U>>i&1)
        right=sum(start[A][i] for i in range(n) if U>>i&1)
        P+=right*H[full^A]
        M+=left*H[full^A]
        E+=right*sum(end[full^A][i] for i in range(n) if U>>i&1)
    S=sum(end[-1][i] for i in range(n) if U>>i&1)
    T=sum(start[-1][i] for i in range(n) if not(U>>i&1))
    root=1<<n
    K=paths(decorate(rows,U))[0][-1]
    assert h+P-K==E+S and M>=2**U.bit_count()-1
    return h,P,E,M,S,T

def join(C,D):
    a=len(C); b=len(D)
    return tuple(r|(((1<<b)-1)<<a) for r in C)+tuple(r<<a for r in D)

def trans(word):
    n=len(word)
    return tuple(((1<<n)-1)^((2<<i)-1) for i in range(n)),sum((c=='1')<<i for i,c in enumerate(word))

def cases():
    for n in range(1,4):
        pairs=list(combinations(range(n),2))
        for bits in product((0,1), repeat=len(pairs)):
            r=[0]*n
            for (i,j),b in zip(pairs,bits):
                if b: r[i]|=1<<j
                else: r[j]|=1<<i
            for U in range(1<<n):
                yield tuple(r),U,profile(r,U)

if __name__=='__main__':
    A=list(cases())
    for (C,U,(h,p,e,m,s,t)),(D,V,(H,P,E,M,S,T)) in product(A,repeat=2):
        j=profile(join(C,D),U|(V<<len(C)))
        assert j[2]==H*e+h*E+p*M+s*(H-T)
        assert j[2]+j[4]>=p*(2**V.bit_count()-1)
        if V==0: assert j[2]+j[4]==H*e
    C,U=trans('100101101011')
    D,V=trans('0')
    assert profile(C,U)[1:5]==(2411,2018,3433,1)
    R=join(C,D); W=U|(V<<len(C))
    h,p,e,m,s,t=profile(R,W)
    assert (h+p,(h+p)-(e+s))==(4823,2805)
    for n in range(1,8):
        C,U=trans('1'+'0'*n);D,V=trans('10')
        h,p,e,m,s,t=profile(join(C,D),U|(V<<len(C)))
        assert e+s==2**n
    print('PASS:',len(A)**2,'literal joins, sharp boundary, and 14-vertex escape')
