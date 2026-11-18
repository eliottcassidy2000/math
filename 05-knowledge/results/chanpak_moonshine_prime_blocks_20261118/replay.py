"""Independent subset checks and exact constructors.  No third-party imports."""
import argparse
from itertools import combinations
from math import isqrt

# row masks D, ports I, coefficient of f(2), coefficient of H.
# The two positive D readers use four ordinal singles on their input.
R = {
    (1,1): ((4,1,2),3,6,1),
    (1,7): ((24,17,3,6,12),5,42,7),
    (-1,1): ((0,1,3,7,15,31,63),3,96,-1),
    (-1,7): ((4,1,42,147,71,27,47,119),28,672,-7)
}
B0 = 32852
CODES=(43134392453,43125956932,43134792502,43122305802,
       35669685332,43135326994,15133096756,43096694850)
CENTER=96*B0

def decode(code, n=9):
    A=[0]*n
    k=n*(n-1)//2-1
    for i,j in combinations(range(n),2):
        if code>>k&1:A[i]|=1<<j
        else:A[j]|=1<<i
        k-=1
    return tuple(A)

def transpose(A):
    return tuple(sum(1<<i for i,r in enumerate(A) if r>>j&1) for j in range(len(A)))

def dp(A):
    n=len(A);f=(1<<n)-1
    assert all(r&~f==0 and not r>>i&1 for i,r in enumerate(A))
    assert all((A[i]>>j&1)+(A[j]>>i&1)==1 for i,j in combinations(range(n),2))
    E=[[0]*n for _ in range(1<<n)]
    for i in range(n):E[1<<i][i]=1
    for s in range(1,f+1):
        for i,u in enumerate(E[s]):
            if not u:continue
            avail=A[i]&~s
            while avail:
                j=avail&-avail; avail^=j
                E[s|j][j.bit_length()-1]+=u
    return E

def reads(A,I):
    E=dp(A);S=dp(transpose(A));f=len(E)-1
    d=sum(E[f]);g=q=0
    for a in range(1,f+1):
        z=sum(v for i,v in enumerate(E[a]) if I>>i&1); b=f^a
        if not b:g+=z;q+=z
        else:
            g+=z*sum(S[b]);q+=z*sum(v for i,v in enumerate(S[b]) if not I>>i&1)
    return d,g,q

def profile(A):
    E=dp(A);h=[1]+[sum(e) for e in E[1:]];f=len(E)-1
    return h[f],sum(h[s]*h[f^s] for s in range(f+1))

def pad(A,k):
    n=len(A);f=(1<<(n+k))-1
    return tuple(r|(f^((1<<n)-1)) for r in A)+tuple(f^((1<<(i+1))-1) for i in range(n,n+k))

def wrap(D,I,A):
    n=len(D);e=n+len(A)
    return (tuple(r|((I>>i&1)<<e) for i,r in enumerate(D))
            +tuple((r<<n)|((1<<n)-1) for r in A)+(((1<<e)-1)^I,))

def select(N,scale=1):
    p,r=divmod(N,scale)
    if r or scale not in (1,7) or p%2!=1 or p%3==0: return None
    # The unique unit residue chooses an offset before height is tested.
    for s in (1,-1):
        for j,k in enumerate(CODES):
            h=1945+6*j
            if (p-s*h)%48:continue
            t,r=divmod(p-s*h,CENTER)
            if r or t<=0 or t&(t-1):return None
            D,I,_,_=R[s,scale]
            return wrap(D,I,pad(decode(k),t.bit_length()-1+4*(s==1)))
    raise AssertionError('uncovered unit of Z/48')

def prime(p):
    return p>=2 and all(p%i for i in range(2,isqrt(p)+1))

def circ(n,steps):
    return tuple(sum(1<<((i+s)%n) for s in steps) for i in range(n))

def triangle(A,B):
    a=len(A);b=len(B);c=a+b
    return (tuple(r|(((1<<b)-1)<<a) for r in A)
            +tuple((r<<a)|(1<<c) for r in B)+((1<<a)-1,))

def insert(Q, port, A):
    sizes=[len(A) if i==port else 1 for i in range(len(Q))]
    choices=[]; n=0
    for k in sizes:
        choices.append(tuple(range(n,n+k)));n+=k
    out=[]
    for i,X in enumerate(choices):
        outside=sum(1<<y for j,Y in enumerate(choices) if Q[i]>>j&1 for y in Y)
        out.extend(((A[k]<<X[0]) if i==port else 0)|outside for k in range(len(X)))
    return tuple(out)

def quotient_by_first(A,O,v):
    E=dp(A); S=dp(transpose(A)); f=len(E)-1;w=1<<v
    o=sum(1<<u for u in O);F=f^o
    return sum(E[s|w][v]*S[f^s][v] for s in range(f+1) if s&~F==0)

def verify():
    for (s,q),(D,I,g,beta) in R.items():
        d,c,u=reads(D,I)
        assert (c,d+u-c)==(g,beta)
        assert profile(wrap(D,I,(0,)))[0]==2*g+beta
        print('reader',s,q,(d,g,u),beta)
    for i,k in enumerate(CODES): assert profile(decode(k))==(1945+6*i,B0)
    off={s*(1945+6*i)%48 for s in (1,-1) for i in range(8)}
    assert off=={p for p in range(48) if p%2 and p%3}
    print('eight exact profile heights',B0,'and all 16 unit classes mod 48')
    for k in range(0,12):
        C=CENTER<<k
        for i in range(8):
            for s in (1,-1):
                p=C+s*(1945+6*i)
                for q in (1,7):assert select(q*p,q)
    for p in (3151847,3155773,51671730073):
        assert prime(p) and select(p) and select(7*p,7)
        print('literal exact prime-pair target',p,7*p)
    shapes=set()
    for steps in range(8):
        U=circ(7,[i if steps>>(i-1)&1 else 7-i for i in (1,2,3)])
        shapes.add(profile(U))
    assert shapes=={(175,1988),(189,2016)}
    U=circ(7,(1,2,3));assert profile(U)==(175,1988)
    assert tuple(sum(dp(insert((2,4,9,3),i,U))[-1]) for i in range(4))==(4151,5789,5789,4151)
    assert all(prime(p) for p in (593,827,265007))
    assert 259*2**10-209==265007
    for n in range(1,6):
        A=pad((),n);T=triangle(U,A)
        h,b=profile(A);p=259*b-209*h
        assert sum(dp(T)[-1])==7*p
        assert quotient_by_first(T,tuple(range(7)),0)==p
    print('carousel and Paley profiles:',sorted(shapes),'first-hit quotient replay passed')
    print('strong four-vertex quotient ports: 4151, 5789, 5789, 4151')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('target',nargs='?',type=int)
    parser.add_argument('--scale',type=int,choices=(1,7),default=1)
    a=parser.parse_args()
    if a.target is None:verify()
    else: print(select(a.target,a.scale))
