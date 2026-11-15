from functools import lru_cache
from math import floor, gcd, lgamma, log, sqrt, isqrt

B,Q,T,S,h=2855,5760,4034,3914,45

def req(b,*args):
    if not b:
        raise AssertionError(args)

@lru_cache(None)
def p(x):
    if not x:
        return 0
    a,e=divmod(x,2)
    return p(a)+(a+(a&1) if not e else 0)

def f(x,y):
    return B+Q*x*y+T*x+S*y+h*(p(x)+p(y))

@lru_cache(None)
def w(s,e):
    if s<0:return ()
    if s==0:return () if e else ((0,0),)
    t,r=divmod(s,2)
    if not r and e:
        return tuple((2*a+1,j) for prev in (0,1) for a,j in w(t-1,prev))
    out=[]
    for prev in (0,1):
        ypar=(t-prev)&1
        for a,j in w(t,prev):
            add=(t+prev+ypar if not r else (t-a+ypar if e else a+prev))
            out.append((2*a+e,j+add))
    return tuple(out)

for s in range(1024):
    for e in (0,1):
        req(sorted(w(s,e))==[(x,p(x)+p(s-x)) for x in range(e,s+1,2)],s,e)
print('fixed-sum polynomial recurrence: 1024 complete coefficient sets')

for H in (0,28,56,252,28000,28*2**180+112):
    vals=[f(16*H+x,16*H+9-x) for x in range(8)]
    mod28={v%28 for v in vals}
    req(vals[1]-vals[0]==46380,H)
    req({v%4 for v in vals}=={1,3} and {v%7 for v in vals}==set(range(7)),H)
    misses={9,15,17,19,21,27}
    req(mod28==set(range(1,28,2))-{(r+90*p(H))%28 for r in misses},H)
    for m in misses:
        target=vals[0]+((m+90*p(H)-vals[0])%28)
        req(vals[0]<target<vals[1] and target%28 not in mod28,H,m)
print('large-height CRT failures: 6 heights, 6 missing classes per window')

# Check double stochasticity of the digit Markov chain as an
# integer statement about indegrees, not a numerical mixing claim.
for M in range(1,101):
    indeg=[0]*(2*M*M)
    for x in range(2*M):
        for c in range(M):
            indeg[(2*x%(2*M))*M+((c+(x+(x&1))//2)%M)]+=1
            indeg[((2*x+1)%(2*M))*M+c]+=1
    req(set(indeg)=={2},M)
print('local digit indegree balance: moduli M=1,...,100')

def row_count(X,x,tx,sy):
    A=Q*x+sy
    Y=X-B-tx*x-h*p(x)
    if Y<0:return 0
    lo,hi=Y//(A+h),Y//A
    return lo+1+sum(A*y+h*p(y)<=Y for y in range(lo+1,hi+1))

def pair_count(X):
    H=isqrt(X//Q)
    if H>=4:
        corner=sum(f(x,y)>X for x in range(H-3,H+1) for y in range(H-3,H+1))
    else:
        corner=sum(f(x,y)>X for x in range(H+1) for y in range(H+1))
    return sum(row_count(X,x,T,S)+row_count(X,x,S,T) for x in range(H+1))-(H+1)**2+corner

for X in (5700, 10000, 49999, 215627, 10**6):
    check=sum(f(x,y)<=X for x in range(X//T+1) for y in range(X//(Q*x+S)+1))
    req(pair_count(X)==check,X)
print('independent absolute pair counts: 5 full brute-force cutoffs')
C=((lgamma(S/Q)-lgamma((S+h)/Q)+lgamma(T/Q)-lgamma((T+h)/Q))/h
    -(1+log(Q))/Q)
for k in (18,22,26,30,34):
    X=1<<k
    exact=pair_count(X)
    estimate=X*(log(X)/Q+C)
    print('height',X,'pairs',exact,'main terms',round(estimate,3))
