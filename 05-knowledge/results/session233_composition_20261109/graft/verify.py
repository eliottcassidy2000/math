import json, math, hashlib
from pathlib import Path
P=Path(__file__).resolve().parent

def check(c,msg):
    if not c: raise RuntimeError(msg)

def mul(x,y):
    h,q,A,G,T,S,E=x; j,r,B,F,U,V,W=y
    return (h*j,q*r,A*r+h*B,F*q+j*G,A*F+T*j+h*U+E*V,S*j,h*W)

def single(b): return (1,2,b,b,0,b,b)
def chain(xs):
    xs=iter(xs); out=next(xs)
    for x in xs: out=mul(out,x)
    return out

def word(w): return chain(single(int(b)) for b in w)
def repeat(x,n):
    out=None
    while n:
        if n&1: out=x if out is None else mul(out,x)
        x=mul(x,x); n//=2
    return out

def literal(adj,cut):
    n=len(adj); N=1<<n
    ends=[[0]*n for _ in range(N)]
    starts=[[0]*n for _ in range(N)]
    for v in range(n): ends[1<<v][v]=starts[1<<v][v]=1
    for s in range(1,N):
        for v in range(n):
            if not s>>v&1: continue
            r=s^(1<<v)
            if not r: continue
            ends[s][v]=sum(ends[r][u] for u in range(n) if r>>u&1 and adj[u]>>v&1)
            starts[s][v]=sum(starts[r][u] for u in range(n) if r>>u&1 and adj[v]>>u&1)
    hs=[sum(z) for z in ends]; hs[0]=1
    ss=[sum(z[v] for v in range(n) if cut>>v&1) for z in starts]
    es=[sum(z[v] for v in range(n) if cut>>v&1) for z in ends]
    full=N-1
    q=sum(hs[s]*hs[full^s] for s in range(N))
    A=sum(hs[s]*ss[full^s] for s in range(N))
    G=sum(es[s]*hs[full^s] for s in range(N))
    T=sum(es[s]*ss[full^s] for s in range(N))
    return (hs[-1],q,A,G,T,ss[-1],es[-1])

cores={h:json.loads((P/f'core{h}.json').read_text()) for h in (37,45)}
for h,c in cores.items():
    check(literal(c['adj'],c['mask'])==tuple(c['core']),f'literal core {h}')
    c['Xstate']=word(c['X']); c['Ystate']=word(c['Y'])
    end=0 if h==37 else 1
    left=chain([single(0),c['core'],c['Ystate'],single(end)])
    right=chain([single(0),c['Xstate'],c['core'],single(end)])
    expected=[0]*7; expected[4]=-2 if h==37 else 2
    check(tuple(b-a for a,b in zip(left,right))==tuple(expected),f'local identity {h}')
    check(c['X'][-1]=='0' and int(c['Y'][0])==end,'flanks')

def family(h,k,j):
    c=cores[h]; xs=[single(1),single(0)]
    if j: xs.append(repeat(c['Xstate'],j))
    xs.append(c['core'])
    if k-j: xs.append(repeat(c['Ystate'],k-j))
    xs.append(single(0 if h==37 else 1))
    return chain(xs)

def H(x): return 2*x[2]+x[0]-x[4]
def combined(K,L,i,j): return H(mul(family(37,K,i),family(45,L,j)))

def selector(a,b,K,L,t):
    # t=a*i+b*j. Choose j's residue modulo a, then shift by a.
    r=(t*pow(b,-1,a))%a
    i0=(t-b*r)//a
    ell=max(0,(i0-K+b-1)//b)
    i=i0-b*ell; j=r+a*ell
    check(0<=i<=K and 0<=j<=L and a*i+b*j==t,'selector')
    return i,j

for h in (37,45):
    for k in (0,1,2,4,8,20,100):
        base=family(h,k,0)
        for j in sorted(set([0,k//2,k])):
            curr=family(h,k,j)
            expect=[0]*7; expect[4]=(-2 if h==37 else 2)*j
            check(tuple(b-a for a,b in zip(base,curr))==tuple(expect),'family')

coin_cases=0
for a in range(2,18):
    for b in range(2,18):
        if math.gcd(a,b)!=1: continue
        for K,L in [(b-1,a-1),(b+2,a+3),(2*b,3*a)]:
            C=(a-1)*(b-1); M=a*K+b*L
            attainable={a*i+b*j for i in range(K+1) for j in range(L+1)}
            for t in range(C,M-C+1):
                check(t in attainable,'brute coin coverage')
                selector(a,b,K,L,t)
            coin_cases+=1

records=[]
for K,L in [(36,44),(36,50),(36,100),(37,44),(40,50)]:
    base=combined(K,L,0,L); C=1584; M=45*K+37*L
    for t in range(C,M-C+1):
        i,j=selector(45,37,K,L,t)
        # Full graph's second family index is L-j.
        # Recompute actual composed states at representative targets.
        if t in {C,C+1,(M//2),M-C-1,M-C}:
            check(combined(K,L,i,L-j)==base+2*t,'physical selected count')
    for i,j in [(0,0),(K,0),(0,L),(K,L),(K//2,L//2)]:
        check(combined(K,L,i,L-j)==base+2*(45*i+37*j),'rectangle')
    lo=base+2*C; hi=base+2*(M-C)
    records.append(dict(K=K,L=L,vertices=1188*K+12*L+23,consecutive_odd_values=M-2*C+1,
                        lower_bit_length=lo.bit_length(),upper_bit_length=hi.bit_length(),
                        lower_sha256=hashlib.sha256(lo.to_bytes((lo.bit_length()+7)//8,'big')).hexdigest()))
# Hostiles: conductor-1 is absent, as is its reflected endpoint; identical h's do not interlace.
check(1583 not in {45*i+37*j for i in range(37) for j in range(45)},'Frobenius hostile')
check(all((45*i+45*j)%45==0 for i in range(8) for j in range(8)),'same-core hostile')
print(json.dumps(dict(status='PASS',literal_cores=[37,45],coin_cases=coin_cases,records=records),indent=2))
