"""Session296 exact rooted funnel source compiler + independent literal DP."""
from itertools import product
from functools import lru_cache
import json

SINGLE=(1,1,1,2,1,1)
PROGRAMS={
97:None,293:("1","10001"),881:("11","1","11"),
2645:("111","11","1"),7937:("01","01111101"),23813:("1011","","","","","1"),
71441:("000011010110001",),1691:("1011","101"),3311:("1","","01","","1"),6551:("101","","01101"),
13031:("010101111011",),25991:("011","","0011001"),
51911:("1","","11","","","1","",""),103751:("1","001111","1",""),
207431:("11","","","00011","1"),
}
def require(q):
    if not q: raise ValueError("verification failed")

@lru_cache(None)
def coeff(w,m):
    t=len(w)
    L=(m-1)**t+sum(m**i*(m-1)**(t-i-1) for i,b in enumerate(w) if b=="1")
    D=(m-1)**t+sum((m-1)**i*m**(t-i-1) for i,b in enumerate(w) if b=="0")
    F={(0,0):1}
    for b in w:
        G={}
        for (last,started),v in F.items():
            def add(k,x): G[k]=G.get(k,0)+x
            add((last,started),(m-1)*v)
            add((1 if b=="1" else 2,started),v)
            if started or b=="0": add((last,1),v)
        F=G
    return L,D,sum(v for (last,_),v in F.items() if last!=2)

def advance(S,w):
    h,e,s,p,g,j=S; t=len(w)
    l1,d1,b1=coeff(w,1); l2,d2,b2=coeff(w,2); l3,d3,_=coeff(w,3)
    return (l2*s+b1*h-d1*e,l1*h-(0 if t else 1)*e,2**t*s+d2*h,
            2*(l3*j+b2*p-d2*g),l2*p-g,3**t*j+d3*p)
def tally(program):
    S=SINGLE
    for w in program: S=advance(S,w)
    return S

def birth(T,r):
    n=len(T)
    return tuple(v|((1<<n) if j==r else 0) for j,v in enumerate(T))+(((1<<n)-1)^(1<<r),)
def funnel(T,w):
    n=len(T); r=n+len(w); tail=((1<<r)-(1<<n)); p=n-1
    return tuple(v|tail|((1<<r) if i!=p else 0) for i,v in enumerate(T))+tuple(
        ((1<<r)-(1<<(n+j+1)))|((1<<r) if b=="1" else 0) for j,b in enumerate(w)
    )+((1<<p)|sum(1<<(n+j) for j,b in enumerate(w) if b=="0"),)
def make(program):
    T=(0,)
    for w in program: T=funnel(T,w)
    return birth(T,len(T)-1)

def count(T,root=None,start=None):
    n=len(T); M=(1<<n)-1
    if root is not None:M^=1<<root
    dp=[[0]*n for _ in range(M+1)]
    for j in range(n):
        if (1<<j)&M and (start is None or j==start):dp[1<<j][j]=1
    for S,row in enumerate(dp):
        if S&~M: continue
        for v,w in enumerate(row):
            if not w:continue
            f=T[v]&M&~S
            while f:
                b=f&-f; f-=b; dp[S|b][b.bit_length()-1]+=w
    return sum(dp[M])

def check(T,h,d=None):
    n=len(T); r=n-1; a=r-1
    for i in range(n):
        require(not ((T[i]>>i)&1))
        for j in range(i):require(((T[i]>>j)&1)+((T[j]>>i)&1)==1)
    require([v for v in range(r) if T[v]>>r&1]==[a])
    for i in range(n):
        S=1<<i
        for _ in range(n):
            S|=sum(1<<j for j in range(n) if any((S>>k&1) and (T[k]>>j&1) for k in range(n)))
        require(S==(1<<n)-1)
    v=(count(T),count(T,r),count(T,start=r))
    require(v[0]==h and v[1]-v[2]==1)
    require(d is None or v[1]==d)
    return v

def inverse(target,max_parent):
    bank={1:[(SINGLE,())]}; hit=[]
    for n in range(3,max_parent+1):
        bank[n]=[]
        for m in (1,*range(3,n)):
            t=n-m-1
            words=[""] if t==0 else ["".join(p)+"1" for p in product("01",repeat=t-1)]
            for S,address in bank[m]:
                for w in words:
                    N=advance(S,w); a=address+(w,)
                    bank[n].append((N,a))
                    if N[0]+N[4]-1==target:hit.append((n,N,a))
    return hit,{k:len(v) for k,v in bank.items()}

def main():
    literal={}
    for h,p in PROGRAMS.items():
        if p is None:continue
        T=make(p); S=tally(p); require(S[0]+S[4]-1==h)
        literal[h]={"address":p,"rows":T,"state":check(T,h,S[0]),"core_state":S}
    chain=(293,881,2645,7937,23813,71441)
    maps={}
    for a,b in zip(chain,chain[1:]):
        require(b==3*a+2)
        T=literal[a]["rows"]
        R=birth(birth(T,len(T)-1),len(T)); Q=literal[b]["rows"]
        require(count(R)==b)
        m,n=len(R),len(Q); c=min(m,n)-2
        phi={i:i for i in range(c)}
        phi[m-2]=n-2; phi[m-1]=n-1
        changes=[(i,j) for i in phi for j in phi if i<j and (R[i]>>j&1)!=(Q[phi[i]]>>phi[j]&1)]
        require(all(m-1 not in e for e in changes))
        maps[b]={"raw_rows":R,"preserved_map":phi,"flips":changes,"delete":[i for i in range(m) if i not in phi],
                 "insert":[i for i in range(n) if i not in phi.values()]}
    hit,heads=inverse(881,12)
    require(len(hit)==2 and {n for n,_,_ in hit}=={9,10})
    for n,S,p in hit:
        require(len(make(p))==n+1);check(make(p),881,S[0])
    print(json.dumps({"status":"PASS","literal":literal,"handoffs":maps,"inverse881":hit,"parent_heads":heads},indent=2))

if __name__=="__main__":main()
