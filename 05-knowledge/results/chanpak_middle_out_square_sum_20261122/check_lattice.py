from itertools import combinations
from math import gcd, inf
from heapq import heappush, heappop

def flip(rows,i,j):
    rows=list(rows);rows[i]^=1<<j;rows[j]^=1<<i;return tuple(rows)
def add(rows,S):
    n=len(rows);return tuple(row|((1-((S>>i)&1))<<n) for i,row in enumerate(rows))+(S,)
def endpoint(rows):
    n=len(rows); limit=1<<n
    cv=[0]*n
    for i,t in enumerate(rows):
        for j in range(n):
            if t>>j&1:cv[j]|=1<<i
    def all_end(rs):
        T=[[0]*n for _ in range(limit)]
        for i in range(n):T[1<<i][i]=1
        for S in range(1,limit-1):
            for i,v in enumerate(T[S]):
                if not v:continue
                e=rs[i]&~S
                while e:
                    j=(e&-e).bit_length()-1;e&=e-1
                    T[S|1<<j][j]+=v
        return T
    return all_end(rows),all_end(cv)
def h(rows):
    return sum(endpoint(rows)[0][-1])
def pol(rows):
    n=len(rows)
    E,B=endpoint(rows)
    mat=[[0]*n for _ in range(n)]
    full=(1<<n)-1
    for S in range(1,full):
        for i,v in enumerate(E[S]):
            if not v:continue
            for j,w in enumerate(B[full^S]):
                if w:
                    mat[i][j]+=v*w
    lin=[0]*n
    q={}
    for i,j in combinations(range(n),2):
        w=mat[i][j]+mat[j][i]
        assert w%2==0 and w>0
        w//=2
        lin[i if rows[i]>>j&1 else j]+=w
        q[i,j]=q[j,i]=-w
    return sum(E[-1]),lin,q

def pair(a,b):
    ha,A,Q=pol(a);hb,B,R=pol(b)
    return ha,hb,tuple(zip(A,[y-x for x,y in zip(A,B)])),{e:(w,R[e]-w) for e,w in Q.items()}

def egcd(a,b):
    x,y,x1,y1=1,0,0,1
    while b:
        q=a//b;a,b=b,a-q*b
        x,x1=x1,x-q*x1;y,y1=y1,y-q*y1
    return a,x,y

def slice_(generators,delta):
    # Progression of lower-read residuals, or None.
    g=a0=0
    for a,d in generators:
        gg,x,y=egcd(g,abs(d))
        if gg:
            g,a0=gg,x*a0+y*a*(1 if d>=0 else -1)
    if g==0:
        if delta: return None
        return 0,gcd(*(a for a,d in generators))
    if delta%g:return None
    k=gcd(*(a-(d//g)*a0 for a,d in generators))
    return (delta//g)*a0,k

def clusters_bounds(U,lin,q,la,mu,c):
    I={i:la*lin[i][0]+mu*lin[i][1] for i in U}
    Q={e:la*q[e][0]+mu*q[e][1] for e in combinations(U,2)}
    # An exact two-variable box bound for high-leverage disjoint pairs;
    # remaining interaction terms are bounded by signs independently.
    L=H=la*c[0]+mu*c[1]
    free=set(U); used=set()
    for (i,j),w in sorted(Q.items(), key=lambda e:abs(e[1]), reverse=True):
        if i in free and j in free and w:
            free.remove(i);free.remove(j);used.add((i,j))
            t=(0,I[i],I[j],I[i]+I[j]+w)
            L+=min(t);H+=max(t)
    for i in free:L+=min(0,I[i]);H+=max(0,I[i])
    for e,w in Q.items():
        if e not in used:L+=min(0,w);H+=max(0,w)
    return L,H

def nearest(A,D,K,side=0,defect=0):
    n=len(A); U=tuple(range(n));q=D
    nodes=pruned_lat=pruned_int=0
    waits=[]; seq=0
    def envelope(U,lin,c, S):
        nonlocal seq,pruned_lat,pruned_int
        G=[lin[i] for i in U]+[q[i,j] for i,j in combinations(U,2)]
        prog=slice_(G,defect-c[1])
        if prog is None:pruned_lat+=1;return
        a,k=prog;a+=c[0]
        Aa=sum(abs(a) for a,d in G);Dd=sum(abs(d) for a,d in G)
        ratio=max(1,Aa//max(1,Dd))
        low=-inf;high=inf
        for l,m in [(0,1)]+[(1,z) for z in sorted({0,1,-1,ratio,-ratio,2*ratio,-2*ratio})]:
            L,H=clusters_bounds(U,lin,q,l,m,c)
            if l==0:
                if L>defect or H<defect:pruned_int+=1;return
            else:
                low=max(low,L-m*defect);high=min(high,H-m*defect)
        if side<0:high=min(K,high)
        if side>0:low=max(K,low)
        if k:
            i=(low-a+k-1)//k
            j=(high-a)//k
            if i>j:pruned_lat+=1;return
            best=min((a+k*max(i,min(j,(K-a)//k+d)) for d in (0,1)),key=lambda x:abs(K-x))
        else:
            if not low<=a<=high:pruned_lat+=1;return
            best=a
        seq+=1
        heappush(waits,(abs(best-K),len(U),seq,U,lin,c,S,ratio))
    envelope(U,A,(0,0),0)
    while waits:
        d,_,_,U,lin,c,S,ratio=heappop(waits)
        nodes+=1
        if not U: return c[0],S,(nodes,pruned_lat,pruned_int)
        i=max(U,key=lambda i: abs(lin[i][0])+ratio*abs(lin[i][1])+sum(abs(q[i,j][0])+ratio*abs(q[i,j][1]) for j in U if i!=j))
        V=tuple(j for j in U if i!=j)
        envelope(V,lin,c,S)
        lin=list(lin)
        for j in V:lin[j]=tuple(a+b for a,b in zip(lin[j],q[i,j]))
        envelope(V,tuple(lin),tuple(a+b for a,b in zip(c,lin[i])),S|(1<<i))
    return None,None,(nodes,pruned_lat,pruned_int)

def reads(A,Q,S):
    n=len(A);V=[i for i in range(n) if S>>i&1]
    G=[A[i] for i in V]+[Q[i,j] for i,j in combinations(V,2)]
    return sum(a for a,d in G),sum(d for a,d in G)

seeds=[(3,(30,52,56,50,32,1),1,5),(5,(160,145,67,7,13,30,59,124),0,2),(7,(988,433,74,34,236,645,938,780,60,286),8,9)]
def E(p):return (1<<(p-1))*((1<<p)-1)
for p,a,i,j in seeds:
    b=flip(a,i,j);ha,hb,A,Q=pair(a,b); K=(E(p+1)-E(p))//2
    assert ha==E(p)-1 and hb==E(p)+1
    cert=[(c,S) for S in range(1<<len(a)) for c,d in [reads(A,Q,S)] if not d]
    print(f'p={p}, n={len(a)}, target={E(p+1)}, eligible masks={len(cert)}, total={1<<len(a)}')
    for side in (-1,1):
        C,S,stats=nearest(A,Q,K,side)
        feasible=[(c,T) for c,T in cert if side*(c-K)>=0]
        opt=min((abs(c-K) for c,T in feasible),default=None)
        assert opt==(None if C is None else abs(C-K))
        if C is not None:
            assert h(add(a,S))==ha+2*C and h(add(b,S))==hb+2*C
        print('  ',side,'center',None if C is None else ha+1+2*C,'tag',S,'states',stats)

from pathlib import Path
bank=Path(__file__).with_name('midout-audit-deck.rows')
if bank.exists():
    stats={}
    for line in bank.read_text().splitlines():
        p,i,j,*lst=map(int,line.split());a=tuple(lst); b=flip(a,i,j)
        ha,hb,A,Q=pair(a,b);K=(E(p+1)-E(p))//2
        assert ha==E(p)-1 and hb==E(p)+1
        C=[(c,S) for S in range(1<<len(a)) for c,d in [reads(A,Q,S)] if not d]
        for side in (-1,1):
            approx,S,steps=nearest(A,Q,K,side)
            valid=[(c,T) for c,T in C if side*(c-K)>=0]
            opt=min((abs(c-K) for c,T in valid),default=None)
            assert opt==(None if approx is None else abs(approx-K))
            d=stats.setdefault((p,side),[])
            d.append((approx,S,steps,tuple(a),(i,j),len(C)))
    for (p,side), d in sorted(stats.items()):
        exact=min((z for z in d if z[0] is not None),key=lambda z:abs(E(p)+2*z[0]-E(p+1)),default=None)
        print('bank',p,side,'nseams',len(d),'center',None if exact is None else E(p)+2*exact[0],'witness',None if exact is None else (exact[1],exact[3],exact[4]),'pop',sum(z[2][0] for z in d),'eligible',sum(z[5] for z in d))

# Reversed lift also has gap +2 as an unordered pair, so scan it separately.
reverse=[]
if bank.exists():
    for line in bank.read_text().splitlines():
        p,i,j,*lst=map(int,line.split()); reverse.append((p,tuple(lst),i,j))
reverse+=seeds[-1:]
rr={}
for p,a,i,j in reverse:
    b=flip(a,i,j);ha,hb,A,Q=pair(a,b)
    K=(E(p+1)-E(p))//2+1
    for side in (-1,1):
        C,S,count=nearest(A,Q,K,side,defect=-2)
        targets=[c for S0 in range(1<<len(a)) for c,d in [reads(A,Q,S0)] if d==-2 and side*(c-K)>=0]
        assert min((abs(c-K) for c in targets),default=None)==(None if C is None else abs(C-K))
        if C is not None:
            assert h(add(a,S))==ha+2*C and h(add(b,S))==hb+2*C-4
            rr.setdefault((p,side),[]).append((E(p)-2+2*C,S,a,(i,j)))
for p in (3,5,7):
    for side in (-1,1):
        v=rr.get((p,side),[])
        nearest_reverse=min(v,key=lambda z:abs(E(p+1)-z[0]),default=None)
        print('reverse',p,side,nearest_reverse)

# The 3-apex exhaustive replay uses one marked orbit per 6-vertex class.
if bank.exists():
    from itertools import permutations
    six={}
    for ln in bank.read_text().splitlines():
        p,i,j,*a=map(int,ln.split())
        if p==3:six.setdefault(tuple(a),set()).add((i,j))
    assert len(six)==2
    for a,marks in six.items():
        s=next(iter(marks))
        orbit=set()
        for pi in permutations(range(6)):
            if all((a[u]>>v&1)==(a[pi[u]]>>pi[v]&1) for u in range(6) for v in range(6)):
                orbit.add(tuple(sorted((pi[s[0]],pi[s[1]]))))
        assert orbit==marks
    print('Two marked six-vertex orbits verified directly.')
