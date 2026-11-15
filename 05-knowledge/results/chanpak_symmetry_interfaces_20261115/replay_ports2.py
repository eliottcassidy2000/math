#!/usr/bin/env python3
"""Exact verification of the full-slot transfers and recursive symmetry bank."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb

U=(14,52,56,50,33,1)
V=(6,44,56,49,35,1)
G=(10,28,57,48,33,3)
C3=(2,4,1)

def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def outlet(A,u):
    return tuple(row+((i!=u)<<len(A)) for i,row in enumerate(A))+(1<<u,)

def substitute(A,B):
    # B lists the nonempty modules in quotient order.
    sizes=[len(x) for x in B]; shifts=[sum(sizes[:i]) for i in range(len(B))]
    masks=[((1<<sizes[i])-1)<<shifts[i] for i in range(len(B))]
    result=[]
    for i,C in enumerate(B):
        exterior=sum(masks[j] for j in range(len(B)) if A[i]>>j&1)
        result.extend(exterior+(row<<shifts[i]) for row in C)
    return tuple(result)

def inflate(A,w,W):
    return substitute(A,[W if i==w else (0,) for i in range(len(A))])

def strong(A):
    for u in range(len(A)):
        seen=1<<u
        for _ in A:
            seen_new=seen
            for j,row in enumerate(A):
                if seen>>j&1:seen_new|=row
            if seen_new==seen:break
            seen=seen_new
        if seen!=(1<<len(A))-1:return False
    return True

def indecomposable(A):
    for m in range(1,(1<<len(A))-1):
        if m.bit_count()>=2 and all((r&m) in (0,m) for i,r in enumerate(A) if not m>>i&1):return False
    return True

def aut(A):
    return [p for p in permutations(range(len(A))) if
            all((A[i]>>j&1)==(A[p[i]]>>p[j]&1) for i,j in combinations(range(len(A)),2))]

def all_values(A,steps):
    """Independent DP: actual Hamilton paths, then ordered subset convolution."""
    n=len(A); U=(1<<n)-1
    dp=[[0]*n for _ in range(U+1)]
    for i in range(n):dp[1<<i][i]=1
    for m in range(1,U+1):
        for i,v in enumerate(dp[m]):
            if v:
                a=A[i]&~m
                while a:
                    t=a&-a;a-=t
                    dp[m|t][t.bit_length()-1]+=v
    H=[1]+[sum(row) for row in dp[1:]]
    subs=[]
    for m in range(U+1):
        a=m;b=[m]
        while a:
            a=(a-1)&m;b.append(a)
        subs.append(tuple(b))
    count=[1]+[0]*U; ans=[count[U]]
    for _ in range(steps):
        count=[sum(count[m^s]*H[s] for s in subs[m]) for m in range(U+1)]
        ans.append(count[U])
    return ans

def packing_reader(A,w):
    """Cycle subsets at fixed start, then ordinary odd-packing recursion."""
    n=len(A); U=(1<<n)-1
    cycles=[]
    for v in range(n):
        C=Counter()
        def go(last,m):
            if m.bit_count()>=3 and A[last]>>v&1:C[m]+=1
            a=A[last]&~m
            while a:
                b=a&-a;a-=b
                go(b.bit_length()-1,m|b)
        go(v,1<<v);cycles.append(C)
    @lru_cache(None)
    def pack(m):
        vals=[0]*(n+1)
        if not m:vals[0]=1;return tuple(vals)
        b=m&-m;v=b.bit_length()-1
        S=pack(m^b)
        for i in range(n):vals[i+1]+=S[i]
        for N,cnt in cycles[v].items():
            if N&m==N and N.bit_count()%2:
                S=pack(m^N)
                for i in range(n):vals[i+1]+=2*cnt*S[i]
        return tuple(vals)
    F=pack(U^(1<<w));O=[0]*(n+1);E=[0]*(n+1)
    for m,z in cycles[w].items():
        target=O if m.bit_count()%2 else E
        for i,p in enumerate(pack(U^m)):target[i]+=z*p
    return (
        tuple(F[i]-(2*E[i-1] if i else 0) for i in range(n+1)),
        (0,)+tuple(O[i]+E[i] for i in range(n)),
        (0,)+tuple(E[i]-O[i] for i in range(n)),
    )

def eval_poly(c,k):return sum(v*k**i for i,v in enumerate(c))

def require(c,msg):
    if not c:raise ValueError(msg)

expected=[
    [(0,-2,0,2,0,1),(0,3,7,4,2),(0,3,-7,4,-2)],
    [(0,0,-2,0,4,0,1),(0,4,10,11,5,2),(0,-4,10,-11,5,-2)],
]
checks=0
for i,A in enumerate((G,outlet(G,5))):
    require(strong(A) and indecomposable(A) and aut(A)==[tuple(range(len(A)))],('quotient',A))
    require(all(not (row>>2&1) for j,row in enumerate(A) if j!=1) and A[1]>>2&1,'entry')
    coeffs=packing_reader(A,2)
    for a,b in zip(coeffs,expected[i]):
        require(a==b+(0,)*(len(a)-len(b)),('coeffs',coeffs))
    for W in ((0,),tt(2),tt(3),C3,U,V):
        X=inflate(A,2,W); pv=all_values(W,len(X)+1); actual=all_values(X,len(X))
        predicted=[0]+[sum(eval_poly(c,k)*pv[k+s] for c,s in zip(coeffs,(0,1,-1))) for k in range(1,len(X)+1)]
        require(actual==predicted,('reader',W,actual,predicted))
        checks+=1

p0=(0,0,18,0,14,0,1)
for W,g in ((U,3),(V,1)):
    require(strong(W) and len(aut(W))==g,('reservoir',W))
    require(all_values(W,6)==[eval_poly(p0,k) for k in range(7)],('slot twin',W))
    X=substitute(C3,[(0,),W,(0,)])
    require(strong(X) and len(aut(X))==g,('wrapped group',W))

p7=(0,24,0,86,0,24,0,1)
for W,g in (((78,52,56,50,33,65,30),3),((70,108,56,49,99,65,12),1)):
    require(strong(W) and len(aut(W))==g,('seven group',W))
    require(all_values(W,7)==[eval_poly(p7,k) for k in range(8)],('seven profile',W))
def labeled(n):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        A=[0]*n
        for j,(u,v) in enumerate(edges):
            if mask>>j&1:A[v]|=1<<u
            else:A[u]|=1<<v
        yield A
nsmall=0
for n in range(1,6):
    groups={}
    for X in labeled(n):
        if not strong(X):continue
        nsmall+=1
        a=len(aut(X));p=tuple(all_values(X,n))
        require(a in (1,3,5) and groups.get(p,a)==a,('small collision',X,a,p))
        groups[p]=a

def skeleton(mods,k):
    # Expand slot words by repeatedly choosing directed cycle length and start.
    sizes=tuple(len(M) for M in mods)
    P=[all_values(M,len(M)) for M in mods]
    L=[[sum((-1)**(j-a)*comb(j,a)*p[a] for a in range(j+1)) for j in range(n+1)] for p,n in zip(P,sizes)]
    one=Counter({(0,0,0):1})
    for l in range(1,sum(sizes)+1):
        for start in range(3):
            v=tuple(sum((start+j)%3==i for j in range(l)) for i in range(3))
            if all(a<=b for a,b in zip(v,sizes)):one[v]+=1
    T=Counter({(0,0,0):1})
    for _ in range(k):
        R=Counter()
        for u,a in T.items():
            for v,b in one.items():
                w=tuple(x+y for x,y in zip(u,v))
                if all(x<=y for x,y in zip(w,sizes)):R[w]+=a*b
        T=R
    return sum(a*L[0][u]*L[1][v]*L[2][w] for (u,v,w),a in T.items())

for mods in [(tt(2),C3,tt(2)),(C3,C3,(0,))]:
    H=all_values(substitute(C3,mods),5)
    require(H==[skeleton(mods,k) for k in range(6)],'lex skeleton')

def convolution(a,b):
    c=Counter()
    for i,x in a.items():
        for j,y in b.items():c[i+j]+=x*y
    return c
exp=[0,1]
for t,n_types in enumerate((2,4,24,4624)):
    require(len(exp)==n_types,('types',t,len(exp)))
    E=Counter(exp); M=(3**(t+1)-1)//2
    require(E[M]==1 and max(e for e in E if e<M)==M-t-1 and E[M-t-1]==1,('top gap',t,E))
    if t:require(min(E)==3**(t-1) and E[min(E)]==2,('two minima',t))
    if t==3:break
    new=[]
    for i,j,k in product(range(len(exp)),repeat=3):
        if (i,j,k)!=min((i,j,k),(j,k,i),(k,i,j)):continue
        new.append(exp[i]+exp[j]+exp[k]+(i==j==k))
    cube=convolution(convolution(E,E),E)
    for d,e in E.items():
        cube[3*d]-=e
    require(all(c%3==0 for c in cube.values()),'orbit division')
    recur=Counter({j:c//3 for j,c in cube.items() if c})
    for d,e in E.items():recur[3*d+1]+=e
    require(recur==Counter(new),('ternary recurrence',t))
    exp=new
print({'functional_full_polynomial_controls':checks, 'independent_lex_controls':2,
       'reservoir_group_orders':(3,1), 'small_strong_tournaments':nsmall, 'ternary_class_totals':(2,4,24,4624)})

def rooted(A):
    n=len(A); U=(1<<n)-1
    s=[[0]*n for _ in range(U+1)]; e=[[0]*n for _ in range(U+1)]
    h=[1]+[0]*U; q=[0]*(U+1)
    for j in range(n):s[1<<j][j]=e[1<<j][j]=1
    for m in range(1,U+1):
        for j in range(n):
            if m>>j&1:
                d=m^(1<<j)
                s[m][j]+=sum(s[d][l] for l in range(n) if A[j]>>l&1)
                e[m][j]+=sum(e[d][l] for l in range(n) if A[l]>>j&1)
        h[m]=sum(s[m])
    for m in range(U+1):
        a=m
        while True:
            q[m]+=h[a]*h[m^a]
            if not a:break
            a=(a-1)&m
    def start(j,P,sign):
        return sum(((-1)**m.bit_count() if sign else 1)*s[m][j]*P[U^m] for m in range(1,U+1))
    return [(start(j,h,False),start(j,h,True),start(j,q,False),start(j,q,True),
             s[U][j],sum(e[m][j]*h[U^m] for m in range(1,U+1))) for j in range(n)]

small=(78,52,56,50,33,65,30)
large=(70,108,56,49,99,65,12)
read_controls=0
for W,x,expected in ((small,4,(165,-27,993,-285)),(large,6,(165,-25,997,-273))):
    rx=rooted(W)[x]
    require(rx[:4]==expected,('rooted signatures',rx))
    h,q=all_values(W,2)[1:]
    for r in range(4):
        M=substitute((2,0),(W,tt(r))) if r else W
        left=outlet(M,x)
        right=outlet(substitute(C3,((0,),M,(0,))),x+1)
        Rx=[h+2**r*rx[0]+(rx[1] if not r else 0),2**(r+1)*q+2*3**r*rx[2]+2*rx[3]]
        require(all_values(left,2)[1:]==Rx and strong(left),('ordinal reader',r))
        Ry=h+2**r*q+2**(r+1)*rx[0]+2*3**r*rx[2]-rx[4]-rx[5]
        require(all_values(right,1)[1]==Ry and strong(right),('triangle reader',r))
        read_controls+=2
print({'internal_reader_controls':read_controls})
