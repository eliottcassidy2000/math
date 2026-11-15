#!/usr/bin/env python3
from itertools import combinations

def profile(A):
    n=len(A); m=(1<<n)-1
    S=[[0]*n for _ in range(m+1)]
    T=[[0]*n for _ in range(m+1)]
    for u in range(n):
        S[1<<u][u]=T[1<<u][u]=1
    h=[1]+[0]*m
    for i in range(1,m+1):
        for u in range(n):
            if not i>>u&1:continue
            j=i^(1<<u)
            for v in range(n):
                if A[u]>>v&1:S[i][u]+=S[j][v]
                if A[v]>>u&1:T[i][u]+=T[j][v]
        h[i]=sum(S[i])
    q=sum(x*h[m^i] for i,x in enumerate(h))
    even=[sum(S[i][u]*h[m^i] for i in range(1,m+1) if i.bit_count()%2==0) for u in range(n)]
    odd=[sum(S[i][u]*h[m^i] for i in range(1,m+1) if i.bit_count()%2==1) for u in range(n)]
    assert [o-e for o,e in zip(odd,even)]==T[m]
    return h[m],q,even,T[m]

def tournaments(n):
    edges=list(combinations(range(n),2))
    for e in range(1<<len(edges)):
        A=[0]*n
        for k,(i,j) in enumerate(edges):
            if e>>k&1:A[j]|=1<<i
            else:A[i]|=1<<j
        yield tuple(A)

def transitive(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def delete(A,u):
    return tuple((row&((1<<u)-1))+((row>>(u+1))<<u) for i,row in enumerate(A) if i!=u)

def outlet(A,u):
    n=len(A)
    return tuple(row+((i!=u)<<n) for i,row in enumerate(A))+(1<<u,)

def common_source(A):
    n=len(A); y=n+2; x=n+1
    T=list(transitive(n+4))
    T[0]^=1<<y; T[y]^=1
    for i,row in enumerate(A,1):
        T[i]=(T[i]&~(((1<<n)-1)<<1))+(row<<1)
    return tuple(T),x,y

def strong_components_closure(A):
    reach=[row|(1<<i) for i,row in enumerate(A)]
    for k in range(len(A)):
        for i in range(len(A)):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return reach

matrices=roots=scc_pairs=0
for n in range(2,6):
    for A in tournaments(n):
        matrices+=1
        h,q,e,t=profile(A)
        reach=strong_components_closure(A)
        for u in range(n):
            roots+=1
            d=A[u].bit_count()
            assert (e[u]>0)==(d>0)
            assert profile(outlet(A,u))[0]==h+2*e[u]
            if d:
                dsum=sum(val for j,val in enumerate(profile(delete(A,u))[3])
                         if A[u]>>(j+(j>=u))&1)
                assert d*h<=2*(d+1)*e[u]-(2*d+1)*dsum
                assert h<=4*e[u]
            if e[u]>=2:assert h<=4*e[u]-3
            for v in range(n):
                if reach[u]>>v&1 and not reach[v]>>u&1:
                    scc_pairs+=1
                    assert e[u]-e[v]>=h

for r in range(1,9):
    A=list(transitive(r+3)); x=r+1; u=r+2
    A[0]^=1<<x; A[x]^=1
    A[x]^=1<<u; A[u]^=1<<x
    for i in range(x):
        # u already loses all spine vertices, the orientation we want.
        assert A[i]>>u&1
    h,_,e,_=profile(A)
    assert (h,e[u])==(1+2**(r+1),1+2**(r-1))

reservoirs=0
for W in [B for n in range(1,5) for B in tournaments(n)]+[transitive(n) for n in range(1,9)]:
    h,q,_,_=profile(W)
    T,x,y=common_source(W)
    H,_,e,_=profile(T)
    assert (H,e[x],e[y])==(h+2*q,2*h+2*q,h+2*q)
    for i,goal in [(x,3*H+2*h),(y,3*H)]:
        B=outlet(T,i)
        assert all(reach==(1<<len(B))-1 for reach in strong_components_closure(B))
        assert profile(B)[0]==goal
    reservoirs+=1
print({'labeled_matrices':matrices,'rooted_outlets':roots,
       'ordered_cross_scc_choices':scc_pairs,'sharp_widths':8,'reservoirs':reservoirs})

def blowout(A,u,W):
    n=len(A); labels=[i for i in range(n) if i!=u]+[u]*len(W)
    return tuple(sum(1<<j for j in range(len(labels)) if i!=j and (
        (A[v]>>labels[j])&1 if labels[j]!=v else (W[i-(n-1)]>>(j-(n-1)))&1))
        for i,v in enumerate(labels))

plus=(10,28,57,48,33,3)
assert all(v==63 for v in strong_components_closure(plus))
assert profile(plus)[0]==33 and profile(outlet(plus,5))[0]==67
assert profile(blowout(plus,2,transitive(2)))[0]==65
assert profile(blowout(outlet(plus,5),2,transitive(2)))[0]==131
for W in [*tournaments(3),*tournaments(4),*[transitive(j) for j in range(1,8)]]:
    h,q,_,_=profile(W)
    assert profile(blowout(plus,2,W))[0]==h+16*q
    assert profile(blowout(outlet(plus,5),2,W))[0]==3*h+32*q

for d in range(1,9):
    A=list(transitive(d+2));A[0]^=1<<(d+1);A[-1]^=1
    s=[2**(d+1)-1]+[2**(d-j+1)+2**j*3**(d-j)-(j==1) for j in range(1,d+1)]+[1+2**d]
    # Adding one sink makes its parity split read each old all-start weight.
    T=tuple(v|(1<<len(A)) for v in A)+(0,)
    H,_,e,_=profile(T)
    assert H==1+2**d and e[:-1]==s

print({'binary_reservoirs':79,'palette_widths':8})
