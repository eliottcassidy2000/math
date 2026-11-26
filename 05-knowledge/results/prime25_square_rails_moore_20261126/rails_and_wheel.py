from collections import Counter, deque
from math import gcd,isqrt
import json
from planted_stars import orbit_chart

def primes():
    return [p for p in range(2,100) if all(p%d for d in range(2,isqrt(p)+1))]
def is_square(x):
    return isqrt(x)**2==x
def graph(labels):
    return [[j for j,b in enumerate(labels) if i!=j and is_square(a+b)] for i,a in enumerate(labels)]
def factor(labels,k):
    G=graph(labels);n=len(G);s,t=2*n,2*n+1
    A=[[] for _ in range(2*n+2)]
    def add(u,v,w):
        A[u].append([v,len(A[v]),w]);A[v].append([u,len(A[u])-1,0])
    for i in range(n):
        add(s,i,k);add(n+i,t,k)
        for j in G[i]:add(i,n+j,1)
    flow=0
    while True:
        dist=[-1]*len(A);dist[s]=0;Q=[s]
        for u in Q:
            for v,b,cap in A[u]:
                if cap and dist[v]<0:dist[v]=dist[u]+1;Q.append(v)
        if dist[t]<0:break
        idx=[0]*len(A)
        def dfs(u,f):
            if u==t:return f
            while idx[u]<len(A[u]):
                z=A[u][idx[u]];v,j,c=z
                if c and dist[v]==dist[u]+1:
                    w=dfs(v,min(f,c))
                    if w:z[2]-=w;A[v][j][2]+=w;return w
                idx[u]+=1
            return 0
        while (w:=dfs(s,n*k)):flow+=w
    U=[i for i in range(n) if dist[i]>=0]; V=[i for i in range(n) if dist[n+i]>=0]
    X=sum(j not in V for i in U for j in G[i])
    assert flow==k*(n-len(U)+len(V))+X
    if flow==n*k:
        used=[{j-n for j,_,c in A[i] if n<=j<2*n and c==0} for i in range(n)]
        assert all(len(s)==k for s in used)
        perms=[]
        for _ in range(k):
            own=[-1]*n
            def aug(i,seen):
                for j in sorted(used[i]-seen):
                    seen.add(j)
                    if own[j]<0 or aug(own[j],seen):
                        own[j]=i;return True
                return False
            assert all(aug(i,set()) for i in range(n))
            p=[None]*n
            for j,i in enumerate(own):p[i]=labels[j];used[i].remove(j)
            perms.append(p)
        return flow,{"witness_permutations":perms}
    return flow,{"cut_left":[labels[i] for i in U],"cut_right":[labels[i] for i in V],"cross_edges":X}
def pack(labels):
    H=graph(labels);n=len(H)
    for k in range(1,n):
        flow,r=factor(labels,k)
        if flow!=n*k:
            return {"vertices":n,"edges":sum(map(len,H))//2,"min_degree":min(map(len,H)),
                    "max_block_count":k-1,"first_failed_flow":[flow,k*n],"cut":r,
                    "witness":factor(labels,k-1)[1] if k>1 else None}
def cycle(p):
    a,b=0,5%p;W=[]
    while not W or (a,b)!=(0,5%p):
        W.append(a);a,b=b,(23*b-a+5)%p
        assert len(W)<=p*p
    return W
def even_factor(labels):
    n=len(labels); E=[(u,v) for u,L in enumerate(graph(labels)) for v in L if u<v]
    K=len(E); AT=[0]*n
    for j,(u,v) in enumerate(E):AT[u]|=1<<j;AT[v]|=1<<j
    nodes=0
    def test(P,U):
        A=[[] for _ in range(n)]
        while P:
            e=P&-P;P^=e;u,v=E[e.bit_length()-1]
            A[u].append(v);A[v].append(u)
        seen=set(); cyc=[]
        for i in range(n):
            if i in seen:continue
            Z=[i];seen.add(i)
            for j in Z:
                for t in A[j]:
                    if t not in seen:seen.add(t);Z.append(t)
            if all(len(A[j])==2 for j in Z):
                if len(Z)%2:return False
                cyc.append(len(Z))
        return cyc
    def walk(P,U):
        nonlocal nodes
        nodes+=1
        while True:
            upd=False
            for row in AT:
                need=2-(P&row).bit_count();rest=U&row
                if need<0 or rest.bit_count()<need:return None
                if not rest:continue
                if not need:U^=rest;upd=True
                elif rest.bit_count()==need:P|=rest;U^=rest;upd=True
            if not upd:break
        if test(P,U) is False:return None
        if not U:return P
        g=min((row&U for row in AT if row&U),key=int.bit_count);e=g&-g
        x=walk(P|e,U^e)
        return x if x is not None else walk(P,U^e)
    P=walk(0,(1<<K)-1)
    if P is None:return {"nodes":nodes,"witness":None}
    return {"nodes":nodes,"cycles":sorted(test(P,0)),
        "witness":[[labels[u],labels[v]] for j,(u,v) in enumerate(E) if P>>j&1]}

def equitable_square(g):
    m=56; n=m//g
    # Each block (r,s) chooses the most square-supported relative shift.
    M={}
    for r in range(1,g+1):
        for s in range(1,g+1):
            candidates=[]
            for b in range(n):
                if r==s and not b:continue
                score=sum(is_square(r+g*k+s+g*((k+b)%n)) for k in range(n))
                candidates.append((score,b))
            M[r,s]=max(candidates)
    D={0:(0,[])}
    for r in range(1,g+1):
        E={}
        for mask,(score,A) in D.items():
            for s in range(1,g+1):
                if mask>>(s-1)&1:continue
                q,b=M[r,s]; U=mask|(1<<(s-1))
                if U not in E or q+score>E[U][0]:
                    E[U]=(q+score,A+[(s,b)])
        D=E
    value,w=D[(1<<g)-1]
    assert all(t[0]<n for t in M.values())
    return {"number_of_orbits":g,"orbit_size":n,"square_slots":value,
            "one_optimal_orbit_map":w,
            "best_single_orbit_contact_count":max(t[0] for t in M.values())}
def parity_types(labels,R):
    G=[set() for _ in labels]; off={v:i for i,v in enumerate(labels)}
    for u,v in R["witness"]:
        G[off[u]].add(off[v]);G[off[v]].add(off[u])
    seen=set();M=[{},{}]
    for u in range(len(G)):
        if u in seen:continue
        p=[u];v=min(G[u])
        while v!=u:
            p.append(v);v=next(z for z in G[v] if z!=p[-2])
        for j,v in enumerate(p):
            z=p[(j+1)%len(p)]
            M[j%2][v]=z;M[j%2][z]=v
        seen.update(p)
    assert len(M[0])==len(M[1])==56
    return {"opposite_parity_entries":[sum((i+j)%2 for i,j in A.items()) for A in M],
            "relative_cycles":sorted([p//2 for p in R["cycles"]]*2)}

def main():
    ps=primes(); assert len(ps)==25 and ps[-1]==97
    cases=[("Q25",list(range(1,26))),("P25",ps),("Q56",list(range(1,57))),
           ("Ogg15",[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71])]
    exact={name:pack(L) for name,L in cases}
    periods={}
    for p in ps:
        v=cycle(p); z=[i for i,m in enumerate(v) if not m]
        assert z==([0] if p==5 else [0,len(v)-1])
        periods[p]=len(v)
    v=cycle(14);clock={}
    for j,m in enumerate(v):clock.setdefault(14//gcd(14,m),[]).append(j)
    assert len(v)==21 and set(clock)=={1,2,7,14}
    C,_=orbit_chart("period_fourteen_site.txt",4,((10,45,2),(27,49,0)))
    m=56;a,b=0,28
    route=[list(C[a,b]) if c==a else list(range(m)) if c==b else [C[c,b][C[a,c][i]] for i in range(m)] for c in range(m)]
    assert all(len({p[i] for p in route})==m for i in range(m))
    def fail(k):return [c for c,p in enumerate(route) if not all(p[(i+k)%m]==(p[i]+k)%m for i in range(m))]
    census={}
    for ell,ph in clock.items():
        fails={tuple(fail(4*v[j])) for j in ph};assert len(fails)==1
        f=list(next(iter(fails)));census[ell]={"phases":ph,"commutes":56-len(f),"exceptions":f}
    reflections=[]
    for s in range(m):
        p=[(s-(i+1)-1)%m+1 for i in range(m)]
        if all(y!=x and is_square(x+y) for x,y in enumerate(p,1)):
            reflections.append([s,sorted({x+y for x,y in enumerate(p,1)})])
    assert reflections==[[25,[25,81]]]
    star={}
    for a in (0,28):
        weights=[]
        for c in range(m):
            if c==a:continue
            t=C[a,c]
            weights.append(sum(is_square(i+y+2) for i,y in enumerate(t)))
        assert sum(weights)==2*exact["Q56"]["edges"]
        star[a]={"max_on_one_block":max(weights),"histogram":dict(sorted(Counter(weights).items()))}
    inv=even_factor(list(range(1,57)))
    inv["alternating_relative_profile"]=parity_types(list(range(1,57)),inv)
    print(json.dumps({"packings":exact,"involutory_pair_Q56":inv,
                      "square_equivariant_caps":{g:equitable_square(g) for g in [4,8]},"pell_periods":periods,"rail_transport":census,
                      "full_square_reflections":reflections,"square_contacts_in_chart":star},indent=2))

if __name__=="__main__":main()
