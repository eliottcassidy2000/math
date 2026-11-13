from itertools import combinations, permutations
def table(n):
    E=tuple(combinations(range(n),2))
    A=[[] for _ in E]
    for i,e in enumerate(E):
        for j,f in enumerate(E):
            if len(set(e)&set(f))==1:
                v,=set(e)&set(f)
                A[i].append((j,int(v==e[0])^int(v==f[0])))
    return E,A
def inc(n,t,u,E,A):
    nE=len(E);D={0:0};q=[0]
    def w(code,e,f,k):return ((code>>e)^(code>>f)^k)&1
    while q:
        i=q.pop()
        for j,k in A[i]:
            d=D[i]^w(t,i,j,k)^w(u,i,j,k)
            if j in D:
                if D[j]!=d:raise ValueError('inconsistent')
            else:D[j]=d;q.append(j)
    valid=set()
    for p in (0,1):
        S={j for j in D if D[j]==p}
        if len(S)!=n-1:continue
        deg=[0]*n
        for j in S:
            for v in E[j]:deg[v]+=1
        if sorted(deg)!=[1,1]+[2]*(n-2):continue
        todo=[min(S)];got=set(todo)
        while todo:
            i=todo.pop()
            for j,k in A[i]:
                if j in S:
                    if not w(t,i,j,k):break
                    if j not in got:got.add(j);todo.append(j)
            else:continue
            break
        else:
            if got==S:valid.add(frozenset(S))
    return valid
for n in (3,4):
    E,A=table(n);idx={e:i for i,e in enumerate(E)};M=(1<<len(E))-1
    for t in range(M+1):
        actual={}
        for p in permutations(range(n)):
            if not all(bool(t>>idx[min(a,b),max(a,b)]&1)==(a<b) for a,b in zip(p,p[1:])):continue
            S=frozenset(idx[min(a,b),max(a,b)]for a,b in zip(p,p[1:]))
            u=t^M^sum(1<<i for i in S)
            for v in (u,u^M):actual.setdefault(v,set()).add(S)
        for u in range(M+1):
            if inc(n,t,u,E,A)!=actual.get(u,set()):raise ValueError((n,t,u))
print('all aligned pairs through n=4, including the double-spine boundary')
