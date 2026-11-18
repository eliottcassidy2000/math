from collections import deque
from math import factorial
def step(X, letter, m):
    X=(0,)+X+(0,)
    Y=[]
    for j in range(1,len(X)-1):
        g,h,b=X[j],X[j+1],X[j-1]
        v=j**3*g
        if letter=='P': v+=j*h-j*j*(j-1)*b
        elif letter=='Q': v+=(j+1)*h+j*(1-j)*(j-1)*b
        else: v+=(2*j+1)*h+j*(1-2*j)*(j-1)*b-2*j*g
        Y.append(v%m)
    return tuple(Y)

for m in (3,9,27,81):
    k=next(j for j in range(1,60) if factorial(j)%m==0)
    seed=(0,tuple(j*j*factorial(j-1)%m for j in range(1,k+1)),
           tuple(j*j*factorial(j-1)%5 for j in range(1,6)))
    D={seed:0};back={};adj={};from_edge={};q=deque([seed])
    while q:
        p,x,y=s=q.popleft();adj[s]=[]
        for c in ('P' if not p else 'Q','R'):
            u,v=step(x,c,m),step(y,c,5)
            if not v[0] or u[0]%3==0:continue
            t=(1-p,u,v)
            back.setdefault(t,set()).add(s);adj[s].append(t)
            if t not in D:
                D[t]=D[s]+1;from_edge[t]=(c,s);q.append(t)
    seen=set();end=[]
    for root in D:
        if root in seen:continue
        dfs=[(root,False)]
        while dfs:
            s,post=dfs.pop()
            if post:end.append(s);continue
            if s in seen:continue
            seen.add(s)
            dfs.append((s,True))
            dfs.extend((t,False) for t in adj[s] if t not in seen)
    done=set();biggest=set();ncomp=0
    for root in reversed(end):
        if root in done:continue
        ncomp+=1;dfs=[root];done.add(root);cur={root}
        while dfs:
            for s in back.get(dfs.pop(),()):
                if s in done:continue
                cur.add(s);done.add(s);dfs.append(s)
        if len(cur)>len(biggest):biggest=cur
    assert done==D.keys()
    C={seed};q=deque([seed])
    while q:
        for s in back.get(q.popleft(),()):
            if s not in C:C.add(s);q.append(s)
    wanted={(h%m,h%5) for h in range(1,10*m,2) if h%3 and h%5}
    got={(s[1][0],s[2][0]) for s in C}
    live={(s[1][0],s[2][0]) for s in biggest}
    assert live==wanted,(m,wanted-live)
    print("mod",10*m,"cutoff",k,"reachable",len(D),"return component",len(C),"largest SCC",len(biggest),"largest SCC projected units",len(live),"reachable projected units",len({(s[1][0],s[2][0]) for s in D}),"max shortest from seed",max(D.values()))
    t=min((s for s in D if s in biggest),key=D.__getitem__)
    word=''
    while t!=seed:
        c,t=from_edge[t];word=c+word
    print('entrance:', word or '(empty)')
