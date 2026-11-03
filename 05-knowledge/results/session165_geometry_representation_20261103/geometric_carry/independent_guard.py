from itertools import combinations, permutations

def need(ok, label):
    if not ok: raise RuntimeError(label)

def paths(adj, vertices):
    return [p for p in permutations(vertices) if all(b in adj[a] for a,b in zip(p,p[1:]))]

def state(adj,q):
    n=len(adj); vertices=tuple(range(n)); incoming=[v for v in vertices if q in adj[v]]
    need(len(incoming)==1, 'unique guard')
    allpaths=paths(adj, vertices)
    s=sum(p[0]==q for p in allpaths)
    r=incoming[0]
    d=sum(p[0]==r for p in paths(adj, tuple(v for v in vertices if v!=q)))
    return (s,len(allpaths)-s,d)

def attach(adj,q):
    n=len(adj); out=[set(x) for x in adj]+[set(range(n))-{q}]
    out[q].add(n)
    return out

def strong(adj):
    for v in range(len(adj)):
        seen={v}; todo=[v]
        while todo:
            for w in adj[todo.pop()]-seen:
                seen.add(w);todo.append(w)
        if len(seen)!=len(adj):return False
    return True

count=0; strong_count=0
for n in range(3,6):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        adj=[set() for _ in range(n)]
        for k,(a,b) in enumerate(edges):
            if mask>>k&1:adj[a].add(b)
            else:adj[b].add(a)
        for q in range(n):
            if sum(q in x for x in adj)!=1:continue
            s,u,d=state(adj,q); out=attach(adj,q)
            need(state(out,n)==(u,s+u+d,s),'C formula')
            need(state(out,q)==(s+d,2*u,0),'D formula')
            count+=1
            if strong(adj):
                need(strong(out),'strength preservation');strong_count+=1
returns=[]
for j in range(1,5):
    adj=[{1},{2},{0}]; q=0
    for letter in 'D'*(j-1)+'CD':
        adj=attach(adj,q)
        if letter=='C':q=len(adj)-1
    t=1+2**j
    need(state(adj,q)==(t,2*t,0),'return state')
    need(strong(adj),'return strength')
    need(len(adj)==j+4,'return order')
    returns.append((j,len(adj),sum(state(adj,q)[:2])))
print({'exhaustive_guarded_tournaments':count,'strong_guarded_tournaments':strong_count,'small_returns':returns,'status':'PASS'})
