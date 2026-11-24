from compact import field_host,decode_field_host
def max_own_overlap(G):
    v=sorted(G);K=[(x,y) for x in v for y in v]
    yes=lambda x,y: x==y or y in G[x]
    def ok(a,b):
        x,y=a; X,Y=b
        return x!=X and y!=Y and yes(x,Y)!=yes(X,y)
    adj=[{j for j in range(i+1,len(K)) if ok(a,K[j])} for i,a in enumerate(K)]
    ans=[]
    def extend(p,c):
        nonlocal ans
        if len(p)>len(ans):ans=p[:]
        while len(p)+len(c)>len(ans):
            i=min(c);c.remove(i)
            extend(p+[i],c&adj[i])
    extend([],set(range(len(K))))
    return [K[i] for i in ans]
C={i:{(i-1)%5,(i+1)%5} for i in range(5)}
for d,G in [(2,C),(3,decode_field_host(field_host(1)[0],1))]:
    u=max_own_overlap(G)
    print(d,len(u),u)
