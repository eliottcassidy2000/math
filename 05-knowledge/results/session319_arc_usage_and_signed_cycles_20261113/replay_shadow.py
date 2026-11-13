from itertools import permutations,combinations
from collections import defaultdict

def need(ok):
    if not ok: raise RuntimeError('control failed')

def census(n):
    E=list(combinations(range(n),2)); I={e:i for i,e in enumerate(E)}; P=list(permutations(range(n)))
    acts=[[(I[min(p[u],p[v]),max(p[u],p[v])],p[u]>p[v]) for u,v in E] for p in P]
    tests=[]
    for p in P:
        x=y=0; js=[]
        for u,v in zip(p,p[1:]):
            j=I[min(u,v),max(u,v)]; js.append(j)
            if u<v: x|=1<<j
            else: y|=1<<j
        tests.append((x,y,js,p[0],p[-1]))
    def data(t):
        w=[0]*len(E);s=[0]*n;e=[0]*n; h=0
        for x,y,js,a,b in tests:
            if t&x==x and not t&y:
                h+=1;s[a]+=1;e[b]+=1
                for j in js:w[j]+=1
        return h,tuple(w),s,e
    def orbit(t):
        return {sum((((t>>j)&1)^f)<<i for i,(j,f) in enumerate(A)) for A in acts}
    def wcanon(w):
        return min(tuple(w[j] for j,_ in A) for A in acts)
    def isstrong(t):
        out=[[] for _ in range(n)]
        for i,(u,v) in enumerate(E):
            if not (t>>i)&1:u,v=v,u
            out[u].append(v)
        for u in range(n):
            seen={u}; stk=[u]
            while stk:
                for v in out[stk.pop()]:
                    if v not in seen:seen.add(v);stk.append(v)
            if len(seen)<n:return False
        return True
    byW=defaultdict(list); unused=[]; seen=set()
    for t in range(1<<len(E)):
        if t in seen: continue
        o=orbit(t); seen|=o
        r=min(o|{q^((1<<len(E))-1) for q in o})
        if t!=r:continue
        h,w,s,e=data(t); strong=isstrong(t)
        need(sum(w)==h*(n-1))
        for u in range(n):need(sum(w[j] for j,q in enumerate(E) if u in q)==2*h-s[u]-e[u])
        byW[wcanon(w)].append((t,h,strong))
        if strong and 0 in w: unused.append((t,h,[q for j,q in enumerate(E) if w[j]==0],w))
    coll=[(w,ts) for w,ts in byW.items() if len(ts)>1]
    print(n,'usage images',len(byW),'collisions',len(coll),'zero strong',len(unused))
    if coll: print('first collision',coll[0])
    if unused: print('first zero strong',unused[0])
    return coll,unused

if __name__=='__main__':
    for n in range(2,8):census(n)
