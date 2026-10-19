from itertools import product
from collections import deque

def run(F):
    init=[tuple(w)+(3,) for w in product(range(2),repeat=F-1)]
    ids={}; states=[]; adj=[]
    def add(s):
        if s[0]>=2:return -1
        if max(s)<=1:return -2
        if s not in ids:
            ids[s]=len(states);states.append(s)
        return ids[s]
    start=[add(s) for s in init]
    j=0
    while j<len(states):
        s=states[j];base=tuple(abs(x-y) for x,y in zip(s,s[1:]))
        adj.append(tuple(add(base+(abs(s[-1]-u),)) for u in [0,1]));j+=1
    n=len(states);pred=[[]for _ in states]
    for i,aa in enumerate(adj):
        for v in aa:
            if v>=0:pred[v].append(i)
    # Reachability to fail: states with no such path are probability0.
    canfail=set(i for i,a in enumerate(adj)if -1 in a);q=deque(canfail)
    while q:
        for i in pred[q.popleft()]:
            if i not in canfail:canfail.add(i);q.append(i)
    # Can reach the p0 closed class or binary absorbing state.
    cansafe=set(i for i,a in enumerate(adj) if -2 in a or i not in canfail);q=deque(cansafe)
    while q:
        for i in pred[q.popleft()]:
            if i not in cansafe:cansafe.add(i);q.append(i)
    mixed=canfail&cansafe
    # Is induced mixed graph acyclic? If yes exact dyadics; longest paths.
    out={i:sum(v in mixed for v in adj[i])for i in mixed};q=deque(i for i in mixed if out[i]==0)
    depth={i:1 for i in q};p={i:0 for i in range(n)if i not in canfail};p.update({i:1 for i in range(n) if i not in cansafe});p[-1]=1;p[-2]=0
    from fractions import Fraction
    order=[]
    while q:
        i=q.popleft();order.append(i);p[i]=sum((p[v]for v in adj[i]),Fraction(0))/2
        depth[i]=1+max((depth[v]for v in adj[i]if v in mixed),default=0)
        for z in pred[i]:
            if z in mixed:
                out[z]-=1
                if out[z]==0:q.append(z)
    print('F',F,'states',n,'mixed',len(mixed),'acyclic',len(order)==len(mixed),'maxdepth',max(depth.values(),default=0),flush=True)
    if len(order)==len(mixed):
        val=sum((p[i]for i in start),Fraction(0))/len(start)
        print('prob',str(val),'scaled_excess',str((val-F*Fraction(2,2**F))*2**F), 'initial_depth',max((depth.get(i,0)for i in start)),flush=True)
    else:
        cycl=[states[i]for i,v in out.items()if v];print('cyclic sample',cycl[:5])
    return len(order)==len(mixed)
if __name__=='__main__':
    import sys
    for f in range(3,int(sys.argv[1])+1):run(f)
