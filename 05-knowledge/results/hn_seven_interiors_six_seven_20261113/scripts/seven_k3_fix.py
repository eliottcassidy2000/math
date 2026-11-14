# Re-check the 10-vertex / 3-group / 7-interior case with group point-counts (4,2,2):
# one colour class holds two port pairs.  WLOG each contact is made by exactly one port point
# (all filters are monotone: deleting contacts never creates a violation).
import itertools
from seven_interiors import critical_graphs, is_clique, popc, circle_ok
def realizable_sizes(adj,P,I,sizes):
    n=len(adj); Imask=sum(1<<v for v in I)
    for Q in itertools.combinations(I,4):
        if is_clique(adj,Q): return False
    for a,b in itertools.combinations(I,2):
        if popc(adj[a]&adj[b]&Imask)>=3: return False
    groups=[]
    for g,sz in zip(P,sizes):
        S=[x for x in I if adj[g]>>x&1]; opts=set()
        for assign in itertools.product(range(sz),repeat=len(S)):
            Ts=[0]*sz
            for x,a in zip(S,assign): Ts[a]|=1<<x
            if all(circle_ok(adj,T,Imask) for T in Ts) and all(popc(adj[x]&Imask&T)<=2 for T in Ts for x in I):
                opts.add(tuple(sorted(Ts)))
        if not opts: return False
        groups.append(sorted(opts))
    order=sorted(range(len(groups)),key=lambda i:len(groups[i]))
    chosen=[]
    def bt(i):
        if i==len(order):
            for a,b in itertools.combinations(I,2):
                if popc(adj[a]&adj[b]&Imask)+sum(1 for T in chosen if T>>a&1 and T>>b&1)>2: return False
            return True
        for Ts in groups[order[i]]:
            if any(popc(T&U)>2 for T in Ts for U in chosen): continue
            chosen.extend(Ts)
            if bt(i+1): return True
            del chosen[-len(Ts):]
        return False
    return bt(0)
crit=critical_graphs(10,3)
surv=0
for gi,adj in enumerate(crit):
    for P in itertools.combinations(range(10),3):
        if not is_clique(adj,P): continue
        I=[v for v in range(10) if v not in P]
        for big in range(3):
            sizes=[4 if t==big else 2 for t in range(3)]
            if realizable_sizes(adj,P,I,sizes): surv+=1; print("survivor",gi,P,sizes)
print("10-vertex critical, 3 colour groups of sizes (4,2,2), 7 interiors: survivors =",surv)
