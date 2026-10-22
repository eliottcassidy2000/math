from common import *
from itertools import combinations

def sig(Q):
    h,l,c=coefficients(Q);k=len(Q)
    return (h,tuple(l),tuple(c[u][v]+c[v][u] for u in range(k) for v in range(u+1,k)))
def cyclic_module(Q,M):
    k=len(Q);M=set(M)
    return all(sum(Q[u][v] for v in M)==1 for u in M) and all(len({Q[x][u] for u in M})==1 for x in range(k) if x not in M)
def flip_inside(Q,M):
    R=[r[:] for r in Q]
    for u in M:
        for v in M:
            if u!=v:R[u][v]=1-Q[u][v]
    return R
for k in [4,5]:
    groups=defaultdict(list)
    for q in range(1<<(k*(k-1)//2)):groups[sig(core(k,q))].append(q)
    col=[g for g in groups.values() if len(g)>1];counts=defaultdict(int);other=[]
    for g in col:
        if len(g)!=2:raise RuntimeError('not pairs')
        Q,R=[core(k,q) for q in g]
        mods=[M for M in combinations(range(k),3) if cyclic_module(Q,M) and flip_inside(Q,M)==R]
        regular=all(sum(r)==(k-1)//2 for r in Q) and all(sum(r)==(k-1)//2 for r in R)
        converse=all(Q[u][v]==R[v][u] for u in range(k) for v in range(k))
        if mods:counts['cyclic_module_flip']+=1
        elif regular and converse:counts['regular_converse']+=1
        else:other.append(g)
    print(k,dict(counts),'unclassified',other)
# Criterion vs direct full one-mask response for every labelled Q through k5.
checks=0;vt_examples=[]
for k in range(1,6):
 for mask in range(1<<(k*(k-1)//2)):
    Q=core(k,mask);R=[list(x) for x in zip(*Q)];h,l,c=coefficients(Q);_,lr,cr=coefficients(R)
    s=[0]*k;e=[0]*k
    for p in permutations(range(k)):
        if all(Q[p[j-1]][p[j]] for j in range(1,k)):s[p[0]]+=1;e[p[-1]]+=1
    criterion=all(2*(e[u]-s[u])+sum(c[u][v]-c[v][u] for v in range(k))==0 for u in range(k))
    if criterion!=(sig(Q)==sig(R)):raise RuntimeError('criterion mismatch')
    checks+=1
print('converse criterion checks:',checks)
