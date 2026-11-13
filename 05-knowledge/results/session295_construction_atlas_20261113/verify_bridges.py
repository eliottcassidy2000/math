from itertools import combinations, permutations
from functools import cache

def profile(arcs, verts):
    verts = tuple(verts)
    n = len(verts)
    @cache
    def walk(S, u, v):
        return sum(all((x,y) in arcs for x,y in zip(p,p[1:]))
                   for w in permutations(i for i in S if i not in (u,v))
                   for p in [(u,*w,v)])
    @cache
    def packing(S):
        if not S: return (1,)
        a = min(S)
        rest = tuple(x for x in S if x != a)
        out = [0,*packing(rest)]
        for k in range(3,len(S)+1,2):
            for W in combinations(rest,k-1):
                Wset = set(W)
                tail = tuple(x for x in rest if x not in Wset)
                for b in W:
                    if (b,a) in arcs:
                        for j,v in enumerate(packing(tail)):
                            out[j+1] += 2*v*walk((a,*W),a,b)
        return tuple(out)
    return packing(verts),walk,packing

def direct_h(arcs,n):
    return sum(all((x,y) in arcs for x,y in zip(p,p[1:]))
               for p in permutations(range(n)))

states=set()
cases=0
for n in range(1,5):
    V=tuple(range(n))
    P=list(combinations(V,2))
    for mask in range(1<<len(P)):
        G={e if mask>>i&1 else e[::-1] for i,e in enumerate(P)}
        F,walk,pack=profile(G,V)
        H=direct_h(G,n)
        if sum(F)!=H: raise ValueError('H identity')
        D2=sum(a*2**i for i,a in enumerate(F))-2*H
        if D2%3 != (((-1)**n-2)*H)%3: raise ValueError('mod 3 identity')
        states.add((H%3,D2%3))
        for zmask in range(1<<n):
            B={v for v in V if zmask>>v&1}; A=set(V)-B
            Z=G|{(n,v) for v in B}|{(v,n) for v in A}
            F1,*_=profile(Z,range(n+1))
            out=[0,*F]
            for k in range(2,n+1,2):
                for S in combinations(V,k):
                    s=set(S)
                    p=sum(walk(S,b,a) for b in B&s for a in A&s)
                    for j,v in enumerate(pack(tuple(x for x in V if x not in s))):
                        out[j+1]+=2*v*p
            if tuple(out)!=F1: raise ValueError('selective profile')
            cases+=1
R={(i,j) if (i,j)!=(0,3) else (j,i) for i,j in combinations(range(4),2)}
R|={(i,4) for i in range(4)}
RF,*_=profile(R,range(5))
RH=direct_h(R,5)
if RH!=5: raise ValueError('odd representative')
states.add((RH%3,(sum(a*2**i for i,a in enumerate(RF))-2*RH)%3))
if len(states)!=5: raise ValueError(states)
print('states mod 3:',sorted(states))
print('selective attachments:',cases)

cycle={(0,1),(1,2),(2,0)}
def attach(C): return C|{(4,0),(1,4),(2,4),(3,4)}
P=cycle|{(i,3) for i in range(3)}
Q=cycle|{(3,i) for i in range(3)}
if profile(P,range(4))[0]!=profile(Q,range(4))[0]:
    raise ValueError('not uniform equivalents')
if (direct_h(attach(P),5),direct_h(attach(Q),5))!=(9,5):
    raise ValueError('hostile')
print('nonisomorphic equal uniform profiles, selective counts:',9,5)
