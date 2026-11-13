"""Exhaustive control of the exact connector predicate on the finite height window."""
from itertools import combinations
from math import isqrt

def ck(ok,msg):
    if not ok: raise AssertionError(msg)

def sq(x): return isqrt(x)**2==x

def graph(P): return {p:{q for q in P-{p} if sq(p+q)} for p in P}

def components(G):
    todo=set(G)
    out=[]
    while todo:
        root=min(todo)
        comp={root}
        queue=[root]
        while queue:
            u=queue.pop()
            ns=G[u]-comp
            comp|=ns
            queue.extend(ns)
        todo-=comp
        out.append(comp)
    return out

def run():
    P={p for p in range(2,72) if all(p%d for d in range(2,isqrt(p)+1))}
    S={p for p in P if p<=29}
    B={2,3,5,7,11,13,17,19,23,31,47}
    O=B|{29,41,59,71}
    ck(sorted(map(len,components(graph(S))))==[2,2,6],"first ten")
    ck(sorted(map(len,components(graph(B-{47}))))==[2,3,5],"first ten Baby")
    ck(sorted(map(len,components(graph(B))))==[3,8],"Baby")
    ck(len(components(graph(O)))==1,"Ogg")
    A=sorted(P-S)
    good=[]
    for i in range(len(A)+1):
        for tup in combinations(A,i):
            U=set(tup)
            G=graph(S|U)
            C=components(G)
            formula=47 in U and (53 in U or {41,59}<=U)
            ck(any(S<=c for c in C)==formula,(U,C))
            ck((len(C)==1)==(formula and not ({37,43,67}&U)),(U,C))
            leaves={p for p in G if len(G[p])==1}
            ck(len(leaves)>=3,(U,leaves))
            if formula: good.append(U)
    ck({frozenset(u) for u in good if not any(v<u for v in good)}
       =={frozenset((47,53)),frozenset((41,47,59))},"minimal choices")
    print("optional subsets:",2**len(A),"; connecting:",len(good),"; at least three leaves in all")
    print("minimal connecting choices: 47,53 and 41,47,59")
if __name__=='__main__':run()
