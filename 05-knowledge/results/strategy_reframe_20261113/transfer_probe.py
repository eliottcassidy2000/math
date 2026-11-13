"""Independent edge-level controls for context composition and network shears."""
from pathlib import Path
from itertools import product,combinations
import json,sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_context_recursion_20261113'))
from typed_glue import glue,decode

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def network_counts(n,edges,s,t):
    counts=[]
    for components in (1,2):
        count=0
        for indices in combinations(range(len(edges)),n-components):
            parent=list(range(n))
            def find(v):
                while parent[v]!=v:v=parent[v]
                return v
            good=True
            for i in indices:
                a,b=map(find,edges[i])
                if a==b:good=False;break
                parent[a]=b
            if good and (components==1 or find(s)!=find(t)):count+=1
        counts.append(count)
    return tuple(counts)

def run():
    contexts=[(decode(n,c),r,u) for n in (2,3) for c in range(1<<(n*(n-1)//2)) for r in range(n) for u in range(n)]
    B=(2,4,9,3);b=3;composition=0
    for (A1,r1,u1),(A2,r2,u2) in product(contexts,repeat=2):
        left=glue(A2,r2,glue(A1,r1,B,b),u1)
        C=glue(A2,r2,A1,u1)
        rC=r2 if r1==u1 else len(A2)+[v for v in range(len(A1)) if v!=u1].index(r1)
        right=glue(C,rC,B,b)
        need(left==right,'edge-level context composition, including labels')
        composition+=1
    # L adds a parallel edge; R appends a series edge and moves the terminal.
    frontier=[('',2,[(0,1)],0,1,1,1)];networks=0;states=set()
    for depth in range(8):
        nxt=[]
        for word,n,edges,s,t,tau,f in frontier:
            need(network_counts(n,edges,s,t)==(tau,f),'network tree/forest response')
            need((tau,f) not in states,'Euclidean tree uniqueness')
            states.add((tau,f));networks+=1
            nxt.append((word+'L',n,edges+[(s,t)],s,t,tau+f,f))
            nxt.append((word+'R',n+1,edges+[(t,n)],s,n,tau,tau+f))
        frontier=nxt
    # Arithmetic claims used in strategy triage; primality independently by trial division.
    def prime(n):return n>=2 and all(n%d for d in range(2,__import__('math').isqrt(n)+1))
    need(43*1234517==53084231 and prime(43) and prime(1234517),'composite obstruction factorization')
    result=dict(status='FINITE-EXACT controls; composition and network laws proved in notes',composition_controls=composition,network_controls=networks,network_depth=7,composite_obstruction=dict(N=53084231,prime_factors=[43,1234517]))
    suffix='' if __debug__ else '.optimized'
    (HERE/('transfer'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)

if __name__=='__main__':run()
