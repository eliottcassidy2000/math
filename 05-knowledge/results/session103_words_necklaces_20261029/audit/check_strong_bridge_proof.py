"""Independent all-subset bridge test, not cycle-interval enumeration."""
from itertools import combinations
from pathlib import Path
import json

def need(p,m):
    if not p:raise RuntimeError(m)

def strong(A,V):
    # Standard definition, deliberately including singleton strong components.
    for start in V:
        seen={start}
        while True:
            nxt=seen|{v for u in seen for v in V if A[u][v]}
            if nxt==seen:break
            seen=nxt
        if seen!=set(V):return False
    return bool(V)

def cuts(A,min_size=3):
    n=len(A);out=[]
    for k in range(min_size,n-min_size+1):
        for R in combinations(range(n),k):
            Q=tuple(v for v in range(n) if v not in R)
            if sum(A[r][q] for r in R for q in Q)==1 and strong(A,R) and strong(A,Q):out.append((R,Q))
    return out

def children(A,C):
    n=len(A);maxima=[v for v in range(n) if sum(A[u][v] for u in range(n))==1]
    if len(maxima)==1:
        a=maxima[0];opts=[('a',C[C.index(a)-1]),('b',a)]
    else:
        a,b=maxima
        if A[b][a]:a,b=b,a
        opts=[('c',C[C.index(a)-1]),('d',a),('e',b)]
    for letter,p in opts:
        B=[r+[int(i==p)] for i,r in enumerate(A)]+[[int(i!=p) for i in range(n)]+[0]]
        D=C[:];D.insert(D.index(p)+1,n)
        yield letter,B,D

def main():
    A=[[0]*4 for _ in range(4)]
    for i,j in ((0,1),(1,2),(2,3),(3,0),(0,2),(1,3)):A[i][j]=1
    seed=cuts(A,1);need(((2,),(0,1,3)) in seed,'singleton hostile')
    stack=[('',A,list(range(4)),False)];rows={};births=0
    while stack:
        w,A,C,old_bad=stack.pop();n=len(A);bs=cuts(A)
        need(bool(bs)==('be' in w or 'ee' in w),('equivalence',w))
        rows[n]=rows.get(n,0)+1
        if bs and not old_bad:
            need(w[-2:] in ('be','ee'),('birth ending',w));births+=1
            for R,Q in bs:
                need(n-1 in Q and len(Q)==3,('birth Q size',w,R,Q))
        for R,Q in bs:
            need(any(set((C+C)[o:o+len(R)])==set(R) for o in range(n)),('interval',w))
        if n<9:
            for l,B,D in children(A,C):stack.append((w+l,B,D,bool(bs)))
    out={'scope':'all legal histories through order9; all vertex subsets; blocks explicitly size>=3','words_by_order':rows,'first_appearance_words':births,'standard_singleton_seed_counterexample':{'R':[2],'Q':[0,1,3],'bridge':[2,3]}}
    Path(__file__).with_name('strong_bridge_audit_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
