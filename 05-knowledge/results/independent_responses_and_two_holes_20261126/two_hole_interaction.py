from collections import Counter
from itertools import combinations
from pathlib import Path
import json, sys

def parse(path):
    it=iter(map(int, Path(path).read_text().split()))
    n=next(it); rows=[]
    for _ in range(n):
        r,k=next(it),next(it)
        rows.append((r,[tuple(next(it) for _ in range(4)) for _ in range(k)]))
    assert next(it,None) is None and n==53
    marker=rows[0][1][0][1:]
    return [(r, [x for x in row if all(a!=b for a,b in zip(x[1:],marker))]) for r,row in rows[1:]]

def replay(path):
    rows=parse(path); holes=[]; used=[set() for _ in range(3)]
    for r,row in rows:
        if r<0:
            holes.append(row);continue
        x=next(y[1:] for y in row if y[0]==r)
        for k,s in enumerate(x): 
            assert s not in used[k]; used[k].add(s)
    A,B=holes;assert (len(A),len(B))==(32,33)
    labels=[h for h in range(56) if h not in {0,1,2,28}]
    slot_owners=[{} for _ in range(3)]
    by_label={}
    for h,(r,row) in zip(labels,rows):
        by_label[h]=(r,row)
        if r>=0:
            x=next(a for a in row if a[0]==r)
            for k,s in enumerate(x[1:]): slot_owners[k][s]=h
    U=[tuple(int(s in used[k]) for k,s in enumerate(x[1:])) for x in A]
    V=[tuple(int(s in used[k]) for k,s in enumerate(x[1:])) for x in B]
    W=[[tuple(int(a==b) for a,b in zip(x[1:],y[1:])) for y in B] for x in A]
    P=[[tuple(U[i][k]+V[j][k]+W[i][j][k] for k in range(3)) for j in range(len(B))] for i in range(len(A))]
    S=[[sum(x) for x in row] for row in P]
    m=min(map(min,S))
    T=Counter()
    independent=nonconstant=bothUnit=toMin=0; first=None
    for a,c in combinations(range(len(A)),2):
        for b,d in combinations(range(len(B)),2):
            z=tuple(W[a][b][k]+W[c][d][k]-W[a][d][k]-W[c][b][k] for k in range(3))
            T[sum(z)]+=1
            if any(z): continue
            independent+=1
            s,x,y,t=S[a][b],S[c][b],S[a][d],S[c][d]
            assert s+t==x+y
            if (s!=x and s!=y): nonconstant+=1
            if abs(s-x)==abs(s-y)==1:
                bothUnit+=1
                if min(s,x,y,t)==m:
                    toMin+=1
                    if first is None: 
                        first=[[A[a][0],A[c][0]],[B[b][0],B[d][0]],[P[a][b],P[c][b],P[a][d],P[c][d]]]
    def repair(i,j):
        x,y=A[i],B[j]
        out=[]
        for j0,z in enumerate([x,y]):
            for k,s in enumerate(z[1:]):
                h=slot_owners[k].get(s)
                if h is not None:
                    r,row=by_label[h]
                    opts=[p for p in row if p[0]!=r and
                          all(t!=a and t!=b for t,a,b in zip(p[1:],x[1:],y[1:]))]
                    cut=lambda p:tuple(sorted({slot_owners[k][t] for k,t in enumerate(p[1:]) if t in slot_owners[k] and slot_owners[k][t]!=h}))
                    least=min(map(lambda p:len(cut(p)),opts),default=None)
                    candidates=[(p[0],cut(p)) for p in opts if len(cut(p))==least]
                    out.append([13 if j0==0 else 33,k,s,h,len(opts),least,candidates])
        return out
    return dict(owner_repairs=[[A[i][0],B[j][0],repair(i,j)] for i,z in enumerate(S) for j,x in enumerate(z) if x==m],
                pair_interaction=dict(Counter(sum(p) for row in W for p in row)),
                minimum=m, minima=[[A[i][0],B[j][0],P[i][j]] for i,z in enumerate(S) for j,x in enumerate(z) if x==m],
                separate_minima=[min(map(sum,U)),min(map(sum,V))],
                mixed_total=dict(T),independent=independent,nonconstant=nonconstant,unit=bothUnit,unit_to_floor=toMin,example=first)

if __name__=='__main__':
    print(json.dumps(replay(sys.argv[1]),indent=2))
