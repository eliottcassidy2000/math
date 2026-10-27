#!/usr/bin/env python3
import itertools,json
from collections import Counter
from fractions import Fraction

def need(x,msg):
    if not x:raise RuntimeError(msg)
def tournament(n,D):
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            x,y=(j,i) if (i,j) in D else (i,j);out[x]|=1<<y
    return out
def fas(out):
    n=len(out);d=[0]*(1<<n)
    for S in range(1,1<<n):
        z=S;best=n*n
        while z:
            v=(z&-z).bit_length()-1;z&=z-1;R=S^(1<<v)
            best=min(best,d[R]+(out[v]&R).bit_count())
        d[S]=best
    return d[-1]
def cost(out,p):return sum(out[v]>>u&1 for i,u in enumerate(p) for v in p[i+1:])
rows=[]
for n in range(5,41):
    e=2*((n-1)//2);source=tournament(n,{(1,e)});U=tournament(n,{(i,i+2) for i in range(n-2)})
    P=list(range(0,n,2))+list(range(1,n,2));pe={tuple(sorted(z)) for z in zip(P,P[1:])}
    need(all(source[u]>>v&1 for u,v in zip(P,P[1:])),'illegal P')
    for i,j in itertools.combinations(range(n),2):
        target_arc=bool(source[i]>>j&1) if (i,j) in pe else not bool(source[i]>>j&1)
        need(target_arc==bool(U[j]>>i&1),'literal complement mismatch')
    tris=[]
    for a,b,c in itertools.combinations(range(n),3):
        if all(sum(U[u]>>v&1 for v in (a,b,c) if v!=u)==1 for u in (a,b,c)):tris.append((a,b,c))
    expected=[(i,i+1,i+2) for i in range(n-2)]+[(i,i+2,i+4) for i in range(n-4)]
    need(set(tris)==set(expected),'triangle classification')
    loads=Counter();total=Fraction(0)
    for t in tris:
        w=Fraction(1,2) if t[1]-t[0]==1 else Fraction(1,4);total+=w
        for a,b in itertools.combinations(t,2):loads[a,b]+=w
    need(max(loads.values())<=1 and total==Fraction(3*n-8,4),'triangle packing')
    p=[]
    for a in range(0,n//4*4,4):p.extend([a+2,a,a+3,a+1])
    p.extend(range(n//4*4,n));f=n-2-n//4
    need(cost(U,p)==f,'rank upper')
    row={'n':n,'FAS_U_proved':f,'source_FAS':1,'complement_difference':f-1,'triangle_weight':str(total),'max_edge_triangle_load':str(max(loads.values())),'subset_DP':n<=18}
    if n<=18:need(fas(U)==f and fas(source)==1,'subset DP')
    rows.append(row)
need(next(r for r in rows if r['n']==10)['complement_difference']==5,'n10 hostile')
print(json.dumps({'status':'FINITE-EXACT controls for all-n proof','rows':rows,'n10_path':[0,2,4,6,8,1,3,5,7,9],'n10_source_exception':[1,8]},indent=2,sort_keys=True))
