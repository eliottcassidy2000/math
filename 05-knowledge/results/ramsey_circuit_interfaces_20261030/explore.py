"""Initial exact weighted-tournament probe, not a standalone proof certificate."""
from itertools import combinations
from fractions import Fraction as F
from collections import Counter

def solve(a,b):
    n=len(b);m=[[F(x) for x in row]+[F(y)] for row,y in zip(a,b)]
    for j in range(n):
        p=next((i for i in range(j,n) if m[i][j]),None)
        if p is None:return None
        m[j],m[p]=m[p],m[j];c=m[j][j];m[j]=[x/c for x in m[j]]
        for i in range(n):
            if i!=j:
                c=m[i][j];m[i]=[x-c*y for x,y in zip(m[i],m[j])]
    return tuple(row[-1] for row in m)

def tournament(n,steps):
    return [[int((j-i)%n in steps) for j in range(n)] for i in range(n)]

def transitive(a,s):
    return all(not (a[i][j] and a[j][k] and a[k][i]) and not(a[j][i] and a[i][k] and a[k][j])
               for i,j,k in combinations(s,3))

def vertices(a):
    n=len(a);sets=[tuple(i for i in range(n) if mask>>i&1) for mask in range(1,1<<n)]
    ts=[s for s in sets if transitive(a,s)]
    maximal=[s for s in ts if not any(set(s)<set(t) for t in ts)]
    rows=[tuple(int(i in s) for i in range(n)) for s in maximal]
    rhs=[1]*len(rows)
    rows += [tuple(-int(i==j) for i in range(n)) for j in range(n)];rhs += [0]*n
    vs=set()
    for sub in combinations(range(len(rows)),n):
        x=solve([rows[i] for i in sub],[rhs[i] for i in sub])
        if x is not None and all(sum(F(u)*v for u,v in zip(row,x))<=b for row,b in zip(rows,rhs)):vs.add(x)
    return maximal,vs

if __name__=='__main__':
    a=tournament(5,{1,2});ms,vs=vertices(a)
    print('maximal',ms,'vertices',len(vs))
    print(Counter(tuple(sorted(v)) for v in vs))
    cap=24;f=[0,1];c=[0,1];plans=[None,None]
    for k in range(2,cap+1):
        best=0;arg=None
        for a in range(k):
            for b in range(min(k-a,k-1)+1):
                for cc in range(min(k-a-b,k-1)+1):
                    for d in range(min(k-b-cc,k-a,k-cc,k-1)+1):
                        e=min(k-a-b,k-a-d,k-cc-d,k-1)
                        val=f[a]+f[b]+f[cc]+f[d]+f[e]
                        if val>best:best=val;arg=(a,b,cc,d,e)
        f.append(best);plans.append(arg)
        cv=max(c[a]+c[b]+c[min(k-a,k-b,k-1)] for a in range(k) for b in range(k-a+1) if b<k)
        c.append(cv)
    print('T5 capacity',f)
    print('C3 capacity',c)
    print('strict improvements',[(k,f[k],c[k],plans[k]) for k in range(1,cap+1) if f[k]>c[k]])
