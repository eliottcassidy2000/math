from collections import Counter
from itertools import combinations
from pathlib import Path
from nested_core_check import graph,orbit_chart
import json

def cycle_type(p):
    seen=set(); lens=[]
    for i in range(len(p)):
        if i in seen:continue
        c=0;j=i
        while j not in seen:seen.add(j);c+=1;j=p[j]
        lens.append(c)
    return tuple(sorted(lens))

def route(f,g):
    return [g[x] for x in f]

def invert(q):
    p=[0]*len(q)
    for i,x in enumerate(q):p[x]=i
    return p

def moments():
    ans=[]
    # Near-layer quotient spectra are complements of C5 and Petersen.
    # The C5 quotient is itself C5; its trace powers are computed below.
    C5=[{(i-1)%5,(i+1)%5} for i in range(5)]
    def trp(A,k):
        out=0
        for s in range(len(A)):
            z=[0]*len(A);z[s]=1
            for _ in range(k):z=[sum(z[t] for t in A[j]) for j in range(len(A))]
            out+=z[s]
        return out
    def qtrace(d,p):
        from math import comb
        if d==2:return (-1)**p*sum(comb(p,j)*trp(C5,j) for j in range(p+1))
        return (-7)**p+5*(1**p)+4*((-2)**p)
    for d,m,n,X,Y in [(2,5,275,2,2970),(3,10,540,6,2700)]:
        dy=57-m;w=n-m; rem=Y-1-w
        # Transfer trace on W is -w because tr(X|W)=0.
        num8=8*rem+w-dy
        assert num8%15==0
        p=num8//15;q=rem-p
        hX=lambda i:{0:n,1:0,2:n*X,3:0,4:n*X*(2*X-1),5:0}[i]
        hY=lambda i:{0:Y,1:0,2:Y*dy,3:0,4:Y*dy*(2*dy-1)}[i]
        from math import comb
        def total(i):
            return dy**i+p*7**i+q*(-8)**i+(-1)**i*sum(comb(i,j)*hX(j) for j in range(i+1))-qtrace(d,i)
        assert [total(i) for i in range(5)]==[hY(i) for i in range(5)]
        assert total(5)%10==0
        ans.append([d,n,Y,p,q,total(5)//10])
    return ans

def main():
    C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
    a,b=0,28; m=56
    T=[]
    for c in range(m):
        T.append(list(C[a,b]) if c==a else list(range(m)) if c==b else route(C[a,c],C[c,b]))
    assert all(len(set(t))==m for t in T)
    assert all(len({t[j] for t in T})==m for j in range(m))
    shifts={k:sum(all(row[(j+k)%m]==(row[j]+k)%m for j in range(m)) for row in T)
            for k in [2,4,8,14,28]}
    rel=Counter()
    intercalates=0
    for p,q in combinations(T,2):
        ts=cycle_type(route(p,invert(q)))
        assert 1 not in ts
        intercalates+=ts.count(2)
        rel[ts]+=1
    assert intercalates==1950 and shifts[4]==52 and shifts[8]==54
    print(json.dumps({"route_latin_order":m,"intercalates":intercalates,
         "commuting_rows":shifts,"relative_types":type_counts(rel),
         "outer_shell_moments":moments()},indent=2))
def type_counts(rel):return {",".join(map(str,k)):n for k,n in sorted(rel.items())}
if __name__=="__main__":main()
