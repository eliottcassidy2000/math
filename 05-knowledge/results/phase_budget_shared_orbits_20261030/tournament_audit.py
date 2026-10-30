#!/usr/bin/env python3
"""Signed-DSU, independent F2 elimination, and literal common-bank controls."""
from itertools import combinations,product
from pathlib import Path
import json

def need(b,msg):
    if not b:raise RuntimeError(msg)

class DSU:
    def __init__(self,n):self.p=list(range(n));self.d=[0]*n;self.sz=[1]*n;self.ok=True
    def find(self,x):
        if self.p[x]!=x:
            p,d=self.find(self.p[x]);self.d[x]^=d;self.p[x]=p
        return self.p[x],self.d[x]
    def join(self,a,b,c):
        ra,da=self.find(a);rb,db=self.find(b)
        if ra==rb:self.ok &= da^db==c;return
        if self.sz[ra]<self.sz[rb]:ra,rb=rb,ra
        self.p[rb]=ra;self.d[rb]=da^db^c;self.sz[ra]+=self.sz[rb]

def solve(n,perms,path,weights=None):
    edges=list(combinations(range(n),2));M=len(edges);idx={e:i for i,e in enumerate(edges)}
    rows=[];ds=DSU(M+1)
    def add(a,b,c):
        ds.join(a,b,c)
        row=(1<<a) if a<M else 0
        if b<M:row^=1<<b
        rows.append((row,c))
    for g in perms:
        for i,(u,v) in enumerate(edges):
            gu,gv=g[u],g[v];j=idx[tuple(sorted((gu,gv)))];add(i,j,int(gu>gv))
    for u,v in zip(path,path[1:]):add(idx[tuple(sorted((u,v)))],M,int(u<v))
    piv={};ok=True
    for row,rhs in rows:
        while row:
            k=row.bit_length()-1
            if k not in piv:piv[k]=(row,rhs);break
            rr,bb=piv[k];row^=rr;rhs^=bb
        if not row and rhs:ok=False
    need(ok==ds.ok,'independent consistency')
    if not ok:return None
    dim=len({ds.find(i)[0] for i in range(M+1)})-1
    need(dim==M-len(piv),'independent rank')
    if weights is None:return dim
    accum={};zr,zd=ds.find(M)
    for e,(w0,w1) in enumerate(weights):
        root,eta=ds.find(e)
        if root not in accum:accum[root]=[1,1]
        accum[root][0]*=(w0,w1)[eta]
        accum[root][1]*=(w0,w1)[1^eta]
    Z=1
    for root,pair in accum.items():Z*=pair[zd] if root==zr else sum(pair)
    return dim,Z

def walecki(t):
    raw=[];lo=0;hi=t-1
    while lo<=hi:
        raw.append(lo);lo+=1
        if lo<=hi:raw.append(hi);hi-=1
    return [2*x%t for x in raw]

def common_arc(u,v,r,m,outer,inner):
    i,j=u%r,u//r;k,l=v%r,v//r
    return ((k-i)%r in outer) if i!=k else ((l-j)%m in inner)

reports=[];weighted_checks=0
for r,m in ((3,3),(3,5),(3,9),(5,5),(5,7),(7,7),(3,19),(11,3)):
    n=r*m;a=[(x+1)%n for x in range(n)];v=[(x+r)%n if x%r==0 else x for x in range(n)]
    vi=[0]*n
    for x,y in enumerate(v):vi[y]=x
    b=[v[a[vi[x]]] for x in range(n)]
    natural=list(range(n));selected=[i+r*j for i in walecki(r) for j in walecki(m)]
    need(len(set(selected))==n and selected[0]==0 and selected[-1]==n-1,'retained endpoint ports')
    bank=solve(n,(a,b),[]);old=solve(n,(a,b),natural);new=solve(n,(a,b),selected)
    need((bank,old,new)==((r+m-2)//2,(r+m-4)//2,0),'all orbit rank formulas')
    outer={4*k%r for k in range(1,(r+1)//2)};inner={4*k%m for k in range(1,(m+1)//2)}
    arc=lambda u,v:common_arc(u,v,r,m,outer,inner)
    need(all(arc(u,v) for u,v in zip(selected,selected[1:])),'explicit selector legality')
    need(all(arc(u,v)==arc(g[u],g[v]) for g in (a,b) for u,v in combinations(range(n),2)),
         'actual common bank, not just difference variables')
    # Literal enumerate ONLY the proved common-bank universe; no hidden quotient.
    natural_count=selected_count=0;weighted=[0,0,0]
    weights=[(e%3+1,(-1)**e*(e%5+1)) for e in range(n*(n-1)//2)] if n<=15 else None
    for bits in product((0,1),repeat=bank):
        ro=(r-1)//2
        O={k if not bits[k-1] else r-k for k in range(1,ro+1)}
        I={k if not bits[ro+k-1] else m-k for k in range(1,(m+1)//2)}
        def legal(p):return all(common_arc(u,v,r,m,O,I) for u,v in zip(p,p[1:]))
        ln,ls=legal(natural),legal(selected)
        natural_count+=ln;selected_count+=ls
        if weights is not None:
            mass=1
            for (u,v),ww in zip(combinations(range(n),2),weights):
                mass*=ww[int(common_arc(u,v,r,m,O,I))]
            weighted[0]+=mass;weighted[1]+=mass*ln;weighted[2]+=mass*ls
    need(natural_count==2**old and selected_count==1,'literal bank census')
    if weights is not None:
        for path,Z in zip(([],natural,selected),weighted):
            need(solve(n,(a,b),path,weights)[1]==Z,'exact signed source-weight sum')
            weighted_checks+=1
    reports.append(dict(r=r,m=m,n=n,unrestricted_dimension=bank,natural_dimension=old,
                        selected_dimension=new,selected_path=selected,
                        literal_banks=2**bank,natural_intersection=natural_count,selected_intersection=selected_count))
# Odd directed-cycle signs and contradictory retained pins must be detected.
d=DSU(3);d.join(0,1,0);d.join(1,2,0);d.join(0,2,1);need(not d.ok,'signed cycle hostile')
need(solve(3,([1,0,2],),[]) is None,'arc-reversing permutation hostile')
out=dict(status='PROVED formulas, independently controlled here',families=reports,
         weighted_intersections=weighted_checks,
         methods=['signed union-find','binary Gaussian elimination','literal common banks'])
Path(__file__).with_name('tournament_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([{k:v for k,v in x.items() if k!='selected_path'} for x in reports],indent=2))
