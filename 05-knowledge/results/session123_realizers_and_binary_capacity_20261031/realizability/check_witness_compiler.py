"""Exact bounded controls for the witness compiler; not a Lenstra implementation."""
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def compile_ilp(v,p,q):
    need(v and all(x>0 for x in v) and p>0 and 2*p<=q,'domain')
    n=len(v);rows=[]
    for i,x in enumerate(v):
        a=[0]*n;a[i]=1;rows.append((a,x-1))
        a=[0]*n;a[i]=-1;rows.append((a,0))
    for i,x in enumerate(v):
        for j,y in enumerate(v):
            a=[0]*n;a[i]+=q*y;a[j]-=q*x
            rows.append((a,(q-p)*x-p*y))
    return rows

def interval(v,z,p,q):
    lo=max(F(q*k+p,q*x) for x,k in zip(v,z));hi=min(F(q*k+q-p,q*x) for x,k in zip(v,z))
    return (lo,hi) if lo<=hi else None

def comb(v,p,q):return [(F(q*k+p,q*v),F(q*k+q-p,q*v)) for k in range(v)]
def intersect(a,b):
    out=[];i=j=0
    while i<len(a) and j<len(b):
        lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
        if lo<=hi:out.append((lo,hi))
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return out

def literal(v,p,q):
    out=[(F(0),F(1))]
    for x in v:out=intersect(out,comb(x,p,q))
    return out

def main():
    instances=assignments=feasible=pointonly=matrix=0
    for n in range(1,5):
      for v in combinations_with_replacement(range(1,6),n):
       for q in range(2,8):
        for p in range(1,q//2+1):
         # Coprime fractions avoid redundant radii, not speed instances.
         import math
         if math.gcd(p,q)>1:continue
         instances+=1;found=[];rows=compile_ilp(v,p,q)
         for z in product(*(range(x) for x in v)):
            assignments+=1
            good=all((q*z[i]+p)*v[j]<=(q*z[j]+q-p)*v[i] for i in range(n) for j in range(n))
            ab=interval(v,z,p,q)
            need(good==(ab is not None),'fiber equivalence')
            if assignments%101==0:
                need(good==all(sum(a*t for a,t in zip(co,z))<=rhs for co,rhs in rows),'compiled matrix');matrix+=1
            if ab is not None:
                found.append(ab);t=ab[0];need(all(F(p,q)<=x*t%1<=1-F(p,q) for x in v),'literal output witness')
                need(t.denominator<=q*max(v),'short denominator')
                need(tuple((x*t).numerator//(x*t).denominator for x in v)==z,'restored floors')
         expected=literal(v,p,q);found.sort()
         need(found==expected,'entire interval fiber equality')
         need(literal(tuple(reversed(v)),p,q)==expected,'label reversal')
         feasible+=bool(found);pointonly+=bool(found) and all(a==b for a,b in found)
    need(all(literal(v,1,3) for v in [(1,2),(1,3),(2,3)]),'pairwise hostile positive sides')
    need(not literal((1,2,3),1,3),'joint hostile failure')
    need(literal((1,3),1,2)==[(F(1,2),F(1,2))],'point-only equality')
    scaling=0
    for v in [(1,2),(1,3,4),(2,5,7),(1,2,3)]:
      for d in (2,3,7):
       out=literal(v,1,4);scaled=literal(tuple(d*x for x in v),1,4)
       want=sorted((F(a+j,d),F(b+j,d)) for j in range(d) for a,b in out)
       need(scaled==want,'periodic scaling');scaling+=1
    return dict(status='PASS',instances=instances,integer_assignments=assignments,feasible_instances=feasible,point_only_instances=pointonly,matrix_controls=matrix,scaling_controls=scaling,jointness_hostile='all three pairs feasible; triple infeasible',algorithm_scope='exact small exhaustive reference, not polynomial-time ILP implementation')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
