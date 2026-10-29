"""Direct exact determinant and inverse tests for the minimal port carrier."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import json

def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return list(map(list,zip(*a)))
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def invdet(a):
 n=len(a);b=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)];d=F(1)
 for j in range(n):
  r=next((r for r in range(j,n) if b[r][j]),None)
  if r is None:return None,F(0)
  if r!=j:b[r],b[j]=b[j],b[r];d=-d
  v=b[j][j];d*=v;b[j]=[x/v for x in b[j]]
  for r in range(n):
   if r!=j:
    v=b[r][j];b[r]=[x-v*y for x,y in zip(b[r],b[j])]
 return [r[n:] for r in b],d
def diag(values):
 n=2*len(values);a=[[0]*n for _ in range(n)]
 for i,v in enumerate(values):a[2*i][2*i+1]=v;a[2*i+1][2*i]=-v
 return a
def main():
 gates=0
 m=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
 for ports in ([[1,0],[0,1],[1,1],[2,-1]],[[1,0,2],[0,1,0],[1,1,2],[2,-1,4]]):
  a,d=invdet(m);g=mm(tr(ports),mm(a,ports));n=len(g)
  for i,j in combinations(range(n),2):
   vals=[]
   for t in (1,-1):
    e=[[0]*n for _ in range(n)];e[i][j]=t;e[j][i]=-t
    k=add(m,mm(ports,mm(e,tr(ports))));ki,kd=invdet(k);vals.append(kd)
    if kd:
     f=add([[int(i==j) for j in range(n)] for i in range(n)],mm(g,e))
     fi,fd=invdet(f)
     if mm(fi,g)!=mm(tr(ports),mm(ki,ports)) or kd!=d*fd:raise RuntimeError('resolvent')
     gates+=1
   if (vals[1]-vals[0])/(4*d)!=g[i][j]:raise RuntimeError('separation')
   gates+=1
 p=[[int(i==j) for j in range(2)] for i in range(8)]
 carriers=[]
 for model in (diag([1,9,1,1]),diag([1,3,3,1])):
  a,d=invdet(model);carriers.append((d,mm(tr(p),mm(a,p))))
 if carriers[0]!=carriers[1]:raise RuntimeError('hidden integral group hostile')
 out={'gates':gates+1,'independent_and_dependent_port_families':2,'same_response_distinct_group_hostile':True}
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
