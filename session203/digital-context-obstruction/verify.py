#!/usr/bin/env python3
import itertools,json
from fractions import Fraction as F
from counting import nine,labelled,relabel
from compiler import move,SEED

def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def inverse(M):
 n=len(M);a=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
 for c in range(n):
  j=next(i for i in range(c,n)if a[i][c]);a[c],a[j]=a[j],a[c];q=a[c][c];a[c]=[x/q for x in a[c]]
  for i in range(n):
   if i!=c:
    q=a[i][c];a[i]=[x-q*y for x,y in zip(a[i],a[c])]
 return [r[n:]for r in a]

def det(M):
 n=len(M);a=[list(map(F,row))for row in M];d=F(1)
 for c in range(n):
  j=next(i for i in range(c,n)if a[i][c])
  if j!=c:a[c],a[j]=a[j],a[c];d=-d
  q=a[c][c];d*=q;a[c]=[x/q for x in a[c]]
  for i in range(c+1,n):
   q=a[i][c];a[i]=[x-q*y for x,y in zip(a[i],a[c])]
 return d

def solve(I,y):return tuple(sum(x*z for x,z in zip(row,y))for row in I)
def glue(a,g,r):
 n=len(a);m=len(g);rows=list(a)+[0]*(m-1)
 for i in range(m):
  ii=0 if i==0 else n+i-1
  for j in range(m):
   if g[i]>>j&1:rows[ii]|=1<<(0 if j==0 else n+j-1)
 for i in range(1,n):
  for j in range(1,m):rows[i]|=1<<(n+j-1)
 rr=0 if r==0 else n+r-1
 return relabel(tuple(rows),[rr]+[i for i in range(n+m-1)if i!=rr])
def word(w):
 p=SEED
 for x in w:p=move(p,x)
 return p

def main():
 cols=(0,1,2,3,4,6,7,8);words=((),('O',),('Z',),('O','O'),('O','Z'),('Z','O'),('Z','Z'),('O','O','Z'))
 states=[word(w)for w in words];D=det([[p[i]for i in cols]for p in states]);require(D==-1152,('span determinant',D))
 intervals=[];layer={SEED}
 for m in range(9):
  require({p[4]for p in layer}==set(range(2**m+1,2**(m+1)+1)),('R interval',m))
  require(len(layer)==2**m,('R uniqueness',m))
  require(all(p[1]==p[5]==1 and p[7]==2**(m+1) and p[2]==p[4]%2 for p in layer),('register invariants',m))
  intervals.append(dict(length=m,states=len(layer),minimum_R=min(p[4]for p in layer),maximum_R=max(p[4]for p in layer)))
  layer={move(p,o)for p in layer for o in ('O','Z')}
 bases=((0,),(2,0),(4,1,2),(6,0,2));p=[nine(a)for a in bases]
 I1=inverse([(x[0],x[2],x[4],x[5])for x in p]);I2=inverse([(x[6],x[3],x[8],x[7])for x in p])
 contexts=marked=0
 for m in range(1,5):
  for g in labelled(m):
   contexts+=1;core=tuple(x>>1 for x in g[1:]);Hcore=nine(core)[0]if core else 1
   for r in range(m):
    marked+=1;out=[nine(glue(a,g,r))for a in bases];hr=solve(I1,[x[0]for x in out]);lr=solve(I2,[x[3]for x in out]);c=out[1][1]
    require(hr[0]==Hcore,('unavoidable h coefficient',g,r,hr,Hcore))
    require(all(q[1]==c*x[1]for q,x in zip(out[1:],p[1:])),('start multiplier',g,r,c,out))
    if r==0:require(lr[0]==0 and lr[2]==0 and lr[1]>=0 and lr[3]>=0,('shared-root L transfer',g,lr))
 return dict(status='PASS',scope='Finite exact controls for universal one-port readout obstruction; not F coverage or counterexample.',span_determinant=str(D),span_words=[list(w)for w in words],register_intervals=intervals,unmarked_contexts=contexts,marked_contexts=marked,context_order_limit=4)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
