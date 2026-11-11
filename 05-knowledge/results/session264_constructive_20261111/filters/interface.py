from itertools import permutations
from collections import Counter
from math import comb
from explore import wrap,cycles

def wrapper_core(r):
 rows=[];side=[]
 for t in range(r):
  k=len(rows);u=k;v=k+1
  rows=[row|(1<<v) for row in rows]
  rows.extend([(1<<k)-1,1<<u]);side.extend([1,0])
 return rows,side

def skeletons(rows,side):
 k=len(rows);out=Counter()
 for perm in permutations(range(k)):
  for cut in range(k+1):
   parts=[perm[:cut],perm[cut:]];qs=[];ms=[];valid=True
   for path in parts:
    q=m=0
    for gap in range(len(path)+1):
     L=path[gap-1] if gap else None;R=path[gap] if gap<len(path) else None
     nonempty=(L is None or side[L]) and (R is None or not side[R])
     empty=(L is None or R is None or rows[L]>>R&1)
     if not (nonempty or empty):valid=False;break
     q+=bool(nonempty);m+=not empty
    qs.append(q);ms.append(m)
    if not valid:break
   if valid:out[(cut,*qs,*ms)]+=1
 return out

def alphabet(n,cs,a,b):
 p=[comb(n+1,i)*a**(n+1-i)*b**i for i in range(n+2)]
 for j,c in cs.items():
  for i in range(n-j+1):
   v=2*c*comb(n-j,i)*a**(n-j-i)*b**i
   p[i]+=a*v;p[i+j+1]+=b*v
 return p

def compile(n,cs,sks,k):
 p=[0]*(n+k+2)
 for (cut,q1,q0,m1,m0),mul in sks.items():
  for i in range(m1+1):
   for j in range(m0+1):
    v=alphabet(n,cs,q0-j,q1-i);coeff=mul*(-1)**(i+j)*comb(m1,i)*comb(m0,j)
    for l,x in enumerate(v):p[l+cut]+=coeff*x
 return p

def rem(p,m=16):return [sum(v*(-1)**(i//(m//2)) for i,v in enumerate(p) if i%(m//2)==j) for j in range(m//2)]

def inputs(r,lo,hi,path):
 rows,side=wrapper_core(r);sks=skeletons(rows,side);print('core',r,rows,side,'categories',len(sks),'skeletons',sum(sks.values()),flush=True)
 with open(path,'w') as f:
  for n in range(lo,hi+1):
   b=compile(n,{},sks,2*r);rb=rem(b);f.write(str(n)+' '+' '.join(map(str,rb))+'\n')
   for d in range(1,n):
    v=compile(n,{j:comb(d-1,j-2) for j in range(2,n+1,2) if d-1>=j-2},sks,2*r)
    f.write(' '.join(str(x-y) for x,y in zip(rem(v),rb))+'\n')
if __name__=='__main__':
 import sys
 if len(sys.argv)>1:inputs(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4])
 else:
  sks=skeletons(*wrapper_core(1))
  for n in range(1,10):
   for c in [{},cycles('10'*(n//2)+'1'*(n%2))]:assert compile(n,c,sks,2)==wrap(n,c)
  print('one-wrap formula equality passes')
