from math import gcd
from collections import Counter

def fib(n):
 a,b=0,1
 for _ in range(n):a,b=b,a+b
 return a

def words(n):
 def rec(i,w,prev):
  if i==n:
   if not ((w&1) and prev):yield w
   return
  yield from rec(i+1,w,0)
  if not prev:yield from rec(i+1,w|(1<<i),1)
 yield from rec(0,0,0)

def setup(n):
 a,b=fib(n-1)-1,fib(n);L=abs(a*a+a*b-b*b);g=gcd(a,b);M=L//(g*g);r=(-a//g*pow(b//g,-1,M))%M if M>1 else 0
 pw=[pow(r,i,M) for i in range(n)]
 def key(w):
  u=sum(pw[i] for i in range(n) if w>>i&1)%M; t=0
  if g==2:
   for i in range(n):
    if w>>i&1:t^=(1,2,3)[i%3]
  return t,u
 ws=list(words(n));keys=[key(w) for w in ws];idx={v:i for i,v in enumerate(keys)}
 if len(idx)!=L:raise ValueError('reader')
 R=[idx[key(sum(((w>>i)&1)<<(n-i if i else 0) for i in range(n)))] for w in ws]
 T=[idx[(t^1 if g==2 else 0,(u+1)%M)] for t,u in keys]
 U=[idx[(t^2 if g==2 else 0,(u+r)%M)] for t,u in keys]
 return L,g,M,r,keys,R,T,U

def compress(vals):
 d={};out=[]
 for x in vals:
  if x not in d:d[x]=len(d)
  out.append(d[x])
 return out

def require(test, message):
 if not test: raise RuntimeError(message)

def order(r,m):
 x=1
 for k in range(1,m*m+1):
  x=x*r%m
  if x==1:return k
 raise RuntimeError("order search")

def run(n):
 L,g,M,r,keys,R,T,U=setup(n)
 idx={v:i for i,v in enumerate(keys)}
 C=[idx[((0,2,3,1)[t] if g==2 else 0,u*r%M)] for t,u in keys]
 require(len(keys)==L,"size")
 require(all(R[R[i]]==i for i in range(L)),"reflection involution")
 require(all(R[C[R[C[i]]]]==i for i in range(L)),"dihedral relation")
 if g==2:
  require(all(keys[R[i]][0]==(0,1,3,2)[keys[i][0]] for i in range(L)),"F4 Frobenius")
 factors=[];z=M;p=3
 while p*p<=z:
  q=1
  while z%p==0:z//=p;q*=p
  if q>1:factors.append(q)
  p+=2
 if z>1:factors.append(z)
 results=[]
 for m in sorted(set(factors+[M])):
  obs=[u%m for t,u in keys]
  double=[(obs[i],obs[R[i]]) for i in range(L)]
  inv=pow(r,-1,m) if m>1 else 0
  require(all(double[R[i]]==double[i][::-1] for i in range(L)),"double reflection")
  require(all(double[C[i]]==(r*double[i][0]%m,inv*double[i][1]%m) for i in range(L)),"double rotation")
  c=compress(obs);seq=[len(set(c))]
  while True:
   d=compress([(c[i],c[R[i]],c[T[i]]) for i in range(L)]);num=len(set(d))
   if num==seq[-1]:break
   seq.append(num);c=d
  require(seq[-1]==L,"observed full TR closure")
  prof=compress([(obs[i],obs[R[i]],obs[R[T[i]]]) for i in range(L)])
  seen={};witness=None
  for i,pair in enumerate(double):
   if pair in seen and double[T[i]]!=double[T[seen[pair]]]:
    j=seen[pair];witness={"states":[keys[j],keys[i]],"double":pair,"after_add1":[double[T[j]],double[T[i]]]};break
   seen[pair]=i
  results.append({"modulus":m,"root":r%m,"double_size":len(set(double)),"TR_refinement_sizes":seq,"triple_profile_size":len(set(prof)),"triple_max_fibre":max(Counter(prof).values()),"double_addition_hostile":witness})
 return {"period":n,"size":L,"g":g,"odd_modulus":M,"odd_root":r,"charts":results}

if __name__=='__main__':
 import json,sys
 normal=[]
 for n in range(1,202,2):
  a,b=fib(n-1)-1,fib(n);g=gcd(a,b);L=abs(a*a+a*b-b*b)
  require(g==(2 if n%3==0 else 1),"gcd normal form")
  M=L//(g*g)
  require(M%2==1,"odd modulus")
  if M>1:
   r=(-a//g*pow(b//g,-1,M))%M
   require((r*r-r-1)%M==0 and pow(r,n,M)==1,"normal-form root")
  normal.append(n)
 require((37*37-37-1)%121==0,"root121")
 require(order(37,121)==55 and order(4,11)==5,"prime-power clock")
 require(order(84,121)==110,"conjugate clock")
 require(11!=0 and 11*11%121==0,"nilpotent")
 receipt={"normal_form_odd_periods":normal,"nonreduced":{"modulus":121,"phi":37,"order":55,"conjugate":84,"conjugate_order":110,"reduced_order":5},"census":[run(n) for n in (9,15,21,25)]}
 print(json.dumps(receipt,indent=2))
