from math import gcd,lcm,isqrt
def lists(A, bound):
 P=[i for i in range(3,bound+1,2) if all(i%j for j in range(3,isqrt(i)+1,2))]
 L=[]
 for p in P:
  x=1;m=0;b=None
  while True:
   if (A*x+23)%p==0:b=m
   x=x*2%p;m+=1
   if x==1:break
  if b is not None:L.append((p,b,m))
 return L
def crt(a,m,b,n):
 g=gcd(m,n)
 u,v=m//g,n//g
 if (b-a)%g:return
 return (a+m*((b-a)//g*pow(u,-1,v)%v),m*v)
def prune(C):
 C=set(C)
 return sorted([(a,m) for a,m in C if not any(m%n==0 and a%n==b and n<m for b,n in C)],key=lambda x:(x[1],x[0]))
def tab(C):
 T={(0,1):1}
 for b,n in prune(C):
  # Multiply protected indicators 1-[k=b mod n].
  U=T.copy()
  for (a,m),c in T.items():
   t=crt(a,m,b,n)
   if t is not None:
    U[t]=U.get(t,0)-c
    if not U[t]: del U[t]
  T=U
 return T
def count(C,M):
 return sum(c*(M//m) for (a,m),c in tab(C).items())
def minimize(C,period):
 p=2; P=[]
 while p*p<=period:
  while period%p==0:
   period//=p
   if p not in P:P.append(p)
  p+=1
 if period>1:P.append(period)
 M=lcm(*(m for a,m in C))
 for p in P:
  while M%p==0:
   d=M//p
   E=lcm(*(m for a,m in C))
   if count(C,E)==count(C+[((a+d)%m,m) for a,m in C],E):M=d
   else:break
 return M
def go(bound):
 A=[(b,m) for p,b,m in lists(99,bound)]
 B=[(b,m) for p,b,m in lists(100,bound)]
 M=lcm(*(m for a,m in A+B))
 for C in (A,B,A+B):
  E=lcm(*(m for a,m in C)); D=minimize(C,E); n=count(C,E)
  print(bound,len(C),len(prune(C)),len(tab(C)),E,D,n//(E//D),flush=True)
def need(ok, label):
 if not ok:raise AssertionError(label)
def run():
 from fractions import Fraction as F
 from collections import Counter
 A=[(r,m) for p,r,m in lists(99,100)]
 B=[(r,m) for p,r,m in lists(100,100)]
 witnesses=[
  {2:5,3:21,5:36,11:62,13:29,29:149,41:69},
  {2:3,3:1,5:7,7:27,11:131,13:33,23:61,41:267},
  {2:69,3:69,5:141,7:4397,11:533,13:149,23:1469,29:3309,41:509}]
 for C,M,n,W in zip((A,B,A+B),(61209720,169909740,9854764920),(32353278,39935280,184192008),witnesses):
  def good(k):return all(k%m!=r for r,m in C)
  E=lcm(*(m for r,m in C))
  need(minimize(C,E)==M and count(C,E)*M==n*E,(C,M,n))
  need(all(good(k) and not good(k+M//p) for p,k in W.items()),W)
  print(M,n,F(n,M))
 C=A+B;low=[(r,m) for r,m in C if 360%m==0]
 S=[k for k in range(360) if all(k%m!=r for r,m in low)]
 tags=Counter((k%3==2,k%5==2) for k in S)
 need(tags=={(False,False):3,(False,True):1,(True,False):4,(True,True):2},tags)
 factor=F(22,23)*F(28,29)*F(39,41)*F(11,13)
 D=factor*sum(F(n,360)*(F(10,11) if a else 1)*(F(6,7) if b else 1) for (a,b),n in tags.items())
 need(D==F(187,10005),D)
 print('360 kernel:',S)
 coarse=[k for k in range(60) if all((a*pow(2,k,p)+23)%p for a in (99,100) for p in (3,5,7,11,13))]
 need(coarse==[9,17,21,29,33,41,53,57],coarse)
 lifts={k:next(j for j in range(k,2760,60) if j%8==1 and j%23==2) for k in coarse}
 for j in lifts.values():
  need((99*pow(2,j,17)+23)%17==0 and (100*pow(2,j,47)+23)%47==0,j)
 print('13# common residues:',coarse,'fixed composite lifts:',lifts)
 def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
 P=[p for p in range(29980,30050) if prime(p)]
 need(P==[29983,29989,30011,30013,30029,30047],P)
 need(99*2**9+23==17*2983 and 100*2**9+23==181*283,'first above 30000')
 print('nearby primes:',P)
if __name__=="__main__":run()
