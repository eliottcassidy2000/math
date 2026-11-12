from math import gcd
h,q,u,v,e,f=3357,50004,25915,20732,1865,1492
M=1<<19
g=[0]*(M+2)
for a in range(1,len(g)):
 x=a//2;g[a]=g[x]+(0 if a&1 else x+(x&1))
print('gcd',gcd(h,q),'horizontal even residual',(u+e)%9,'vertical',(v+f)%9)
for t in range(-4000,4001):
 rhs=4-q*t-(u-v-e+f)
 if rhs%h or t%2:continue
 D=rhs//h;hits=[]
 for A in range(max(1,-1-t),min(M-1,M-1-t),2):
  if A%2==0:continue
  B=A+1+t
  if B<=0 or B%2:continue
  if g[A+1]-g[A]+g[B-1]-g[B]==D:
   hits.append((A,B,D))
   if len(hits)==3:break
 print(t,D,hits,flush=True)
