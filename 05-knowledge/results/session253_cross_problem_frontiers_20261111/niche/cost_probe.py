from math import comb,isqrt
import json

def factor(n):
 out={};p=2
 while p*p<=n:
  while n%p==0:out[p]=out.get(p,0)+1;n//=p
  p+=1
 if n>1:out[n]=1
 return out

def period(m,k):
 P=m
 for p in factor(m):
  q=p
  while q<=k:P*=p;q*=p
 return P

def costs(m,ks,starts):
 INF=10**100;c=[0]+[INF]*(m-1)
 for k,b in zip(ks,starts):
  a=[INF]*m
  for n in range(b,b+period(m,k)):
   v=comb(n,k) if n>=k else 0
   a[v%m]=min(a[v%m],v)
  c=[min(c[s]+a[(r-s)%m] for s in range(m)) for r in range(m)]
 return c

def direct(N,ks,starts):
 S={0}
 for k,b in zip(ks,starts):
  vals=[];n=b
  while (comb(n,k) if n>=k else 0)<=N:
   vals.append(comb(n,k) if n>=k else 0);n+=1
  S={a+b for a in S for b in vals if a+b<=N}
 return S
receipts=[]
for ks,starts in [([2,4,6,8],[2,3,5,7]),([4,6],[3,5])]:
 for m in range(2,61):
  c=costs(m,ks,starts);S=direct(m-1,ks,starts)
  for N in range(m):
   if (c[N]==N)!=(N in S):raise Exception((ks,m,N,c[N],S))
  receipts.append((ks,m,sum(x<10**100 for x in c)))
print(json.dumps({'costed_fiber_moduli':len(receipts),'all_pass':True,'restricted_hostile':{'ranks':[4,6],'modulus':4,'target':3,'cost':costs(4,[4,6],[3,5])[3]},'sun_modular_images_full':all(n==m for ks,m,n in receipts if len(ks)==4)},indent=2))
