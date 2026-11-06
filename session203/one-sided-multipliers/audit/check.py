#!/usr/bin/env python3
from itertools import permutations
import json

def require(x,m):
 if not x:raise RuntimeError(m)

def rows(n,k):
 a=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j)if k>>b&1 else(j,i);a[u]|=1<<v;b+=1
 return a

def data(a):
 n=len(a);N=1<<n;d=[[0]*n for _ in range(N)];s=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for S in range(1,N):
  for v in range(n):
   if not(S>>v&1):continue
   Z=S^(1<<v)
   if Z==0:d[S][v]=s[S][v]=1
   else:
    d[S][v]=sum(d[Z][w]for w in range(n)if Z>>w&1 and a[w]>>v&1)
    s[S][v]=sum(s[Z][w]for w in range(n)if Z>>w&1 and a[v]>>w&1)
  H[S]=sum(d[S])
 L=sum(d[S][0]*H[(N-1)^S]for S in range(1,N)if S&1)
 R=sum(s[S][0]*H[(N-1)^S]for S in range(1,N)if S&1)
 return H[-1],d[-1][0],s[-1][0],L,R,H[(N-1)^1]

def edge_sector(a):
 if a==1:return[2,0]
 bits=[int(c)for c in bin(a-1)[2:]];r=[2,0]+[0]*len(bits)
 for i,b in enumerate(bits,2):
  r[i]|=1
  if b:r[1]|=1<<i
  else:r[i]|=2
  for j in range(2,i):r[j]|=1<<i
 return r

def product(G,U):
 g=len(G);r=G[:]+[0]*(len(U)-1)
 for j in range(1,len(U)):
  y=g+j-1
  if U[0]>>j&1:r[0]|=1<<y
  else:r[y]|=1
  for x in range(1,g):r[x]|=1<<y
  for i in range(1,j):
   x=g+i-1
   if U[i]>>j&1:r[x]|=1<<y
   else:r[y]|=1<<x
 return r

def det(a):
 out=0;n=len(a)
 for p in permutations(range(n)):
  inv=sum(p[i]>p[j]for i in range(n)for j in range(i+1,n));v=(-1)**inv
  for i in range(n):v*=a[i][p[i]]
  out+=v
 return out

def main():
 factors=[1,3,5,9,17,33];profiles=[];checks=0
 for a in factors:
  G=edge_sector(a);h,e,s,L,R,d=data(G)
  require((h,e,R,d)==(2*a-1,a-1,2*a,a),('core',a,h,e,R,d));profiles.append({'a':a,'rows':G,'profile_h_e_R_d':[h,e,R,d]})
  for n in range(1,5):
   for code in range(1<<(n*(n-1)//2)):
    U=rows(n,code);hu,eu,su,Lu,Ru,du=data(U)
    observed=data(product(G,U))[0];parent=data(product([2,0],U))[0]
    require(observed==a*(hu+Lu)-su==a*parent+(a-1)*su,('product',a,n,code));checks+=1
 # Independent profiles for the four asserted coordinate-basis test inputs.
 tests=[[0],[0,1],[0,1,3],[2,4,1]];basis=[]
 for U in tests:
  h,e,s,L,R,d=data(U);basis.append([h,L,s,d])
 require(basis==[[1,1,1,1],[1,2,0,1],[1,4,0,1],[3,3,1,1]],'basis values')
 determinant=det(basis);require(determinant!=0,'basis rank')
 slacks=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   G=rows(n,code);h,e,s,L,R,d=data(G);f=R+e-h-d
   require(f>=0 and ((f==0)==(G[0].bit_count()<=1)),('zero slack control',n,code));slacks+=1
 print(json.dumps({'status':'PASS','literal_product_checks':checks,'rooted_inputs_per_factor':75,'factors':factors,'core_profiles':profiles,'basis':basis,'basis_determinant':determinant,'independent_zero_slack_controls':slacks,'scope':'One-sided rooted cores only; all-order classification relies on the reviewed proof and inherited universal zero-slack theorem'},indent=2,sort_keys=True))
if __name__=='__main__':main()
