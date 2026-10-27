#!/usr/bin/env python3
from fractions import Fraction as F
from math import factorial
from itertools import combinations
import json

def check(x,msg):
 if not x:raise RuntimeError(msg)

def shortcut_sign(n):
 clock=0
 while n!=1:
  n=n//2 if not n%2 else (3*n+1)//2;clock+=1
  check(clock<10000,'shortcut test cap')
 return (-1)**clock,clock

def W(n):
 if n==1:return F(1)
 t=K=0
 while n!=1:
  n=3*n+1;a=0
  while not n%2:n//=2;a+=1
  K+=(a-1)//2;t+=1
  check(t<10000,'odd test cap')
 return F(2*factorial(K)*factorial(t),factorial(K+t+1))

def folded_by_literal_doubling(N,m):
 coefficient=0;x=m;sgn=1
 while x<=4*N:
  if x>N:coefficient+=(3*(x%3==2)-1)*sgn
  x*=2;sgn=-sgn
 return coefficient

def spatial():
 signs={n:shortcut_sign(n)[0] for n in range(1,1025)}
 weights={n:W(n) for n in range(1,1025,2)}
 rows=[];attainers=0
 for N in range(1,257):
  direct=sum((3*(n%3==2)-1)*signs[n] for n in range(N+1,4*N+1))
  K={m:folded_by_literal_doubling(N,m) for m in range(1,4*N+1,2)}
  check(K[1]==-3,'root kernel')
  check(sum(v*signs[m] for m,v in K.items())==direct,'doubling fold')
  norm=sum((F(v*v)/weights[m] for m,v in K.items() if m>1),F())
  value=sum((F(v)/weights[m]*v for m,v in K.items() if m>1),F())
  check(norm>0 and value*value/norm==norm,'row norm attainer J=K')
  rows.append(direct);attainers+=1
 windows=[]
 for X in [4,8,16,32,64,128]:
  col=[sum(F(folded_by_literal_doubling(N,m)**2)/weights[m] for N in range(X,2*X)) for m in range(3,8*X,2)]
  windows.append({'X':X,'actual_centered_spatial_energy':sum((rows[N-1]+3)**2 for N in range(X,2*X)),'source_operator_squared_norm_lower':str(max(col)),'hilbert_schmidt_squared':str(sum(col)),'Cauchy_upper':str(F(13,3)*sum(col))})
 return {'literal_shortcut_sources':1024,'folded_rows':256,'row_attainers':attainers,'windows':windows}

def fas_first_vertex(T):
 n=len(T);dp=[0]*(1<<n)
 for mask in range(1,1<<n):
  dp[mask]=min(dp[mask^(1<<v)]+sum(T[u][v] for u in range(n) if u!=v and mask>>u&1) for v in range(n) if mask>>v&1)
 return dp[-1]

def scaling():
 rows=[]
 for n in range(5,15):
  T=[[False if i==j else (i<j)!=(abs(i-j)==2) for j in range(n)] for i in range(n)]
  e=2*((n-1)//2)
  S=[[False if i==j else (i<j)!=(set([i,j])==set([1,e])) for j in range(n)] for i in range(n)]
  P=list(range(0,n,2))+list(range(1,n,2));pe={tuple(sorted(t)) for t in zip(P,P[1:])}
  check(all(S[a][b] for a,b in zip(P,P[1:])),'scaling path')
  for i,j in combinations(range(n),2):check((S[i][j] if (i,j) in pe else not S[i][j])==T[j][i],'scaling exact complement')
  f=fas_first_vertex(T);f0=fas_first_vertex(S)
  check(f==n-2-n//4 and f0==1,'independent first-vertex DP')
  rows.append({'n':n,'source_fas':f0,'target_fas':f,'difference':f-f0})
 return rows

x=spatial();result={'status':'PASS','arithmetic':'exact integers/Fractions','spatial':x,'scaling':scaling()}
print(json.dumps(result,indent=2,sort_keys=True))
