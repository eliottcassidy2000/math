#!/usr/bin/env python3
"""Focused literal-DP replay; explicit errors survive python -O."""
import json,sys
from pathlib import Path
from itertools import permutations
from math import gcd
from central import endpoints,circulant,apex,data,cut
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'session275_exact_selection_20261111/switches'))
from carrier import profile,matrix,mul,from_matrix
COUNT=0

def check(a,b,label):
 global COUNT
 COUNT+=1
 if a!=b:raise RuntimeError((label,a,b))
def f(a):
 if not a:return 1
 q,r=divmod(a,2)
 return f(q)+(q+(q&1) if not r else 0)
def L(a):
 bits=[int(x) for x in bin(a)[2:]] if a else []
 n=len(bits);T=[sum(x<<(i+1) for i,x in enumerate(bits))]
 for i,x in enumerate(bits):T.append((1-x)|sum(1<<(j+1) for j in range(i+1,n)))
 return tuple(T)
def conv(T):return tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T)))
def fusion(T,S):
 n=len(T);m=len(S);V=[0]*(n+m-1)
 for X,mp in [(T,list(range(n))),(S,[0]+list(range(n,n+m-1)))]:
  for i in range(len(X)):
   for j in range(len(X)):
    if X[i]>>j&1:V[mp[i]]|=1<<mp[j]
 for i in range(1,n):V[i]|=sum(1<<j for j in range(n,n+m-1))
 return tuple(V)
def H(T):return sum(endpoints(T)[-1])
def F(p,A,B):
 h=p['h'];return p['H']+h*(f(A)+f(B)-2)+(p['U']-h)*B+(p['V']-h)*A+p['Q']*A*B-(p['d']-h)*(B&1)-(p['e']-h)*(A&1)
def inverse(N,p,width):
 h=p['H'];a=p['U']-p['h'];b=p['V']-p['h'];Q=p['Q'];delta=N-h
 if delta<0:return []
 ans=[]
 for A in range(delta//b+1):
  bound=delta//a if not A else min(delta//a,delta//(Q*A))
  for B in range(bound+1):
   rem=N-F(p,A,B)
   if rem>=0 and rem%2==0 and rem//2<=width:ans.append((A,B,rem//2))
 return ans

def run():
 C=(56,13,25,48,34,6);Cp=(52,13,24,49,34,6);p=profile(C);pp=profile(Cp)
 check(p,dict(H=41,h=15,Q=120,k=0,U=78,V=57,d=6,e=9),'C41 profile')
 check(pp,{**p,'H':43},'C43 profile')
 for T in [C,Cp]:check(H(T),sum(all(T[a]>>b&1 for a,b in zip(w,w[1:])) for w in permutations(range(6))),'independent permutations')
 for a in range(256):
  T=L(a);check(H(T),f(a),'digital f');q=profile(T);check((q['U'],q['d']),(a+1,1-(a&1)),'digital interface')
 for A in range(8):
  for B in range(8):
   T=fusion(fusion(L(A),C),conv(L(B)));Tp=fusion(fusion(L(A),Cp),conv(L(B)))
   check(H(T),F(p,A,B),'sandwich absolute target');check(H(Tp),F(p,A,B)+2,'central translation')
   check((F(p,A,B)-41)%6,0,'C5 content')
 for k in range(2,50):
  check(f((1<<k)-2),(1<<(k-1))+1,'carry before');check(f((1<<k)-1),1,'carry after')
  for B in [0,1,17]:check(F(p,(1<<k)-1,B)-F(p,(1<<k)-2,B),p['Q']*B+p['V']-p['h']-p['h']*(1<<(k-1))-(p['e']-p['h']),'carry jump')
 core=circulant(9,1);D=data(core);menus={cut(D,m):m for m in range(512) if m.bit_count()==4}
 C9=apex(core,menus[15381]);C9p=apex(core,menus[15383]);p9=profile(C9);p9p=profile(C9p)
 check(p9,dict(H=15381,h=3357,Q=50004,k=0,U=24089,V=29272,d=1865,e=1492),'C9 profile');check(p9p,{**p9,'H':15383},'C9 twin');check(gcd(3357,5183),1,'content removed')
 for A,B in [(0,0),(1,0),(0,1),(1,1),(2,1),(3,2)]:check(H(fusion(fusion(L(A),C9),conv(L(B)))),F(p9,A,B),'C9 sandwich')
 for t in range(3357):check((F(p9,0,2*t)-15381)%(2*3357),(2*20732*t)%(2*3357),'primitive local residues')
 check(len({F(p9,0,2*t)%(2*3357) for t in range(3357)}),3357,'all odd local residues')
 c7=circulant(7,0);d7=data(c7);m7={cut(d7,m):m for m in range(128) if m.bit_count()==2};p7=profile(apex(c7,m7[559]));check({(F(p7,A,B)+j)%7 for A in range(14) for B in range(14) for j in [0,2]},{0,1,2,5,6},'primitive but local obstruction')
 for q in [p,p7,p9]:
  h=q['h'];a=q['U']-h;b=q['V']-h;c=q['d']-h;e=q['e']-h
  formula=gcd(gcd(gcd(2*h,2*a),gcd(2*b,a-c)),gcd(b-e,q['Q']))
  actual=0
  for A in range(4):
   for B in range(4):actual=gcd(actual,F(q,A,B)-q['H'])
  check(formula,actual,'exact gcd generators')
 targets=[F(p,A,B)+j for A,B in [(0,0),(1,0),(0,1),(2,3),(5,2)] for j in [0,2]]
 for N in targets:
  ans=inverse(N,p,1)
  if not ans:raise RuntimeError(('decoder false negative',N))
  for A,B,j in ans:check(F(p,A,B)+2*j,N,'decoder')
 for N in [45,51,57,63,69]:check(inverse(N,p,1),[],'content rejection')
 result=dict(status='PASS',checks=COUNT,C41=C,C43=Cp,C9_15381=C9,C9_15383=C9p,integer_target_examples=[dict(N=N,solutions=inverse(N,p,1)) for N in targets],scope='literal sandwich tests; all-order assertions use proofs, not finite extrapolation')
 print(json.dumps(result,indent=2))
if __name__=='__main__':run()
