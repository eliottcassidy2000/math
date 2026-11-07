from fractions import Fraction as F
from math import comb
import json,pathlib
D=pathlib.Path(__file__).resolve().parent

def pairf(z):
 s=[i for i in range(z.bit_length()) if z>>i&1]
 return sum((1<<(j-i-1))+int(j==i+1) for i in s for j in s if i<j)
def v(z):
 b=bin(z)[2:];return F(int(b[::-1],2),1<<len(b)) if z else F(0)
def runs(z,n):
 s=f'{z:0{n}b}' if n else '';return sum(c=='1' and (i==0 or s[i-1]=='0') for i,c in enumerate(s))
words=0
for n in range(13):
 cnt=[0]*(n+1)
 for z in range(1<<n):
  r=runs(z,n);cnt[r]+=1;words+=1
  delta=z*v(z)-pairf(z)
  assert 0<=delta<=F(3*r,4)
  if z:assert delta>0
  st=bin(z)[2:];i=0;prefix=0;runexpr=F(0)
  while i<len(st):
   if st[i]=='0':prefix*=2;i+=1;continue
   j=i
   while j<len(st) and st[j]=='1':j+=1
   u=v(prefix);runexpr+=u+(1-u)/(1<<(j-i))
   prefix=(prefix<<(j-i))+((1<<(j-i))-1);i=j
  assert delta==runexpr
 for r,x in enumerate(cnt):assert x==(comb(n+1,2*r) if 2*r<=n+1 else 0)
for m in range(1,101):
 z=((1<<(2*m))-1)//3
 assert z*v(z)-pairf(z)==F(2*m,3)-F(2,9)+F(2,9*(4**m))
states=[list(map(int,l.split())) for l in (D/'CORE_STATES.txt').read_text().splitlines()]
products=0
for h,q,A,G,T,S,E in states:
 for a in range(5):
  for b in range(5):
   for p in range(1<<a):
    for r in range(1<<b):
     tp=int(f'{p:0{a}b}'[::-1],2) if a else 0
     sr=int(f'{r:0{b}b}'[::-1],2) if b else 0
     alpha=(1<<b)*(q*p+A)+h*r;gamma=(1<<a)*(q*sr+G)+h*tp;Q=q*(1<<(a+b))
     tau=p*G+h*pairf(p)+T+(p&1)*S+(q*p+A)*sr+h*pairf(r)+E*(r>>(b-1) if b else 0)
     B=F(A*G,q)+F(h*A,q)*v(p)+F(h*G*r,q*(1<<b))+F(h*h*r,q*(1<<b))*v(p)-T-(p&1)*S-E*(r>>(b-1) if b else 0)
     assert F(alpha*gamma,Q)-tau==h*(p*v(p)-pairf(p)+r*v(r)-pairf(r))+B
     products+=1
print(json.dumps({'status':'PASS','binary_words':words,'core_product_cases':products,'alternating_runs_checked':100,'max_word_length':12,'run_count_formula':'exact histogram agreement'},indent=2))
