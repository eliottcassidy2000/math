"""Exact Session184 return-tree 7-lane extraction. Standard library only."""
from math import isqrt
import json
from pathlib import Path

def need(b,s):
 if not b: raise RuntimeError(s)
def sig(word):
 h,p,f,g,j,e,s=1,1,0,0,0,0,0
 for index,(typ,a) in enumerate(word):
  H,P,F,G,J,E,S=(1,2,a,a,0,a,a) if typ==1 else (3,12,3*a,3*a,a*(a-1),a,a)
  if index==0:
   h,p,f,g,j,e,s=H,P,F,G,J,E,S;continue
  h,p,f,g,j,e,s=h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H
 return h,p,f,g,j,e,s
STATES=[('',[(1,1)]*3),('2',[(1,1),(1,1),(3,1)]),('4',[(1,1),(1,0),(1,1),(1,0),(1,1)]),('22',[(1,1),(1,0),(1,1),(3,2)]),('42',[(1,1),(1,0),(3,1),(1,1)]),('422',[(1,1),(1,0),(3,1),(3,1)]),('4222',[(3,1),(1,1),(1,1),(1,1),(1,1),(1,0),(1,1)])]
def graph(word,t):
 rows=[sum(1<<j for j in range(i+1,t)) for i in range(t)]
 ports=[]
 for typ,a in word:
  n=len(rows); rows=[x|(((1<<typ)-1)<<n) for x in rows]
  rows += [0] if typ==1 else [2<<n,4<<n,1<<n]
  ports+=list(range(n,n+a))
 n=len(rows);rows.append((1<<t)-1)
 for i in range(t,n):
  if i in ports:rows[i]|=1<<n
  else:rows[n]|=1<<i
 return rows

def strong(rows):
 for i in range(len(rows)):
  s=1<<i
  while True:
   q=s
   for j,r in enumerate(rows):
    if s>>j&1:q|=r
   if q==s:break
   s=q
  if s!=(1<<len(rows))-1:return False
 return True

def hamilton(rows):
 n=len(rows);incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
 dp=[None]*(1<<n);dp[0]=[0]*n
 for S in range(1,1<<n):
  d=[0]*n;b=S
  while b:
   z=b&-b;v=z.bit_length()-1;b-=z;T=S^z
   if not T:d[v]=1;continue
   a=incoming[v]&T;x=dp[T];total=0
   while a:
    q=a&-a;a-=q;total+=x[q.bit_length()-1]
   d[v]=total
  dp[S]=d
 return sum(dp[-1])
def prime(p):
 if p<2:return False
 return all(p%d for d in range(2,isqrt(p)+1))
def val(n):
 k=0
 while n%7==0:k+=1;n//=7
 return k

def run():
 tree=[]
 for address,w in STATES:
  h,p,f,g,j,e,s=sig(w);c=h-j
  tree.append({'address':address,'F':f,'c':c,'F_mod7':f%7,'c_mod7':c%7})
 need([x['address'] for x in tree if x['c_mod7']==0]==['422'],'unique7family')
 need(all(x['F_mod7']==0 for x in tree),'slope7')
 w=STATES[5][1];need(sig(w)==(9,576,189,333,100,3,9),'literal signature')
 physical=[]
 for t in (1,2,3,4,5,8):
  rows=graph(w,t);H=hamilton(rows)
  need(strong(rows),'strong');need(H==189*2**t-91,'literal DP')
  physical.append({'t':t,'n':len(rows),'rows':rows,'H':H,'cofactor':H//7})
 roots=[];a=0;L=3
 for j in range(1,13):
  need((27*pow(2,a,7**j)-13)%7**j==0,'root')
  choices=[a+k*L for k in range(7) if (27*pow(2,a+k*L,7**(j+1))-13)%7**(j+1)==0]
  need(len(choices)==1,'unique lift')
  roots.append({'j':j,'modulus':L,'root':a,'next_root':choices[0]})
  a=choices[0];L*=7
 for t in range(1,1001):
  H=7*(27*2**t-13);v=val(H)
  need((v==1)==(t%3!=0),'exact valuation1')
  for r in roots[:4]:need((v>=r['j']+1)==((t-r['root'])%r['modulus']==0),'valuation root')
 witnesses=[]
 for t in range(1,36):
  q=27*2**t-13
  if prime(q):witnesses.append({'t':t,'n':t+9,'p':q,'H':7*q})
 need([x['t'] for x in witnesses]==[1,4,8,13,19,28,32,35],'prime list')
 out={'status':'PASS','tree':tree,'direct_graph_DP':physical,'seven_adic_roots':roots,'trial_division_prime_witnesses':witnesses,'scope':'Unique 7-divisible ray among inherited seven states, all-clock valuations; no infinite-prime or spectrum-completeness claim.'}
 Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':run()
