"""Exact finite controls for all-order return-collapse proof; stdlib only."""
from itertools import product
import hashlib,json
from pathlib import Path
T=(1,0,2,1);b=(3,1,3,1);z=(1,0,2,0)
def need(x,m):
 if not x:raise RuntimeError(m)
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def state(w,s):
 for c in w:s=step(s,c)
 return s
def scale(a,s):return tuple(a*x for x in s)
def balanced(s):return s[0]==s[2]and s[1]==s[3]
def onray(s):return s==scale(s[1],b)
def K(v,w):
 D,E,U,d=v;A,B,C,e=w
 return D*e+d*A+E*B+U*C-E*e-d*B-U*e-d*C+d*e

def count_graph(word):
 rows=[2,0];root=0
 for c in word:
  n=len(rows);rows.append(0)
  for v in range(n):
   if c!='O'and v==root:rows[n]|=1<<v
   else:rows[v]|=1<<n
  if c=='J':root=n
 n=len(rows);full=(1<<n)-1
 def run(starts,cycle):
  dp=[{}for _ in range(1<<n)]
  for x in starts:dp[1<<x][x]=1
  for mask in range(1,full+1):
   for v,c in dp[mask].items():
    avail=rows[v]&~mask&full
    while avail:
     bit=avail&-avail;avail-=bit;x=bit.bit_length()-1
     row=dp[mask|bit];row[x]=row.get(x,0)+c
  return sum(c for v,c in dp[full].items()if not cycle or rows[v]&1)
 return n,run(range(n),False),run([0],True)

def main():
 basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
 plane=[(1,0,1,0),(0,1,0,1)]
 returns=[];words=0
 for l in range(1,11):
  for x in product('OZJ',repeat=l):
   w=''.join(x);v=state(w,b);words+=1
   if not onray(v):continue
   lam=v[1];need(lam>=3 and lam%2==1,'physical return multiplier')
   need(all(balanced(state(w,e))for e in basis),'whole-space image in L')
   need(all(onray(state(w,e))for e in plane),'plane image in b ray')
   returns.append((w,lam))
 need(returns,'nonvacuous return universe')
 A='OJOJ';B='ZZOJ'
 need(state(A,b)==state(B,b)==scale(9,b),'equal physical returns')
 need(state(A,z)==scale(5,b)and state(B,z)==(23,7,23,7),'different first sidecars')
 counts={}
 for r in range(1,11):
  states={state(''.join(x),z)for x in product([A,B],repeat=r)}
  if r>=2:need(states=={scale(mu*9**(r-2),b)for mu in [45,69,71]},'exact three-state family')
  counts[r]=len(states)
 # Test arbitrary return compositions from separately declared short menu.
 short=[x for x in returns if len(x[0])<=5]
 for (a,la),(c,lc) in product(short,repeat=2):
  qq=state(a+c,z);need(onray(qq),'two arbitrary returns collapse')
  mu=qq[1];need(la*lc<2*mu<2*la*lc,'joint ratio bounds')
  for d,ld in short:
   need(state(a+c+d,z)==scale(mu*ld,b),'later choice changes only scale')
 C='ZOJ';w=state(C,b);q=state(C,z)
 need(w==scale(5,b)and q==(11,3,11,3),'single return boundary')
 need((K(T,w),K(T,q))==(25,19),'single return joint values')
 need(all(139%d for d in range(2,12)),'139 prime')
 literal=count_graph('JZOOZZJ');need(literal==(9,139,1),'literal prime target')
 return dict(status='PASS',word_length_cap=10,complete_words=words,return_words=len(returns),short_return_menu=len(short),
  pair_tests=len(short)**2,triple_tests=len(short)**3,AB_state_counts=counts,single_return_prime=dict(word='JZOOZZJ',vertices=literal[0],H=literal[1],C=literal[2]),
  verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Exact finite controls for reviewed all-order obstruction candidate; no universal target coverage.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
