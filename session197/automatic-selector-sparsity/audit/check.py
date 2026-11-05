"""Independent exact binary formulas, canonical DFA fixtures and finite bands."""
from itertools import product
from math import gcd,comb,isqrt
from collections import deque
import json,random,hashlib

def need(b,s):
 if not b:raise RuntimeError(s)
def fword(word):
 # Unfolded zero-position sum; independent of the author's halving recurrence.
 return 1+sum(int(word[:i]or'0',2)+(int(word[:i]or'0',2)&1)for i,b in enumerate(word)if b=='0')
def f(n):return fword(bin(n)[2:])
def valuation(n):
 v=0
 while n%2==0:n//=2;v+=1
 return v

def odd(n):return n>>valuation(n)
def pair_encoding(a,b):
 x,y=bin(a)[2:],bin(b)[2:];k=max(len(x),len(y));return x.zfill(k),y.zfill(k)

formula_count=padding_count=empty_suffix=0;digest=hashlib.sha256()
for P in range(1,18):
 p=bin(P)[2:]
 for t in range(6):
  for S in range(1<<t):
   s=format(S,f'0{t}b')if t else''
   q=int(''.join('1'if c=='0'else'0'for c in s)[::-1]or'0',2)
   for h,e,R in((1,1,1),(3,1,3),(5,1,7),(13,2,18)):
    K=h*q+R*(1<<t);A=(h+K)*(P+1);B=h*(P+1);C=h+K;D=h*(f(P)+f(S)+1)+R*(S+1)-e*((S+1)%2)
    for a,b in product(range(1,9),repeat=2):
     word=p+'1'*a+'0'*b+s;n=int(word,2);J=h*fword(word)+R*(n+1)-e*((n+1)%2);Y=A<<(a+b);Z=(B<<a)+(C<<b)
     need(J==Y-Z+D,'direct digit objective');need(len(bin(Y)[2:])-len(word)==A.bit_length()-len(p)-len(s),'length-controlled auxiliary')
     need(B>0 and C>0 and Z>0,'positive band coefficients')
     # Add an auxiliary track longer than either retained track. Project,
     # trim common zero columns, and check exact canonical pair recovery.
     k=max(n.bit_length(),Z.bit_length(),Y.bit_length(),J.bit_length())+2
     x,y=bin(n)[2:].zfill(k),bin(Z)[2:].zfill(k)
     while len(x)>1 and x[0]==y[0]=='0':x,y=x[1:],y[1:]
     need((x,y)==pair_encoding(n,Z),'projection/common-zero normalization')
     padding_count+=1;formula_count+=1;empty_suffix+=not t
     digest.update(f'{n},{J},{Z}\n'.encode())

bands=distinguish=0
pairs=[(b,c)for b in range(1,25)for c in range(1,25)]+[(3<<20,5),(5,3<<20),(7<<40,9<<37)]
for B,C in pairs:
 for delta0 in(-11,-3,0,2,13):
  H=odd((B<<max(delta0,0))+(C<<max(-delta0,0)));ell=H.bit_length()-1
  lo=-max(0,valuation(B)-valuation(C)+ell);hi=max(0,valuation(C)-valuation(B)+ell)
  allowed=[]
  for delta in range(min(lo,-80)-3,max(hi,80)+4):
   good=odd((B<<max(delta,0))+(C<<max(-delta,0)))<=H
   if good:need(lo<=delta<=hi,'finite odd-part support');allowed.append(delta)
  need(delta0 in allowed,'nonempty selected diagonal')
  for T in(1,2,5):
   M=1+max(abs(d-delta0)for d in allowed)//T
   for i,j in product(range(6),repeat=2):
    need((delta0+T*M*(i-j)in allowed)==(i==j),'distinguishing suffixes');distinguish+=1
  bands+=1

def run(delta,state,word):
 for b in word:state=delta[state][int(b)]
 return state

def shortest(delta,start,targets,within=None):
 q=deque([(start,'')]);seen={start}
 while q:
  x,w=q.popleft()
  if x in targets:return w
  for b,y in enumerate(delta[x]):
   if y not in seen and(within is None or y in within):seen.add(y);q.append((y,w+str(b)))
 raise RuntimeError('path missing')

def cycle_phase(delta,state,bit):
 seen={};word='';x=state
 while x not in seen:seen[x]=len(word);word+=str(bit);x=delta[x][bit]
 mu=seen[x];period=len(word)-mu;exponent=mu or period
 return exponent,period,run(delta,state,str(bit)*exponent)

fixtures=empty=closed=leaky=rectpoints=escapechecks=0
for n in range(1,4):
 for flat in product(range(n),repeat=2*n):
  delta=[flat[2*i:2*i+2]for i in range(n)]
  reach=[[i==j or j in delta[i]for j in range(n)]for i in range(n)]
  for k in range(n):
   for i in range(n):
    for j in range(n):reach[i][j]|=reach[i][k]and reach[k][j]
  for accmask in range(1<<n):
   fixtures+=1;accept={i for i in range(n)if accmask>>i&1};live={i for i in range(n)if reach[0][i]and any(reach[i][a]for a in accept)}
   if not live:empty+=1;continue
   parts={tuple(j for j in live if reach[i][j]and reach[j][i])for i in live}
   complete=[set(S)for S in parts if all(y in S for x in S for y in delta[x])]
   if complete:
    closed+=1;S=complete[0];r=min(S);prefix='1'+shortest(delta,0,{r});a,t1,r1=cycle_phase(delta,r,1);b,t0,r2=cycle_phase(delta,r1,0);T=t1*t0//gcd(t1,t0);suffix=shortest(delta,r2,accept,S)
    # The explicit canonical DFA starts with one new leading1 edge, then
    # uses delta; any leading0 goes to a new dead state.
    for i,j in product(range(4),repeat=2):
     word=prefix+'1'*(a+T*i)+'0'*(b+T*j)+suffix
     need(word[0]=='1'and run(delta,0,word[1:])in accept,'canonical complete-component rectangle');rectpoints+=1
   else:
    leaky+=1;s=len(live)
    for v in live:
     escaped=shortest(delta,v,set(range(n))-live)
     need(1<=len(escaped)<=s,'shortest live escape bound')
     remaining=sum(run(delta,v,''.join(bits))in live for bits in product('01',repeat=s))
     need(remaining<1<<s,'uniform escaping block');escapechecks+=1
    counts=[int(i==0)for i in range(n)]
    for length in range(20):
     q,r=divmod(length,s);need(sum(counts[i]for i in accept)<=((1<<s)-1)**q*(1<<r),'full word-count leakage bound')
     nxt=[0]*n
     for i,v in enumerate(counts):
      for j in delta[i]:nxt[j]+=v
     counts=nxt

prime_checks=0
for n in range(1,301):
 primes=[p for p in range(2,2*n+1)if all(p%d for d in range(2,isqrt(p)+1))]
 central=comb(2*n,n);need(central*(2*n+1)>=4**n,'central lower bound');need(central<=(2*n)**len(primes),'prime-count upper product bound');prime_checks+=1
print(json.dumps(dict(status='PASS',direct_rectangle_values=formula_count,empty_suffix_values=empty_suffix,forced_padding_controls=padding_count,formula_digest=digest.hexdigest(),odd_part_bands=bands,distinguishability_controls=distinguish,complete_DFA_fixtures=fixtures,empty_fixtures=empty,complete_live_fixtures=closed,leaky_fixtures=leaky,canonical_rectangle_points=rectpoints,live_escape_controls=escapechecks,prime_lower_bound_controls=prime_checks),sort_keys=True,indent=2))
