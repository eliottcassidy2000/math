"""Standalone exact checks for field-invariant physical word encodings."""
from itertools import product
from math import isqrt
from pathlib import Path
import json

T=(1,0,2,1)
B=(3,1,3,1)
MAT={
 'I':((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),
 'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0)),
}
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))

def need(ok,msg):
 if not ok: raise RuntimeError(msg)

def mv(A,s): return tuple(sum(a*b for a,b in zip(r,s)) for r in A)
def step(s,c):
 D,E,U,d=s
 if c=='I': return s
 if c=='O': return D,0,2*U,d
 if c=='Z': return D-E+U,d,2*U-d,d
 if c=='J': return D-E+U,D-E,D+U,D
 raise ValueError(c)
def state(w,seed=T):
 for c in w: seed=step(seed,c)
 return seed

def K(v,w): return sum(v[i]*Q[i][j]*w[j] for i in range(4) for j in range(4))
def gate(kind,v,w):
 D,E,U,d=v;g,s,A,b=w
 if kind=='K': return K(v,w)
 if kind=='N': return (D-E+U)*g-(D-E+d)*s+(D+3*U-d)*A-(U+E)*b
 if kind=='M': return (D+2*U)*g-(U+d)*s+(2*D-E+6*U-d)*A-(U+E)*b
 raise ValueError(kind)

def cap_length(B):
 need(B>=2,'B at least two')
 power=5;L=-1
 while power<=B**3: L+=1;power*=5
 return L

def word_bank(L):
 bank={T:''}
 frontier={T:''}
 for _ in range(L):
  nxt={}
  for s,w in frontier.items():
   for c in 'OZJ':
    z=step(s,c)
    if z not in nxt: nxt[z]=w+c
    if z not in bank: bank[z]=w+c
  frontier=nxt
 return bank

def height_bank(B):
 if B<2:return {}
 bank={T:''};queue=[T]
 for s in queue:
  for c in 'OZJ':
   z=step(s,c)
   if z[2]<=B and z not in bank:bank[z]=bank[s]+c;queue.append(z)
 return bank

def selectors_ok(e,p=None,orthogonal=True):
 vals=[x*x-x for x in e]+[sum(e)-1]
 if orthogonal: vals +=[e[i]*e[j] for i in range(len(e)) for j in range(i+1,len(e))]
 return all(v==0 if p is None else v%p==0 for v in vals)

def equations(word,L,p=None):
 need(len(word)<=L,'padding bound')
 s=T;res=[]
 for c in word+'I'*(L-len(word)):
  e=tuple(int(c==a) for a in 'IOZJ')
  res += [x*x-x for x in e]
  res += [e[i]*e[j] for i in range(4) for j in range(i+1,4)]
  res += [sum(e)-1]
  z=step(s,c)
  res +=[z[j]-sum(e[i]*mv(MAT[a],s)[j] for i,a in enumerate('IOZJ')) for j in range(4)]
  s=z
 need(all(x==0 if p is None else x%p==0 for x in res),'polynomial recurrence')
 need(len(res)==15*L,'equation count excluding initial state')
 return s

def source_graph(word):
 rows=[2,0];root=0
 for c in word:
  if c=='I':continue
  n=len(rows);rows.append(0)
  for u in range(n):
   if c!='O' and u==root:rows[n]|=1<<u
   else:rows[u]|=1<<n
  if c=='J':root=n
 return rows,root

def counts(rows,root=0,profile=False):
 n=len(rows);size=1<<n
 dp=[[0]*n for _ in range(size)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,size):
  for v in range(n):
   if not(mask>>v&1):continue
   total=dp[mask][v]
   if not total:continue
   todo=rows[v]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;u=bit.bit_length()-1
    dp[mask|bit][u]+=total
 h=[sum(r) for r in dp];h[0]=1
 sp=[[0]*n for _ in range(size)];sp[1<<root][root]=1
 for mask in range(1,size):
  if not(mask>>root&1):continue
  for v in range(n):
   total=sp[mask][v]
   if not total:continue
   todo=rows[v]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;u=bit.bit_length()-1
    sp[mask|bit][u]+=total
 cyc=sum(sp[-1][v] for v in range(n) if rows[v]>>root&1)
 if not profile:return h[-1],cyc
 U=sum(sum(sp[mask])*h[(size-1)^mask] for mask in range(1,size) if mask>>root&1)
 return (h[-1],dp[-1][root],U,h[(size-1)^(1<<root)]),sum(sp[-1])

def c1_graph(alpha,beta):
 sigma={'O':'Z','Z':'O','J':'J','I':'I'}
 w=beta+''.join(sigma[c] for c in reversed(alpha))+'J'
 return source_graph(w)[0]

def c2_graph(kind,alpha,beta):
 V,c=source_graph(alpha);W,a=source_graph(beta)
 rows=V[:]
 def birth():
  nonlocal rows,c
  n=len(rows);rows.append(1<<c)
  for u in range(n):
   if u!=c:rows[u]|=1<<n
  return n
 old=c;x=birth()
 if kind=='M':c=x;x=birth()
 entries={old,x}
 n=len(rows);m=len(W)
 G=[((1<<m)-1)^(1<<v)^W[v] for v in range(m)]
 rows +=[r<<n for r in G]
 for u in range(n):
  for v in range(m):
   if v==a and u in entries:rows[n+v]|=1<<u
   else:rows[u]|=1<<(n+v)
 return rows

def isprime(p):
 if p<2:return False
 return all(p%d for d in range(2,isqrt(p)+1))
def nextprime(n):
 p=n+1
 while not isprime(p):p+=1
 return p

def main():
 # Complete small selector fibers, including the characteristic-2/3 hostile.
 selector={}
 for p in [2,3,5,7]:
  full=[e for e in product(range(p),repeat=4) if selectors_ok(e,p)]
  naive=[e for e in product(range(p),repeat=4) if selectors_ok(e,p,False)]
  need(len(full)==4,'exact one-hot over finite fields')
  selector[p]={'full':len(full),'naive':len(naive)}
 need(selector[2]['naive']>4 and selector[3]['naive']>4,'small-characteristic hostile')
 # Exact all-word transitions and growth, deliberately separate implementations.
 nwords=0
 for l in range(8):
  for letters in product('OZJ',repeat=l):
   w=''.join(letters);s=T
   for c in w:
    need(step(s,c)==mv(MAT[c],s),'matrix versus direct recurrence');s=step(s,c)
   need(5**(l+1)<=s[2]**3,'source growth witness')
   need(s[2]<=2**(l+1),'source upper height')
   need(equations(w,l+1)==s,'padding identity')
   nwords+=1
 # Independent literal source and full graft counters.
 bank=word_bank(2);physical=0
 for s,w in word_bank(3).items():
  rows,r=source_graph(w);actual,starts=counts(rows,r,True)
  need(actual==s and starts==1,'literal source profile');physical+=1
 graph_pairs=0
 for v,alpha in bank.items():
  for w,beta in bank.items():
   need(counts(c1_graph(alpha,beta))==(K(v,w),1),'literal C1 compiler')
   for kind in 'NM':
    need(counts(c2_graph(kind,alpha,beta))==(gate(kind,v,w),2),'literal C2 compiler')
   graph_pairs+=1
 # Polynomial output ranges and exact target equivalence with inherited cutoffs.
 records=[];supports={}
 for kind in ['K','C2']:
  for n in [7,13,15,21,25,41,43,67,101,127]:
   if kind=='C2' and n<13:
    records.append({'kind':kind,'target':n,'accepted':False,'reason':'proved minimum13'});continue
   Bcap=isqrt(4*(n-1)) if kind=='K' else (n-1)//6
   L=cap_length(Bcap);height=height_bank(Bcap);words=word_bank(L)
   need(set(height)<=set(words),'word cap contains exact height bank')
   gates=['K'] if kind=='K' else ['N','M']
   support={15} if kind=='C2' else set()
   small={15} if kind=='C2' else set()
   for G in gates:
    for v in words:
     for w in words:support.add(gate(G,v,w))
    for v in height:
     for w in height:small.add(gate(G,v,w))
   accepted=n in support
   need(accepted==(n in small),'full length-bank versus complete bounded inverse')
   M=(1 if kind=='K' else 11)*4**(L+1)
   if kind=='C2':M=max(M,15)
   need(max(support)<=M,'all exact scalar range bounds')
   p=nextprime(max(n,M));modaccept=n%p in {h%p for h in support}
   need(modaccept==accepted,'large-prime exactness')
   if n in [7,21,67]:need(not accepted,'known C1/C2 exclusion')
   records.append({'kind':kind,'target':n,'B':Bcap,'L':L,'profiles':len(words),'accepted':accepted,'prime':p,'variables':16*L+8,'C1_equations':30*L+9})
   supports[kind,n]=support
 # Explicit positive equation and graph witness, chosen by exhaustive exact level.
 n=101;L=cap_length(isqrt(4*(n-1)));bank=word_bank(L);witness=None
 for v,a in bank.items():
  for w,b in bank.items():
   if K(v,w)==n:witness=(a,b);break
  if witness:break
 need(witness is not None,'positive target101')
 alpha,beta=witness
 v=equations(alpha,L);w=equations(beta,L)
 need(K(v,w)==n,'positive terminal equation')
 graph=c1_graph(alpha,beta)
 need(counts(graph)==(n,1),'positive decoded graph')
 result={'status':'PASS','source_words_tested':nwords,'literal_sources':physical,'literal_pairs_three_gates':graph_pairs,'finite_field_selectors':selector,'target_controls':records,'positive_C1':{'target':n,'alpha':alpha,'beta':beta,'left':v,'right':w,'rows':graph,'literal_counts':counts(graph)},'scope':'Finite exact controls; all-order completeness uses frozen physical source theorems.'}
 print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
