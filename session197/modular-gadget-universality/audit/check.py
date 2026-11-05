from pathlib import Path
from math import gcd
from functools import lru_cache
import json,hashlib

def need(x,s):
 if not x:raise RuntimeError(s)
def count_word(word):
 # Direct zero-position contribution, evaluated by the prefix's integer value.
 value=1
 for i,bit in enumerate(word):
  if bit=='0':
   n=int(word[:i]or'0',2);value+=n+(n&1)
 return value
def transition(word,M,n=0,f=1):
 for bit in word:
  if bit=='0':f=(f+n+(n&1))%M;n=(2*n)%M
  else:n=(2*n+1)%M
 return n,f

def loop_data(M):
 q=M;k=0
 while q%2==0:q//=2;k+=1
 a=2
 while pow(2,a,q)!=1%q:a+=1
 delta=((1<<a)-1)*((1<<k)-1)+1
 return k,a,delta

def general_word(M,u,z):
 k,a,d=loop_data(M);r=(u-1)%M;width=(M-1).bit_length();suffix=format(r,'0%db'%width)
 change=count_word(suffix)
 # Independent exhaustive residue search rather than a modular inverse.
 reps=next(j for j in range(M//2)if(change+j*d)%M==z)+M//2
 return ('1'*a+'0'*k)*reps+suffix

def dyadic(k,u,z,neutrals=1):
 M=1<<k;r=(u-1)%M;suffix=format(r,'0%db'%k);j=((count_word(suffix)-z)//2)%(M//2)
 terms=[]
 for i in range(k-1):
  if j>>i&1:terms.extend([2]if i==0 else[i+1,2])
 return ('1'+'0'*k)*neutrals+''.join('1'*a+'0'*k for a in terms)+suffix

def graph(L,seed=False):
 m=(L-1).bit_length();a=[0]*(m+1)
 for i in range(m):
  for j in range(i+1,m):a[i]|=1<<j
  if(L-1)>>i&1:a[i]|=1<<m
  else:a[m]|=1<<i
 if seed:
  # Add a universal source, then reverse only the old root-to-source edge.
  a=[(row<<1)for row in a];a[-1]|=1;a=[((1<<(m+1))-1)<<1]+a;a[0]^=1<<(m+1)
 return a

def literal(a):
 n=len(a);full=(1<<n)-1;h=cycles=0
 def paths(S,v,start):
  nonlocal h,cycles
  if S==full:
   h+=1;cycles+=start==0 and bool(a[v]&1);return
  out=a[v]&~S
  while out:
   b=out&-out;out-=b;paths(S|b,b.bit_length()-1,start)
 for v in range(n):paths(1<<v,v,v)
 return h,cycles

general=loop_states=suffix_states=0;digest=hashlib.sha256()
for M in range(2,66,2):
 k,a,d=loop_data(M);W='1'*a+'0'*k
 need(gcd(d,M)==2,'exact loop gcd')
 for start_n in range(0,5*M,M):
  for start_f in range(1,M,2):
   n,v=transition(W,M,start_n,start_f)
   need(n==0 and(v-start_f)%M==d%M,'arbitrary-state loop');loop_states+=1
 for u in range(M):
  r=(u-1)%M;word=format(r,'0%db'%((M-1).bit_length()))
  for start_f in (1,3,17):
   n,v=transition(word,M,0,start_f%M)
   need(n==r and(v-start_f)%M==(count_word(word)-1)%M,'suffix translation');suffix_states+=1
  for z in range(1,M,2):
   word=general_word(M,u,z);n,f=transition(word,M)
   need(word[0]=='1'and(n+1)%M==u and f==z,'general joint residues')
   need(len(word)<=(M-1)*(M+k)+(M-1).bit_length(),'effective length bound')
   need(count_word(word)%M==z,'direct zero-position sum')
   digest.update(f'{M},{u},{z},{word}\n'.encode());general+=1
short=padding=prefix=0
for k in range(1,8):
 M=1<<k
 for u in range(M):
  for z in range(1,M,2):
   word=dyadic(k,u,z);n,v=transition(word,M)
   need((n+1)%M==u and v==z,'dyadic correction')
   need(len(word)<=4*k*k+2*k+1,'quadratic bound')
   bigger=dyadic(k,u,z,3);nn,vv=transition(bigger,M)
   need((nn+1)%M==u and vv==z and int(bigger,2)>int(word,2),'arbitrarily large neutral loops');padding+=1;short+=1
for k in (16,32,64,128):
 M=1<<k
 for u,z in ((0,M-1),(1,1),(M-1,7),(M//2,15)):
  word=dyadic(k,u,z);need((int(word,2)+1)%M==u and count_word(word)%M==z,'large bit-length construction');need(len(word)<=4*k*k+2*k+1,'large word bound');short+=1
for k in range(1,8):
 M=1<<k
 for h,e,R in ((1,1,1),(3,1,3),(5,1,7),(13,2,18)):
  for u in range(M):
   for N in range(1,M,2):
    z=(pow(h,-1,M)*(N-R*u+e*(u%2)))%M
    need(z%2==1,'required gadget residue odd')
    word=dyadic(k,u,z);n,f=transition(word,M)
    need((h*f+R*((n+1)%M)-e*((n+1)%2))%M==N,'fixed-prefix residue');prefix+=1
Ls=set(range(2,129))
for k in (1,2):
 for u in range(1<<k):
  for z in range(1,1<<k,2):Ls.add(int(dyadic(k,u,z),2)+1)
physical=0
for L in sorted(Ls):
 G=graph(L);T=graph(L,True);g,c=literal(G);h,C=literal(T)
 need(g==count_word(bin(L-1)[2:]),'literal gadget path count')
 need(h==g+L-L%2 and C==1,'literal strong singleton-prefix count and cycle')
 need(all(T[i]>>j&1 for i in range(len(T)-1)for j in range(i+1,len(T)-1)),'transitive root deletion');physical+=1
result=dict(status='PASS',general_joint_pairs=general,general_word_digest=digest.hexdigest(),loop_state_controls=loop_states,suffix_translation_controls=suffix_states,dyadic_joint_pairs=short,neutral_padding_controls=padding,fixed_prefix_pairs=prefix,literal_graph_pairs=physical,scope='Constructive congruence coverage only; no exact target or primality implication.')
print(json.dumps(result,sort_keys=True,indent=2))
