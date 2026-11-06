"""Independent fixed-depth inverse, physical compiler, and prime-chain audit.
No author-code imports. All guards remain enabled under python -O.
"""
from pathlib import Path
from math import gcd,isqrt
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json,hashlib
ROOT=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def ceil(a,b):return -((-a)//b)
T=(1,0,2,1)
M={'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
P=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1));UU=P[2]
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j]for k in range(4))for j in range(4))for i in range(4))
def act(c,v):return tuple(sum(x*y for x,y in zip(row,v))for row in M[c])
def state(w,s=T):
 for c in w:s=act(c,s)
 return s
def matrix(w):
 a=tuple(tuple(int(i==j)for j in range(4))for i in range(4))
 for c in w:a=mul(M[c],a)
 return mul(Q,a)
def binary(u):
 need(u>=2,'physical binary source domain')
 return ''.join('O'if b=='1'else'Z'for b in bin(u-1)[3:])
@lru_cache(None)
def g(u):
 need(u>=1,'scalar domain')
 return 1 if u==1 else state(binary(u))[0]
def tau(u):return(g(u),u%2,u,1)
def dagger(w):return w[::-1].translate(str.maketrans('OZ','ZO'))
def row(A,u):
 t=tau(u);r=tuple(sum(t[i]*A[i][j]for i in range(4))for j in range(4))
 return r[0],r[2],r[1],r[3]
def read(r,h):a,b,c,d=r;return a*g(h)+b*h+c*(h%2)+d
def valid(r):a,b,c,d=r;return a>0 and a%2 and b>=a and -a<=c<=0 and a-b<=d<=0

def inv(r,n,vmin=2,vmax=None,deep_checks=False):
 # Uniform-depth tree, rather than the author's adaptive leaf stopping.
 # The complete binary level K is searched; modular inverses are recomputed locally.
 need(valid(r),'initial numerical row');a,b,c,d=r;d0=gcd(a,b);m=a//d0
 if n<a:return [],dict(nodes=0,budget=1,singleton_leaves=0,candidates=0)
 vmax=min(1+(n-a)//b,vmax)if vmax is not None else 1+(n-a)//b
 if vmax<vmin:return [],dict(nodes=0,budget=1,singleton_leaves=0,candidates=0)
 K=0
 while 4**K*b*b<n*d0:K+=1
 budget=2**(K+1)-1;level=[(r,1,0)];nodes=singleton=candidates=0;found=[]
 for depth in range(K+1):
  nextlevel=[]
  for rr,s,off in level:
   nodes+=1;aa,bb,cc,dd=rr
   need(valid(rr) and aa==a and gcd(a,bb)==d0,'preserved row/gcd')
   need(0<=off<s and s==1<<depth,'index map range')
   need((bb-(1<<depth)*b)%a==0,'fixed-depth coefficient congruence')
   lo=max(1,ceil(n,a+bb),ceil(vmin+off,s));hi=min(1+(n-a)//bb,(vmax+off)//s)
   need(Fraction(1+(n-a)//bb-ceil(n,a+bb))<1+Fraction(n*a,bb*bb),'strict un-intersected width')
   if deep_checks:
    for h in(1,2,3,7):need(read(rr,h)==read(r,s*h-off),'exact original index map including singleton')
   if lo>hi:continue
   if depth<K:
    nextlevel.extend([((a,2*bb,0,dd),2*s,off),((a,a+2*bb,-a,dd+cc-bb),2*s,s+off)])
    continue
   need(bb*bb>=n*d0 and hi-lo<2*m,'uniform-depth leaf width')
   for parity in(0,1):
    rhs=n-dd-(bb+cc)*parity
    if rhs%d0:continue
    il,ih=ceil(lo-parity,2),(hi-parity)//2
    if il>ih:continue
    if m==1:residue=0
    else:
     local=pow(2*bb//d0,-1,m)
     prepared=pow(2*b//d0,-1,m)*pow(pow(2,-1,m),depth,m)%m
     need(local==prepared,'precomputed inverse at fixed depth')
     residue=local*(rhs//d0)%m
    t=il+(residue-il)%m
    need(t+m>ih,'at most one candidate per parity')
    if t>ih:continue
    h=2*t+parity;candidates+=1;singleton+=int(h==1)
    if read(rr,h)==n:
     v=s*h-off;need(vmin<=v<=vmax and read(r,v)==n,'restored exact original index')
     found.append(v)
  level=nextlevel
 need(nodes<=budget and candidates<=2**(K+1) and len(found)==len(set(found)),'full-tree budgets and uniqueness')
 return sorted(found),dict(nodes=nodes,budget=budget,singleton_leaves=singleton,candidates=candidates)

# Complete finite cone controls, including interior gamma/delta and arbitrary windows.
rows=set()
for a in(1,3,5,7,9,15,21):
 for b in(a,a+1,a+3,2*a,3*a):
  for c in(-a,-a//2,0):
   for d in(a-b,(a-b)//2,0):rows.add((a,b,c,d))
need(all(valid(r)for r in rows),'test rows valid');numerical_row_count=len(rows)
controls=singletons=nodes=0
for r in sorted(rows):
 expected={}
 for v in range(2,304):
  n=read(r,v)
  if n<=301:expected.setdefault(n,[]).append(v)
 for n in range(1,302):
  got,stats=inv(r,n,deep_checks=True);need(got==expected.get(n,[]),'complete numerical row control')
  singletons+=stats['singleton_leaves'];nodes+=stats['nodes'];controls+=1
 for n in range(1,302,13):
  for left,right in((2,2),(3,19),(5,5),(7,103),(40,11)):
   got,stats=inv(r,n,left,right,True)
   need(got==[v for v in expected.get(n,[])if left<=v<=right],'restricted original-index window')
   controls+=1
need(singletons>0,'singleton arithmetic leaf exercised')

# Rebuild full transpose-closed length-five menu, not the twelve chosen representatives.
raw=['','J']+['J'+''.join(z)+'J'for k in range(2,6)for z in product('OZJ',repeat=k-2)]
forms={}
for w in sorted(raw,key=lambda z:(len(z),z)):
 A=matrix(w)
 if A in forms:continue
 C=mul(mul(tuple(zip(*P)),A),P)
 bound=min(Fraction(4*C[i][j],UU[i]*UU[j])for i in range(4)for j in range(4))
 need(bound.denominator==1 and bound>=2,'complete product bound')
 forms[A]=(w,int(bound))
need(len(forms)==35 and all(tuple(zip(*A))in forms for A in forms),'35-matrix transpose closure')
for u in range(2,129):
 for A,(w,bound)in forms.items():
  r=row(A,u);s=state(dagger(w),tau(u));D,E,U,d=s
  need(r==(d,U-d,E-d,D-E-U+d),'physical adjoint row')
  need(valid(r) and r[1]>=u-1,'physical row cone and complexity lower bound')

def target(n):
 out=[];stats=dict(rows=0,nodes=0)
 for A,(w,bound)in forms.items():
  for u in range(2,isqrt(4*n//bound)+1):
   vs,ss=inv(row(A,u),n,u,4*n//(bound*u));stats['rows']+=1;stats['nodes']+=ss['nodes']
   out.extend((w,u,v)for v in vs)
 return sorted(out,key=lambda x:(len(x[0]),x[0],x[1],x[2])),stats
CAP=501;direct={}
for A,(w,bound)in forms.items():
 for u in range(2,isqrt(4*CAP//bound)+1):
  r=row(A,u)
  for v in range(u,4*CAP//(bound*u)+1):
   n=read(r,v)
   if n<=CAP:direct.setdefault(n,[]).append((w,u,v))
for n in range(1,CAP+1):
 hit,stats=target(n);need(hit==sorted(direct.get(n,[]),key=lambda x:(len(x[0]),x[0],x[1],x[2])),'complete target control')
critical,critical_stats=target(2134117);need(critical==[],'exhaustive finite-menu hostile')

# Independently verify the recursive Lucas certificates using complete n-1 factors.
cert=json.loads((ROOT/'inputs/PRIME-CERTIFICATES.json').read_text());verified={};modchecks=0
for key in sorted(cert,key=int):
 p=int(key);c=cert[key];need(c['prime']==p,'certificate label')
 if p==2:
  need(c['factorization']=={}and c['base']is None,'Lucas base case');verified[p]=True;continue
 powers={int(q):e for q,e in c['factorization'].items()};prod=1
 for q,e in powers.items():
  need(q in verified and isinstance(e,int)and e>0,'prior certified prime')
  prod*=q**e
 need(prod==p-1,'complete factorization')
 a=c['base'];need(isinstance(a,int)and 1<a<p and pow(a,p-1,p)==1,'full Fermat equation');modchecks+=1
 for q in powers:
  need(gcd(pow(a,(p-1)//q,p)-1,p)==1,'full Lucas order condition');modchecks+=1
 verified[p]=True
need(len(verified)==26,'prime certificate total')

def graph(w):
 rows=[2,0];root=0
 for c in w:
  k=len(rows);rows.append(0)
  for v in range(k):
   if c!='O'and v==root:rows[k]|=1<<v
   else:rows[v]|=1<<k
  if c=='J':root=k
 return rows
menu=json.loads((ROOT/'inputs/MENUS.json').read_text())['5'];allowed={r['word']for r in menu}
need(all(matrix(r['word'])==tuple(map(tuple,r['matrix']))for r in menu),'compact-menu identity')
receipts=[]
for item in json.loads((ROOT/'inputs/TARGETED-WITNESSES.json').read_text()):
 n,u,v,w=item['target'],item['left_index'],item['right_index'],item['word']
 need(item['prime_triple']==[n,3*n-2,5*n-4]and all(p in verified for p in item['prime_triple']),'target prime triple')
 need(w in allowed and min(u,v)>=2 and read(row(matrix(w),u),v)==n,'actual twelve-form indices')
 need([g(u),g(v)]==[item['g_left'],item['g_right']],'literal g coordinates')
 for ww,uu,vv in item['normalization_trace']:need(min(uu,vv)>=2 and read(row(matrix(ww),uu),vv)==n,'actual normalization trace')
 word=binary(v)+w+dagger(binary(u))+'J';s=state(word[:-1]);rows=graph(word);N=len(rows)
 need(word==item['graph_word'] and list(s)==item['source_profile']and s[0]-s[1]+s[2]==n,'full compiled source profile')
 need(N==item['vertices']==2+len(word) and N<=n.bit_length()-1+6,'actual adjacency order and bound')
 need(all(not(rows[i]>>i&1)and rows[i]<(1<<N)for i in range(N))and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(N)for j in range(i)),'literal tournament orientation')
 receipts.append(dict(target=n,template=w,indices=[u,v],source_profile=s,vertices=N,rows=rows))
# Strong normalization check for the only nontrivial boundary absorption used here.
for v in(2,3,7,191317):need(read(row(matrix('JOZOJ'),4),v)==read(row(matrix('JOZJ'),8),v),'boundary absorption control')
need(row(matrix('JOZOJ'),4)==row(matrix('JOZJ'),8),'all-index boundary absorption row identity')
print(json.dumps(dict(status='PASS',numerical_rows=numerical_row_count,
 fixed_row_controls=controls,fixed_row_nodes=nodes,scalar_singleton_candidates=singletons,
 full_matrix_count=len(forms),complete_small_targets=CAP,exhaustive_hostile=dict(target=2134117,witnesses=critical,statistics=critical_stats),
 Lucas_certified_primes=len(verified),Lucas_modular_checks=modchecks,target_graph_certificates=receipts,
 literal_large_path_enumeration=False,large_first_witness_search_statistics='author-only, not replayed',
 scope='Complete finite-template inverse and exact targeted certificates; universal occupancy remains OPEN.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
