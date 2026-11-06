"""Independent compact-central quotient audit. No author-code imports."""
import hashlib,itertools,json,math
from pathlib import Path
from fractions import Fraction

T=(1,0,2,1);B=(3,1,3,1);AUX=(1,0,2,0)
M={
 'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
P=(T,(1,0,4,1),B,(5,2,6,3))
def need(p,m):
 if not p:raise RuntimeError(m)
def act(c,v):return tuple(sum(x*y for x,y in zip(r,v))for r in M[c])
def word(s,w):
 for c in w:s=act(c,s)
 return s
def dagger(w):return w[::-1].translate(str.maketrans('OZJ','ZOJ'))
def cov(v):return tuple(sum(v[i]*Q[i][j]for i in range(4))for j in range(4))
def dot(v,w):return sum(x*y for x,y in zip(v,w))
def K(v,w):return dot(cov(v),w)
def positive(v):
 d,e,u,z=v;return d-e-z,u-d+e-z,e,z-e

def rank2(w,q):
 # Exact rational elimination, rather than the author's minor guard.
 i=next(j for j,x in enumerate(w)if x)
 r=Fraction(q[i],w[i])
 return any(Fraction(y)!=r*x for x,y in zip(w,q))
def determinant(A):
 A=[[Fraction(x)for x in row]for row in A];sgn=1
 for j in range(len(A)):
  i=next((i for i in range(j,len(A))if A[i][j]),None)
  if i is None:return 0
  if i!=j:A[i],A[j]=A[j],A[i];sgn=-sgn
  for k in range(j+1,len(A)):
   s=A[k][j]/A[j][j]
   A[k]=[a-s*b for a,b in zip(A[k],A[j])]
 return sgn*math.prod(A[i][i]for i in range(len(A)))

def bank(H,paired=False,prune=False):
 initial=(B,AUX)if paired else T
 height=lambda v:v[0][2]if paired else v[2]
 out={initial:(0,'')};buckets=[set()for _ in range(H+1)];buckets[height(initial)].add(initial)
 for u in range(H+1):
  for v in sorted(buckets[u]):
   cost,w=out[v]
   for c in 'JZO':
    z=tuple(act(c,x)for x in v)if paired else act(c,v);h=height(z)
    need(h>u,'strict monotone physical cover')
    if h>H or(prune and not rank2(*z)):continue
    recipe=(cost+int(c=='J'),w+c)
    if z not in out:buckets[h].add(z);out[z]=recipe
    elif recipe<out[z]:out[z]=recipe
 return out

def prime(n):return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def put(D,k,c):D[k]=min(D.get(k,10**9),c)
def images(S,R,N):
 f={};g={};tests=0
 for v,(a,_)in S.items():
  vq=cov(v)
  for(w,q),(b,_)in R.items():
   m=dot(vq,w);h=dot(vq,q)
   need(3<=h<m<2*h and m%2==h%2==1,'joint physicality and parity')
   if m>N:continue
   tests+=1
   if prime(m):put(f,m,a+b);put(g,m,a+b)
   power=2
   while m+2*(power-1)*h<=N:
    n=m+2*(power-1)*h
    if prime(n):put(f,n,a+b)
    power*=2
 return f,g,tests

def graph(w):
 arcs={(0,1)};n=2;r=0
 for c in w:
  for v in range(n):arcs.add((n,v)if c!='O'and v==r else(v,n))
  if c=='J':r=n
  n+=1
 return n,arcs

def literal(w):
 n,arcs=graph(w);hp=0;hc=0;steps=0
 def visit(start,last,used):
  nonlocal hp,hc,steps
  steps+=1
  if len(used)==n:
   hp+=1
   if start==0 and(last,0)in arcs:hc+=1
   return
  for x in range(n):
   if x not in used and(last,x)in arcs:visit(start,x,used|{x})
 for s in range(n):visit(s,s,{s})
 return dict(word=w,vertices=n,paths=hp,cycles=hc,partial_directed_paths=steps)


import re
ROOT=Path(__file__).resolve().parent.parent
if not(ROOT/'MENUS.json').is_file():ROOT=ROOT/'session201-compact-central'
PP=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1));UU=(1,1,2,2)
I=tuple(tuple(int(i==j)for j in range(4))for i in range(4))
def transpose(A):return tuple(zip(*A))
def mul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j]for k in range(4))for j in range(4))for i in range(4))
def matrix(w):
 a=I
 for c in w:a=mul(M[c],a)
 return a

def central(limit):
 out=['']
 for k in range(1,limit+1):
  for chars in itertools.product('JZO',repeat=k):
   if chars[0]==chars[-1]=='J':out.append(''.join(chars))
 return out

def unitword(u):
 need(u>=2,'unit source index range')
 if u==2:return ''
 if u%2:return unitword((u+1)//2)+'Z'
 return unitword(u//2)+'O'
def gtable(N):
 g=[0]*(N+1);g[1]=1
 for n in range(2,N+1):
  if n%2==0:g[n]=g[n//2]
  else:
   q=(n+1)//2;g[n]=g[q]+q-q%2
 return g
GG=gtable(20001)
def tau(u):return(GG[u],u%2,u,1)
def target(A,u,v):return dot(tau(u),tuple(dot(r,tau(v))for r in A))


def component(w):
 seen={w};todo=[w]
 for s in todo:
  for j in range(len(s)-1):
   if s[j:j+2]in('JZ','OJ'):
    t=s[:j]+('OJ'if s[j:j+2]=='JZ'else'JZ')+s[j+2:]
    if t not in seen:seen.add(t);todo.append(t)
 return seen
raw=central(8);seen=set();normal=[]
for w in raw:
 if w in seen:continue
 co=component(w);seen|=co
 need(all(matrix(s)==matrix(w)for s in co),'rewrite matrix equality')
 if not any(s and(s[0]!='J'or s[-1]!='J')for s in co):normal.append(min(co))
need(len(normal)==235,'retained component inventory')
author_menus=json.loads((ROOT/'MENUS.json').read_text());quotients={};allpairs={};index_checks=0
for limit in(5,6,8):
 retained={w for w in normal if len(w)<=limit};done=set();items=[];pairs=[]
 for w in sorted(retained,key=lambda w:(len(w),w)):
  p=min(component(dagger(w)))
  need(p in retained and min(component(dagger(p)))==w,'dagger component involution')
  need(mul(Q,matrix(w))==transpose(mul(Q,matrix(p))),'all-index transpose identity')
  for u,v in((2,3),(3,2),(3,5),(11,23),(55,3)):
   need(target(mul(Q,matrix(w)),u,v)==target(mul(Q,matrix(p)),v,u),'physical swapped indices')
   index_checks+=1
  if w in done:continue
  done|={w,p};items.append(w);pairs.append((w,p))
 need(len(done)==len(retained),'quotient covers every retained component')
 need(len(items)=={5:12,6:26,8:137}[limit],'compact quotient size')
 # Author may choose the other member; compare component-pair sets exactly.
 expected={frozenset((w,p))for w,p in pairs}
 actual={frozenset((r['word'],r['dagger_partner']))for r in author_menus[str(limit)]}
 need(actual==expected and len(actual)==len(author_menus[str(limit)]),'author menus exactly one representative per orbit')
 for r in author_menus[str(limit)]:need(tuple(tuple(row)for row in r['matrix'])==mul(Q,matrix(r['word'])),'exact expanded author matrix')
 quotients[limit]=items;allpairs[limit]=pairs
explicit=['','J','JJ','JJJ','JOZJ','JJZJ','JJJJ','JOOZJ','JJZZJ','JOZJJ','JJZJJ','JJJJJ']
need({frozenset((w,min(component(dagger(w)))))for w in explicit}=={frozenset(x)for x in allpairs[5]},'twelve explicitly printed templates')
forms={};minima={}
for w in raw:
 A=mul(Q,matrix(w));CC=mul(mul(transpose(PP),A),PP)
 b=min(Fraction(4*CC[i][j],UU[i]*UU[j])for i in range(4)for j in range(4))
 need(b.denominator==1 and b>=2,'exact product bound')
 if A not in forms or len(w)<len(forms[A][0]):forms[A]=(w,int(b))
 minima[len(w)]=min(minima.get(len(w),10**9),int(b))
need(len({A for A,(w,b)in forms.items()if len(w)<=5})==35,'transpose-closed length-five matrix count')

repairs=json.loads((ROOT/'REPAIRS.json').read_text());rawrep=json.loads((ROOT/'RAW-REPAIRS.json').read_text())
need(rawrep['cap']==2134117 and repairs['gaps5']==rawrep['gaps5']and len(rawrep['gaps5'])==65,'fixed finite-prefix target list')
old=json.loads((ROOT/'dependencies/REPAIR_GRAPHS.json').read_text());oldgaps={r['H']for r in old}|{7,21}
need(sorted(set(rawrep['gaps5'])-oldgaps)==[r['target']for r in repairs['repairs']],'all 36 additional exceptions covered')
certificate_receipts=[];graph_header=[]
for r in repairs['repairs']:
 n,u,v,w=r['target'],r['left_index'],r['right_index'],r['center']
 need(target(mul(Q,matrix(w)),u,v)==n,'direct arithmetic certificate')
 ww=unitword(v)+w+dagger(unitword(u))+'J';s=word(T,ww[:-1]);N,arcs=graph(ww)
 need(ww==r['graph_word']and list(s)==r['source_profile']and s[0]-s[1]+s[2]==n,'complete literal source profile certificate')
 need(N==r['vertices']==1+(u-1).bit_length()+(v-1).bit_length()+len(w),'exact certificate graph order')
 need(N<=n.bit_length()-1+6,'conditional repair size bound')
 rows=[sum(1<<j for j in range(N)if(i,j)in arcs)for i in range(N)]
 graph_header.append((n,N,rows));certificate_receipts.append(dict(target=n,vertices=N,word=ww,profile=s))
need(len(graph_header)==36 and max(N for n,N,rows in graph_header)==24,'all 36 compiled repair graphs')
for r in old:need(r['order']==len(r['rows'])<=r['H'].bit_length()-1+6,'all 27 inherited graph size bounds')
need(next(r for r in repairs['repairs']if r['target']==2134117)['center']=='JOOZZJ','largest finite hostile repair')
need(target(mul(Q,matrix('JOOZZJ')),5,14256)==2134117,'explicit largest endpoint witness')
for length in range(7):
 b=minima[length];j=max(j for j in range(-1,10)if Fraction(b,4)>=Fraction(2)**j)
 need(length+3-j<=6,'all relevant menu lengths imply logarithmic size bound')
# Independent C++ inputs: twelve reconstructed forms, both index orders;
# and all repair adjacencies generated by this checker from literal words.
head=['// Generated by independent.py from independently reconstructed words.']
head+=['const int NF=12;','const long long AF[12][4][4]={'+','.join('{'+','.join('{'+','.join(map(str,row))+'}'for row in mul(Q,matrix(w)))+'}'for w in explicit)+'};']
head+=['const int BF[12]={'+','.join(str(forms[mul(Q,matrix(w))][1])for w in explicit)+'};']
head+=['const int NG=36;','const int EXPECTED[36]={'+','.join(str(n)for n,N,rows in graph_header)+'};',
 'const int ORDERS[36]={'+','.join(str(N)for n,N,rows in graph_header)+'};',
 'const uint32_t ROWS[36][24]={'+','.join('{'+','.join(map(str,rows))+'}'for n,N,rows in graph_header)+'};']
(Path(__file__).parent/'independent_data.h').write_text('\n'.join(head)+'\n')
print(json.dumps(dict(status='PASS',quotient_sizes={k:len(v)for k,v in quotients.items()},retained_sizes={k:sum(len(w)<=k for w in normal)for k in(5,6,8)},
 daggers_checked=index_checks,explicit_twelve=explicit,generated_repair_profiles=certificate_receipts,inherited_order_checks=len(old),
 strict_scope='All-index dagger quotient and physical certificates verified; raw finite census is separately replayed; every uniform occupancy lemma remains OPEN.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
