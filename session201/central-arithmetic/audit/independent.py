"""Independent bounded-central-word audit. No author-code imports."""
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
if not(ROOT/'central_forms.h').is_file():ROOT=ROOT/'session202-central-arithmetic'
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

words=central(8);forms={};minima=[None]*9
for w in words:
 A=mul(Q,matrix(w));CC=mul(mul(transpose(PP),A),PP)
 coeff=min(Fraction(4*CC[i][j],UU[i]*UU[j])for i in range(4)for j in range(4))
 need(coeff.denominator==1 and coeff>=2,'exact integer product coefficient')
 b=int(coeff);minima[len(w)]=b if minima[len(w)]is None else min(b,minima[len(w)])
 left=mul(transpose(PP),A)
 for j in(0,2):need(all(left[i][j]>=0 for i in range(4))and left[2][j]>0 and left[3][j]>0,'all-index positive inverse coefficients')
 if A not in forms or(len(w),w)<(len(forms[A][0]),forms[A][0]):forms[A]=(w,b)
need(len(words)==1095 and len(forms)==592,'word and matrix census')
need(minima==[2,4,4,8,12,20,36,60,100],'all-word coefficient minima by actual length')
need(all(transpose(A)in forms for A in forms),'full matrix set transpose closure')
for u in range(2,4097):
 need(word(T,unitword(u))==tau(u),'direct recursive physical unit word')
 need(len(unitword(u))==(u-1).bit_length()-1,'exact unit source order')
 need(act('O',tau(u))==tau(2*u)and act('Z',tau(u))==tau(2*u-1),'two source-index maps')

# Independently certify every generated header entry, not merely a header hash.
header=(ROOT/'central_forms.h').read_text()
hforms=json.loads(re.search(r'forms\[\]\[4\]\[4\]=(.*);',header).group(1).replace('{','[').replace('}',']'))
hlens=list(map(int,re.search(r'lengths\[\]=\{(.*?)\}',header).group(1).split(',')))
hbounds=list(map(int,re.search(r'lower4\[\]=\{(.*?)\}',header).group(1).split(',')))
need(len(hforms)==len(hlens)==len(hbounds)==len(forms),'header shape')
for A,k,b in zip(hforms,hlens,hbounds):
 key=tuple(tuple(r)for r in A)
 need(key in forms and len(forms[key][0])==k and forms[key][1]==b,'header exact symbolic matrix/length/bound')
AMAX=max(abs(x)for A in hforms for row in A for x in row)
need(18*AMAX*100000001<2**63,'documented-range signed-int64 expression bound')

# Rewrite components and exact boundary absorption, retaining both indices.
def comp(w):
 seen={w};todo=[w]
 for x in todo:
  for j in range(len(x)-1):
   s=x[j:j+2]
   if s in('JZ','OJ'):
    y=x[:j]+('OJ'if s=='JZ'else'JZ')+x[j+2:]
    if y not in seen:seen.add(y);todo.append(y)
 return seen
seen=set();kept=[];certs=[];absorption=0
for w in words:
 if w in seen:continue
 group=comp(w);seen|=group
 need(all(matrix(x)==matrix(w)for x in group),'global exact local-rewrite equality')
 loose=sorted(x for x in group if x and(x[0]!='J'or x[-1]!='J'))
 if not loose:kept.append(min(group));continue
 x=loose[0];start=next(j for j,c in enumerate(x)if c=='J');end=max(j for j,c in enumerate(x)if c=='J')+1
 lead,core,tail=x[:start],x[start:end],x[end:]
 need(len(core)<len(w)and core in words,'strict shorter central core')
 for u,v in((2,2),(3,5),(5,3),(11,17),(17,11)):
  uu=word(tau(u),dagger(tail))[2];vv=word(tau(v),lead)[2]
  need(word(tau(u),dagger(tail))==tau(uu)and word(tau(v),lead)==tau(vv),'absorbed endpoints remain indexed physical unit sources')
  need(target(mul(Q,matrix(w)),u,v)==target(mul(Q,matrix(core)),uu,vv),'directed-index boundary absorption equality')
  absorption+=1
 certs.append((w,x,core))
need(len(kept)==235 and len(certs)==376,'reduced menu and discarded component counts')
need([sum(len(w)==k for w in kept)for k in range(9)]==[1,1,1,1,3,8,21,55,144],'reduced menu by length')

N=5001
def scan(formitems,N,ordered):
 out={};tests=0
 for A,(w,b)in formitems:
  utop=math.isqrt(4*N//b)if ordered else 2*N//b
  for u in range(2,utop+1):
   row=tuple(sum(tau(u)[i]*A[i][j]for i in range(4))for j in range(4))
   for v in range(u if ordered else 2,4*N//(b*u)+1):
    n=dot(row,tau(v));tests+=1
    need(n>=3 and n%2==1,'physical positive odd target')
    if n<=N:out[n]=min(out.get(n,99),len(w))
 return out,tests
full,tests=scan(forms.items(),N,True)
redforms=[(mul(Q,matrix(w)),(w,forms[mul(Q,matrix(w))][1]))for w in kept]
reduced,rtests=scan(redforms,N,False)
reducedwrong,_=scan(redforms,N,True)
need(full==reduced,'reduced menu preserves exact minimum central lengths with both index orders')
source=bank(N-1);direct={s[0]-s[1]+s[2]for s in source if s[0]-s[1]+s[2]<=N}
need(set(full)==direct,'complete physical source comparison at finite target cap')
gaps=[n for n in range(3,N+1,2)if n not in full]
expected=[7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647]
need(gaps==expected,'complete small gap list')

# Complete direct target scan and a separately solved congruence progression.
def egcd(a,b):
 if b==0:return a,1,0
 g,x,y=egcd(b,a%b);return g,y,x-(a//b)*y

def invert(n,progression):
 found=set()
 for A,(w,b)in forms.items():
  for u in range(2,math.isqrt(4*n//b)+1):
   row=tuple(sum(tau(u)[i]*A[i][j]for i in range(4))for j in range(4));alpha,gamma,beta,delta=row
   upper=4*n//(b*u)
   if not progression:
    vs=range(u,upper+1)
   else:
    vs=[]
    for q in(0,1):
     h,x,y=egcd(2*beta,alpha);rhs=n-delta-(beta+gamma)*q
     if rhs%h:continue
     modulus=alpha//h;res=(x*(rhs//h))%modulus
     low=(u+1-q)//2;high=(upper-q)//2
     first=low+(res-low)%modulus
     vs.extend(2*t+q for t in range(first,high+1,modulus))
   for v in vs:
    if dot(row,tau(v))==n:found.add((A,u,v))
 return found
inverse_checks=0;critical={}
for n in list(range(3,301,2))+expected+[1877,17383]:
 ans=invert(n,True)
 need(bool(ans)==(n not in expected),'exact progression membership')
 if n in(1877,17383)or n<=101:
  need(ans==invert(n,False),'progression candidates equal full direct inverse domain')
 if n in(1877,17383):
  need(min(len(forms[A][0])for A,u,v in ans)==8,'complete length-seven nonmembership')
  critical[n]=[dict(word=forms[A][0],u=u,v=v)for A,u,v in sorted(ans)]
 inverse_checks+=1

# Literal graph enumeration, independently using DFS of directed paths.
def literal_arcs(n,arcs):
 need(all((i,i)not in arcs for i in range(n))and all(((i,j)in arcs)+((j,i)in arcs)==1 for i in range(n)for j in range(i)),'literal tournament arcs')
 hp=0;hc=0;steps=0
 adjacent=[tuple(j for j in range(n)if(i,j)in arcs)for i in range(n)]
 def visit(start,last,mask):
  nonlocal hp,hc,steps
  steps+=1
  if mask==(1<<n)-1:
   hp+=1
   if start==0 and(last,0)in arcs:hc+=1
   return
  for j in adjacent[last]:
   if not(mask>>j&1):visit(start,j,mask|1<<j)
 for j in range(n):visit(j,j,1<<j)
 return hp,hc,steps
literal_critical=[]
for n,v in((1877,5),(17383,55)):
 w='JJJJJOJJ';u=3;ww=unitword(v)+w+dagger(unitword(u))+'J';count,arcs=graph(ww)
 need(target(mul(Q,matrix(w)),u,v)==n,'named critical arithmetic witness')
 hp,hc,steps=literal_arcs(count,arcs)
 need((hp,hc)==(n,1)and count==1+(u-1).bit_length()+(v-1).bit_length()+len(w),'literal critical graph and exact size')
 literal_critical.append(dict(target=n,vertices=count,word=ww,paths=hp,cycles=hc,partial_paths=steps))
repairs=json.loads((ROOT/'dependencies/REPAIR_GRAPHS.json').read_text());literal_repairs=[]
need(sorted(r['H']for r in repairs)==[n for n in expected if n not in(7,21)],'exact finite repair target set')
for r in repairs:
 rows=r['rows'];n=len(rows);need(all(0<=z<1<<n for z in rows),'bounded row bits')
 arcs={(i,j)for i in range(n)for j in range(n)if rows[i]>>j&1}
 hp,hc,steps=literal_arcs(n,arcs);need((hp,hc)==(r['H'],r['C'])and hc in(2,3),'literal repair paths and cycles')
 literal_repairs.append(dict(target=hp,vertices=n,cycles=hc,partial_paths=steps))
# Verify the exact size coefficient from each minimum, without floating logs.
for k,b in enumerate(minima):
 # If b>=4*2^j then uv<=n/2^j, so size <=floor(log2 n)+k+3-j.
 j=max(x for x in range(-1,10)if Fraction(b,4)>=Fraction(2)**x)
 need(k+3-j<=7,'all nine lengths imply logarithmic graph bound')
print(json.dumps(dict(status='PASS',raw_words=len(words),distinct_matrices=len(forms),all_word_lower4_minima=minima,normal_words=len(kept),discarded_components=len(certs),absorption_map_checks=absorption,
 target_cap=N,full_arithmetic_pair_tests=tests,reduced_directed_pair_tests=rtests,finite_targets=len(full),complete_physical_sources=len(source),finite_gaps=gaps,
 reduced_wrong_index_order_loses=sorted(set(full)-set(reducedwrong)),progression_controls=inverse_checks,critical_inverse=critical,literal_critical=literal_critical,literal_repairs=literal_repairs,
 header_matrices_checked=len(hforms),maximum_matrix_entry=AMAX,int64_absolute_expression_bound_at_max_cap=18*AMAX*100000001,
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='All-order compiler and reductions reviewed; finite exact family comparison and literal graphs only. Universal occupancy remains OPEN.'),indent=2,sort_keys=True))
