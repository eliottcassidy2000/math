from itertools import product,combinations
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
def trim(f):
 while len(f)>1 and not f[-1]:f.pop()
 return f
def add(f,g):
 r=[0]*max(len(f),len(g))
 for i,x in enumerate(f):r[i]+=x
 for i,x in enumerate(g):r[i]+=x
 return trim(r)
def mul(f,g):
 r=[0]*(len(f)+len(g)-1)
 for i,x in enumerate(f):
  for j,y in enumerate(g):r[i+j]+=x*y
 return trim(r)
def basis(k,p):
 f=[0]*p+[1]
 for j in range(1,k+1):f=mul(f,[-j*j,0,1])
 return f
def step_poly(f,c):
 if c!='D':return [0]*(1 if c=='X'else 2)+f
 r=[0,0]+f
 for k,x in enumerate(f):
  for j in range(1,k+1,2):r[k-j+1]+=2*x*comb(k,j)
 return trim(r)
def step_vec(v,c,p):
 if c=='X'and p==1:return v[:],2
 r=[0]*(len(v)+1)
 for k,x in enumerate(v):
  r[k]+=x*(k*k+6*k+(3 if p==1 else 5)if c=='D'else(k+1)**2)
  r[k+1]+=x
  if c=='D'and k:r[k-1]+=x*2*k*(k if p==1 else k+1)*(2*k+1)
 return trim(r),3-p if c=='X'else p
def apply_poly(f,w):
 for c in w:f=step_poly(f,c)
 return f
def apply_vec(v,w,p):
 for c in w:v,p=step_vec(v,c,p)
 return v,p
def det(A):
 A=[list(map(Fraction,r))for r in A];n=len(A);z=Fraction(1)
 for k in range(n):
  s=next((s for s in range(k,n)if A[s][k]),None)
  if s is None:return 0
  if s!=k:A[s],A[k]=A[k],A[s];z=-z
  d=A[k][k];z*=d
  for i in range(k+1,n):
   f=A[i][k]/d
   for j in range(k+1,n):A[i][j]-=f*A[k][j]
 return z
def minors(A,K):
 n=len(A);m=len(A[0]);count=0
 for k in range(1,K+1):
  for I in combinations(range(n),k):
   for J in combinations(range(m),k):need(det([[A[i][j]for j in J]for i in I])>=0,'TN minor');count+=1
 return count
words=sorted(set([''.join(w)for n in range(7)for w in product('PD',repeat=n)]+[''.join(w)for n in range(7)for w in product('DX',repeat=n)if w.count('X')%2==0]))
columns=series=prodmin=hmin=0;examples=[]
for p in (1,2):
 for w in words:
  for k in range(5):
   f=apply_poly(basis(k,p),w);v,q=apply_vec([0]*k+[1],w,p);need(q==p,'return parity');g=[0]
   for j,x in enumerate(v):g=add(g,[x*t for t in basis(j,p)])
   need(f==g,'full virtual basis operator');columns+=1
  f=[0]*p+[1];v=[1];mu=[1]
  for _ in range(7):
   f=apply_poly(f,w);v,q=apply_vec(v,w,p);need(sum(f)==v[0] and q==p,'independent moments');mu.append(sum(f))
  series+=1
  if len(w)<=3:
   T=[[0]*4 for _ in range(4)]
   for k in range(4):
    z,_=apply_vec([0]*k+[1],w,p)
    for j in range(min(4,len(z))):T[j][k]=z[j]
   prodmin+=minors(T,4);hmin+=minors([[mu[i+j]for j in range(4)]for i in range(4)],4)
  if w in ('PD','DP','DXDX'):examples.append({'p':p,'word':w,'moments':mu[:5]})
# Positive pivots and exact all-order comparison polynomials.
pivots=0
for p in (1,2):
 ell=Fraction(3 if p==1 else 5)
 for n in range(100):
  if n:ell=n*n+6*n+(3 if p==1 else 5)-Fraction(2*n*(n if p==1 else n+1)*(2*n+1),ell)
  need(ell>0,'positive pivot')
  if p==1:need(ell>=(n+1)**2,'odd pivot bound')
  elif n>=8:need(ell>=n*n+2*n-4,'even pivot bound')
  pivots+=1
# Literal directed-word counts on actual constructed graphs.
def trans(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def join(a,b):
 n=len(a);m=len(b);return [x+(((1<<m)-1)<<n)for x in a]+[x<<n for x in b]
def wrap(a):
 n=len(a);return [x|(1<<n)for x in a]+[1<<(n+1),(1<<n)-1]
def hp(a):
 N=(1<<len(a))-1
 def go(S,v):
  if S==N:return 1
  z=0;U=a[v]&~S
  while U:
   b=U&-U;U-=b;z+=go(S|b,b.bit_length()-1)
  return z
 return sum(go(1<<v,v)for v in range(len(a)))
direct=0
for p in (1,2):
 for w in words:
  if len(w)>3:continue
  for reps in (1,2):
   n=p+reps*(2*(len(w)-w.count('X'))+w.count('X'))
   if n>10:continue
   a=trans(p);f=[0]*p+[1]
   for _ in range(reps):
    for c in w:a=wrap(a)if c=='D'else join(a,trans(1 if c=='X'else 2))
    f=apply_poly(f,w)
   need(hp(a)==sum(f),'literal actual word count');direct+=1
# Independently named-vertex mirror equality, including nontransitive inputs.
mirrors=arcs=0
for a in ([0],[2,0],[2,4,1],[6,12,24,17,3]):
 n=len(a)
 for t in range(1,7):
  for j in range(t):
   N=n+2*t;rows=a[:]+[0]*(2*t);O=lambda k:n+k;I=lambda k:n+t+k
   for x in range(t):
    for y in range(x+1,t):rows[O(x)]|=1<<O(y);rows[I(y)]|=1<<I(x)
    for v in range(n):rows[I(x)]|=1<<v;rows[v]|=1<<O(x)
    for y in range(t):rows[O(x)if x==j and y==j else I(y)]|=1<<(I(y)if x==j and y==j else O(x))
   # Target block labels: initial I>j, central [I_j, module, O_j], final O>j.
   pre=[I(k)for k in range(t-1,j,-1)];middle=[I(k)for k in range(j-1,-1,-1)]+list(range(n))+[O(k)for k in range(j)];post=[O(k)for k in range(j+1,t)]
   order=pre+middle+[O(j),I(j)]+post
   expected=join(join(trans(t-1-j),wrap(join(join(trans(j),a),trans(j)))),trans(t-1-j))
   actual=[sum(1<<k for k,v in enumerate(order)if rows[u]>>v&1)for u in order]
   need(actual==expected,'literal physical mirror decomposition');mirrors+=1;arcs+=N*(N-1)//2
bad=[0]*4+[1];vals=[1]
for _ in range(2):bad=step_poly(bad,'D');vals.append(sum(bad))
need(vals==[1,17,241]and vals[0]*vals[2]-vals[1]**2==-48,'prefix hostile')
print(json.dumps({'status':'PASS','distinct_words_through6':len(words),'full_basis_columns':columns,'moment_series':series,'production_minors':prodmin,'Hankel_minors':hmin,'positive_pivots':pivots,'literal_graph_counts':direct,'literal_mirror_graphs':mirrors,'mirror_arcs':arcs,'prefix_hostile':vals,'examples':examples},sort_keys=True,indent=2))
