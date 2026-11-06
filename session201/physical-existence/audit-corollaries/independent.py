"""Independent asymmetric-budget and J-crossing audit. No author-code imports."""
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


def cuts(w,seed):
 l=[T]
 for c in w:l.append(act(c,l[-1]))
 r=[seed]
 for c in w[::-1]:r.append(act(dagger(c),r[-1]))
 return l,r[::-1]
cone=((1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1))
C=[[K(act('J',v),w)for w in cone]for v in cone];A=[act('J',v)[2]for v in cone]
residual=[[16*C[i][j]-7*A[i]*A[j]for j in range(4)]for i in range(4)]
need(C==[[2,1,4,3],[1,1,2,2],[4,2,7,6],[3,2,6,5]],'exact J pairing matrix')
need(A==[2,1,4,3]and residual==[[4,2,8,6],[2,9,4,11],[8,4,0,12],[6,11,12,17]],'exact nonnegative J residual')
words=rectangles=noninteger=Jcross=joint_checks=0;floor_product_failures=0
for le in range(8):
 for letters in itertools.product('OZJ',repeat=le):
  ww=''.join(letters)
  for seed in(T,B):
   l,r=cuts(ww,seed);n=K(l[0],r[0]);words+=1
   # Probe every possible left-height threshold and a noninteger point just
   # below each next left cover, rather than using one fixed L menu.
   choices={Fraction(2),Fraction(5,2)}
   for x in l:
    choices.add(Fraction(x[2]));
    if x[2]>2:choices.add(Fraction(2*x[2]-1,2));choices.add(Fraction(3*x[2]-1,3))
   for L in sorted(choices):
    R=max(Fraction(seed[2]),Fraction(3*(n-1),1)/L)
    need(L>=2 and R>=seed[2]and L*R>=3*(n-1),'admissible exact real budgets')
    i=max(i for i,x in enumerate(l)if x[2]<=L)
    need(l[i][2]<=L and r[i][2]<=R,'last-left cut simultaneously meets both budgets')
    need(l[i][2]<=math.floor(L)and r[i][2]<=math.floor(R),'rounding exact after real-budget test')
    need(ww[:i].count('J')+dagger(ww[i:]).count('J')==ww.count('J'),'rectangle preserves J cost')
    if seed==B:
     _,q=cuts(ww,AUX);actual=word(T,ww)
     need(K(l[i],q[i])==2*actual[2]-actual[3]and all(0<=x<=y for x,y in zip(q[i],r[i])),'rectangle preserves synchronized sidecar')
     joint_checks+=1
    rectangles+=1;noninteger+=L.denominator>1 or R.denominator>1
    floor_product_failures+=math.floor(L)*math.floor(R)<3*(n-1)
   cross=next((i for i,(v,w)in enumerate(zip(l,r))if v[2]>=w[2]),None)
   if cross not in(None,0)and ww[cross-1]=='J':
    candidate=min(max(l[cross-1][2],r[cross-1][2]),max(l[cross][2],r[cross][2]))
    need(candidate<=math.isqrt((16*n)//7),'conditional J-crossing cutoff')
    Jcross+=1
sharp=[];graphs=[]
for k in range(1,33):
 d=2**k+1;v=word(T,'Z'*k+'J');n=K(act('J',v),v)
 need(v==(2*d-1,d-1,2*d,d)and n==7*d*d-2*d,'sharp physical family k>=1')
 height=act('J',v)[2];need(height==4*d-1,'sharp J height')
 ww='Z'*k+'JJJ'+'O'*k;l,r=cuts(ww,T)
 need(min(max(v[2],w[2])for v,w in zip(l,r))==height,'sharp best cut across specified word')
 need(Fraction(height**2,n)<Fraction(16,7),'approach sharp coefficient from below')
 sharp.append(dict(k=k,target=n,best_height=height,coefficient=str(Fraction(height**2,n))))
 if k<=2:
  g=literal(ww+'J');need((g['paths'],g['cycles'])==(n,1),'literal J-sharp graph DFS')
  graphs.append(g)
need(word(T,'J')==(3,1,3,1)and word(T,'J')!=(3,1,4,2),'empty-Z boundary regression')
print(json.dumps(dict(status='PASS',J_pairing_matrix=C,J_residual=residual,complete_word_length_cap=7,words=words,rectangle_checks=rectangles,
 noninteger_rectangle_checks=noninteger,real_valid_but_floor_product_small=floor_product_failures,synchronized_rectangle_checks=joint_checks,J_crossings=Jcross,
 sharp_family=sharp,literal_graph_controls=graphs,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 scope='Exact finite controls for independently reviewed all-order rectangle and conditional J-crossing results; no universal occupancy or globally optimal target-representation claim.'),indent=2,sort_keys=True))
