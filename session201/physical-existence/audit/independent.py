"""Independent adjacent-cut and even-diagonal audit. No author-code imports."""
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


# Matrix proof audit through an independent raw basis calculation.
cone=((1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1))
expected={
'O':((2,2,3,3),(2,2,0,0),(4,1,3,3),(4,1,0,0)),
'Z':((2,2,4,4),(2,2,1,1),(3,0,3,0),(3,0,3,0)),
'J':((2,1,4,3),(1,2,2,3),(4,2,2,3),(3,3,3,3))}
residuals={};rows={}
for c in 'OZJ':
 s=dagger(c)
 residual=tuple(tuple(3*K(act(c,v),w)-act(c,v)[2]*act(s,w)[2]-3*v[3]*w[3]for w in cone)for v in cone)
 need(residual==expected[c]and min(x for row in residual for x in row)>=0,'nonnegative adjacent-cut residual')
 residuals[c]=residual;rows[c]=dict(A=[act(c,v)[2]for v in cone],B=[act(s,v)[2]for v in cone])

def C(n):return math.isqrt(3*(n-1))
def F(v):return v[0]-v[1]+v[2]
def G(v):return v[0]+2*v[2]
def Hjoint(v):return 2*v[2]-v[3]
def cuts(w,seed):
 left=[T]
 for c in w:left.append(act(c,left[-1]))
 right=[seed]
 for c in w[::-1]:right.append(act(dagger(c),right[-1]))
 right.reverse()
 return left,right
cut_count=0;joint_cut_count=0;boundaries={'full_initial_equal':0,'balanced_no_crossing':0}
for length in range(9):
 for ww in itertools.product('OZJ',repeat=length):
  w=''.join(ww);actual=word(T,w)
  for base,label in ((T,'full'),(B,'balanced')):
   left,right=cuts(w,base);n=K(left[0],right[0]);maxima=[max(v[2],z[2])for v,z in zip(left,right)]
   need(all(K(v,z)==n for v,z in zip(left,right)),'every cut preserves same target')
   cross=next((i for i,(v,z)in enumerate(zip(left,right))if v[2]>=z[2]),None)
   if cross is None:
    need(base==B and not w and maxima==[3]and n==5,'balanced empty endpoint case')
    choose=0;boundaries['balanced_no_crossing']+=1
   elif cross==0:
    need(base==T and not w and maxima==[2]and n==3,'full empty endpoint case')
    choose=0;boundaries['full_initial_equal']+=1
   else:
    i=cross;c=w[i-1];v=left[i-1];z=right[i]
    aa=act(c,v)[2];bb=act(dagger(c),z)[2]
    need((maxima[i-1],maxima[i])==(bb,aa),'orientation of two adjacent candidate maxima')
    need(aa*bb<=3*(n-v[3]*z[3]),'actual adjacent inequality')
    choose=i if aa<=bb else i-1
   need(maxima[choose]<=C(n),'chosen physical cut meets tighter cutoff')
   for i in range(len(w)+1):need(w[:i].count('J')+dagger(w[i:]).count('J')==w.count('J'),'cut cost preservation')
   if base==B:
    _,aux=cuts(w,AUX)
    need(n==G(actual),'joint primary target')
    need(all(K(v,q)==Hjoint(actual)and all(0<=a<=b for a,b in zip(q,z))for v,z,q in zip(left,right,aux)),'common-cut sidecar and boundedness')
    joint_cut_count+=1
   cut_count+=1

N=5001;cap=C(N);small=bank(cap);sync=bank(cap,True);direct=bank(N-1)
df={};dg={3:0};dj={}
for s,(cost,ww)in direct.items():
 f,g,h=F(s),G(s),Hjoint(s)
 need(f>=s[2]+1,'complete direct forward cutoff')
 if f<=N:put(df,f,cost)
 if g<=N:put(dg,g,cost);put(dj,(g,h),cost)
fimage={};bimage={3:0};jimage={};mixed={3:0};pair_tests=0
for v,(cv,av)in small.items():
 co=cov(v)
 for w,(cw,aw)in small.items():
  n=dot(co,w)
  if n<=N and max(v[2],w[2])<=C(n):put(fimage,n,cv+cw)
 for(w,q),(cw,beta)in sync.items():
  m,h=dot(co,w),dot(co,q)
  if m>N or max(v[2],w[2])>C(m):continue
  pair_tests+=1;cost=cv+cw
  need(3<=h<m<2*h,'strict joint cone')
  put(bimage,m,cost);put(jimage,(m,h),cost);put(mixed,m,cost)
  power=2
  while m+2*(power-1)*h<=N:
   n=m+2*(power-1)*h
   need(max(v[2],w[2])<=C(n//power)and 6*power<=n+1,'branch cutoff and finite exponent bound')
   put(mixed,n,cost);power*=2
need(fimage==df and bimage==dg and jimage==dj and mixed==df,'complete smaller inverse and exact minimum J costs')
for target in(443,617):
 need(target in df and target not in dg and not any(m+2*h==target for m,h in dj),'unchanged no-k1 hostiles')
need((15,13)not in dj and sorted(h for m,h in dj if m==15)==[9,11],'unsynchronized false-positive rejection')

sharp=[]
for k in range(33):
 d=2**(k+1)+1;v=word(T,'Z'+'O'*k);w=word(T,'Z'*(k+1)+'J');n=3*d*d
 need(v==(3,int(k==0),3*2**k,1)and w==(2*d-1,d-1,2*d,d),'sharp source profiles including k0 endpoint')
 A=act('O',v)[2];BB=act('Z',w)[2]
 need(K(act('O',v),w)==n and(A,BB)==(3*(d-1),3*d),'sharp adjacent values')
 need(Fraction(A*BB,n-v[3]*w[3])==Fraction(9*(d-1),3*d-1),'sharp one-letter quotient')
 ww='Z'+'O'*(k+1)+'J'+'O'*(k+1);l,r=cuts(ww,T)
 best=min(max(a[2],b[2])for a,b in zip(l,r))
 l2,r2=cuts(ww[:-1],B)
 need(F(word(T,ww))==G(word(T,ww[:-1]))==n and best==3*(d-1),'sharp full word')
 need(l2==l[:-1]and r2==r[:-1],'balanced boundary is exactly shifted full cuts')
 need(min(max(a[2],b[2])for a,b in zip(l2,r2))==best,'balanced sharp-word minimum')
 sharp.append(dict(k=k,n=n,best=best,coefficient=str(Fraction(best**2,n-1))))

family=0;profiles_checked=0;literal_graphs=[];first_profiles={}
for p in range(13):
 for q in range(1,14):
  x=word(T,'O'*p+'Z'*q+'O'*p)
  need(G(x)==3+2**q*(2**(2*p+2)-1),'all-form palindromic balanced identity')
  family+=1
for a in range(2,42,2):
 p=(a-2)//2;q=a+2;r=1<<p;n=((1<<(a+1))-1)**2+2
 recipes=['O'*a+'Z'*a,'O'*p+'Z'*q+'O'*(p+1),'O'*p+'Z'*q+'O'*p+'J','Z'*(p+1)+'O'*q+'Z'*p]
 states=[word(T,w)for w in recipes];x=word(T,'O'*p+'Z'*q+'O'*p)
 need(len(set(states))==4 and all(F(s)==n for s in states),'four distinct numerical profiles with same F')
 need(dagger(recipes[1])==recipes[3]and G(x)==n,'dagger mate and balanced precursor')
 need([states[i][2]for i in(0,1,3)]==[32*r**4-4*r**2+1,64*r**4-32*r**3+2*r,32*r**4+16*r**3-r+1],'explicit U separation formulas')
 need(states[2][3]==x[0]>1 and states[0][3]==states[1][3]==states[3][3]==1,'deletion separates fourth profile')
 need(states[0][1]==1 and states[1][1]==0,'endpoint separates original and twin')
 if a<=4:
  first_profiles[str(a)]=states
  for idx,ww in enumerate(recipes):
   receipt=literal(ww+'J');need((receipt['paths'],receipt['cycles'])==(n,1),'literal profile closing target')
   receipt.update(a=a,recipe=idx+1);literal_graphs.append(receipt)
  balanced=literal('O'*p+'Z'*q+'O'*p+'OJ')
  need((balanced['paths'],balanced['cycles'],balanced['vertices'])==(n,1,2*a+4),'literal balanced graph size and target')
  balanced.update(a=a,recipe='balanced');literal_graphs.append(balanced)
 profiles_checked+=4
print(json.dumps(dict(status='PASS',residual_matrices=residuals,height_rows=rows,word_cut_checks=cut_count,joint_word_cut_checks=joint_cut_count,boundaries=boundaries,
 target_cap=N,forward_height_cap=N-1,forward_profiles=len(direct),smaller_height_cap=cap,smaller_source_profiles=len(small),synchronized_states=len(sync),
 pair_tests=pair_tests,distinct_joint_values=len(jimage),C1_targets=len(df),balanced_targets=len(dg),minimum_J_costs_preserved=True,
 sharp_family=sharp,palindromic_family_checks=family,even_diagonal_profiles_checked=profiles_checked,first_diagonal_profiles=first_profiles,literal_graph_controls=literal_graphs,
 scope='Exact finite comparison and examples; all-order inequalities and family proofs reviewed separately; coefficient3 sharp only for specified-word cuts, not globally.',
 verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
