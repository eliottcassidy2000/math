from pathlib import Path
import itertools,json
P=Path(__file__).parent
# Independent contracted-gap subset recurrence, copied from prior audited implementation.
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def coeff(a,v):
 ex=[u for u in range(len(a)) if u!=v];n=len(ex);N=1<<n
 if not n:return [1,1]
 dp=[[[] for j in ex] for S in range(N)]
 for j,u in enumerate(ex):dp[1<<j][j]=[1,int(bool(a[v]>>u&1))]
 for S in range(1,N):
  for j,u in enumerate(ex):
   if not dp[S][j]:continue
   for l,w in enumerate(ex):
    if not S>>l&1:
     t=[int(bool(a[u]>>w&1)),int(bool(a[u]>>v&1) and bool(a[v]>>w&1))]
     dp[S|1<<l][l]=add(dp[S|1<<l][l],mul(dp[S][j],t))
 out=[]
 for j,u in enumerate(ex):out=add(out,mul(dp[-1][j],[1,int(bool(a[u]>>v&1))]))
 return (out+[0]*(len(a)+1))[:len(a)+1]
def tournaments(n):
 edges=list(itertools.combinations(range(n),2))
 for bits in itertools.product((0,1),repeat=len(edges)):
  a=[0]*n
  for (i,j),b in zip(edges,bits):
   if b:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a

def matching(a,v):
 # Vertex-disjoint matching size via exhaustive edge subsets (small controls).
 A=[j for j in range(len(a)) if a[v]>>j&1];B=[j for j in range(len(a)) if a[j]>>v&1]
 edges=[(i,j) for i in A for j in B if a[i]>>j&1];best=0
 for subset in range(1<<len(edges)):
  picked=[edges[j] for j in range(len(edges)) if subset>>j&1]
  if len({i for i,j in picked})==len(picked)==len({j for i,j in picked}):best=max(best,len(picked))
 return best
hist={};num=0
for n in range(1,6):
 for a in tournaments(n):
  for v in range(n):
   nu=matching(a,v);c=coeff(a,v)
   if [j for j,x in enumerate(c) if x]!=list(range(nu+2)):raise RuntimeError(('support',a,v,nu,c))
   hist[nu]=hist.get(nu,0)+1;num+=1
# Exact endpoint/cover tables via literal path enumeration, for small L/R and interval factors.
def paths(a):
 out=[]
 for S in range(1<<len(a)):
  vs=[j for j in range(len(a)) if S>>j&1]
  out.append([p for p in itertools.permutations(vs) if all(a[i]>>j&1 for i,j in zip(p,p[1:]))])
 return out

def side(a,cut):
 pp=paths(a);F=len(pp)-1;alpha=tau=0
 for S in range(F+1):
  for p in pp[S]:
   for q in pp[F^S]:
    if q and cut>>q[0]&1:
     alpha+=1
     if p and cut>>p[-1]&1:tau+=1
 e=sum(bool(p) and bool(cut>>p[-1]&1) for p in pp[F])
 return len(pp[F]),alpha,tau,e

def beta(a,cut):
 pp=paths(a);F=len(pp)-1
 return sum((not p or not(cut>>p[-1]&1)) for S in range(F+1) for p in pp[S] for q in pp[F^S])
def build(L,R,cl,cr,U):
 # Layout x, U, R, L, independent of producer's L,U,R,x convention.
 m,l,r=len(U),len(L),len(R);V=[('x',0)]+[('u',i) for i in range(m)]+[('r',i) for i in range(r)]+[('l',i) for i in range(l)]
 a=[0]*len(V);groups={'l':0,'u':1,'r':2}
 for i,(t,z) in enumerate(V):
  for j,(s,w) in enumerate(V):
   if i==j:continue
   if t==s:b=({'l':L,'r':R,'u':U}[t][z]>>w&1)
   elif t=='x':b=(cl>>w&1) if s=='l' else (cr>>w&1) if s=='r' else 0
   elif s=='x':b=not((cl>>z&1) if t=='l' else (cr>>z&1) if t=='r' else 0)
   else:b=groups[t]<groups[s]
   if b:a[i]|=1<<j
 return a

def induced(a,vertices):return [sum((a[u]>>v&1)<<j for j,v in enumerate(vertices)) for u in vertices]
small=[a for n in range(4) for a in tournaments(n)];graphs=[];factorchecks=0
for L in small:
 for cl in range(1<<len(L)):
  h,al,ta,e=side(L,cl)
  for R in small:
   hr=len(paths(R)[-1])
   for cr in range(1<<len(R)):
    be=beta(R,cr)
    rx=build([],R,0,cr,[]);hrx=len(paths(rx)[-1]);offset=h*hrx-hr*ta
    # Include a nontransitive module in addition to the two producer controls.
    for U,hu,qu in [([0],1,2),([2,0],1,4),([2,4,1],3,12)]:graphs.append((al*be*qu+offset*hu,build(L,R,cl,cr,U),False))
    factorchecks+=1
# Verify listed interval decompositions independently.
W=json.loads((P/'dependencies/interval_witnesses.json').read_text());factors=json.loads((P/'frozen/INTERVAL_FACTORS.json').read_text())
for w,f in zip(W,factors):
 a=w['rows'];v=w['root'];n=len(a)
 if f['dual']:a=[sum((a[j]>>i&1)<<j for j in range(n)) for i in range(n)]
 x=f['x'];L=f['L'];R=f['R'];A=[j for j in range(n) if a[v]>>j&1];B=[j for j in range(n) if a[j]>>v&1]
 if set(B)!=set(L) or set(A)-{x}!=set(R) or x not in A:raise RuntimeError('star partition')
 if any(not(a[i]>>j&1) for i in L for j in R):raise RuntimeError('cross cut')
 la=induced(a,L);ra=induced(a,R);cl=sum((a[x]>>i&1)<<j for j,i in enumerate(L));cr=sum((a[x]>>i&1)<<j for j,i in enumerate(R))
 h,al,ta,e=side(la,cl);hr=len(paths(ra)[-1]);be=beta(ra,cr);hrx=len(paths(induced(a,R+[x]))[-1])
 got={'alpha':al,'beta':be,'hL':h,'hR':hr,'hRx':hrx,'tau':ta,'d':h*hrx-hr*ta}
 if any(f[k]!=z for k,z in got.items()) or al*be!=36:raise RuntimeError(('factor record',got,f))
# Transitive closed forms including empty boundary.
transitive_controls=0
for n in range(8):
 a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 for cut in range(1<<n):
  h,al,ta,e=side(a,cut)
  x=sum(1<<(n-i-1) for i in range(n) if cut>>i&1)
  y=sum((1<<(j-i-1))+int(j==i+1) for i in range(n) for j in range(i+1,n) if cut>>i&1 and cut>>j&1)
  if (al,ta)!=(x,y):raise RuntimeError('binary formulas')
  transitive_controls+=1
(P/'GRAPHS.txt').write_text(''.join(f'{h} {int(strong)} {len(a)} '+ ' '.join(map(str,a))+'\n' for h,a,strong in graphs))
print(json.dumps({'degree_checks':num,'matching_histogram':hist,'factor_contexts':factorchecks,'literal_composites_prepared':len(graphs),'interval_factor_records':len(factors),'transitive_formula_controls':transitive_controls},indent=2))
