import itertools,json,sympy as sp
x=sp.symbols('x');checks=0
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def rank2(rows,n):
 rows=rows[:];r=0
 for j in range(n):
  k=next((i for i in range(r,len(rows)) if rows[i]>>j&1),None)
  if k is None:continue
  rows[r],rows[k]=rows[k],rows[r]
  for i in range(len(rows)):
   if i!=r and rows[i]>>j&1:rows[i]^=rows[r]
  r+=1
 return r
def compose(A,B):
 return [__import__('functools').reduce(int.__xor__,[A[j] for j in range(len(A)) if b>>j&1],0) for b in B]
def power(A,m):
 R=[1<<j for j in range(len(A))]
 while m:
  if m&1:R=compose(A,R)
  A=compose(A,A);m//=2
 return R
hostiles=[]
for f in [(x**3+x+1)**2*(x**3+x*x+1)**4,(x**3+x+1)**3*(x**3+x*x+1)**3]:
 p=sp.Poly(f,x,modulus=2);n=p.degree();low=sum((int(p.nth(j))%2)<<j for j in range(n));A=[1<<(j+1) for j in range(n-1)]+[low]
 fixes={m:2**(n-rank2([v^(1<<j) for j,v in enumerate(power(A,m))],n)) for m in [1,2,4,7,14,28]}
 cycles={m:sum(int(sp.mobius(m//d))*fixes[d] for d in sp.divisors(m))//m for m in fixes}
 req(cycles=={1:1,2:0,4:0,7:9,14:288,28:9216},'hostile graph counts')
 hostiles.append({'polynomial':str(p.as_expr()),'fixed_counts':fixes,'cycles':cycles})
def key(A,perm):return tuple(A[perm[i]][perm[j]] for i in range(len(A)) for j in range(i+1,len(A)))
def canon(A):return min(key(A,p) for p in itertools.permutations(range(len(A))))
classes={}
for n in range(5):
 C={}
 for bits in itertools.product([0,1],repeat=n*(n-1)//2):
  A=[[0]*n for _ in range(n)]
  for b,(i,j) in zip(bits,itertools.combinations(range(n),2)):A[i][j]=b;A[j][i]=1-b
  C.setdefault(canon(A),A)
 classes[n]=C
req([len(classes[n]) for n in range(5)]==[1,1,1,2,4],'canonical classes')
basis=[(n,k) for n in classes for k in sorted(classes[n])];idx={b:i for i,b in enumerate(basis)};cols=[]
for j,(n,k) in enumerate(basis):
 target=0
 if n<4:
  A=classes[n][k];B=[row+[0] for row in A]+[[1]*n+[0]];target=1<<idx[(n+1,canon(B))]
 cols.append((1<<j)^target)
seen=set();hist={}
for w in range(1<<len(basis)):
 if w in seen:continue
 v=w;size=0
 while v not in seen:
  seen.add(v);size+=1;v=__import__('functools').reduce(int.__xor__,[cols[j] for j in range(len(cols)) if v>>j&1],0)
 req(v==w,'dual permutation cycle');hist[size]=hist.get(size,0)+1
req(hist=={1:16,2:24,4:48,8:32},'tournament cutoff4')

color_checks=0
for w in range(1<<len(basis)):
 v=__import__('functools').reduce(int.__xor__,[cols[j] for j in range(len(cols)) if w>>j&1],0)
 if v==w:continue
 def color(z):
  j=next(j for j in range(len(basis)) if z>>j&1)
  target=cols[j]^(1<<j)
  req(target!=0,'nonfixed leading source exists')
  return (z & target).bit_count()%2
 req(color(w)!=color(v),'proper source-coordinate coloring')
 color_checks+=1
print(json.dumps({'status':'PASS','cutoff':4,'state_count':1<<len(basis),'directed_nonfixed_edges_checked':color_checks,'all_checks':checks},indent=2))
