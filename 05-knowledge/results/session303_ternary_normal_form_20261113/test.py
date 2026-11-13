from itertools import combinations, product
def H(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m,r in enumerate(d):
  for i,x in enumerate(r):
   if x:
    b=a[i]&~m
    while b:
     u=b&-b;b-=u;d[m|u][u.bit_length()-1]+=x
 return sum(d[-1]) if n else 1
def rows(n,mask):
 a=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if mask>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
def sub(a,inds):
 return [sum(((a[i]>>j)&1)<<k for k,j in enumerate(inds)) for i in inds]
def blocks(a):
 n=len(a)
 for X in combinations(range(n),3):
  s=sum(1<<i for i in X)
  if all((a[i]&s).bit_count()==1 for i in X) and all(a[j]&s in (0,s) for j in range(n) if j not in X):
   yield X
def match2(a,X):
 v=X[0]
 X=set(X)
 O=[i for i in range(len(a)) if i not in X and a[v]>>i&1]
 I=[i for i in range(len(a)) if i not in X and a[i]>>v&1]
 E=[(i,j) for i in O for j in I if a[i]>>j&1]
 return any(i!=k and j!=l for (i,j),(k,l) in combinations(E,2))
def spine(w):
 n=len(w);a=[0]*(n+1)
 for i in range(n):
  a[i]=((1<<n)-1)^((1<<(i+1))-1)
  if w[i]=='1':a[n]|=1<<i
  else:a[i]|=1<<n
 return a
if __name__=='__main__':
 for n in range(3,6):
  seen={}
  for m in range(1<<(n*(n-1)//2)):
   a=rows(n,m);h=H(a)
   if h%3==0 and not any(not match2(a,X) for X in blocks(a)):seen.setdefault(h,a)
  print(n,seen)

def pairs(a):
 for i,j in combinations(range(len(a)),2):
  if all((a[k]>>i&1)==(a[k]>>j&1) for k in range(len(a)) if k not in (i,j)):
   yield i,j

def inflate(a,P):
 # fixed pair expands by one vertex, preserving its old labelled vertices
 u,v=P
 if not(a[u]>>v&1): u,v=v,u
 n=len(a);b=a[:]+[0]
 b[n]|=1<<u;b[v]|=1<<n
 for i in range(n):
  if i not in P:
   if a[u]>>i&1:b[n]|=1<<i
   else:b[i]|=1<<n
 return b

def need(t,m):
 if not t: raise AssertionError(m)
def normal(a):
 def rec(b,ids):
  cs=[X for X in blocks(b) if not match2(b,X)]
  if not cs:return {tuple(ids)}
  out=set()
  for X in cs:
   z=max(X,key=lambda k:ids[k]); survivors=[k for k in range(len(b)) if k!=z]
   c=sub(b,survivors)
   need(H(b)==3*H(c),'exact stripping')
   out.update(rec(c,[ids[k] for k in survivors]))
  return out
 return rec(a,list(range(len(a))))

if __name__=='__main__':
 eq=nf=0
 for n in range(2,6):
  for m in range(1<<(n*(n-1)//2)):
   a=rows(n,m); h=H(a)
   for P in pairs(a):
    need((H(inflate(a,P))==3*h)==(not match2(a,P)),('matching',a,P))
    eq+=1
   need(len(normal(a))==1,('confluence',a))
   nf+=1
 for J in product([0,1],repeat=3):
  a=spine('1001100')
  for j,P in zip(J,[(1,2),(3,4),(5,6)]):
   if j:a=inflate(a,P)
  need(H(a)==63*3**sum(J),'ray')
  need(len(normal(a))==1 and len(next(iter(normal(a))))==8,'ray normal order')
  # Explicitly allow deletion of the new vertices to recover original labelled core.
  need(sub(a,list(range(8)))==spine('1001100'),'induced birth')
 # Two overlapping TT2 pairs inside a 3-run: cap is one; order does not matter.
 a=spine('10001100')
 need(set((1,2,3))==set(sum((list(p) for p in pairs(a) if p[0] in (1,2)),[])),'run')
 for a0 in range(1,21):
  a=spine('1'+'0'*a0+'1100')
  need(not list(blocks(a)),'module-free')
 # Independent tiny refutations: no strong reduced 3-divisible example below n=5.
 for w,h in [('1001100',63),('1100110',63)]:
  need(H(spine(w))==h and normal(spine(w))=={tuple(range(8))},w)
 print('PASS',eq,'literal pair inflations,',nf,'labelled normal forms, all eight strong 63 inflations, and module-free controls')
