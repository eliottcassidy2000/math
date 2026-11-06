import itertools,json

def graphs(n):
 es=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(es)):
  a=[0]*n
  for k,(i,j) in enumerate(es):
   u,v=(i,j) if mask>>k&1 else (j,i);a[u]|=1<<v
  yield a

def paths(a):
 s=[0]*len(a)
 for p in itertools.permutations(range(len(a))):
  if all(a[u]>>v&1 for u,v in zip(p,p[1:])):s[p[0]]+=1
 return s

def reach(a,i):
 z=1<<i
 while True:
  w=z
  for j in range(len(a)):
   if z>>j&1:w|=a[j]
  if w==z:return z
  z=w

def strong(a):return all(reach(a,i)==(1<<len(a))-1 for i in range(len(a)))
def sub(a,ids):return [sum(1<<j for j,v in enumerate(ids) if a[u]>>v&1) for u in ids]
def add(a,A):
 n=len(a);return [r|((not(A>>i&1))<<n) for i,r in enumerate(a)]+[A]
counts=dict(strong=0,subsets=0,mixed_subsets=0,pairs=0,scc_roots=0)
for n in range(1,6):
 for a in graphs(n):
  if not strong(a):continue
  counts['strong']+=1;s=paths(a)
  for mask in range(1,1<<n):
   m=mask.bit_count();assert sum(s[i] for i in range(n) if mask>>i&1)>=m*m//4+1;counts['subsets']+=1
  if n>=3:
   for root in range(n):
    ids=[i for i in range(n) if i!=root];b=sub(a,ids);rs=[reach(b,i) for i in range(n-1)];left=set(range(n-1));comps=[]
    while left:
     i=min(left);C={j for j in left if rs[i]>>j&1 and rs[j]>>i&1};left-=C;comps.append(sorted(C))
    allcomps=comps[:];comps.sort(key=lambda C:-sum(bool(b[C[0]]>>D[0]&1) for D in allcomps if C!=D))
    C=comps[0];cs=paths(sub(b,C));entry=[j for j,i in enumerate(C) if a[root]>>ids[i]&1];q=1
    for D in comps[1:]:q*=sum(paths(sub(b,D)))
    assert s[root]==q*sum(cs[j] for j in entry)
    assert entry and q<=s[root] and s[root]%q==0
    elig=sum(x<=s[root]//q for x in cs);assert elig<=4*(s[root]//q)-1
    counts['scc_roots']+=1
  if 3<=n<=4:
   for A in range(1,(1<<n)-1):
    b=add(a,A);assert strong(b);t=paths(b);p=sum(s[i] for i in range(n) if A>>i&1);assert t[-1]==p
    assert all(t[i]>=s[i] for i in range(n))
    for mask in range(1<<n):
     m=mask.bit_count();inc=p+sum(t[i]-s[i] for i in range(n) if mask>>i&1);assert inc>=(m+1)//2;counts['mixed_subsets']+=1
   b=[r|(1<<n) for r in a]+[1<<(n+1),(1<<n)-1];t=paths(b);assert strong(b) and t[-2:]==[sum(s)]*2;counts['pairs']+=1
print(json.dumps(dict(status='PASS',method='Independent permutation enumeration, literal subsets and SCCs',scope='All strong labeled inputs through n5; mixed and pair additions through n4',counts=counts),indent=2))
