import itertools as it, json, sys, random
sys.path.insert(0,'frozen')
import compiler as c

def need(b,m):
 if not b: raise RuntimeError(m)

def profile(a):
 n=len(a);h=s=e=L=R=C=T=0
 def valid(p): return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
 # A permutation and one or two nondecreasing cuts uniquely represent labeled paths.
 for p in it.permutations(range(n)):
  if valid(p): h+=1;s+=p[0]==0;e+=p[-1]==0
  for i in range(n+1):
   if valid(p[:i]) and valid(p[i:]):
    C+=1
    if i: L+=p[i-1]==0;R+=p[0]==0
   if i and p[0]==0 and valid(p[:i]):
    for j in range(i,n+1):T+=valid(p[i:j]) and valid(p[j:])
 d=b=0
 for p in it.permutations(range(1,n)):
  d+=valid(p)
  for i in range(n): b+=valid(p[:i]) and valid(p[i:])
 return(h,s,e,L,R,d,C,b,T)

def child(a,op):
 n=len(a);edges={(u,v) for u in range(n) for v in range(n) if a[u]>>v&1}
 edges.update((u,n) for u in range(n) if u!=0 or op=='O')
 if op!='O':edges.add((n,0))
 order=list(range(n+1)) if op!='J' else [n]+list(range(n))
 return [sum(1<<j for j,v in enumerate(order) if (u,v) in edges) for u in order]

def graph(n,code):
 a=[0]*n
 for k,(u,v) in enumerate(it.combinations(range(n),2)):
  if code>>k&1:a[u]|=1<<v
  else:a[v]|=1<<u
 return a

checks=0;counts={};rng=random.Random(203)
for n in range(1,7):
 codes=range(1<<(n*(n-1)//2)) if n<=4 else sorted(rng.sample(range(1<<(n*(n-1)//2)),24 if n==5 else 6))
 counts[n]=len(codes)
 for code in codes:
  a=graph(n,code);p=profile(a)
  for op in 'OZJ':
   actual=profile(child(a,op));expected=c.move(p,op)
   if n==1 and op=='Z': expected=expected[:1]+(0,)+expected[2:]
   need(actual==expected,('transfer',n,code,op,p,actual,expected));checks+=1
# All words, no state merging or B pruning. Independent of BFS control flow.
support={1};front=[c.SEED];wordcounts=[1]
for depth in range(1,10):
 front=[c.move(p,op) for p in front for op in 'OZJ'];wordcounts.append(len(front))
 support.update(c.base(p) for p in front if c.base(p)<=997)
info,w=c.inventory(997)
need(info['complete'] and support==set(w) and len(support)==330,'unpruned support')
# The 139 literal recipe checked with the permutation counter, not DP.
r,w139=c.inventory(139);word=w139[139];p=profile(c.literal(word));need(p[1]==1 and c.base(p)==139,'139 literal')
# Resource truncation cannot claim completeness, including zero and negative limits.
guards=[]
for limit in [-1,0,1,2,10]:
 r,ww=c.inventory(137,limit);need(not r['complete'] and r['reason']=='RESOURCE_LIMIT','guard');guards.append(r)
print(json.dumps(dict(status='PASS',method='Permutation enumeration with nondecreasing cuts; no subset DP',input_counts=counts,child_checks=checks,unpruned_word_counts_by_depth=wordcounts,support997=sorted(support),support997_count=len(support),target137_absent=137 not in support,target139_word=word,target139_profile=p,resource_checks=guards),indent=2,sort_keys=True))
