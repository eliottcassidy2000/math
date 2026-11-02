"""Exact tests; no imports from the inherited repository; Python stdlib only.
Run: python3 -B audit.py; python3 -B -O audit.py.
"""
from explore import tournament, extend, H, strong, literal, endpoint_Q, endpoint_Q_dp, c3_blowup
from itertools import combinations
from functools import reduce
from math import gcd
from collections import Counter
import hashlib,json

ledger=[]
def check(x,label):
 if not x:raise RuntimeError(label)
def emit(row):
 ledger.append(row);print(json.dumps(row,sort_keys=True,separators=(',',':')))
def induced(a,vs):
 return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def currents(a):
 s,e,q=endpoint_Q_dp(a);n=len(a)
 b=[s[i]+e[i]+sum(q[j][i] for j in range(n)) for i in range(n)]
 beta=[[q[i][j]-q[j][i] for j in range(n)] for i in range(n)]
 return b,beta

def seams():
 bases=cases=literal_checks=0
 for m in range(2,5):
  for code in range(1<<(m*(m-1)//2)):
   a=tournament(m,code);bases+=1
   s,e,q=endpoint_Q(a)
   check((s,e,q)==endpoint_Q_dp(a),'independent exposed-word and DP boundary')
   b,B=currents(a);g=reduce(gcd,b+[B[i][j] for i in range(m) for j in range(i+1,m)])
   actualg=0
   for A in range(1<<m):
    u=extend(a,A);cap=[b[i]+sum(B[i][j] for j in range(m) if A>>j&1) for i in range(m)]
    for C in range(1<<m):
     d1=extend(u,C);d0=extend(u,C|(1<<m));h1=H(d1);h0=H(d0)
     predicted=sum(cap[i]*(((C>>i)&1)-((A>>i)&1)) for i in range(m))
     check(h1-h0==predicted,'seam formula');cases+=1
     w=extend(a,C)
     check(h1==H(extend(w,A|(1<<m))) and h0==H(extend(w,A)),'common-child reparenting')
     if m<=3 or (code==4 and A in [1,5] and C in [2,10]):
      check(literal(d1)==h1 and literal(d0)==h0,'literal child counts');literal_checks+=2
     if 0<A<(1<<m)-1 and 0<C<(1<<m)-1:actualg=gcd(actualg,predicted)
     if A&~C==0 and A!=C:
      lower=2*sum(H(induced(a,[j for j in range(m) if j!=i])) for i in range(m) if C>>i&1 and not(A>>i&1))
      check(predicted>=lower,'nested seam positivity floor')
   if m>=4:check(actualg==g,'exact seam lattice')
 emit({'test':'seam identities','bases':bases,'two-ear completions':cases,'literal child checks':literal_checks})

def clone_paths():
 for n in [4,5,6]:
  count=steps=0;maxsteps=0
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code)
   if not strong(a):continue
   count+=1;vs=list(range(n));removed=[]
   while len(vs)>4:
    v=next(v for v in vs if strong(induced(a,[w for w in vs if w!=v])))
    removed.append(v);vs.remove(v)
   pair=next((r,t) for r,t in combinations(vs,2) if all(bool(a[r]>>w&1)==bool(a[t]>>w&1) for w in vs if w not in [r,t]))
   r,t=pair;cur=a
   for v in removed:
    core=[w for w in range(n) if w!=v];check(strong(induced(cur,core)),'old strong common core')
    nxt=list(cur)
    for w in core:
     out=(w!=r and bool(cur[r]>>w&1))
     if out:nxt[v]|=1<<w;nxt[w]&=~(1<<v)
     else:nxt[v]&=~(1<<w);nxt[w]|=1<<v
    nxt=tuple(nxt);check(strong(nxt),'clone parent strong');check(induced(cur,core)==induced(nxt,core),'shared core literal')
    # The common core is strong; both deleted vertices have nonconstant attachment.
    check(0<sum((1<<i) for i,w in enumerate(core) if cur[v]>>w&1)<(1<<(n-1))-1,'old nonconstant ear')
    check(0<sum((1<<i) for i,w in enumerate(core) if nxt[v]>>w&1)<(1<<(n-1))-1,'new nonconstant ear')
    cur=nxt;steps+=1
   outs=[v for v in vs if v not in pair];module=[v for v in range(n) if v not in outs]
   check(H(induced(cur,module))==1,'transitive final module')
   check(all(bool(cur[x]>>w&1)==bool(cur[r]>>w&1) for x in module for w in outs),'homogeneous final module')
   check(H(cur)==1+(1<<(n-2)),'hub count')
   maxsteps=max(maxsteps,len(removed))
  emit({'test':'constructive parent overlap','order':n,'strong labelled parents':count,'replacement steps':steps,'maximum radius bound':maxsteps})

def unit_atlas():
 for n in [5,6]:
  pairs=list(combinations(range(n),2));idx={p:k for k,p in enumerate(pairs)};rows={}
  for code in range(1<<len(pairs)):
   a=tournament(n,code)
   if strong(a):rows[code]=(H(a),a)
  p={c:c for c in rows}
  def root(c):
   while p[c]!=c:p[c]=p[p[c]];c=p[c]
   return c
  counts=Counter();witness={}
  def edge(c,d,kind):
   if d not in rows or d<=c:return
   h,k=rows[c][0],rows[d][0]
   if abs(h-k)<=2:
    p[root(c)]=root(d);counts[kind]+=1
    witness.setdefault((min(h,k),max(h,k)),(c,d,kind))
  for c,(h,a) in rows.items():
   for r in range(n):
    if not strong(induced(a,[v for v in range(n) if v!=r])):continue
    es=[(v,idx[tuple(sorted((r,v)))]) for v in range(n) if v!=r]
    for v,k in es:edge(c,c^(1<<k),'flip')
    for (v,k),(w,l) in combinations(es,2):
     if ((a[r]>>v)^(a[r]>>w))&1:edge(c,c^(1<<k)^(1<<l),'exchange')
  components={}
  for c,(h,a) in rows.items():components.setdefault(root(c),[]).append((c,h))
  summary=sorted([{'vertices':len(v),'image':sorted({h for c,h in v})} for v in components.values()],key=lambda v:v['image'])
  expected=[list(range(9,16,2))] if n==5 else [[15,17,19],[23,25],[27,29,31,33],[37],[41,43,45]]
  check([s['image'] for s in summary]==expected,'full unit component images')
  if n==6:check((25,27) not in witness,'25/27 forbidden physical seam')
  emit({'test':'full labelled unit atlas','child order':n,'strong children':len(rows),'components':summary,'root-marked undirected edge counts':dict(counts)})

def missing_seam_repair():
 n=6;pairs=list(combinations(range(n),2));rows={}
 for code in range(1<<len(pairs)):
  a=tournament(n,code)
  if strong(a):rows[code]=(H(a),a)
 count=0;first=None
 for code,(h,a) in rows.items():
  if h!=25:continue
  for k in range(len(pairs)):
   target=code^(1<<k)
   check(target not in rows or rows[target][0]!=27,'no H25/H27 one-arc bridge')
  for k,l in combinations(range(len(pairs)),2):
   target=code^(1<<k)^(1<<l)
   if target not in rows or rows[target][0]!=27:continue
   common=set(pairs[k])&set(pairs[l]);check(len(common)==1,'every 25/27 two-flip bridge shares pivot')
   r=next(iter(common));v=next(v for v in pairs[k] if v!=r);w=next(v for v in pairs[l] if v!=r)
   check(not strong(induced(a,[v for v in range(n) if v!=r])),'every shortest bridge has nonstrong parent')
   check(((a[r]>>v)^(a[r]>>w))&1==0,'every shortest bridge changes two equal-membership bits')
   count+=1
   if first is None:first=(code,target,r)
 check(count==1440,'complete 25/27 distance-two count')
 base=(4,1,2,7,15);diamond={s:H(extend(base,s)) for s in [16,17,24,25]}
 check(diamond=={16:27,17:23,24:33,25:25},'nonstrong-parent cut diamond')
 check(all(strong(extend(base,s)) for s in diamond),'all diamond children strong')
 emit({'test':'missing scalar seam repair','labelled H25-to-H27 distance-two bridges':count,'all pivots':'nonstrong deletion parent, same-direction two-bit change','first child codes and pivot':first,'explicit base masks':base,'explicit base H':H(base),'strong-child diamond':diamond})

def divisibility():
 for k,quot in [(3,tournament(3,2)),(4,tournament(4,4))]:
  check(strong(quot),'strong quotient');a=c3_blowup(quot);n=len(a);b,B=currents(a)
  g=reduce(gcd,b+[B[i][j] for i in range(n) for j in range(i+1,n)])
  check(g%(2*3**(k-2))==0,'block-rotation seam divisor')
  row={'test':'C3 substitution seam obstruction','quotient order':k,'base order':n,'base masks':a,'H':H(a),'seam gcd':g,'forced divisor':2*3**(k-2)}
  if k==3:
   deltas=set();best=10**30;witness=None
   for A in range(1,(1<<n)-1):
    cap=[b[i]+sum(B[i][j] for j in range(n) if A>>j&1) for i in range(n)]
    ss=[0]*(1<<n)
    for C in range(1,1<<n):
     bit=C&-C;ss[C]=ss[C^bit]+cap[bit.bit_length()-1]
    for C in range(1,(1<<n)-1):
     d=ss[C]-ss[A];deltas.add(d)
     if d and abs(d)<best:best=abs(d);witness=(A,C,d)
   check(best==30 and g==6,'saturated seam lattice misses minimal step')
   A,C,d=witness;u=extend(a,A);h1=H(extend(u,C));h0=H(extend(u,C|(1<<n)))
   check(h1-h0==d,'literal large two-ear witness via child DP')
   row.update({'nonconstant attachment pairs':((1<<n)-2)**2,'distinct seam gaps':len(deltas),'minimum nonzero seam':best,'witness cuts':(A,C),'witness child counts x_to_y,y_to_x':(h1,h0)})
  emit(row)

seams();clone_paths();unit_atlas();missing_seam_repair();divisibility()
print('PASS semantic_sha256='+hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest())
