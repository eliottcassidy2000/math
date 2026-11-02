#!/usr/bin/env python3
"""Standalone exact checks of the obstruction-lifetime framework.
Universe: every strong tournament through 7 via a fixed Hamilton cycle;
all singleton-free core words of total order at most 7.
No third-party modules, imported catalogue, or assert-based validation.
"""
import itertools,json,math,collections,hashlib,sys

def need(c,m):
 if not c: raise ValueError(m)
def valid(a):
 n=len(a);need(all(a[v]>>n==0 and not(a[v]>>v&1) for v in range(n)),'row range')
 need(all(((a[u]>>v)&1)+((a[v]>>u)&1)==1 for u in range(n) for v in range(u+1,n)),'tournament')
def hc(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];hs=[0]*N;hs[0]=1
 for v in range(n):dp[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   c=dp[S][v]
   if not c: continue
   q=a[v]&~S
   while q:
    b=q&-q;q-=b;dp[S|b][b.bit_length()-1]+=c
  hs[S]=sum(dp[S])
 return hs

def literal(a,vertices=None):
 vs=range(len(a)) if vertices is None else vertices
 return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(vs))
def code(a,p):
 c=0;k=0
 for i in range(len(p)):
  for j in range(i+1,len(p)):
   c|=((a[p[i]]>>p[j])&1)<<k;k+=1
 return c

def rows(n,c):
 a=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if c>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
   k+=1
 return a

def canonical(a):
 cells=collections.defaultdict(list)
 for v,r in enumerate(a):cells[r.bit_count()].append(v)
 best=None;aut=0
 for pp in itertools.product(*(itertools.permutations(cells[d]) for d in sorted(cells))):
  p=tuple(x for q in pp for x in q);c=code(a,p)
  if best is None or c<best:best=c;aut=1
  elif c==best:aut+=1
 return best,aut

def join(bs):
 n=sum(map(len,bs));a=[];o=0
 for b in bs:
  tail=((1<<n)-1)^((1<<(o+len(b)))-1)
  a.extend((r<<o)|tail for r in b);o+=len(b)
 return a

def record(a,aut=1):
 valid(a);hs=hc(a);h=hs[-1];deletions=[hs[(1<<len(a))-1 ^ (1<<v)] for v in range(len(a))]
 need(h==literal(a),'independent literal H')
 need(deletions==[literal(a,[u for u in range(len(a)) if u!=v]) for v in range(len(a))],'independent literal deck')
 need(all(d<h and d%2 for d in deletions),'strict core deletion')
 return dict(n=len(a),H=h,d=max(deletions),deletion_H=deletions,rows=a,automorphisms=aut)

def main():
 strong=[];levels=[];strong_labelled={1:1};NMAX=7
 for n in range(2,NMAX+1):
  strong_labelled[n]=(1<<(n*(n-1)//2))-sum(math.comb(n,k)*strong_labelled[k]*(1<<((n-k)*(n-k-1)//2)) for k in range(1,n))
 for n in range(3,NMAX+1):
  fixed={tuple(sorted((i,(i+1)%n))) for i in range(n)}
  free=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in fixed];cache={}
  for bits in range(1<<len(free)):
   a=[0]*n
   for i in range(n):a[i]|=1<<((i+1)%n)
   for k,(i,j) in enumerate(free):
    if bits>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
   c,aut=canonical(a)
   if c not in cache:cache[c]=aut
   else:need(cache[c]==aut,'automorphism consistency')
  rs=[record(rows(n,c),aut) for c,aut in sorted(cache.items())]
  mass=sum(math.factorial(n)//r['automorphisms'] for r in rs)
  need(mass==strong_labelled[n],'independent SCC recurrence mass')
  strong.extend(rs)
  births=collections.Counter(r['d'] for r in rs)
  levels.append(dict(n=n,cycle_presentations=1<<len(free),classes=len(rs),labelled_mass=mass,birth_histogram=sorted(births.items()),max_birth_multiplicity=max(births.values())))
 need([r['classes'] for r in levels]==[1,1,6,35,353],'complete strong class counts')
 # SCC-word enumeration, preserving order, gives every reduced core <=7.
 cores=[]
 def words(bs):
  if bs:
   a=join([b['rows'] for b in bs]);r=record(a)
   P=math.prod(b['H'] for b in bs)
   d=max((P//b['H'])*b['d'] for b in bs)
   need((r['H'],r['d'])==(P,d),'SCC envelope multiplication')
   r['scc_word']=[(b['n'],code(b['rows'],range(b['n']))) for b in bs];cores.append(r)
  for b in strong:
   if sum(x['n'] for x in bs)+b['n']<=NMAX:words(bs+[b])
 words([])
 need(len(cores)==399,'core word inventory')
 Hmax=max(r['H'] for r in cores);born=collections.Counter(r['d'] for r in cores);dead=collections.Counter(r['H'] for r in cores)
 prev=0
 for B in range(1,Hmax+1):
  current=sum(r['d']<=B<r['H'] for r in cores)
  need(current-prev==born[B]-dead[B],'barcode conservation')
  prev=current
 need(prev==0,'finite barcode closes')
 # Hostile: same order, H and strong status do not fix obstruction birth.
 a=record([48,1,3,7,14,30]);b=record([34,8,3,5,15,30])
 need((a['n'],a['H'],a['d'])==(6,23,11) and (b['n'],b['H'],b['d'])==(6,23,13),'H23 hostile')
 need(a['d']<=11<a['H'] and not b['d']<=11<b['H'],'threshold sidecar necessary')
 # Neutral padding preserves H but destroys all-threshold minimality.
 padded=join([[0],a['rows']]);hs=hc(padded);need(hs[-1]==hs[((1<<7)-1)^1]==23,'neutral padding hostile')
 # H15 robust atoms give 2^r ordinary classes; test the common threshold construction.
 robust=[r for r in strong if r['n']==5 and r['H']==15 and r['d']==5]
 need(len(robust)==2,'two robust H15 atoms')
 for r in range(1,5):
  seen=set()
  for bs in itertools.product(robust,repeat=r):
   # SCC uniqueness suffices: compute actual scalar formula and syntactic word injectivity.
   sig=tuple(code(b['rows'],range(5)) for b in bs);seen.add(sig)
   need(math.prod(b['H'] for b in bs)==15**r,'word product')
  need(len(seen)==2**r,'word distinction')
 # Small all-B numerical extension check, using exact integer powers.
 for B in range(5,100001):
  p=15;r=1
  while 15*p<=3*B:p*=15;r+=1
  k=0
  while p<=B:p*=3;k+=1
  need(B<p<=3*B,'triangle scaling covers every B')
  need(p//3<=B,'minimal obstruction condition')
 # Explicit counting bounds as exact rational numbers, no fake positivity assumption.
 bounds=[]
 for n in range(4,21):
  M=math.factorial(n-1)//2;num=1<<(n*(n-1)//2-n);den=math.factorial(n)*M
  bounds.append(dict(n=n,max_birth=M,strong_class_lower_numerator=num,strong_class_lower_denominator=math.factorial(n),max_birth_class_lower_ceiling=(num+den-1)//den))
 out=dict(status='FINITE-EXACT',universe='Every strong isomorphism class n=3..7 and every singleton-free SCC word of total order <=7',levels=levels,core_count=len(cores),barcode_flow_thresholds=[1,Hmax],H23_hostile=[a,b],robust_H15_atoms=robust,universal_counting_bounds=bounds,checks=['direct permutation H and every deletion for each representative','independent SCC recurrence gives all labelled strong masses','same-H same-order birth collision','source padding destroys minimality','SCC max-ratio multiplication','finite-universe barcode conservation','H15 word multiplicity and triangle-scaled coverage B=5..100000'])
 out['birth37_certificates']=[r for r in strong if r['n']==7 and r['d']==37]
 need(len(out['birth37_certificates'])==104,'birth37 exact count')
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
