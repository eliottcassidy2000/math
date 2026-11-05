from itertools import permutations,combinations
from functools import lru_cache
import json

def need(x,s):
 if not x:raise RuntimeError(s)
def graphs(n):
 E=list(combinations(range(n),2))
 for z in range(1<<len(E)):
  r=[0]*n
  for k,(a,b) in enumerate(E):r[a if z>>k&1 else b]|=1<<(b if z>>k&1 else a)
  yield tuple(r)
def path(r,p):return all(r[a]>>b&1 for a,b in zip(p,p[1:]))
@lru_cache(None)
def inv(r):
 n=len(r);h=0;e=[0]*n;s=[0]*n;R=[0]*n;K=[[0]*n for _ in r];C=0
 for p in permutations(range(n)):
  if path(r,p):
   h+=1;s[p[0]]+=1;e[p[-1]]+=1
   if n>=3 and p[0]==0 and r[p[-1]]&1:C+=1
  for cut in range(n):
   if path(r,p[:cut]) and path(r,p[cut:]):
    R[p[cut]]+=1
    if cut:K[p[cut-1]][p[cut]]+=1
 return h,e,s,R,K,C
def gadget(L):
 w=[] if L==1 else [1-int(b)for b in bin(L-1)[2:][::-1]];G=[0]*(len(w)+1)
 for i,b in enumerate(w,1):
  G[0 if b else i]|=1<<(i if b else 0)
  for j in range(i+1,len(G)):G[i]|=1<<j
 return tuple(G)
def join(T,entries,G=(0,)):
 n=len(T);r=list(T)+[x<<n for x in G]
 for a in range(n):
  for b in range(len(G)):r[n if b==0 and a in entries else a]|=1<<(a if b==0 and a in entries else n+b)
 return tuple(r)
def profile(r,a):
 h,e,s,R,K,C=inv(r);return h,e[a],R[a]
def permute(r,f):
 out=[0]*len(r)
 for a in range(len(r)):
  for b in range(len(r)):
   if r[a]>>b&1:out[f[a]]|=1<<f[b]
 return tuple(out)
from pathlib import Path
import sys
P=Path(__file__).parent
sys.path.insert(0,str(P/'source'))
from verify_weak_start_one import oracle
@lru_cache(None)
def coefficients(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=coefficients(L//2);return g,Q,2*B
 m=(L+1)//2;g,Q,B=coefficients(m);return g+m-m%2,Q+B,2*B
bound=101;queue=[(1,1,1,1)];seen=set(queue)
for h,e,R,d in queue:
 for L in range(1,bound//R+1):
  g,Q,B=coefficients(L);t=(g*h-(L%2)*e+L*R,h-int(L==1)*e,Q*h+B*R,h)
  if t[0]<=bound and t not in seen:seen.add(t);queue.append(t)
expected={};core_tests=0
for D,E,U,d in queue:
 if (D,E,U)==(1,0,2):continue
 need(set(oracle(D,E,U))=={d},'independent strong quadruple readout');core_tests+=1
 for k in range(1,5):expected.setdefault((D,(1<<k)*U),set()).add(d)
controls=0
for D in range(1,bound+1):
 for U in range(1,16*D+1):
  need(set(oracle(D,0,U))==expected.get((D,U),set()),'complete independent weak image');controls+=1
literal=0
for n in range(1,6):
 for T in graphs(n):
  h,e,starts,R,K,C=inv(T)
  for r in range(n):
   if starts[r]!=1:continue
   others=[v for v in range(n)if v!=r];d=sum(path(T,p)for p in permutations(others))
   need(d in oracle(h,e[r],R[r]),'literal source deletion response');literal+=1
large=[]
for k in (20,80,200):
 got=oracle(9,0,9*(1<<k));need(set(got)=={1,3} and all(v[1]==k for v in got.values()),'large exact dyadic ambiguity')
 need(oracle(9,0,9*(1<<k)+1)=={},'large inexact dyadic rejection')
 need(oracle(1,0,1<<k)=={1:((),k)},'large singleton tail')
 large.append(k)
need(oracle(1,0,1)=={}and oracle(1,1,1)=={1:((),0)}and oracle(1,1,2)=={},'singleton boundary')
print(json.dumps(dict(status='PASS',independent_forward_quadruples=len(queue),strong_core_tests=core_tests,weak_numeric_controls=controls,literal_start_one_roots=literal,large_tail_exponents=large,collision_deletions=sorted(oracle(9,0,18)),method='Independent forward quadruple closure plus literal permutation/cut inventories and exact dyadic edge controls'),sort_keys=True,indent=2))
