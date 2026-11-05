"""Independent graph-DP verification plus bounded arithmetic coverage.
Source bank enumeration is input provenance, not an independent completeness proof.
Run from this directory: python3 independent_check.py [LOGX=20].
"""
import sys,json,math,itertools,time
from functools import lru_cache
from pathlib import Path
from probe_h127 import stats,apex,join,g
def require(c,m):
 if not c: raise RuntimeError(m)
rows=[list(map(int,s.split())) for s in open('source_bank7.graphs.txt')]
@lru_cache(None)
def triples(out):
 n=len(out);full=(1<<n)-1;hs,en=stats(out)
 rev=tuple(sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n))
 _,start=stats(rev)
 return tuple((hs[full],en[full][v],sum(hs[S]*start[full^S][v] for S in range(full+1) if not S>>v&1)) for v in range(n))
for n,r,h,e,R,*out in rows:
 require(len(out)==n,'order mismatch')
 require(all(not(out[i]>>i&1) for i in range(n)),'loop')
 require(all(((out[i]>>j&1)+(out[j]>>i&1))==1 for i in range(n) for j in range(i)),'not tournament')
 require(triples(tuple(out))[r]==(h,e,R),('bank triple',n,r,h,e,R))
bank=sorted({tuple(row[2:5]) for row in rows})
require(len(bank)==2732,'bank size')
# Every labelled tournament of order at most five: separate permutation counting.
perm_tests=0
for n in range(1,6):
 pairs=list(itertools.combinations(range(n),2))
 for code in range(1<<len(pairs)):
  out=[0]*n
  for k,(i,j) in enumerate(pairs):
   a,b=(i,j) if code>>k&1 else (j,i);out[a]|=1<<b
  hc=sum(all(out[p[k]]>>p[k+1]&1 for k in range(n-1)) for p in itertools.permutations(range(n)))
  require(stats(tuple(out))[0][-1]==hc,'permutation-v-DP')
  perm_tests+=1
# Independent one-apex table and direct graph checks.
for v in range(128):
 out,a=apex(v);require(triples(tuple(out))[a]==(g(v),int(v%2==0),v+1),'apex triple')
for u in range(16):
 S,r=apex(u)
 for v in range(16):
  N=(u+1)*(v+1)+g(u)*g(v)-(u%2==0 and v%2==0)
  require(stats(join(v,S,r))[0][-1]==N,'F1 direct')
logx=int(sys.argv[1]) if len(sys.argv)>1 else 20;X=1<<logx
G=[1]*(X+1)
for v in range(1,X+1):
 u=v//2;G[v]=G[u]+(u+(u&1) if v%2==0 else 0)
f1=bytearray(X+1);f0=bytearray(X+1);allcover=bytearray(X+1)
for v in range(X+1):
 if G[v]<=X:f0[G[v]]=1
for h,e,R in bank:
 if h<=X:f0[h]=1
for u in range(math.isqrt(X)+1):
 for v in range(u,X//(u+1)):
  N=(u+1)*(v+1)+G[u]*G[v]-int(u%2==0 and v%2==0)
  if N<=X:f1[N]=1
for v in range(X+1):
 if G[v]<=X: require(bool(f1[G[v]]),'plain one-apex subset F1')
for N in range(1,X+1,2):allcover[N]=f1[N] or f0[N]
residual_before=[N for N in range(1,X+1,2) if not allcover[N]]
witness={}
for h,e,R in bank:
 for v in range(X//R):
  N=h*G[v]-(e if v%2==0 else 0)+(v+1)*R
  if N<=X:
   if not allcover[N] and N%2: witness.setdefault(N,dict(v=v,bank=[h,e,R]))
   allcover[N]=1
miss=[N for N in range(1,X+1,2) if not allcover[N]]
require(miss==[7,21],('coverage',miss))
isprime=bytearray(b'\1')*(X+1);isprime[:2]=b'\0\0'
for p in range(2,math.isqrt(X)+1):
 if isprime[p]:isprime[p*p:X+1:p]=b'\0'*(((X-p*p)//p)+1)
primes=[p for p in range(3,X+1,2) if isprime[p] and p!=7]
sevenp=[7*p for p in range(5,X//7+1,2) if isprime[p]]
need_prime=[p for p in primes if not(f0[p] or f1[p])]
need_7p=[n for n in sevenp if not(f0[n] or f1[n])]
# Exact coefficient obstruction: odd v implies gcd(h,R)|N.
for N,w in witness.items():
 if w['v']%2:
  h,e,R=w['bank'];require(N%math.gcd(h,R)==0,'odd branch gcd')
report=dict(status='FINITE-EXACT bounded independent audit',X=X,graph_bank_rows=len(rows),distinct_graphs=len({tuple(row[5:]) for row in rows}),distinct_bank_triples=len(bank),bank_completeness='Source class generation, not independently proved here; every listed graph and tuple is independently validated.',permutation_controls=perm_tests,apex_graph_controls=128,F1_graph_controls=256,F0_apex_subset_F1_controls=X+1,all_misses=miss,F1_only_misses=sum(not f1[N] for N in range(1,X+1,2)),F0_F1_misses=len(residual_before),prime_atom_count=len(primes),prime_atoms_needing_F2=len(need_prime),seven_prime_atom_count=len(sevenp),seven_prime_atoms_needing_F2=len(need_7p),first_prime_witnesses={str(p):witness[p] for p in need_prime[:12]},first_7p_witnesses={str(n):witness[n] for n in need_7p[:12]},distinct_first_completion_triples=len({tuple(w['bank']) for w in witness.values()}),bank_gcd_hR={str(d):sum(math.gcd(h,R)==d for h,e,R in bank) for d in sorted({math.gcd(h,R) for h,e,R in bank})})
Path('independent_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
