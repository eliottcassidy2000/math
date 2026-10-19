#!/usr/bin/env python3
"""Exact three-tree Gold cover for one full ten-vertex root core.
Standard library only. Reconstructs every two-chain atom and counts each
comparison event directly by ideal DP, not the gap/profile counting formula.
The tiny separation certificate is literal integer data below.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import product
import hashlib,json,sys
Q=(0,0,0,6,0,19,2,18,86,151)
ROOTS=(0,3);K=7;N=10;DEN=200
TREES=(((3,7),(0,2)),((5,8),(0,6)),((6,7),(5,8)))
# (tree number, cell number), zero-based; cells TT,TF,FT,FF.
SAFE=((0,0),(0,1),(1,1),(1,2),(1,3),(2,2),(2,3))
# (positive coefficient, tree, cell).
CERTS=(
 ((1,0,2),(1,2,0)),
 ((1,0,2),(1,1,0),(2,2,1)),
 ((1,0,3),(1,2,0)),
 ((1,0,3),(1,1,0),(2,2,1)),
)
checks=0
def require(ok,message):
 global checks
 checks+=1
 if not ok:raise ValueError(message)
def inc(p,a,b):return a!=b and not((p[a]>>b&1) or (p[b]>>a&1))
DELTA=tuple(K-sum(inc(Q,q,r) for r in range(N)) for q in range(N))
require(DELTA==(0,4,1,0,2,1,0,1,2,3),'slack vector')
def explicit_profiles(side):
 if side==0:
  return sorted((0,a,b,0,c,0,0,d,0,0) for a in range(5) for b in range(2) for c in range(3) for d in range(min(1,a,c)+1))
 return sorted((0,0,0,0,c,e,0,f,g,h) for c in range(2) for e in range(c,2) for f in range(c,2) for g in range(c,3) for h in range(f,4))
def general_profiles(side):
 eligible=[q for q in range(N) if q not in ROOTS and (not any(Q[q]>>r&1 for r in ROOTS) if side==0 else not any(Q[r]>>q&1 for r in ROOTS))]
 out=[]
 for v in product(*(range(DELTA[q]+1) for q in eligible)):
  t=[0]*N
  for q,z in zip(eligible,v):t[q]=z
  if sum(z>0 for z in t)>K:continue
  if any((t[a]<t[b] if side==0 else t[a]>t[b]) for a in range(N) for b in range(N) if Q[b]>>a&1):continue
  out.append(tuple(t))
 return sorted(out)
LP=explicit_profiles(0);UP=explicit_profiles(1)
require(LP==general_profiles(0),'complete lower profile formula')
require(UP==general_profiles(1),'complete upper profile formula')
require((len(LP),len(UP))==(46,48),'one-sided universe sizes')
def model(a,b):
 h=max(a);j=max(b);n=N+h+j;p=list(Q)+[0]*(h+j)
 for i in range(h):
  c=N+i
  for z in range(i):p[c]|=1<<(N+z)
  for q in range(N):
   if i<h-a[q]:p[q]|=1<<c
 for i in range(j):
  d=N+h+i
  for z in range(h+i):p[d]|=1<<(N+z)
  for q in range(N):
   if i>=b[q]:p[d]|=1<<q
 for _ in range(n):
  old=p.copy()
  for q in range(n):
   for r in range(n):
    if p[q]>>r&1:p[q]|=p[r]
  if p==old:break
 require(all(not(p[q]>>q&1) for q in range(n)),'acyclic model')
 require(tuple(z&((1<<N)-1) for z in p[:N])==Q,'induced core unchanged')
 require(max(sum(inc(p,q,r) for r in range(n)) for q in range(n))<=K,'seven-thin model')
 require({q for q in range(n) if q in ROOTS or any(inc(p,q,r) for r in ROOTS)}==set(range(N)),'full root core')
 x,y=ROOTS
 require(not(p[x]&~p[y]),'critical predecessor inclusion')
 require(all(not(p[q]>>y&1) or p[q]>>x&1 for q in range(n)),'critical successor inclusion')
 return tuple(p)
def count(p,events=()):
 p=list(p)
 for a,b in events:p[b]|=1<<a
 full=(1<<len(p))-1
 @lru_cache(None)
 def dp(s):
  if s==full:return 1
  return sum(dp(s|1<<q) for q,pred in enumerate(p) if not(s>>q&1) and not(pred&~s))
 return dp(0)
SINGLES=tuple(sorted(set(q for t in TREES for q in t)))
QUESTIONS=((),)+tuple((q,) for q in SINGLES)+tuple(tuple(sorted(t)) for t in TREES)
require(len(QUESTIONS)==len(set(QUESTIONS))==9,'exactly nine direct counts per atom')
min_safe=[None]*len(SAFE);min_certs=[None]*len(CERTS);witness_safe=[None]*len(SAFE);witness_certs=[None]*len(CERTS)
rows=[];atoms=0;digest=hashlib.sha256();min_order=100;max_order=0
for ia,a in enumerate(LP):
 for ib,b in enumerate(UP):
  if any(x+y>d for x,y,d in zip(a,b,DELTA)):continue
  p=model(a,b);counts={q:count(p,q) for q in QUESTIONS};E=counts[()]
  require(E>0,'nonempty extension set');cells=[];M=[]
  for first,second in TREES:
   tt=counts[tuple(sorted((first,second)))];f=counts[(first,)];s=counts[(second,)];cc=(tt,f-tt,s-tt,E-f-s+tt)
   require(min(cc)>=0 and sum(cc)==E,'four comparison cells partition')
   mm=tuple(E-cc[2*(r//2)]-cc[2*(r//2)+1]-cc[r] for r in range(4))
   cells.append(cc);M.append(mm)
  for z,(t,r) in enumerate(SAFE):
   v=M[t][r];require(DEN*v>E,'universally safe margin exceeds delta')
   f=Fraction(v,E)
   if min_safe[z] is None or f<min_safe[z]:min_safe[z]=f;witness_safe[z]=[ia,ib,E,v]
  for z,cert in enumerate(CERTS):
   c=sum(x for x,t,r in cert);v=sum(x*M[t][r] for x,t,r in cert)
   require(DEN*v>c*E,'positive certificate exceeds delta times coefficient sum')
   f=Fraction(v,c*E)
   if min_certs[z] is None or f<min_certs[z]:min_certs[z]=f;witness_certs[z]=[ia,ib,E,v]
  atoms+=1;min_order=min(min_order,len(p));max_order=max(max_order,len(p))
  line=json.dumps([list(a),list(b),E,cells],separators=(',',':'))+'\n';digest.update(line.encode())
  if '--table' in sys.argv:rows.append(line)
require(atoms==1764,'complete compatible pair count')
# Independently check the finite logical elimination: every possible chosen
# bad path from each tree is either SAFE or contains an entire CERT support.
patterns=0;safe_patterns=0;covered_patterns=0
for bad in product(range(4),repeat=3):
 patterns+=1
 if any(bad[t]==r for t,r in SAFE):safe_patterns+=1;continue
 require(any(all(bad[t]==r for c,t,r in cert) for cert in CERTS),'certificate eliminates failure pattern')
 covered_patterns+=1
require((patterns,safe_patterns,covered_patterns)==(64,60,4),'all three-tree failure patterns covered')
result={'status':'PASS','guards':checks,'core':Q,'roots':ROOTS,'slack':DELTA,'lower_profiles':len(LP),'upper_profiles':len(UP),'two_sided_atoms':atoms,'ambient_order_range':[min_order,max_order],'trees':TREES,'delta':f'1/{DEN}','safe_cells':[list(x) for x in SAFE],'safe_minima':list(map(str,min_safe)),'certificate_minima':list(map(str,min_certs)),'safe_witnesses':witness_safe,'certificate_witnesses':witness_certs,'certificate_max_coefficient':max(c for cert in CERTS for c,t,r in cert),'all_failure_patterns':patterns,'nonautomatic_failure_patterns':covered_patterns,'exact_count_table_sha256':digest.hexdigest()}
if '--table' in sys.argv:sys.stdout.write(''.join(rows))
else:print(json.dumps(result,indent=2))
