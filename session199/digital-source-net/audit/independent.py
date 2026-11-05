"""Independent exact source-net identities and all-corner discrepancy controls."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
from bisect import insort
from hashlib import sha256
import random,re,json
ROOT=Path(__file__).resolve().parent.parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def f(n):
 ans=1
 while n:
  r=n//2
  if not n%2:ans+=r+r%2
  n=r
 return ans

def identity(n):
 s=bin(n)[2:];k=len(s);B=1<<k
 C=B-1-int(s[::-1],2);runs=list(re.finditer('0+',s))
 correction=B*(1+len(runs))
 for run in runs:
  ell=run.end()-run.start();i=k-run.end();m=n% (1<<i)
  correction-=m*((B>>i)-(B>>(i+ell)))
 need(B*f(n)-n*C==correction,'exact digital identity')
 need(B<=correction<=B*(1+k//2),'correction range')
 need((correction==B)==('0'not in s),'strict equality boundary')
 need(2*correction<=(k+2)*B,'shell perturbation upper bound')
 return correction
checks=0
for n in range(1,200001):identity(n);checks+=1
rng=random.Random(199998)
for k in range(1,1025):
 n=(1<<(k-1))+rng.getrandbits(k-1);identity(n);checks+=1
# At each step the entire prefix count is known. Check both limits at every
# first-coordinate jump and both limits at every second-coordinate jump.
# These corners control the exact supremum over all real anchored thresholds.
cornerchecks=0;extrema=[]
for k in range(1,11):
 B=1<<k;N=B//2;ordered=[];worst=Q(0)
 for j in range(N+1):
  if j:insort(ordered,Q(f(B//2+j-1),B//2+j-1))
  widths=[j] if j==0 else [j-1,j]
  levels=sorted({Q(0),Q(1)}|{z for z in ordered if 0<=z<=1})
  cursor=0
  for y in levels:
   while cursor<len(ordered)and ordered[cursor]<y:cursor+=1
   lo=cursor;hi=cursor
   while hi<len(ordered)and ordered[hi]==y:hi+=1
   counts=[lo] if y==1 else [lo,hi]
   for width in widths:
    for c in counts:
     err=abs(c-width*y);worst=max(worst,err);cornerchecks+=1
     need(err<=2*k+4,'all critical anchored corners')
  # The reflected net baseline is independently tested on the same prefix.
 extrema.append({'k':k,'exact_anchored_supremum':str(worst),'bound':2*k+4})
# Every aligned first-coordinate bit cylinder gives a shifted uniform grid.
netcells=0
for k in range(1,11):
 B=1<<k
 for h in range(k+1):
  size=1<<(k-h)
  for prefix in range(1<<h):
   zs=sorted(B-1-int(f'{n:0{k}b}'[::-1],2)for n in range(prefix*size,(prefix+1)*size))
   if len(zs)>1:need(all(zs[i+1]-zs[i]==1<<h for i in range(len(zs)-1)),'exact shifted grid')
   need(zs[0]<(1<<h),'grid shift range')
   netcells+=1

def source_graph(U):
 rows=[{1},set()];root=0
 for bit in bin(U-1)[3:]:
  n=len(rows);rows.append(set())
  for i in range(n):
   if bit=='0'and i==root:rows[n].add(i)
   else:rows[i].add(n)
 return rows
literal=0
for U in range(2,65):
 rows=source_graph(U);H=ends=starts=deleted=0
 for p in permutations(range(len(rows))):
  if all(p[i+1]in rows[p[i]]for i in range(len(p)-1)):
   H+=1;ends+=(p[-1]==0);starts+=(p[0]==0)
 for p in permutations(range(1,len(rows))):
  deleted+=all(p[i+1]in rows[p[i]]for i in range(len(p)-1))
 need((H,ends,starts,deleted)==(f(U-1),U%2,1,1),'literal physical net source')
 literal+=1
report={'status':'PASS','proof_sha256':sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),
 'candidate_verifier_sha256':sha256((ROOT/'verify.py').read_bytes()).hexdigest(),
 'independent_verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
 'exact_identity_controls':checks,'all_real_anchor_corner_checks':cornerchecks,
 'exact_small_shell_suprema':extrema,'dyadic_net_cells':netcells,
 'literal_permutation_sources':literal,
 'scope':'Theorems1-2 PASS; the separate exploratory fixed-bank census is not certified by this audit; no exact coverage claim.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
