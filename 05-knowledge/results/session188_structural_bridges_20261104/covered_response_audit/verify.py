"""Independent finite coefficient audit for equal-profile cyclic quotient query."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
from pathlib import Path
import json,hashlib,sys
# This helper belongs to a separately written packet, not the producer under audit.
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'covered_quotient'))
from verify import cycle,rows,hp,extend,profile,need


def coefficients(m,converse=False):
 # States track actual consecutive run words, with x/y used at most once.
 # Old block runs are distinguishable by position and later filled by ordered covers.
 R=[0]*5
 for i in range(3):R[i]|=1<<((i+(-1 if converse else 1))%3)
 # a=whole block 0, c=whole block 1; x=3 -> y=4.
 R[3]|=1<<4
 for i in range(3):
  if i==0:R[i]|=1<<3
  else:R[3]|=1<<i
  if i==1:R[i]|=1<<4
  else:R[4]|=1<<i
 D=defaultdict(int)
 for i in range(5):
  c=[0,0,0];p=0
  if i<3:c[i]=1
  else:p=1<<(i-3)
  D[tuple(c),p,i]+=1
 out=defaultdict(int)
 for length in range(1,3*m+3):
  N=defaultdict(int)
  for (c,p,u),v in D.items():
   if p==3 and min(c)>0:out[c]+=v
   for w in range(5):
    if not R[u]>>w&1:continue
    a=list(c);b=p
    if w<3:
     if a[w]==m:continue
     a[w]+=1
    else:
     bit=1<<(w-3)
     if p&bit:continue
     b|=bit
    N[tuple(a),b,w]+=v
  D=N
 return out


def equal_profile_identity(m):
 A=coefficients(m);B=coefficients(m,True);left=defaultdict(int)
 for c in set(A)|set(B):left[tuple(sorted(c))]+=A[c]-B[c]
 right=defaultdict(int)
 def add(c,v):
  if min(c)>0 and max(c)<=m:right[tuple(sorted(c))]+=v
 for k in range(m):
  z=k+1
  add((k+2,k+1,k+1),-z);add((k+1,k+1,k+1),-2*z)
  add((k,k+2,k+2),z);add((k,k+1,k+2),2*z)
 need({c:v for c,v in left.items() if v}=={c:v for c,v in right.items() if v},'equal-profile coefficient identity')
 return {'max_runs_per_block':m,'nonzero_cubic_profile_monomials':sum(v!=0 for v in left.values())}


def gap(q):
 q=list(q)+[0,0]
 return -sum((k+1)*(q[k+2]+2*q[k+1])*(q[k+1]**2-q[k]*q[k+2]) for k in range(len(q)-2))


def full_table_checks():
 root=Path(__file__).resolve().parent.parent/'covered_response'
 raw=json.loads((root/'full_response.json').read_text())
 need(raw['basis_coordinates']==2116,'complete named basis')
 table={(r[0],r[1]):r[2:] for r in raw['coordinate_data']}
 need(len(table)==2116,'distinct retained coordinates')
 n=9;checks=0;hist=defaultdict(int)
 for (a,c),row in table.items():
  As=[u for u in range(1<<n) if u&a==u];Cs=[u for u in range(1<<n) if u&c==u]
  for source in range(3):
   z=sum((-1)**((a.bit_count()-u.bit_count())+(c.bit_count()-v.bit_count()))*table[u,v][source] for u in As for v in Cs)
   need(z==row[source+3],'retained Mobius coefficient');checks+=1
  d=a.bit_count()+c.bit_count()
  if d==4:need(row[3]==row[5],'whole-converse quartic control')
  if row[3]!=row[4]:hist[d]+=1
 literals=[]
 for a,c in ((7,56),(341,170),(511,511),(17,258)):
  values=[hp(extend(T,a,c)) for T in raw['rows']]
  for j,value in enumerate(values):
   reconstructed=sum(row[j+3] for (u,v),row in table.items() if u&a==u and v&c==v)
   need(value==reconstructed,'independent literal full response')
  literals.append([a,c]+values)
 need(hist=={2:54,3:432,4:540},'complete difference degrees')
 return {'coefficient_checks':checks,'difference_nonzero_degrees':dict(hist),'literal_named_queries':literals,
         'source_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py','full.cpp','full_response.json')}}


def main():
 symbolic=[equal_profile_identity(m) for m in range(1,8)]
 cases=[]
 # Identical ordered-cover profiles, but distinct physical orientations allowed.
 for B in ([0],[2,0],[6,4,0],[2,4,1]):
  C=[((1<<len(B))-1)^r^(1<<i) for i,r in enumerate(B)]
  blocks=[B,C,B];T,offsets=cycle(blocks);U,_=cycle(blocks,True);a=(1<<len(B))-1;c=a<<len(B)
  x=hp(extend(T,a,c));y=hp(extend(U,a,c));q=profile(B)
  need(x-y==gap(q),'literal actual equal-profile query')
  cases.append({'children':blocks,'profile':q,'mask_a':a,'mask_c':c,'T_count':x,'U_count':y,'difference':x-y})
 root=Path(__file__).resolve().parent;source=root.parent/'covered_response'/'PROOF.md'
 out={'status':'PASS','independent_run_coefficient_checks':symbolic,'literal_cases':cases,'retained_full_table':full_table_checks(),'independent_helper_sha256':hashlib.sha256((root.parent/'covered_quotient'/'verify.py').read_bytes()).hexdigest(),'method':'No producer imports. Explicit bounded run words, collected as independent cubic profile variables, and actual full named attachments.','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 if source.exists():out['audited_proof_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
