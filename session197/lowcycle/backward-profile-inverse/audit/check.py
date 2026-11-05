from functools import lru_cache
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent

def need(v,s):
 if not v:raise RuntimeError(s)
@lru_cache(None)
def constants(L):
 m=(L-1).bit_length()
 if not m:return 1,1,1
 w=[1-((L-1)>>i&1)for i in range(m)]
 # Physical Hamilton words split the increasing base before/after the apex.
 g=int(w[0]==1)+int(w[-1]==0)
 for i in range(m):
  if w[i]:continue
  for j in range(m):
   if not w[j]:continue
   if j==i+1:g+=1
   elif i>j:g+=1<<(i-j-1)
 # Q reserves the apex as first vertex of its path; choose the first base vertex.
 Q=1+sum(1<<(m-i-1)for i in range(m)if w[i])
 return g,Q,1<<m

def step(old,L):
 h,e,R=old;g,Q,B=constants(L)
 return g*h-(L%2)*e+L*R,h-(L==1)*e,Q*h+B*R
@lru_cache(None)
def solve(h,e,R):
 if h==1:return ()if(e,R)==(1,1)else(1,)if(e,R)==(0,2)else None
 if h==3:return (2,)if(e,R)==(1,3)else None
 if h<3 or h%2==0 or not 1<=e<=h or not h<=R<=2*(h-e):return None
 # Deliberately try BOTH last-step types, without the new numerical branch theorem.
 old=(R-h+e,R-h,h-e)
 if 1<=old[0]<h:
  z=solve(*old)
  if z is not None:
   need(step(old,1)==(h,e,R),'singleton reconstruction');return z+(1,)
 if not 1<=e<h:return None
 for L in range(2,R//e):
  g,Q,B=constants(L);rem=R-Q*e
  if rem<0 or rem%B:continue
  R0=rem//B
  if R0<e:continue
  if L%2:ends=[g*e+L*R0-h]
  elif h==g*e+L*R0:ends=range(e+1)
  else:continue
  for e0 in ends:
   old=(e,e0,R0);z=solve(*old)
   if z is not None:
    need(step(old,L)==(h,e,R),'nonsingleton reconstruction');return z+(L,)
 return None

def physical(recipe):
 a=[0];r=0
 for L in recipe:
  m=(L-1).bit_length();g=[0]*(m+1)
  for i in range(m):
   w=1-((L-1)>>i&1)
   g[0 if w else i+1]|=1<<(i+1 if w else 0)
   for j in range(i+1,m):g[i+1]|=1<<(j+1)
  n=len(a);b=a+[row<<n for row in g]
  for i in range(n):
   for j in range(len(g)):
    if i==r and j==0:b[n]|=1<<r
    else:b[i]|=1<<(n+j)
  a=b;r=n
 return a,r

known={tuple(map(int,s.split(',')[:3])):int(s.split(',')[3])for s in (P/'physical-profiles.csv').read_text().splitlines()}
accepted={};trials=0
for h in range(1,102,2):
 for e in range(h+1):
  for R in range(1,2*h+1):
   t=(h,e,R);recipe=solve(*t);trials+=1
   need((recipe is not None)==(t in known),'complete physical bank versus generic inverse')
   if recipe is not None:accepted[t]=recipe
rows=[]
for t,recipe in sorted(accepted.items()):
 old=(1,1,1)
 for L in recipe:old=step(old,L)
 need(old==t,'profile replay')
 h,e,R=t;d=1 if h==1 else R-h+e if 3*e>h or 3*e==h and R>h else e
 need(known[t]==d,'literal deletion profile law')
 if h>3:need((recipe[-1]==1)==(3*e>h or 3*e==h and R>h),'last step necessity')
 a,r=physical(recipe);rows.append(','.join(map(str,[len(a),r,*t,d,*a]))+'\n')
(P/'recipes.csv').write_text(''.join(rows))
# The source's large returned recipes are checked with the independent physical
# constants, without claiming a direct adjacency count on the large graphs.
source=json.loads((P/'source-normal.json').read_text());large=[]
for c in source['large_full_profile_controls']:
 old=(1,1,1)
 for L in c['inverse_recipe']:old=step(old,L)
 need(list(old)==c['profile'],'large independent-constant replay')
 large.append(old)
print(json.dumps({'status':'PASS','tested_triples':trials,'physical_distinct_profiles':len(known),'last_step_and_deletion_checks':len(accepted),'generic_inverse_does_not_use_last_step_partition':True,'large_profile_replays':large,'recipes_sha256':hashlib.sha256((P/'recipes.csv').read_bytes()).hexdigest()},sort_keys=True,indent=2))
