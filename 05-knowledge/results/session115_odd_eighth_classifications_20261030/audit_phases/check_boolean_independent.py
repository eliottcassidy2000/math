from itertools import product
from collections import Counter
from math import isqrt
import json
checks=0
def gate(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
def fft(mask,n):
 a=[1-2*((mask>>i)&1) for i in range(n)];h=1
 while h<n:
  for k in range(0,n,2*h):
   for j in range(k,k+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
  h*=2
 return a
def linear_mask(u,d):return sum(((u&x).bit_count()%2)<<x for x in range(1<<d))
def normal(mask):
 for j in range(3):
  if mask>>(1<<j)&1:mask^=linear_mask(1<<j,3)
 return mask
classes={normal(mask) for mask in range(0,256,2)}
gate(len(classes)==16,'16 classes')
profiles=Counter()
for m in sorted(classes):
 W=fft(m,8);profile=tuple(sorted(Counter(W).items()));profiles[profile]+=1
 A=m>>3&1;B=m>>5&1;D=m>>6&1;C=D^(m>>7&1);k=A^B^C
 gate(k==m.bit_count()%2,'kappa')
 expected=Counter({8:1,0:7}) if m==0 else (Counter({4:3,-4:1,0:4}) if not k else (Counter({-6:1,2:7}) if (A,B,C,D)==(1,1,1,1) else Counter({6:1,2:4,-2:3})))
 gate(Counter(W)==expected,'profile')
 z=(A,B,C,D,1^B^C^D,1^A^C^D,1^A^B^D)
 primes=(19,2,59,29,5,89,7);n=1;sig=1
 for p,b in zip(primes,z):
  if b:n*=p;sig*=p+1
 gate(sig%30==0 and isqrt(sig//30)**2==sig//30,'arithmetic')
 if (A,B,C,D)==(1,1,1,1):gate((n,sig)==(65018,108000),'unique negative peak source')
GL=[]
for cols in product(range(1,8),repeat=3):
 perm=[]
 for x in range(8):
  y=0
  for j in range(3):
   if x>>j&1:y^=cols[j]
  perm.append(y)
 if len(set(perm))==8:GL.append(perm)
gate(len(GL)==168,'GL count')
orbits=set()
for m in classes:
 orb=frozenset(normal(sum(((m>>g[x])&1)<<x for x in range(8))) for g in GL)
 orbits.add(orb)
gate(sorted(map(len,orbits))==[1,1,7,7],'four orbits')
rows=[]
for d in (3,4):
 n=1<<d;full=(1<<n)-1
 defects={linear_mask(u,d)^(full if eps else 0)^(1<<v) for u in range(n) for eps in (0,1) for v in range(n)}
 gate(len(defects)==2*n*n,'unique defects')
 selected=0
 for m in range(1<<n):
  W=fft(m,n);isdef=max(map(abs,W))==n-2;selected+=isdef
  gate(isdef==(m in defects),'extremal iff one defect')
  if d==3:gate(isdef==bool(m.bit_count()%2),'odd exhaustion only d3')
 rows.append({'dimension':d,'all_functions':1<<n,'defects':len(defects),'extremal_functions':selected})
gate(max(map(abs,fft(7,16)))==10,'d4 odd weight hostile')
print(json.dumps({'status':'PASS','gates':checks,'normalized_classes':len(classes),'profiles':{str(k):v for k,v in profiles.items()},'relative_GL_orbits':sorted(map(len,orbits)),'exhaustive_boolean_rows':rows},indent=2))
