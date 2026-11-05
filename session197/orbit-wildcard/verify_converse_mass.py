from verify import *
from collections import Counter

def key(r,p=None):
 n=len(r);p=p or tuple(range(n));v=0
 for k,(i,j) in enumerate(it.combinations(range(n),2)):
  if r[p[i]]>>p[j]&1:v|=1<<k
 return v
def canon(r):return min(key(r,p) for p in it.permutations(range(len(r))))
def norm_universe(n):
 edges=list(it.combinations(range(n),2));free=[e for e in edges if e[1]>e[0]+1];orbit={};groups=[]
 for e in free:
  if e in orbit:continue
  f=(n-1-e[1],n-1-e[0]);g={e,f};groups.append(g)
  for x in g:orbit[x]=len(groups)-1
 gate(len(groups)==(n-1)**2//4,'normalized orbit exponent')
 for mask in range(1<<len(groups)):
  r=[0]*n
  for i,j in edges:
   win=j==i+1 or (mask>>orbit[i,j]&1 if j>i+1 else True)
   if win:r[i]|=1<<j
   else:r[j]|=1<<i
  yield r

def class_profile(r):
 n=len(r);base=key(r);opposite=((1<<(n*(n-1)//2))-1)^base
 ps=list(it.permutations(range(n)));G=[p for p in ps if key(r,p)==base]
 antis=[p for p in ps if key(r,p)==opposite and all(p[p[i]]==i for i in range(n))]
 gate(antis,'self converse anti-involution')
 words=[p for p in ps if all(r[a]>>b&1 for a,b in zip(p,p[1:]))]
 values=[]
 for tau in antis:
  cent=[g for g in G if all(g[tau[i]]==tau[g[i]] for i in range(n))]
  fixed=[w for w in words if w==tuple(tau[x] for x in reversed(w))]
  gate(len(fixed)%len(cent)==0 and len(fixed)%2==1,'centralizer free odd fibre')
  values.append((len(fixed)//len(cent),len(fixed),len(cent)))
 gate(len(set(values))==1,'anti-involution conjugacy data')
 return values[0][0],len(G),len(words),len(antis)

receipt=[];individual=0
for n in range(0,7):
 freq=Counter();reps={}
 for r in norm_universe(n):
  validate(r);gate(all(r[i]>>(i+1)&1 for i in range(n-1)),'forced Hamilton word')
  rev=tuple(reversed(range(n)))
  gate(all(((r[i]>>j)&1)==((r[rev[j]]>>rev[i])&1) for i in range(n) for j in range(n)),'reflection anti')
  c=canon(r);freq[c]+=1;reps[c]=r
 for c,r in reps.items():
  fibre,ag,h,anti_count=class_profile(r)
  gate(freq[c]==fibre and fibre%2==1,'individual normalized fibre')
  gate(h%ag==0 and (h//ag)%2==1,'ordinary free action');individual+=1
 if n<=5:
  full=set()
  for r in graphs(n):
   c=canon(r);F=(1<<n)-1;rev=[F^(1<<i)^r[i] for i in range(n)]
   if c==canon(rev):full.add(c)
  gate(full==set(reps),'independent complete self-converse class inventory')
 gate(sum(freq.values())==1<<((n-1)**2//4),'mass total')
 if n>=3:gate(len(reps)%2==0,'class parity')
 receipt.append({'n':n,'normalized_presentations':sum(freq.values()),'self_converse_classes':len(reps),'fibre_histogram':dict(sorted(Counter(freq.values()).items()))})
print(json.dumps({'status':'PASS','individual_class_fibres':individual,'orders':receipt},sort_keys=True,indent=2))
