import sys,itertools,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent))
from inherited import direct_signature,strong,tables

def canon(r):
 n=len(r);return min(tuple(sum(1<<j for j in range(n) if r[p[i]]>>p[j]&1) for i in range(n)) for p in itertools.permutations(range(n)))
reps=[]
for n in (1,3,4,5):
 seen=set();pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for j,(u,v) in enumerate(pairs):
   if bits>>j&1:r[u]|=1<<v
   else:r[v]|=1<<u
  if not strong(r):continue
  c=canon(r)
  if c not in seen:seen.add(c);reps.append(c)
print('reps',reps,flush=True)
bank=[]
for i,r in enumerate(reps):
 for m in range(1<<len(r)):
  z=direct_signature(r,{j for j in range(len(r)) if m>>j&1})
  bank.append((i,m,z))
Path(__file__).with_name('bank.json').write_text(json.dumps({'rows':reps,'bank':bank}))

def mul(x,y):
 h,p,f,g,j,e,s=x;H,P,F,G,J,E,S=y
 return h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H

def search(F,cs):
 hits=[];n=0;seen=set()
 def rec(w,z):
  nonlocal n
  n+=1
  if z in seen:return
  seen.add(z)
  for i,m,v in bank:
   zz=v if z is None else mul(z,v)
   if zz[2]>F:continue
   ww=w+[(i,m)]
   if zz[2]==F and m and zz[0]-zz[4] in cs:hits.append((ww,zz))
   if zz[1]<=F:rec(ww,zz)
 rec([],None)
 print('search',F,cs,'states',n,'hits',hits[:30],len(hits),flush=True)
 return hits
for F,cs in [(63,[-35,-29]),(126,[-35,-29]),(189,[-109,-107]),(378,[-109,-107]),(567,[-269]),(1134,[-269])]:search(F,cs)
