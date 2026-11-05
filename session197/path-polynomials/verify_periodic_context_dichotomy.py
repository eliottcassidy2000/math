#!/usr/bin/env python3
from verify_extension import *

def substitute(q,r,t):
 outside=[i for i in range(len(q)) if i!=r];n=len(outside);out=[0]*(n+len(t))
 for i,u in enumerate(outside):
  for j,v in enumerate(outside):
   if q[u]>>v&1:out[i]|=1<<j
  for j in range(len(t)):
   if q[u]>>r&1:out[i]|=1<<(n+j)
   else:out[n+j]|=1<<i
 for i,row in enumerate(t):out[n+i]|=row<<n
 validate(out);return tuple(out)

def induced(q,vs):return tuple(sum(((q[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs)
def active(q,r):return any((q[r]>>o&1) and(q[o]>>i&1) and(q[i]>>r&1) for o in range(len(q)) for i in range(len(q)))
def run():
 contexts=columns=passives=0
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):
   q=tournament(n,mask)
   for r in range(n):
    op=skeleton_operator(q,r);N=n-1;act=active(q,r);contexts+=1
    O=[i for i in range(n) if q[r]>>i&1];I=[i for i in range(n) if q[i]>>r&1]
    if not act:
     exterior=mul(bin_polynomial(induced(q,I)),bin_polynomial(induced(q,O)))
    for par,k in product((0,1),range(6)):
     f=basis(k,par);got=apply_skeleton(op,f)
     lower=mul(mon(N-2),D(f)) if act else mul(mon(N),f)
     cs=coordinates(add(got,scale(lower,-1)),(N+par)%2)
     require(all(c>=0 for c in cs),'positive sub-operator');columns+=1
     if not act:require(got==mul(exterior,f),'passive full-profile multiplication');passives+=1
 # Period transitions through all short D/X words, allowing odd parity changes.
 transitions=0
 for length in range(1,6):
  for letters in product('DX',repeat=length):
   if 'D' not in letters:continue
   for par in (0,1):
    for k in range(7):
     p=basis(k,par)
     for c in letters:p=D(p) if c=='D' else mul(mon(1),p)
     cs=coordinates(p,(par+letters.count('X'))%2)
     require(cs[k]>=(k+1)**2 and cs[k+1]>=1,'period stay/up')
     if k:require(cs[k-1]>=1,'period down')
     transitions+=1
 tri=(2,4,1);pad=(2,0);oddact=join(tri,(0,));passcyc=join(tri,(0,));ident=(0,)
 # In oddact choose a triangle root; in passcyc choose the exterior sink root.
 words=[[(tri,0)],[(oddact,0)],[(pad,0),(tri,1)],[(passcyc,3)],[(passcyc,3),(tri,2)],[(ident,0)],[]]
 controls=[]
 for seed in ((0,),transitive(2),tri):
  for word in words:
   t=seed;h=h_subsets(t)[-1];counts=[h];ops=[skeleton_operator(q,r) for q,r in word];p=bin_polynomial(seed)
   activeword=any(active(q,r) for q,r in word);c=1
   if not activeword:
    for q,r in word:c*=h_subsets(substitute(q,r,(0,)))[-1]
   for m in range(1,4):
    for (q,r),op in zip(word,ops):t=substitute(q,r,t);p=apply_skeleton(op,p)
    require(len(t)<=18,'bounded physical direct check')
    count=h_subsets(t)[-1];require(count==value(p,1),'literal context count')
    if activeword:
     for k in range(m//2+1):require(count>=h*(k+1)**(2*m-4*k),'quantitative growth path bound')
    else:require(count==h*c**m,'exact passive geometric family')
    counts.append(count)
   controls.append({'seed':list(seed),'word':[(list(q),r) for q,r in word],'active':activeword,'counts':counts})
 # Longer profile-only checks of the quantitative bound for both initial parities.
 longchecks=0
 for word in ('D','DX','XD','XXD','DXDX','DXXD','XXDDX'):
  for par in (1,2):
   p=mon(par)
   for m in range(1,21):
    for c in word:p=D(p) if c=='D' else mul(mon(1),p)
    h=value(p,1);require(h>=(1+m//4)**m,'all-period quantitative lower bound');longchecks+=1
 print(json.dumps({'status':'PASS','contexts':contexts,'virtual_columns':columns,'passive_columns':passives,'period_transition_controls':transitions,'long_lower_bounds':longchecks,'physical_context_words':controls},sort_keys=True,indent=2))
if __name__=='__main__':run()
