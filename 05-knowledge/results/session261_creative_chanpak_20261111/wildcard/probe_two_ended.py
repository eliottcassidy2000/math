from check_automaton import *
from itertools import product
from functools import lru_cache
A=[[False,True],[False,False]]
init,acc,tr=compile_automaton(2,A)
@lru_cache(None)
def value(w):
 v=init
 for c in w:v=step(v,tr[c])
 return sum(x for x,a in zip(v,acc) if a)
def op(w,c):return (c%4,)+w if c>=4 else w+(c,)
# A finite-depth success is only a candidate, never a uniform constructor.
for n in range(7):
 cand=[];one=0
 for w in product(range(4),repeat=n):
  h=value(w)
  if h<3:continue
  pp=[c for c in range(8) if value(op(w,c))==2*h+1]
  qq=[c for c in range(8) if value(op(w,c))==2*h-1]
  for p in pp:
   for q in qq:
    one+=1;front=[(w,h)];good=True
    for _ in range(4):
     nxt=[]
     for u,hh in front:
      for c,e in [(p,1),(q,-1)]:
       v=op(u,c)
       if value(v)!=2*hh+e:good=False;break
       nxt.append((v,value(v)))
      if not good:break
     if not good:break
     front=nxt
    if good:cand.append((w,h,p,q))
 print({'seed_length':n,'immediate_pairs':one,'depth4_candidates':cand[:10]},flush=True)
