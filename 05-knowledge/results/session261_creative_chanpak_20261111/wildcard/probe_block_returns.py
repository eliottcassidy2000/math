from check_automaton import *
from itertools import product
A=[[False,True],[False,False]]
init,acc,tr=compile_automaton(2,A)
def out(v):return sum(x for x,a in zip(v,acc) if a)
def move(v,b):
 for c in b:v=step(v,tr[c])
 return v
blocks=[b for n in range(1,4) for b in product(range(4),repeat=n)]
wave=[(init,())];candidates=0;survivors=[]
for n in range(7):
 for v,w in wave:
  h=out(v)
  if h<3:continue
  plus=[];minus=[]
  for b in blocks:
   hh=out(move(v,b))
   if hh==2*h+1:plus.append(b)
   if hh==2*h-1:minus.append(b)
  for p in plus:
   for q in minus:
    candidates+=1; frontier=[(v,h)];good=True
    for depth in range(4):
     nxt=[]
     for vv,hh in frontier:
      for b,eps in [(p,1),(q,-1)]:
       vvv=move(vv,b);hhh=out(vvv)
       if hhh!=2*hh+eps:good=False;break
       nxt.append((vvv,hhh))
      if not good:break
     if not good:break
     frontier=nxt
    if good:survivors.append((w,h,p,q))
 print({'seed_length':n,'seeds':len(wave),'immediate_pairs_cumulative':candidates,'depth4_survivors':survivors[:10]},flush=True)
 if n<6:wave=[(step(v,tr[c]),w+(c,)) for v,w in wave for c in range(4)]
