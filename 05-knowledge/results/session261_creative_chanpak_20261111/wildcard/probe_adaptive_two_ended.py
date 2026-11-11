from check_automaton import *
from functools import lru_cache
from itertools import product
A=[[False,True],[False,False]]
init,acc,tr=compile_automaton(2,A)
@lru_cache(None)
def value(w):
 v=init
 for c in w:v=step(v,tr[c])
 return sum(x for x,a in zip(v,acc) if a)
def op(w,c):return (c%4,)+w if c>=4 else w+(c,)
@lru_cache(None)
def depth(w,d):
 if not d:return True
 h=value(w)
 for e in [-1,1]:
  if not any(value(op(w,c))==2*h+e and depth(op(w,c),d-1) for c in range(8)):return False
 return True
for n in range(7):
 bank=[w for w in product(range(4),repeat=n) if value(w)>=3]
 print({'length':n,'surviving_seeds_by_adaptive_depth':{d:sum(depth(w,d) for w in bank) for d in range(1,5)}},flush=True)
