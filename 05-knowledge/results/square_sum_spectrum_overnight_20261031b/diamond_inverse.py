#!/usr/bin/env python3
"""Well-founded exact inverse membership and bounded target-entry search."""
import argparse
from functools import lru_cache
import json
from diamond_profiles import profiles,apply,single
from audit_feedback import require

BLOCKS=profiles()


def g_values(h):
 return tuple(k*(h//5) for k in range(6)) if h%5==0 else (0,h)


@lru_cache(None)
def unary_F(q):
 if q==0:return 1
 parent=q//2
 return unary_F(parent)+(parent+parent%2 if q%2==0 else 0)


@lru_cache(None)
def inverse(F,q,g,h):
 if min(F,h)<1 or q<0 or F%2!=1 or F<h or g not in g_values(h):return None
 # Actual marked cuts have nonnegative integral interior weight.
 if q<F-h+g or (q-F+h-g)%2:return None
 power=h
 while power%5==0:power//=5
 if power!=1:return None
 if h==1:
  if g==q%2 and F==unary_F(q):return tuple(('S',int(b)) for b in bin(q)[2:]) if q else ()
  return None
 candidates=[]
 if g==0 and q%2==0:
  previous_q=q//2
  for previous_g in g_values(h):candidates.append(((F-previous_q-previous_g,previous_q,previous_g,h),('S',0)))
 if g==h and q>=h and (q-h)%2==0:
  for previous_g in g_values(h):candidates.append(((F,(q-h)//2,previous_g,h),('S',1)))
 previous_h=h//5
 digit=(q*pow(previous_h,-1,32))%32
 for b in BLOCKS:
  if b['q']!=digit or g!=(5-b['t'])*previous_h:continue
  numerator=q-digit*previous_h
  if numerator<0 or numerator%32:continue
  previous_q=numerator//32
  for previous_g in g_values(previous_h):
   numerator=F-b['p']*previous_q-(5-b['s'])*previous_g-(b['F']-5)*previous_h
   if numerator%5==0:candidates.append(((numerator//5,previous_q,previous_g,previous_h),('D',b['mask'])))
 for state,letter in candidates:
  if state==(F,q,g,h):continue
  if state[0]<1 or state[1]<0:continue
  require(state[0]+state[1]+state[3]<F+q+h,'strict inverse rank')
  word=inverse(*state)
  if word is not None:return word+(letter,)
 return None


def target(F):
 require(isinstance(F,int) and F>=1,'target must be a positive integer')
 original=F;tail=0
 while True:
  h=1
  while h<=F:
   for q in range((32*F)//5+1):
    for g in g_values(h):
     if g==h:continue # no trailing p=0 generator in the normalized core
     if q<F-h+g or (q-F+h-g)%2:continue
     word=inverse(F,q,g,h)
     if word is not None:
      word=word+(('D',15),)*tail
      state=(1,0,0,1)
      for kind,a in word:state=single(state,a) if kind=='S' else apply(state,BLOCKS[a])
      require(state[0]==original,'reconstructed target count')
      return {'target':original,'word':word,'state':state,'vertices':2+sum(1 if k=='S' else 4 for k,a in word)}
   h*=5
  if F%5:break
  F//=5;tail+=1
 return {'target':original,'word':None,'scope':'absent from the declared singleton/strong-four-block grammar at every word length'}


def main():
 p=argparse.ArgumentParser();p.add_argument('targets',nargs='*',type=int);a=p.parse_args()
 targets=a.targets or [7,21,43,93,173,247,249,251,1001]
 # This finite profile inequality proves the canonical entry bound q<=32F/5.
 for b in BLOCKS:
  if b['p']:
   require(5*32<=32*b['p'] and 5*b['q']<=32*b['F'],'entry coefficient bound')
 results=[target(F) for F in targets]
 if not a.targets:
  require(all(r['word'] is None for r in results[:6]),'known exact grammar obstructions')
  require(all(r['word'] is not None for r in results[6:]),'positive inverse controls')
 print(json.dumps({'status':'passed','results':results,'memo_states':inverse.cache_info().currsize,
                   'rank':'F+q+h','canonical_q_bound':'floor(32F/5)',
                   'interior_guard':'q>=F-h+g and q-F+h-g even',
                   'scope':'Exact inverse certificates or exhaustion for the stated targets; universal acceptance above247 remains open.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
