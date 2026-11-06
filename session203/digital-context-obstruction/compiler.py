#!/usr/bin/env python3
"""Complete start-two compiler with the sharper proved order cap; not an occupancy oracle."""
import argparse,heapq,itertools,json
from counting import literal,check_tournament
SEED=(1,1,0,1,2,1,4,2,3)
COST={'O':1,'Z':1,'J':1,'K1':2,'K2':3}
def base(p):return p[0]+p[3]-p[1]
def move(p,op):
 h,s,e,L,R,d,C,b,T=p
 if op=='O':return(h,s,0,L,2*R,d,2*C,2*b,3*T)
 if op=='Z':return(h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b)
 if op=='J':return(h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C)
 if op=='K1':return(h-e+3*R-2*d,2*s,R-d,C-2*L+2*T-b,2*h-e+4*R-d,h-e+R,4*C-6*L+10*T-4*b,2*C-2*L+2*T,3*C-2*L+6*T-b)
 if op=='K2':return(3*h-2*e+6*R-2*d,2*s,h+R-d,3*C-2*L+4*T-2*b,3*h-e+8*R-d,h+2*R,12*C-12*L+28*T-8*b,4*C-2*L+6*T,7*C-4*L+18*T-2*b)
 raise ValueError(op)
def inventory(bound,max_states=None,stop_at_target=False):
 if not isinstance(bound,int)or bound<1:raise ValueError('positive integer bound required')
 if max_states is not None and max_states<1:raise ValueError('state budget must be positive')
 cap=max(2,bound.bit_length()+1);serial=itertools.count()
 heap=[(2,next(serial),SEED)];seen={SEED:(2,())};witness={};processed=0
 while heap:
  n,_,p=heapq.heappop(heap)
  if n!=seen[p][0]:continue
  word=seen[p][1];processed+=1
  if p[1]==2:
   witness.setdefault(base(p),(word,p))
   if stop_at_target and base(p)==bound:
    return dict(complete=False,reason='TARGET_WITNESS_FOUND',bound=bound,order_cap=cap,states=len(seen),processed=processed),witness
  for op in ('O','Z','J','K1','K2')if p[1]==1 else('O','Z','J'):
   nn=n+COST[op]
   if nn>cap:continue
   t=move(p,op)
   if base(t)>bound:continue
   if t in seen and seen[t][0]<=nn:continue
   if t not in seen and max_states is not None and len(seen)>=max_states:
    return dict(complete=False,reason='RESOURCE_LIMIT',bound=bound,order_cap=cap,states=len(seen),processed=processed,state_limit=max_states),witness
   seen[t]=(nn,word+(op,));heapq.heappush(heap,(nn,next(serial),t))
 return dict(complete=True,bound=bound,order_cap=cap,states=len(seen),processed=processed,start_two_base_count=len(witness)),witness

def compile_target(target,max_states=1000000):
 if not isinstance(target,int)or target<1 or target%2==0:raise ValueError('positive odd B required')
 info,w=inventory(target,max_states,True)
 if target in w:
  word,p=w[target];rows=literal(word);check_tournament(rows)
  info.update(status='PHYSICAL_START_TWO_WITNESS',target_B=target,source_word=list(word),profile=list(p),root=0,rows=list(rows),order=len(rows))
 elif info['complete']:info.update(status='ABSENT_START_TWO_PROFILE',target_B=target,scope='Complete all-order absence only of start 2 at this B. No tournament-count exclusion and no conclusion about shifted F profiles.')
 else:info.update(status='INCOMPLETE',target_B=target,scope='The resource limit permits no absence conclusion.')
 return info
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('target',type=int);ap.add_argument('--max-states',type=int,default=1000000);ap.add_argument('--inventory',action='store_true');a=ap.parse_args()
 if a.target<1 or a.target%2==0:raise SystemExit('positive odd B required')
 if a.inventory:
  result,w=inventory(a.target,a.max_states);result['start_two_bases']=sorted(w)
  result['scope']='Exact all-order start-two B support through bound only if complete=true. No eventual-coverage conclusion.'
 else:result=compile_target(a.target,a.max_states)
 print(json.dumps(result,sort_keys=True,indent=2))
