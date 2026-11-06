#!/usr/bin/env python3
"""Exact target-bounded source-grammar compiler. No occupancy oracle."""
import argparse,json
from collections import deque
SEED=(1,1,0,1,2,1,4,2,3)
def move(p,op):
 h,s,e,L,R,d,C,b,T=p
 if op=='O':return(h,s,0,L,2*R,d,2*C,2*b,3*T)
 if op=='Z':return(h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b)
 if op=='J':return(h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C)
 raise ValueError(op)
def base(p):return p[0]+p[3]-p[1]
def literal(word):
 a=[2,0]
 for op in word:
  n=len(a);a=[x|(1<<n)for x in a]+[0]
  if op in'ZJ':a[0]^=1<<n;a[n]=1
  if op=='J':
   order=[n]+list(range(n));a=[sum(1<<j for j,v in enumerate(order)if a[u]>>v&1)for u in order]
 return a
def inventory(bound,max_states=None):
 if bound<1:raise ValueError('positive bound required')
 order_cap=bound.bit_length()+1;depth_cap=max(0,order_cap-2)
 seen={SEED:''};todo=deque([SEED]);witness={1:''};levels=[0]*(depth_cap+1)
 while todo:
  p=todo.popleft();word=seen[p];levels[len(word)]+=1;witness.setdefault(base(p),word)
  if len(word)==depth_cap:continue
  for op in'OZJ':
   t=move(p,op)
   if base(t)>bound or t in seen:continue
   if max_states is not None and len(seen)>=max_states:
    return {'complete':False,'reason':'RESOURCE_LIMIT','bound':bound,'states':len(seen),'state_limit':max_states,'scope':'No absence conclusion is authorized.'},witness
   seen[t]=word+op;todo.append(t)
 return {'complete':True,'bound':bound,'order_cap':order_cap,'word_length_cap':depth_cap,'states':len(seen),'level_counts':levels,'unit_start_base_count':len(witness),'scope':'Exact all-order unit-start B support through bound, using the inherited complete source grammar. This does not decide other starts or F(8,137).'},witness
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('target',type=int);ap.add_argument('--max-states',type=int,default=1000000);ap.add_argument('--inventory',action='store_true');args=ap.parse_args()
 if args.target<1 or args.target%2==0:raise SystemExit('target must be a positive odd B, not a general child H count')
 info,w=inventory(args.target,args.max_states)
 if args.inventory:info['unit_start_bases']=sorted(w)
 elif args.target in w:
  word=w[args.target];p=SEED
  for op in word:p=move(p,op)
  info.update(status='PHYSICAL_UNIT_START_WITNESS',target_B=args.target,source_word=word,profile=list(p),root=0,rows=literal(word))
 elif info['complete']:
  info.update(status='ABSENT_UNIT_START_PROFILE',target_B=args.target,qualification='This excludes only start1 at this B; it does not exclude a Hamilton count, other starts, or a staggered certificate.')
 else:info['status']='INCOMPLETE'
 print(json.dumps(info,sort_keys=True,indent=2))
