#!/usr/bin/env python3
"""Exact controls for the common-defect and arbitrary-prefix actuator."""
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
spec=importlib.util.spec_from_file_location('compiler',ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
C3=[2,4,1];lo=c.signature(C3,0);hi=c.signature(C3,1)
def L(b):return 13*c.g(b)+32*b+6*(b&1)+18
def E(b):return 10*c.g(b)-4*b+6*(b&1)-12

def main():
 for b in range(4096):
  assert L(b)+2-c.count(hi,0,8*b+6)==E(b)
  assert L(b)-c.count(lo,0,16*b+12)==E(b)
 prefixes=[]
 for p in range(32):
  q=4*p+3;m=q-(5*c.g(q)-3)//2
  assert m>=2
  b=16**m*q+9*(16**m-1)//15
  assert E(b)==0
  assert L(b)+2==c.count(hi,0,8*b+6)
  assert L(b)==c.count(lo,0,16*b+12)
  prefixes.append({'p':p,'q':q,'m':m,'B_bitlength':b.bit_length(),'target_bitlength':L(b).bit_length()})
 literal=[]
 b=921
 for o,addr,expected in [(0,16*b+12,L(b)),(1,8*b+6,L(b)+2)]:
  rows=c.construct(C3,o,0,addr)
  actual=c.hp_subsets(rows)[-1]
  assert actual==expected
  literal.append({'B':b,'cut':o,'right_address':addr,'order':len(rows),'HP':actual})
 state_spec=importlib.util.spec_from_file_location('state_probe',HERE/'probe.py')
 sp=importlib.util.module_from_spec(state_spec);state_spec.loader.exec_module(sp)
 rooted=[]
 for b in range(8):
  bits=b.bit_length();rev=int(f'{b:0{bits}b}'[::-1],2) if bits else 0
  for o,addr,increment in [(0,16*b+12,0),(1,8*b+6,2)]:
   actual=sp.state(c.construct(C3,o,0,addr))
   expected=[3*c.g(b)+36*b+30+increment,3,39*(1<<bits)-3*rev,3]
   assert actual==expected
   rooted.append({'B':b,'cut':o,'state':actual})
 response_count=0;response_max=None;response_equality=[]
 for q in range(1,1<<16):
  m=q-(5*c.g(q)-3)//2
  if m<1:continue
  bits=q.bit_length();rev=int(f'{q:0{bits}b}'[::-1],2)
  delta=6*(1<<bits)-10*rev
  assert delta<=-6,(q,m,delta)
  response_count+=1
  response_max=delta if response_max is None else max(delta,response_max)
  if delta==-6:response_equality.append(q)
 receipt={'status':'FINITE-EXACT controls of arbitrary-prefix theorem',
          'defect_identity_universe':'0 <= B < 4096, both C3 replacements',
          'prefix_universe':'Every 0 <= p < 32; no filtered prefixes',
          'prefixes':prefixes,'literal_controls':literal,'rooted_pair_controls':rooted,
          'response_boundary_control':{'universe':'1 <= q < 65536, m(q)>=1','eligible_count':response_count,'maximum_delta':response_max,'equality_prefixes':response_equality},
          'original_target_preserved_during_padding':False,
          'both_padded_targets_preserved_by_rebase':True,
          'full_spectrum_claim':False}
 (HERE/'carry_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'prefixes':len(prefixes),'largest_B_bitlength':max(x['B_bitlength'] for x in prefixes),'literal_controls':literal},indent=2))

if __name__=='__main__':main()
