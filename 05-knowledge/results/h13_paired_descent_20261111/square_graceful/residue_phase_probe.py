#!/usr/bin/env python3
"""Focused controls for the second residue phase and ternary carries."""
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
spec=importlib.util.spec_from_file_location('compiler',ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
C3=[2,4,1];lo=c.signature(C3,3);hi=c.signature(C3,1)
def L(b):return 13*c.g(b)+32*b+6*(b&1)+18
def N(b):return 3*c.g(b)+12*b+10
def D(b):return L(b)-N(3*b)

def main():
 for b in range(1024):
  beta=b&1;G=[c.g(3*b+r) for r in range(3)]
  even=[G[0]+3*b+beta,G[0],G[1]+3*b+2-beta]
  odd=[G[1],G[2]+3*b+2+beta,G[2]]
  for r in range(3):
   assert c.g(6*b+r)==even[r]
   assert c.g(6*b+3+r)==odd[r]
  assert D(2*b)==D(b)+4*beta
 for b in range(128):
  for t in range(7):
   for s in range(1<<t):assert L(b)!=N((1<<t)*b+s)
 rays=[]
 for j in range(65):
  b=1002*(1<<j);cc=3*b;target=45090*(1<<j)-4415
  assert L(b)==N(cc)==target
  assert c.count(lo,0,2*cc+1)==target
  assert c.count(hi,0,4*cc+3)==target+2
  rays.append({'j':j,'address_bitlength':b.bit_length(),'target_bitlength':target.bit_length()})
 literal=[]
 for o,addr,expected in [(3,6013,40675),(1,12027,40677)]:
  rows=c.construct(C3,o,0,addr);actual=c.hp_subsets(rows)[-1]
  assert actual==expected
  literal.append({'cut':o,'right_address':addr,'order':len(rows),'HP':actual})
 state_spec=importlib.util.spec_from_file_location('state_probe',HERE/'probe.py')
 sp=importlib.util.module_from_spec(state_spec);state_spec.loader.exec_module(sp)
 rooted=[]
 for cc in range(8):
  bits=cc.bit_length();rev=int(f'{cc:0{bits}b}'[::-1],2) if bits else 0
  for o,addr,increment in [(3,2*cc+1,0),(1,4*cc+3,2)]:
   actual=sp.state(c.construct(C3,o,0,addr))
   expected=[N(cc)+increment,3,15*(1<<bits)-3*rev,3]
   assert actual==expected
   rooted.append({'C':cc,'cut':o,'state':actual})
 receipt={'status':'FINITE-EXACT controls of arbitrary-quantifier proofs',
          'signatures':{'lower':lo,'upper':hi},
          'seed':{'B':1002,'gB':c.g(1002),'C':3006,'gC':c.g(3006),'defect':D(1002)},
          'ternary_carry_universe':'0 <= B < 1024, all three carries, both input bits',
          'dyadic_no_go_control':'0 <= B < 128, 0 <= t <= 6, every suffix 0 <= s < 2^t',
          'ray_control':rays,'literal_witnesses':literal,'rooted_pair_controls':rooted,'residue_phase':[1,0],
          'all_address_coverage_claim':False}
 (HERE/'residue_phase_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'seed':receipt['seed'],'literal_witnesses':literal,'ray_exponents':len(rays)},indent=2))

if __name__=='__main__':main()
