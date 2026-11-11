#!/usr/bin/env python3
"""Reproduce eight source-rewrite languages; no completeness filter assumed."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
COMPILER=ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('compiler',COMPILER)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
S4=[2,4,9,3];S5b=[2,4,25,3,11]
def pow2(n):return n>0 and n&(n-1)==0
def ray(b,offset,minimum):return b-offset>=minimum and pow2(b-offset)
def source(b,upper):return 13*c.g(b)+32*b+6*(b&1)+18+2*upper
RULES=[
 ('U0',True,S4,8,0,lambda b:4*b+3,lambda b:b>=2 and pow2(b+2),(8,-4,6,-16)),
 ('U1',True,S4,11,0,lambda b:2*b+1,lambda b:b==0 or ray(b,0,2),(8,-8,6,-8)),
 ('U2',True,S4,4,1,lambda b:b,lambda b:b==4 or ray(b,4,16),(8,-7,5,-12)),
 ('L0',False,S4,0,1,lambda b:b,lambda b:b>=1 and pow2(b+1),(8,0,6,-14)),
 ('L1',False,S4,1,0,lambda b:4*b+2,lambda b:b in (0,3) or ray(b,3,8),(8,-2,6,-8)),
 ('L2',False,S4,4,0,lambda b:4*b+2,lambda b:b==2 or ray(b,2,8),(8,-6,6,-12)),
 ('L3',False,S4,2,1,lambda b:b,lambda b:b==6 or ray(b,6,16),(8,-5,5,-10)),
 ('L4',False,S5b,11,0,lambda b:b,lambda b:b>=74 and (b-2)%9==0 and pow2((b-2)//9),(4,-3,-1,10)),
]
def defect(coeff,b):a,d,e,f=coeff;return a*c.g(b)+d*b+e*(b&1)+f

def main():
 results={};covered=set();literal=[]
 for name,upper,rows,o,a,address,predicate,coeff in RULES:
  sig=c.signature(rows,o);hits=[]
  for b in range(1<<16):
   lhs=source(b,upper)-c.count(sig,a,address(b))
   assert lhs==defect(coeff,b),(name,b,'defect')
   assert (lhs==0)==predicate(b),(name,b,'iff')
   if lhs==0:hits.append(b);covered.add(b)
  # First three exact hits give positive actual-tournament controls, not just formulas.
  for b in hits[:3]:
   tournament=c.construct(rows,o,a,address(b));actual=c.hp_subsets(tournament)[-1]
   assert actual==source(b,upper),(name,b,'literal')
   literal.append({'rule':name,'B':b,'address':address(b),'order':len(tournament),'HP':actual})
  for q in range(256):
   aa,bb,cc,dd=coeff
   assert defect(coeff,2*q)==defect((aa,aa+2*bb,aa,dd),q)
   assert defect(coeff,2*q+1)==defect((aa,2*bb,0,dd+bb+cc),q)
  results[name]={'signature':sig,'defect':coeff,'hits':hits}
 receipt={'status':'FINITE-EXACT controls for arbitrary-word iff proofs',
          'universe':'All 0 <= B < 65536, all eight rules, no inherited target filters',
          'compiler_sha256':hashlib.sha256(COMPILER.read_bytes()).hexdigest(),
          'rules':results,'literal_witnesses':literal,'union_count':len(covered),
          'first_uncovered_addresses':[b for b in range(100) if b not in covered][:20],
          'hostile_controls':'Every address outside each stated equality language is checked to have nonzero defect; the union is explicitly not exhaustive.',
          'defect_transitions':'Both digits for every q in [0,255], all eight coefficient states',
          'complete_spectrum_claim':False}
 (HERE/'grammar_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({'rule_hit_counts':{k:len(v['hits']) for k,v in results.items()},'union_count':len(covered),'literal_controls':len(literal),'first_uncovered_addresses':receipt['first_uncovered_addresses']},indent=2))

if __name__=='__main__':main()
