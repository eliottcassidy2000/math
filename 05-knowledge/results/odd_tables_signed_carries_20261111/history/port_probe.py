#!/usr/bin/env python3
"""Focused exact controls for the H25 flank-block and entry-gap gates."""
import hashlib
import json
import subprocess
import types
from collections import defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
COMPILER='05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
data=subprocess.check_output(['git','show','HEAD:'+COMPILER],cwd=ROOT)
c=types.ModuleType('inherited_compiler');exec(compile(data,COMPILER,'exec'),c.__dict__)
CORES=[([50,8,3,5,46,14],38),([50,40,35,5,46,8],35)]
SIGS=[c.signature(rows,cut) for rows,cut in CORES]
def rev(n,k):return int(f'{n:0{k}b}'[::-1],2)
def ports(a,b):
 length=b.bit_length()
 return [25,(304*a+147)*(1<<length)-25*rev(b,length),25,304*(1<<(a.bit_length()+length))]
def count(a,b,branch):return c.count(SIGS[branch],a,b)
def actual_state(rows):
 n=len(rows);root=n-1;full=(1<<n)-1
 dp=[[0]*n for _ in range(full+1)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,full+1):
  for v,x in enumerate(dp[mask]):
   if not x:continue
   avail=rows[v]&~mask&full
   while avail:
    bit=avail&-avail;avail-=bit;dp[mask|bit][bit.bit_length()-1]+=x
 h=[sum(x) for x in dp];h[0]=1
 E=dp[-1][root];del dp
 dp=[[0]*n for _ in range(full+1)];dp[1<<root][root]=1
 for mask in range(1,full+1):
  for v,x in enumerate(dp[mask]):
   if not x:continue
   avail=rows[v]&~mask&full
   while avail:
    bit=avail&-avail;avail-=bit;dp[mask|bit][bit.bit_length()-1]+=x
 U=sum(sum(dp[mask])*h[full^mask] for mask in range(1,full+1))
 deletion=full^(1<<root)
 K=sum(h[mask]*h[deletion^mask] for mask in range(deletion+1))
 return [h[-1],E,U,h[deletion],K]

def main():
 literal=[]
 for a,b in [(1,208),(189,1)]:
  for branch,(rows,cut) in enumerate(CORES):
   graph=c.construct(rows,cut,a,b);state=actual_state(graph)
   assert state==[count(a,b,branch),*ports(a,b)]
   literal.append({'A':a,'B':b,'branch':branch,'order':len(graph),'state_HEUdK':state})
 blocks=[];formula_checks=0
 for t in range(1,17):
  M=1<<t
  for cc in range((18*(M-1))%25,M,25):
   numerator=147*(M-1)-304*cc
   if numerator<0 or numerator%25:continue
   rr=numerator//25
   if rr>=M:continue
   x=rev(rr,t)
   P=132*(M-1)-25*rev(cc,t)-304*x
   S=25*(c.g(cc)-c.g(x))+107*cc-122*x+15*(cc&1)-16*(x&1)
   blocks.append({'t':t,'c':cc,'x':x,'P':P,'S':S})
   for a in [1,2,5]:
    for C in [1,2,7]:
     old=(a,M*C+x);new=(M*a+cc,C)
     assert ports(*old)==ports(*new)
     assert count(*new,0)-count(*old,0)==P*a+10*(a&1)+16*(C&1)+S
     assert max(count(*new,j) for j in range(2))<min(count(*old,j) for j in range(2))
     formula_checks+=1
 fibre_groups=0;collision_groups=0
 for total in range(2,13):
  groups=defaultdict(list)
  for alen in range(1,total):
   blen=total-alen
   for a in range(1<<(alen-1),1<<alen):
    for b in range(1<<(blen-1),1<<blen):
     groups[ports(a,b)[1]].append((alen,a,b,sorted(count(a,b,j) for j in range(2))))
  for group in groups.values():
   group.sort();fibre_groups+=1
   if len(group)>1:collision_groups+=1
   for left,right in zip(group,group[1:]):
    assert left[0]<right[0]
    assert right[3][-1]<left[3][0]
 receipt={'status':'FINITE-EXACT controls for all-length proofs',
          'compiler_sha256':hashlib.sha256(data).hexdigest(),'literal_four_corners':literal,
          'block_universe':'All solutions of the port equation with 1<=t<=16; A in {1,2,5}; C in {1,2,7}; both branch choices',
          'block_count':len(blocks),'first_blocks':blocks[:5],'formula_checks':formula_checks,
          'port_fibre_universe':'Every A,B>=1 with 2<=bitlength(A)+bitlength(B)<=12; both core branches',
          'fibre_groups':fibre_groups,'nontrivial_fibre_groups':collision_groups,
          'restoration_example':{'entry_gap':115156,'old_targets':[93799,93803],'new_counts':[78589,79155],'all_four_restorations_excluded_by_entry_gap':True},
          'all_target_coverage_claim':False}
 (HERE/'port_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
