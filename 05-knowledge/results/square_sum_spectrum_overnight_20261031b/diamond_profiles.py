#!/usr/bin/env python3
"""Literal four-vertex block profiles and reversible (h,q) history."""
from itertools import permutations,product
from functools import lru_cache
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'square_sum_spectrum_overnight_20261030'))
from audit_feedback import subset_dp,require
ROWS=(8,1,3,6)
def ispath(p,rows):return all(rows[a]>>b&1 for a,b in zip(p,p[1:]))
def profiles():
 n=4;paths=[p for p in permutations(range(n)) if ispath(p,ROWS)]
 covers=[(p[:k],p[k:]) for p in permutations(range(n)) for k in range(n+1) if ispath(p[:k],ROWS) and ispath(p[k:],ROWS)]
 result=[]
 for mask in range(16):
  s=sum(mask>>p[0]&1 for p in paths);t=sum(not(mask>>p[-1]&1) for p in paths)
  p=sum(bool(P) and not(mask>>P[-1]&1) for P,Q in covers)
  q=sum(bool(Q) and (mask>>Q[0]&1) for P,Q in covers)
  rows=list(ROWS)+[mask]
  for i in range(n):
   if not(mask>>i&1):rows[i]|=1<<n
  F=subset_dp(rows)
  result.append({'mask':mask,'H':len(paths),'E':len(covers),'s':s,'t':t,'p':p,'q':q,'F':F})
 return result

def apply(state,b):
 F,q,g,h=state;H,E=b['H'],b['E']
 return (H*F+b['p']*q+(H-b['s'])*g+(b['F']-H)*h,
         E*q+b['q']*h,(H-b['t'])*h,H*h)

def single(state,b):
 F,q,g,h=state
 return (F,2*q+h,h,h) if b else (F+q+g,2*q,0,h)

def actual(word):
 # Fixed source: parent vertex0 beats the marked external vertex1.
 rows=[2,0];x=1
 for kind,cut in word:
  block=ROWS if kind=='D' else (0,);n=len(rows);size=len(block)
  rows.extend([0]*size)
  for old in range(n):
   if old!=x:
    for j in range(size):rows[old]|=1<<(n+j)
  for i in range(size):
   for j in range(size):
    if block[i]>>j&1:rows[n+i]|=1<<(n+j)
   if cut>>i&1:rows[x]|=1<<(n+i)
   else:rows[n+i]|=1<<x
 parents=[v for v in range(len(rows)) if v!=x]
 paths=[p for p in permutations(parents) if ispath(p,rows)]
 h=len(paths);t=sum(rows[p[-1]]>>x&1 for p in paths)
 q=0
 for p in permutations(parents):
  for k in range(len(p)):
   P,Q=p[:k],p[k:]
   if ispath(P,rows) and ispath(Q,rows) and rows[x]>>Q[0]&1:q+=1
 return (subset_dp(rows),q,h-t,h),rows

def main():
 bs=profiles();digits=[b['q'] for b in bs]
 require(all(b['H']==5 and b['E']==32 for b in bs),'four-block counts')
 require(len(set(digits))==16 and all(0<=d<32 for d in digits),'injective residue digit')
 decoder={b['q']:b['mask'] for b in bs};count=0
 seed=(1,0,0,1)
 for length in range(4):
  residues=set()
  for word in product(range(16),repeat=length):
   state=seed
   for a in word:state=apply(state,bs[a])
   F,q,g,h=state;residues.add(q%(32**length))
   decoded=[]
   while h>1:
    previous_h=h//5
    digit=(q*pow(previous_h,-1,32))%32
    require(digit in decoder,'terminal digit is admitted')
    decoded.append(decoder[digit]);q=(q-digit*previous_h)//32;h=previous_h
   require(tuple(decoded[::-1])==word and q==seed[1],'exact reverse word')
   count+=1
  require(len(residues)==16**length,'distinct finite residue addresses')
 # A nonconstant cut which really avoids eventual divisibility by five.
 state=seed;repeated=[]
 for r in range(1,21):
  state=apply(state,bs[2]);repeated.append(state[0])
  if r>=3:require(state[0]%5==(3*pow(2,r-2,5))%5!=0,'repeat-cut nonzero mod five')
  require(state[0]==(100*32**r+(216*r+359)*5**r)//675,'closed repeated-cut formula')
  require((100*32**r+(216*r+359)*5**r)%675==0,'closed formula integrality')
  require(state[0]%3==((r-1)*pow(2,r,3))%3,'source-specific mod three filter')
 # Every nonzero XOR displacement occurs: any four fixed linear bit probes
 # have a kernel direction and necessarily identify two admitted digits.
 xor_hostiles=[]
 for shift in range(1,32):
  pair=next(((d,d^shift) for d in digits if d^shift in decoder),None)
  require(pair is not None,'all nonadaptive linear kernels obstructed')
  xor_hostiles.append([shift,*pair])
 @lru_cache(None)
 def decision_tree(candidates,depth):
  if len(candidates)==1:return ('digit',candidates[0])
  if depth==0 or len(candidates)>2**depth:return None
  for mask in range(1,32):
   halves=tuple(tuple(d for d in candidates if (d&mask).bit_count()%2==b) for b in (0,1))
   if min(map(len,halves))==0 or max(map(len,halves))>2**(depth-1):continue
   zero=decision_tree(halves[0],depth-1)
   if zero is None:continue
   one=decision_tree(halves[1],depth-1)
   if one is not None:return ('query',mask,zero,one)
  return None
 tree=decision_tree(tuple(sorted(digits)),4)
 require(tree is not None,'four adaptive linear probes')
 for d in digits:
  node=tree;used=0
  while node[0]=='query':
   node=node[2+(d&node[1]).bit_count()%2];used+=1
  require(node==('digit',d) and used==4,'decision tree literal decoding')
 mixed_controls=0;records=[]
 for mask in range(16):
  for bit in (0,1):
   for word in [(('D',mask),('S',bit)),(('S',bit),('D',mask))]:
    expected=seed
    for kind,a in word:expected=apply(expected,bs[a]) if kind=='D' else single(expected,a)
    found,rows=actual(word);require(found==expected,'independent mixed-block state')
    records.append((word,found,rows));mixed_controls+=1
 first=next(r for r in records if r[0]==(('D',2),('S',0)))
 second=next(r for r in records if r[0]==(('S',0),('D',9)))
 require(first[1][1]==second[1][1] and first[1][3]==second[1][3] and first[1]!=second[1],
         'mixed-generator ambiguity survives equal source and vertex clock')
 print(json.dumps({'status':'passed','rows':ROWS,'profiles':bs,'q_digit_set':sorted(digits),
                   'decoded_words':count,'repeat_mask':2,'repeat_F':repeated,
                   'independent_mixed_state_controls':mixed_controls,
                   'nonadaptive_linear_probe_minimum':5,'xor_kernel_hostiles':xor_hostiles,
                   'adaptive_linear_probe_minimum':4,'adaptive_probe_tree':tree,
                   'mixed_same_h_q_hostile':[{'word':w,'state_Fqgh':state,'rows':rows,'vertices':len(rows)} for w,state,rows in (first,second)],
                   'scope':'Exact cut-profile enumeration and finite decoder controls; all-depth inversion is proved by division and the distinct residues.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
