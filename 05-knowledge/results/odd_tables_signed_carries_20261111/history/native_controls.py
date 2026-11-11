#!/usr/bin/env python3
"""Independent literal controls for the native zero-digit operation."""
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('port_probe',HERE/'port_probe.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
def converse(rows):return [sum(1<<j for j in range(len(rows)) if rows[j]>>i&1) for i in range(len(rows))]
def state(rows):
 h,e,u,d,k=p.actual_state(rows);_,s,v,_,_=p.actual_state(converse(rows))
 return [h,e,s,u,v,d,k]
def birth(rows,pred):
 ans=rows[:];n=len(rows);ans[pred]|=1<<n;ans.append(((1<<n)-1)^(1<<pred));return ans
def inflate(rows,star,k):
 labels=list(range(star))+[star]*k+list(range(star+1,len(rows)));out=[0]*len(labels)
 for i in range(len(labels)):
  for j in range(i+1,len(labels)):
   forward=labels[i]==labels[j] or rows[labels[i]]>>labels[j]&1
   a,b=(i,j) if forward else (j,i);out[a]|=1<<b
 return out
def strong(rows):
 n=len(rows)
 for source in range(n):
  seen=1<<source;todo=[source]
  while todo:
   avail=rows[todo.pop()]&~seen
   while avail:
    bit=avail&-avail;avail-=bit;seen|=bit;todo.append(bit.bit_length()-1)
  if seen!=(1<<n)-1:return False
 return True

def main():
 cases=0;strong_cases=0
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   rows=p.c.tournament(n,code)
   for root in range(n):
    Y=birth(rows,root);a=n;Y=birth(Y,a);b=n+1
    before=state(Y)
    Z=birth(Y,a);cc=n+2;Z=birth(Z,cc)
    after=state(Z)
    assert [after[0],after[4],after[2]]==[3*before[0],3*before[0]+2*before[2],before[0]]
    assert sum(bool(Z[v]>>(n+3)&1) for v in range(len(Z)))==1
    assert sum(bool(Z[v]>>cc&1) for v in range(len(Z)))==1
    if strong(Y):assert strong(Z);strong_cases+=1
    cases+=1
 data=json.loads((HERE/'native_repair_receipt.json').read_text())
 witness=data['all_distance_two_patches'][1];inflations=[]
 for k in [1,2,3]:
  rows=inflate(witness['rows'],7,k);actual=state(rows);D=1<<k
  expected=[2673*D+645,1,891*D+215,169130*3**k+95820*D-6453,
            3267*D+789,1485*D+359,237530*3**k+134800*D-8982]
  assert actual==expected,(k,actual,expected)
  assert strong(rows)
  inflations.append({'k':k,'order':len(rows),'HESUVdK':actual})
 receipt={'status':'FINITE-EXACT independent controls of the proved universal operation',
          'small_universe':'Every labelled tournament of orders1..3 and every marked root, followed by two births and the native zero operation',
          'rooted_cases':cases,'strong_input_cases':strong_cases,
          'literal_inflation_controls':inflations,
          'desired_plus_two_target_claim':False}
 (HERE/'native_controls_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
