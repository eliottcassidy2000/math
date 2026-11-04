#!/usr/bin/env python3
"""All-order H=13 fixed-two-apex one-letter continuation certificate."""
import hashlib,itertools,json
from pathlib import Path
from verify import value,graph,dp,strong_word,require,state,step,Span,INIT

EXPECTED={
 '120':{'zero':[33,23,29,13],'coefficients':[[7,32],[13,14],[13,18]]},
 '212':{'zero':[37,29,25,13],'coefficients':[[9,33],[13,18],[13,15]]},
 '230':{'zero':[37,29,29,13],'coefficients':[[5,36],[13,18],[13,18]]},
 '302':{'zero':[33,29,25,13],'coefficients':[[5,33],[13,18],[13,15]]},
 '310':{'zero':[37,23,37,13],'coefficients':[[3,40],[13,14],[13,26]]},
 '321':{'zero':[43,25,37,13],'coefficients':[[7,41],[13,15],[13,26]]},
}

def main():
 # The source cutoff n<=5 follows from 5^5>13^3 and Busch's cited bound.
 require(5**5>13**3,'exact cutoff inequality')
 one_apex_values=set()
 for m in range(5):
  for bits in itertools.product(range(2),repeat=m):
   rows=[0]*(m+1)
   for i,bit in enumerate(bits,1):
    if bit:rows[0]|=1<<i
    else:rows[i]|=1
    for j in range(1,i):rows[j]|=1<<i
   one_apex_values.add(dp(rows))
 require(one_apex_values=={1,3,5,9,11},'complete one-apex small core spectrum')
 require(all(value(w)!=13 for m in range(3) for w in itertools.product(range(4),repeat=m)), 'no smaller strong core')
 found=[''.join(map(str,w)) for w in itertools.product(range(4),repeat=3) if value(w)==13]
 require(found==list(EXPECTED),'complete five-name core presentations')
 # Certify the inherited clock-polynomial operator identity on an invariant
 # reachable span independently within this script, not from a saved receipt.
 span=Span();span.add(INIT);R=[INIT];i=0
 while i<len(R):
  for c in range(4):
   v=step(R[i],c)
   if span.add(v):R.append(v)
  i+=1
 for v in R:
  v1=step(v,3);v2=step(v1,3);v3=step(v2,3);v4=step(v3,3)
  require(all(a-6*b+11*c-6*d==0 for a,b,c,d in zip(v4,v3,v2,v1)),'tail-clock recurrence')
 checks=0;cases={}
 for word,rec in EXPECTED.items():
  w=tuple(map(int,word));require(strong_word(w),'strong core')
  require(value(w)==dp(graph(w))==13,'core DP')
  require([value(w+(c,)) for c in range(4)]==rec['zero'],'zero-clock row')
  rows=[]
  for q in range(4):
   vals=[]
   for c in range(4):
    h=value(w+(3,)*q+(c,));require(h==dp(graph(w+(3,)*q+(c,))),'independent child DP');checks+=1
    if q>=1 and c<3:
     a,b=rec['coefficients'][c];require(h==a+b*2**q,'positive-clock initial')
    if c==3:require(h==13,'neutral sink')
    vals.append(h)
   rows.append(vals)
  require(not any(39 in row or (41 in row and 43 in row) for row in rows[:2]),'q=0,1 obstruction')
  require(all(a+4*b>43 and b>0 for a,b in rec['coefficients']),'q>=2 strict bound')
  cases[word]={**rec,'q1':rows[1]}
 # The restriction to append at one end really matters.
 parent=(1,2,0,3);lower=parent+(1,);upper=(3,)+parent
 require([value(w) for w in (parent,lower,upper)]==[13,41,43],'two-ended escape')
 require([dp(graph(w)) for w in (parent,lower,upper)]==[13,41,43],'two-ended escape DP')
 result={'status':'PASS','five_name_core_words':found,'cases':cases,'independent_child_DPs':checks,'one_apex_counts_at_order_at_most_five':sorted(one_apex_values),
 'all_order_no_39':True,'all_order_no_common_parent_41_43':True,
 'two_ended_escape':{'parent':'1203','H':13,'append_1':'12031','lower_H':41,'prepend_3':'31203','upper_H':43,
 'parent_strong':strong_word(parent),'lower_strong':strong_word(lower),'upper_strong':strong_word(upper)}}
 raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('append_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','complete_core_presentations':6,'independent_child_DPs':checks,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
