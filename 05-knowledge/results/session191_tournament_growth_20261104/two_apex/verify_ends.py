#!/usr/bin/env python3
"""Complete core/end-extension certificate for the H13 -> H41 dead branch."""
import hashlib,itertools,json
from pathlib import Path
from verify import value,graph,dp,strong_word,require
EXPECTED=['10310','10321','11310','11321','12011','12031','13101','31011','31031']

def extensions(w,bound):
 for births in range(bound+1):
  for left in range(births+1):
   for extra in itertools.product(range(4),repeat=births):
    yield extra[:left]+w+extra[left:],births,left

def main():
 require(5**5>13**3 and 5**7>41**3,'source and middle cutoffs')
 require(5**10>125**3 and 5**10>127**3,'target cutoffs')
 cores=[w for w in itertools.product(range(4),repeat=3) if value(w)==13]
 require(len(cores)==6,'six H13 cores')
 middle={};initial_tests=0
 for core in cores:
  for w,b,l in extensions(core,2):
   h=value(w);require(h==dp(graph(w)),'independent first-stage DP');initial_tests+=1
   if h==41 and strong_word(w):middle[w]={'H13_core':''.join(map(str,core)),'births':b,'left_births':l}
 require([''.join(map(str,w)) for w in sorted(middle)]==EXPECTED,'all nine H41 middle cores')
 # Keep all counts, including nonstrong outputs; no strength filter can hide
 # a target. Cache literal duplicate presentations but count the full universe.
 checked={};hist={};tests=0
 for core in sorted(middle):
  for w,b,l in extensions(core,3):
   if w not in checked:
    h=value(w);require(h==dp(graph(w)),'independent second-stage DP');checked[w]=h
   h=checked[w];tests+=1;hist[h]=hist.get(h,0)+1
   require(h not in (125,127),'forbidden continuation appears')
 require(initial_tests==342 and tests==2817,'complete presentation universes')
 # A positive next-generation control: the other initial child can continue.
 controls=[((3,1,2,0),43),((3,1,2,0,0),131),((3,3,1,2,0),133)]
 for w,h in controls:require(value(w)==dp(graph(w))==h and strong_word(w),'upper positive control')
 result={'status':'PASS','H13_core_words':[''.join(map(str,w)) for w in cores],
 'first_stage_presentations':initial_tests,'H41_middle_words':EXPECTED,
 'middle_provenance':{''.join(map(str,w)):v for w,v in sorted(middle.items())},
 'second_stage_presentations':tests,'independent_distinct_second_stage_DPs':len(checked),
 'second_stage_count_histogram':hist,'missing_targets':[125,127],
 'upper_positive_control':{'parent_word':'3120','H':43,'append0_word':'31200','lower_H':131,'prepend3_word':'33120','upper_H':133}}
 raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('ends_receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','middle_cores':len(middle),'first_stage_presentations':initial_tests,'second_stage_presentations':tests,
 'independent_distinct_second_stage_DPs':len(checked),'missing_targets':[125,127],'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
