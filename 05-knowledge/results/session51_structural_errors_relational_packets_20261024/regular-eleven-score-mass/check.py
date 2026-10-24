#!/usr/bin/env python3
"""Atlas-free regular11 count: score DP and complete nonrigid group mass."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb,factorial
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
CHECKS=Counter()
def need(ok,message):
 CHECKS[message]+=1
 if not ok:raise ValueError(message)

@lru_cache(None)
def grouped(scores):
 n=len(scores)
 if any(d<0 or d>=n for d in scores) or sum(scores)!=n*(n-1)//2:return 0
 if n<=1:return 1
 target=scores[0];groups=tuple(sorted(Counter(scores[1:]).items()));answer=0
 def choose(k,left,remaining,weight):
  nonlocal answer
  if k==len(groups):
   if left==0:answer+=weight*grouped(tuple(sorted(remaining)))
   return
  degree,multiplicity=groups[k]
  for out in range(min(multiplicity,left)+1):
   if degree==0 and out<multiplicity:continue
   if degree==n-1 and out>0:continue
   choose(k+1,left-out,remaining+[degree]*out+[degree-1]*(multiplicity-out),weight*comb(multiplicity,out))
 choose(0,target,[],1)
 return answer

@lru_cache(None)
def subsets(scores):
 # Independently distinguish the LAST vertex and visit every literal subset.
 n=len(scores)
 if any(d<0 or d>=n for d in scores) or sum(scores)!=n*(n-1)//2:return 0
 if n<=1:return 1
 target=scores[-1];answer=0
 for outgoing in combinations(range(n-1),target):
  outgoing=set(outgoing)
  residual=tuple(sorted(scores[j]-(j not in outgoing) for j in range(n-1)))
  answer+=subsets(residual)
 return answer

def run():
 for record in json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']:
  need(hashlib.sha256((HERE/record['archived_path']).read_bytes()).hexdigest()==record['sha256'],'inherited source hash')
 inherited=json.loads((HERE/'sources/nonrigid-receipt.json').read_text())
 classes=inherited['q1_classes']+[r for r in inherited['q3_candidates'] if r['error']]+inherited['q7_classes']+inherited['q9_classes']
 need(len(classes)==18 and len({tuple(r['rows']) for r in classes})==18,'eighteen inherited nonrigid representatives without duplicate control')
 need(all(sum(r['vertex_orbit_sizes'])==11 and r['full_automorphism_order']>1 for r in classes),'nonrigid group and orbit typing')
 need(all(factorial(11)%r['full_automorphism_order']==0 for r in classes),'group orders divide eleven factorial')
 nonrigid_mass=sum(factorial(11)//r['full_automorphism_order'] for r in classes)
 pointed_nonrigid=sum(len(r['vertex_orbit_sizes']) for r in classes)
 defect_nonrigid=sum(r['error'] for r in classes)
 need(pointed_nonrigid==78 and defect_nonrigid==660,'nonrigid pointed and defect inputs')
 controls={}
 for n in (1,3,5,7,9,11):
  scores=((n-1)//2,)*n;g=grouped(scores);s=subsets(scores)
  need(g==s,'grouped score count matches literal subset recurrence')
  controls[n]=g
 need([controls[n] for n in (1,3,5,7)]==[1,2,24,2640],'small independently certified regular controls')
 profile=(4,)*5+(5,)*5;profile_count=grouped(profile)
 need(profile_count==subsets(profile),'almost regular fixed profile independent count')
 labelled=controls[11]
 need(comb(10,5)*profile_count==labelled,'fixed root deletion labelled score profile identity')
 difference=labelled-nonrigid_mass
 need(difference>=0 and difference%factorial(11)==0,'remaining regular labelled mass is rigid classes')
 rigid=difference//factorial(11);regular=rigid+len(classes);pointed=11*rigid+pointed_nonrigid
 error=110*rigid+defect_nonrigid
 need(error==11*(pointed-regular),'regular error from pointed class difference')
 need(error==10*pointed-120,'structural regular correction agrees with mass route')
 need(pointed%11==1 and error%110==0,'pointed and error divisibility guards')
 histogram=Counter(r['error'] for r in classes);histogram[110]+=rigid
 need(sum(histogram.values())==regular and sum(e*c for e,c in histogram.items())==error,'full regular defect histogram arithmetic')
 strong_error=error+9166;total_error=error+19480
 need(total_error-strong_error==10314,'inherited SCC padding')
 return {'status':'PASS','scope':'Atlas-free regular-eleven count from a complete inherited nonrigid classification and labelled score counting. No almost-regular10 or regular11 atlas is read. U11 and unrestricted a11 are not computed.','labelled_regular_counts':controls,'grouped_memo_states':grouped.cache_info().currsize,'literal_subset_memo_states':subsets.cache_info().currsize,'fixed_profile_almost_regular10_labelled_count':profile_count,'all_profiles_factor':comb(10,5),'labelled_regular11':labelled,'nonrigid_regular11_labelled_mass':nonrigid_mass,'rigid_regular11_classes':rigid,'regular11_classes':regular,'almost_regular10_classes':pointed,'regular11_error_histogram':sorted(histogram.items()),'regular11_error':error,'strong11_error':strong_error,'total11_error':total_error,'unrestricted_order11_count_formula':f'(U11-{total_error})/11','checks':dict(sorted(CHECKS.items()))}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
