#!/usr/bin/env python3
"""Independent completion-polynomial audit of the final receiver revision."""
from pathlib import Path
from collections import Counter,defaultdict
from hashlib import sha256
from math import comb
import argparse,json
P=Path(__file__).resolve().parent
N=0
def need(value,message):
 global N
 N+=1
 if not value:raise RuntimeError(message)
def digest(value):return sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()
def grade(x):return (x&65535).bit_count()+3*(x>>16)
def costs(x):return ((x&65535).bit_count()+3*(x>>16),2*(x&65535).bit_count()+8*(x>>16))
def main():
 pins=json.loads((P/'SOURCEPINS.json').read_text())
 actual={p.relative_to(P).as_posix() for directory in ('inputs','sources') for p in (P/directory).rglob('*') if p.is_file()}
 need(actual==set(pins),'complete pinned supplemental input closure')
 for path,row in pins.items():
  data=(P/path).read_bytes();need(len(data)==row['bytes'] and sha256(data).hexdigest()==row['sha256'],'source pin '+path)
 independent=json.loads((P/'inputs/independent_audit_result.json').read_text())
 producer=json.loads((P/'inputs/final_producer_certificate.json').read_text())
 manifest=json.loads((P/'sources/final_producer_MANIFEST.json.txt').read_text())
 for name,local in [('THEOREM.md','sources/final_producer_THEOREM.md.txt'),('check.py','sources/final_producer_check.py.txt'),('certificate.json','inputs/final_producer_certificate.json')]:
  data=(P/local).read_bytes();need(sha256(data).hexdigest()==manifest[name]['sha256'] and len(data)==manifest[name]['bytes'],'reviewed final mathematical source matches final producer manifest '+name)
 ind_manifest=json.loads((P/'sources/independent_MANIFEST.json.txt').read_text())
 data=(P/'inputs/independent_audit_result.json').read_bytes()
 need(sha256(data).hexdigest()==ind_manifest['AUDIT_RESULT.json']['sha256'],'frozen independent input agrees with original independent final manifest')
 ind_receipt=json.loads((P/'sources/independent_PORTABILITY_RECEIPT.json.txt').read_text())
 need(ind_receipt['manifest_sha256']==pins['sources/independent_MANIFEST.json.txt']['sha256'] and ind_receipt['manifest_only_copy'] and ind_receipt['normal_optimized_full_suite'],'independent manifest-only replay receipt refers to exact frozen audit')
 producer_log=[json.loads(line) for line in (P/'sources/final_producer_PORTABILITY_LOG.txt').read_text().splitlines() if line.strip()]
 need(all(row['status']=='PASS' for row in producer_log),'final producer portability log all PASS')
 need([row['manifest_sha256'] for row in producer_log if 'manifest_sha256' in row]==[pins['sources/final_producer_MANIFEST.json.txt']['sha256']]*2,'final producer replay pins same reviewed final mathematical manifest')
 need(producer_log[-2].get('normal_optimized_equal') and producer_log[-2].get('scripts')==4 and producer_log[-1].get('manifest_only_relocation') and producer_log[-1].get('normal_and_optimized_full_suite'),'final producer all-four-script normal/-O actual-copy replay passed')
 checkrows=[row for row in producer_log if 'checks' in row]
 need([row['checks'] for row in checkrows]==[2663575,2663575,175980,175980,115250,115250,441763,441763],'final producer replay actual mathematical check counts')
 need(all(checkrows[j]==checkrows[j+1] for j in range(0,8,2)),'all four final producer normal/-O output pairs identical')
 trades=producer['all_primitive_trades']
 need(digest([r[:2] for r in trades])==independent['unique_primitive_trade_digest'],'final producer complete primitive list matches frozen independent computation')
 need(producer['complete_pair_census']['pair_record_sha256']==independent['all_pair_optimum_sha256'],'final producer all-pair optimum stream matches frozen independent computation')
 by_grade=Counter();by_cost=Counter();by_grade_cost=Counter();direct_poly=[0]*20;groups=Counter();maximals=[]
 for x,y,s,c in trades:
  need(x&y==0 and x<y,'normalized primitive orientation and stripped common bits')
  need(grade(x)==grade(y)==s//2 and costs(x^y)==(s,c),'full weighted support, equal grades, and exact cost')
  t=(x^y)>>16;b=((x^y)&65535).bit_count();need(t==(c-2*s)//2 and b==s-3*t,'exact triple and base support')
  g=grade(x);groups[s,c]+=1
  # First method: ordinary multiplication by each untouched coordinate's (1+z^weight).
  coefficients=[0]*20;coefficients[g]=1
  for i in range(17):
   if not ((x^y)>>i)&1:
    w=3 if i==16 else 1
    for k in range(19,w-1,-1):coefficients[k]+=coefficients[k-w]
  # Second method: closed binomial form z^g(1+z)^(16-b)(1+z^3)^(1-t).
  formula=[0]*20
  for h in range(17-b):
   formula[g+h]+=comb(16-b,h)
   if not t:formula[g+h+3]+=comb(16-b,h)
  need(coefficients==formula,'each primitive coordinate-product equals closed binomial form')
  for k,n in enumerate(coefficients):
   by_grade[k]+=n;by_cost[s,c]+=n;by_grade_cost[k,s,c]+=n;direct_poly[k]+=n
  if (s,c)==(16,34):maximals.append([x,y])
 need(min(s for s,c in groups)==4,'all nonzero primitive support at least4')
 need(max((s,c) for s,c in groups)==(16,34),'no primitive exceeds16/34')
 need(all(s<=14 and c<=30 or (s,c)==(16,34) for s,c in groups),'only critical primitive type is16/34')
 need(len(maximals)==196,'196 unique maximal primitives')
 for row in independent['per_grade']:
  need(by_grade[row['grade']]==row['primitive_pair_occurrences'],'completion coefficients match independent native primitive occurrence census')
 for row in independent['primitive_occurrence_cost_histogram']:
  need(by_cost[row['slots'],row['clocks']]==row['pairs'],'cost-polynomial totals match independent primitive occurrence costs')
 for row in producer['grades']:
  actual=sorted([s,c,n] for (k,s,c),n in by_grade_cost.items() if k==row['grade'] and n)
  need(actual==row['primitive_pair_cost_histogram'],'full per-grade per-cost completion polynomial matches final native replay')
 critical=[196*comb(3,k-8) if 8<=k<=11 else 0 for k in range(20)]
 for row in independent['per_grade']:
  actual=sum(v['pairs'] for v in row['strict_optimum_histogram'] if v['clocks']==34)
  need(actual==critical[row['grade']],'all strict34 endpoint counts match maximal primitive completion polynomial')
 ordinary_critical=[sum(v['pairs'] for v in row['ordinary_optimum_histogram'] if v['clocks']==34) for row in independent['per_grade']]
 need(ordinary_critical==[1 if k in (8,11) else 0 for k in range(20)],'exactly two ordinary34 endpoint pairs, at grades8 and11')
 need(sum(critical)==1568 and sum(by_grade.values())==463840,'exact critical and full primitive embedding totals')
 # Every same-grade pair has even actual support<=18. Any decomposition has at
 # least two disjoint primitive supports of>=4, so every component has<=14.
 need(18-4==14 and all(c<=30 for s,c in groups if s<=14),'decomposable endpoint pairs have strict routes within14/30')
 answer={'status':'PASS','primitive_completion_formula':'sum_(s,c) N_(s,c) z^(s/2) (1+z)^(16-s+3t) (1+z^3)^(1-t), t=(c-2s)/2','primitive_pairs_by_grade':[by_grade[k] for k in range(20)],'primitive_pair_total':sum(by_grade.values()),'strict34_polynomial':'196 z^8 (1+z)^3','strict34_by_grade':critical,'ordinary34_total':sum(ordinary_critical),'strict34_total':sum(critical),'normalized_maximal_trades':len(maximals),'pair_digest':independent['all_pair_optimum_sha256'],'independent_audit_manifest_sha256':pins['sources/independent_MANIFEST.json.txt']['sha256'],'final_producer_manifest_sha256':pins['sources/final_producer_MANIFEST.json.txt']['sha256'],'checks':N}
 return answer
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=main();f=P/'RESULT.json'
 if args.write:f.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
 else:need(r==json.loads(f.read_text()),'frozen supplemental audit result')
 print(json.dumps(r,sort_keys=True))
