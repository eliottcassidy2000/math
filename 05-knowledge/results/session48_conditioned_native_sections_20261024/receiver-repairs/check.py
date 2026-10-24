#!/usr/bin/env python3
"""Fresh 16-by-128 group algebra, exact repair hypergraph and literal controls."""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
from hashlib import sha256
import argparse,json,math
P=Path(__file__).resolve().parent;I=P/'inputs';CHECKS=0;FULL=(1<<128)-1

def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)
def read(s):return json.loads((I/s).read_text())
def plus(a,b):return ((a[0]+b[0])%128,a[1]^b[1])
def minus(a,b):return ((a[0]-b[0])%128,a[1]^b[1])
def pack(a):return 16*a[0]+a[1]
def unpack(q):return(q//16,q%16)
def shifted(row,label):
 x,y=label
 return [((row[j^y]<<x)|(row[j^y]>>(128-x)))&FULL for j in range(16)]
def extend(d,label):
 out=[r[:] for r in d]+[[0]*16]
 for k in range(1,len(out)):
  v=shifted(d[k-1],label);out[k]=[a|b for a,b in zip(out[k],v)]
 return out
def support(labels):
 d=[[1]+[0]*15]
 for label in labels:d=extend(d,label)
 return d
def holes(row):return sorted(16*x+y for y in range(16) for x in range(128) if not row[y]>>x&1)
def size(row):return sum(x.bit_count() for x in row)
def full(row):return all(v==FULL for v in row)
def failure(d,label,grades):
 for k in grades:
  v=[a|b for a,b in zip(d[k],shifted(d[k-1],label))]
  if not full(v):
   q=holes(v)[0];previous=pack(minus(unpack(q),label))
   need(q in holes(d[k]) and previous in holes(d[k-1]),'literal adjacent-hole difference witness')
   return {'grade':k,'hole':q,'previous_hole':previous}
 return None

def budgets(h,l,u,c):
 rc=[2]*(h-1)+[c];W=lambda n:sum(sorted(rc,reverse=True)[:min(n,h)])
 rows=[]
 for e in range(8):
  maxima={s:0 for s in ('entry','adjacent','outside','helper')};supports=dict(maxima)
  for k in range(l,u+e+1):
   r0=max(l,k-e)
   for r in range(h+1):
    if not 0<=k-r<=e:continue
    if r<l:dd=min(r+l,2*h-r-l);de=l-r
    elif r>u:dd=min(r+u,2*h-r-u);de=r-u
    else:continue
    maxima['entry']=max(maxima['entry'],W(dd)+4*de);supports['entry']=max(supports['entry'],dd+de)
   for r in range(r0+1,u+1):
    if 0<=k-r<e:
     dd=min(2*r-1,2*h-2*r+1);maxima['adjacent']=max(maxima['adjacent'],W(dd)+4);supports['adjacent']=max(supports['adjacent'],dd+1)
   dd=2*min(r0,h-r0)
   if 0<k-r0<e:maxima['outside']=max(maxima['outside'],W(dd)+8);supports['outside']=max(supports['outside'],dd+2)
   maxima['helper']=max(maxima['helper'],W(dd)+15);supports['helper']=max(supports['helper'],dd+2)
  limit=49+c if h==18 else 45+c;slimit=20 if h==18 else 18
  need(max(maxima.values())<=limit and max(supports.values())<=slimit,'all exact feasible-grade routing budgets')
  rows.append({'free_outside':e,'repair_cost':c,'interior':[l,u+e],'stage_clock_maxima':maxima,'stage_support_maxima':supports,'clock_limit':limit,'slot_limit':slimit})
 return rows

def main():
 old=read('native_packet_decomposition_20261023/native/domains.json')['sources'];e=read('metagraph_order8_electrical_20261023/receipt.json')['order8'];X=e['higher_layers'][-1]['lifts'][0];Y=e['radical_lifts_mod4']
 labels={r['source']:(X[r['source']]%128,sum((Y[j][r['source']]%2)<<(j-1) for j in range(1,5))) for r in old}
 costs={r['source']:sum(a!=b for a,b in zip(*(s['path'] for s in r['exact_all_shared_clocks_target_domain']))) for r in old}
 cheap=sorted(s for s in labels if costs[s]==2);outside=sorted(set(labels)-set(cheap))
 need(len(cheap)==20 and len(outside)==8,'actual source20+8 decomposition')
 for r in old:need(pack(labels[r['source']])==(r['code'][0]<<4)|sum(v<<j for j,v in enumerate(r['code'][1:])),'old labels reconstructed from original potentials')
 for q in range(2048):
  row=[0]*16;row[q%16]=1<<(q//16)
  for label in ((127,15),(1,1)):
   moved=shifted(row,label);want=pack(plus(unpack(q),label));need(holes([FULL^v for v in moved])==[want],'literal 16-column translation and wrap')
 source=read('external_pin_sections_20261024/erasure/receipt.json');oldrepair=read('external_pin_sections_20261024/theory/repair_labels_receipt.json')
 levels=[];exceptions=[];forbidden=defaultdict(list);triple_all=[];quadrows=[];quadpatterns=Counter();pairfail=[]
 for n in range(5):
  profiles=Counter();universal=set(range(21-n));nocomplete=[]
  for erased in combinations(cheap,n):
   d=support([labels[s] for s in cheap if s not in erased]);profile=tuple(size(v) for v in d);profiles[profile]+=1
   need(profile==profile[::-1],'each fresh receiver complement support counts')
   complete={k for k,v in enumerate(d) if full(v)};universal&=complete
   if not complete:nocomplete.append(list(erased))
   if n==3:
    good=[]
    for s in outside:
     f=failure(d,labels[s],(8,9,10));need(f is None,'each actual repair on every one of1140 triples');good.append(s)
    triple_all.append({'erased':list(erased),'good_sources':good})
    if not full(d[8]) or not full(d[9]):
     hs={k:holes(d[k]) for k in (7,8,9,10)};exceptions.append({'erased':list(erased),'complete_grades':sorted(complete),'holes':hs})
     for k in (8,9,10):
      for q in hs[k]:
       for p in hs[k-1]:forbidden[pack(minus(unpack(q),unpack(p)))].append({'erased':list(erased),'grade':k,'hole':q,'previous_hole':p})
   if n==4:
    good=[];bad={}
    for s in outside:
     f=failure(d,labels[s],(8,9))
     if f is None:good.append(s)
     else:bad[str(s)]=f
    need(len(good)>=2,'every fourth-freeze set has at leasttwo actual repairs');quadpatterns[len(good)]+=1
    quadrows.append({'erased':list(erased),'good_sources':good,'failure_witnesses':bad})
    # Pairs involving a successful singleton are certified by adjacent full grades.
    for p,q in combinations(outside,2):
     if p in good or q in good:continue
     dp=extend(d,labels[p]);f=failure(dp,labels[q],(8,9,10))
     if f is not None:pairfail.append({'erased':list(erased),'pair':[p,q],**f})
  row={'erasure_count':n,'total_pin_sets':sum(profiles.values()),'universal_complete_grades':sorted(universal),'no_complete_grade_pin_sets':nocomplete,'coverage_profiles':[{'coverage':list(k),'pin_sets':v} for k,v in sorted(profiles.items())]};levels.append(row)
  need(sum(profiles.values())==math.comb(20,n),'complete frozen-set universe')
  if n<4:
   previous=source['levels'][n];need(row['coverage_profiles']==previous['profiles'],'all1351 old profiles independently agree')
   need(sorted(universal)==previous['complete_grades_for_every_pin_set'],'old universal band')
 need(len(exceptions)==6 and len(forbidden)==32,'six exceptional sections and exactly32 forbidden labels')
 need(sorted(forbidden)==oldrepair['forbidden_labels'],'literal forbidden-label set matches source certificate')
 for row,oldrow in zip(exceptions,source['adaptive_repair']['exceptional_original_centers']):
  need(row['erased']==oldrow['erased'] and {str(k):v for k,v in row['holes'].items()}==oldrow['original_holes'],'every exceptional adjacent-grade hole independently recovered')
 sets=[set(r['good_sources']) for r in quadrows]
 transversals={n:[list(t) for t in combinations(outside,n) if all(set(t)&s for s in sets)] for n in (1,2)}
 need(not transversals[1] and transversals[2]==[[1923,2087],[1923,2835]],'sharp two-source protected availability certificate')
 need(dict(quadpatterns)=={8:4522,7:241,6:47,5:17,4:9,3:5,2:4},'exactfour-freeze distribution')
 need(len(pairfail)==6,'exact six two-source receiver failures')
 # One explicit mixed relation: actual target(1841)=target(1457), not an indexed-owner equality.
 erase={176,1425,2071};linked=1457;repair=1841;base=[s for s in cheap if s not in erase|{linked}];d=support([labels[s] for s in base]);block=plus(labels[linked],labels[repair]);mh={}
 for k in (8,9,10):mh[str(k)]=holes([a|b for a,b in zip(d[k],shifted(d[k-2],block))])
 need(mh=={'8':[477],'9':[],'10':[1098]},'exact mixed same-target relation obstruction')
 # Separate literal subset enumeration, not the support DP, proves this relational boundary.
 enum=[set() for _ in range(19)];witness={};sums=[(0,0)]*(1<<len(base));grades=[0]*(1<<len(base))
 for mask in range(1<<len(base)):
  if mask:
   bit=(mask&-mask).bit_length()-1;prev=mask&(mask-1);sums[mask]=plus(sums[prev],labels[base[bit]]);grades[mask]=grades[prev]+1
  for bit in (0,1):
   k=grades[mask]+2*bit;q=pack(plus(sums[mask],block)) if bit else pack(sums[mask]);enum[k].add(q)
   if k in (8,9,10):witness.setdefault((k,q),[base[i] for i in range(len(base)) if mask>>i&1]+([linked,repair] if bit else []))
 for k in (8,9,10):need(sorted(set(range(2048))-enum[k])==mh[str(k)],'literal131072 full-record relation assignments match')
 native=json.loads((P/'native_certificate.json').read_text());need(all(native['labels'][str(s)]==pack(labels[s]) for s in labels),'all arithmetic labels join native certificate')
 need(sum(sorted(costs.values(),reverse=True)[:14])==39 and native['helper_direct_clock_bound']==29,'actual boundary direct cost39+29, not56+29')
 for pair in transversals[2]:need(sum(sorted([2]*16+[costs[s] for s in pair],reverse=True)[:14])+29==60,'ten-frozen protected-pair direct clock60')
 for c in (3,4):need(sum(sorted([2]*16+[c],reverse=True)[:14])+29==55+c,'eleven-frozen selected-repair clock58 or59')
 return {'status':'FINITE-EXACT independently rebuilt native receiver arithmetic','source_commit':'c09b642a629b7d12b529b88c2697abaee78e73e7','cheap_sources':cheap,'outside_sources':outside,'labels':{s:pack(labels[s]) for s in labels},'costs':costs,'levels':levels,'triple_universe':triple_all,'triple_exceptional_sections':exceptions,'forbidden_triple_labels':sorted(forbidden),'forbidden_triple_witnesses':dict(forbidden),'quadruple_distribution':dict(quadpatterns),'quadruple_universe':quadrows,'minimum_protected_pairs':transversals[2],'two_repair_slot_failures':pairfail,'mixed_relation':{'frozen':sorted(erase),'equal_target_sources':[repair,linked],'holes':mh,'literal_assignments':131072},'triple_stage_budgets':sum((budgets(18,8,10,c) for c in (3,4)),[]),'quadruple_stage_budgets':sum((budgets(17,8,9,c) for c in (3,4)),[]),'global_clock_bound':68,'core_direct_clock_bound':39,'helper_direct_clock_bound':29,'protected_pair_ten_frozen_clock_bound':60,'single_successful_repair_eleven_frozen_clock_bounds':[58,59],'checks':CHECKS}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();result=json.loads(json.dumps(main()));out=P/'certificate.json'
 if a.write:out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 else:need(result==json.loads(out.read_text()),'frozen exhaustive certificate exact')
 print(json.dumps({'status':'PASS','checks':CHECKS,'triples':1140,'quadruples':4845,'minimum_protected_pairs':result['minimum_protected_pairs'],'global_clocks':68},sort_keys=True))
