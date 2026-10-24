#!/usr/bin/env python3
"""Exact finite controls for full-record relational descent and bridge costs."""
from itertools import combinations,product
from pathlib import Path
import argparse,json
HERE=Path(__file__).resolve().parent;COUNT=0
def need(ok,why):
 global COUNT
 COUNT+=1
 if not ok:raise RuntimeError(why)
def bits(n):return list(product((0,1),repeat=n))
def hamming(a,b):return sum(x!=y for x,y in zip(a,b))
def interval_tests():
 count=0;nonzero=0
 for h in range(1,9):
  for l in range(h+1):
   for u in range(l,h+1):
    for e in range(6):
     for t in range(6):
      for tt in range(6):
       for k in range(h+e+6):
        I=list(range(max(0,k-t-u),min(e,k-t-l)+1));J=list(range(max(0,k-tt-u),min(e,k-tt-l)+1))
        if not I or not J:continue
        delta=max(0,abs(t-tt)-(u-l));actual=min(abs(i-j) for i in I for j in J)
        need(actual==delta,'clipped interval gap is exact')
        i,j=min(product(I,J),key=lambda p:(abs(p[0]-p[1]),p))
        r,s=k-t-i,k-tt-j
        need(l<=r<=u and l<=s<=u and abs(i-j)==delta,'constructive nested-outside bridge grades')
        count+=1;nonzero+=bool(delta)
 return {'cases':count,'positive_outside_gap_cases':nonzero}
def grade_distance():
 pairs=0
 for h in range(7):
  sets=[frozenset(c) for k in range(h+1) for c in combinations(range(h),k)]
  for r in range(h+1):
   for s in range(h+1):
    A=[a for a in sets if len(a)==r];B=[b for b in sets if len(b)==s]
    need(max(len(a^b) for a in A for b in B)==min(r+s,2*h-r-s),'sharp actual grade-only symmetric difference')
    pairs+=len(A)*len(B)
 return {'enumerated_subset_pairs':pairs}
def descent_and_hostiles():
 # Exhaust EVERY feasible subset of the two-bit full-record universe, every
 # nonempty boundary relation and every requested support in C2.
 states=bits(2);cases=0
 for fmask in range(16):
  F=[x for i,x in enumerate(states) if fmask>>i&1]
  C={b:{(x[0]+x[1])%2 for x in F if x[0]==b} for b in (0,1)}
  for qmask in range(4):
   Q={q for q in (0,1) if qmask>>q&1}
   unary=all(Q<=C[b] for b in (0,1))
   relational=True
   for gmask in (1,2,3):
    G={b for b in (0,1) if gmask>>b&1}
    direct={(x[0]+x[1])%2 for x in F if x[0] in G}
    joined=set().union(*(C[b] for b in G))
    need(direct==joined,'full-boundary relational union identity')
    relational &= Q<=direct
   need(unary==relational,'all frozen fibres iff all nonempty boundary relations')
   cases+=1
 need(all({(x+y)%2 for y in (0,1)}=={0,1} for x in (0,1)),'each frozen x is charge onto')
 need({(x+y)%2 for x,y in states if x==y}=={0},'cross-receiver equality kills charge1')
 # Same-charge full representatives may carry DIFFERENT compatible choices.
 A=[('a0',0),('a1',0)];B=[('b0',0),('b1',0)]
 compat=lambda a,b:a[0][-1]!=b[0][-1]
 need(not compat(A[0],B[0]) and any(compat(a,b) for a in A for b in B),'fixed charge representatives lose existential full witness')
 need(all(any(x==y for y in (0,1)) for x in (0,1)),'each equality edge is locally onto')
 need(not any(x==y and y==z and z!=x for x,y,z in bits(3)),'cycle obstruction survives all local marginals')
 NA={(0,0),(1,1)};NB={(0,1),(1,0)}
 need(all({x[i] for x in NA}=={x[i] for x in NB}=={0,1} for i in range(2)),'equal neutral coordinate envelopes')
 need(min(hamming(a,b) for a in NA for b in NB)==1,'relational neutral exact gap differs from envelope zero')
 need(hamming((0,0),(1,1))==2,'neutral diagonal needs2 not1')
 # Two individually nonempty frozen sections; tight event budget disconnects.
 endpoint0=(0,0,0,1,1,1);endpoint1=(1,1,1,0,0,0)
 need(sum(endpoint0)==sum(endpoint1)==3 and hamming(endpoint0,endpoint1)==6,'binary ternary-relation six-slot barrier')
 return {'exhaustive_descent_cases':cases,'hostiles':['omitted receiver relation','noncanonical same-charge witnesses','cyclic separator inconsistency','relational neutral envelope false gap','neutral internal disconnection','ternary fixed-grade six-slot barrier']}
def native_budget_tests():
 # Uniform upper tests use ACTUAL old cost multiset; source identities remain
 # the native checker's responsibility. Frozen records never contribute work.
 base=[2]*20+[3]*5+[4]*3
 need(sum(sorted(base,reverse=True)[:18])==47,'native core top18 cost')
 count=0;interior=0;bridge=0;stage_max={};boundary=0
 bb=bits(3)
 for c in (3,4):
  Rcost=sorted([2]*17+[c],reverse=True)
  W=lambda a:sum(Rcost[:min(a,18)])
  outside_costs=[3]*5+[4]*3;outside_costs.remove(c);outside_costs.sort(reverse=True)
  WE=lambda a:sum(outside_costs[:min(a,7)])
  stage_max[c]={'entry':0,'consecutive':0,'outside':0,'helper':0,'boundary':0}
  for e in range(8):
   for mask in range(1,256):
    G=[b for i,b in enumerate(bb) if mask>>i&1];tmin=min(map(sum,G));tmax=max(map(sum,G));lo,hi=7+tmax,11+e+tmin
    for k in range(22+e):
     if not lo<=k<=hi:
      need(min(k,21+e-k)<=9,'outside common interior direct core differences at most18')
      need(sum(sorted(base,reverse=True)[:min(18,2*min(k,21+e-k))])+29<=76,'literal-derived global boundary cost')
      boundary+=1;continue
     interior+=1
     for b in G:
      t=sum(b);n=k-t;r0=max(7,n-e)
      need(7<=r0<=11 and 0<=n-r0<=e,'within-boundary migration remains in certified grades')
      for r in range(19):
       if not 0<=n-r<=e:continue
       if r<7:
        ss=r+7+(7-r);cost=W(r+7)+WE(7-r)
        # The actual seven-source cost multiset pays the complete entry.
        need(ss<=14 and cost<=40,'entry below band')
        stage_max[c]['entry']=max(stage_max[c]['entry'],cost)
       elif r>11:
        ss=(25-r)+(r-11);cost=W(25-r)+WE(r-11)
        need(ss<=14 and cost<=40,'entry above band')
        stage_max[c]['entry']=max(stage_max[c]['entry'],cost)
      for r in range(r0+1,12):
       if not 0<=n-r<e:continue
       ss=min(2*r-1,37-2*r)+1;cost=W(ss-1)+4
       need(ss<=18 and cost<=36+c,'consecutive grade step')
       stage_max[c]['consecutive']=max(stage_max[c]['consecutive'],cost)
      for tag,extra in [('outside',8),('helper',15)]:
       cost=W(2*min(r0,18-r0))+extra
       need(2*min(r0,18-r0)+2<=20 and cost<=49+c,'normalized skeleton event')
       stage_max[c][tag]=max(stage_max[c][tag],cost)
     # Endpoint tuples suffice: the compiler crosses directly, not along
     # a possibly disconnected one-bit graph of Gamma.
     for b in G:
      for bp in G:
       t,tp=sum(b),sum(bp);a=max(0,k-t-11,k-tp-11);z=min(e,k-t-7,k-tp-7)
       need(a<=z,'common outside count without changing an outside record')
       r,s=k-t-a,k-tp-a;m=hamming(b,bp);d=min(r+s,36-r-s)
       need((r+s)%2==m%2 and d+m<=20,'ternary parity support saving')
       cost=W(d)+2*m
       need(cost<=38+c,'exact ternary bridge cost upper bound')
       stage_max[c]['boundary']=max(stage_max[c]['boundary'],cost);bridge+=1
     count+=1
 return {'all_nonempty_ternary_tables':255,'common_interior_cases':interior,'direct_boundary_cases':boundary,'bridge_cases':bridge,'stage_maxima':stage_max}
def native_target_relation():
 receipt=json.loads((HERE/'native_receipt.json').read_text())['routing']['literal_route']
 domain=json.loads((HERE/'inputs/domains.json').read_text())['sources']
 targets={r['source']:{int(o['target']==2234):o['target'] for o in r['exact_all_shared_clocks_target_domain']} for r in domain}
 need(receipt['allowed_bit_patterns']==[0,7],'retained literal equality relation')
 for state in receipt['states']:
  need(len({targets[s][int(s in state['core_V'])] for s in receipt['boundary_sources']})==1,'literal equality of existing native target classes')
 # A receiver-INTERNAL equal-target relation has a different effect.
 labels={r['source']:(r['code'][0]<<4)|sum(v<<j for j,v in enumerate(r['code'][1:])) for r in domain}
 cheap=[r['source'] for r in domain if sum(a!=b for a,b in zip(*[x['path'] for x in r['exact_all_shared_clocks_target_domain']]))==2]
 add=lambda a,b:((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
 base=[s for s in cheap if s not in {176,1425,2071,1457}]
 need(len(base)==16,'actual sixteen unaffected native cheap records')
 tabs=[set() for _ in range(17)];tabs[0]={0}
 for i,z in enumerate(base):
  for r in range(i+1,0,-1):tabs[r]|={add(q,labels[z]) for q in tabs[r-1]}
 linked=add(labels[1841],labels[1457])
 holes={r:sorted(set(range(2048))-(tabs[r]|{add(linked,q) for q in tabs[r-2]})) for r in (8,9,10)}
 need(holes=={8:[477],9:[],10:[1098]},'actual receiver-internal same-target relation destroys old complete sections')
 return {'boundary_sources':receipt['boundary_sources'],'native_targets':[1965,2234],'checked_arrivals':len(receipt['states']),'new_owner_identifications':0,'receiver_internal_counterexample':{'frozen_cheap_sources':[176,1425,2071],'equal_target_sources':[1841,1457],'holes_by_receiver_grade':holes}}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 out={'interval_bridge':interval_tests(),'sharp_grade_distance':grade_distance(),'descent':descent_and_hostiles(),'native_budget':native_budget_tests(),'native_equal_target_relation':native_target_relation(),'checks':COUNT}
 out=json.loads(json.dumps(out,sort_keys=True));path=HERE/'theory_receipt.json'
 if args.write:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 else:need(out==json.loads(path.read_text()),'frozen theory receipt')
 print('PASS',COUNT,'exact conditions')
if __name__=='__main__':main()
