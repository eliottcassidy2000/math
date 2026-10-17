#!/usr/bin/env python3
"""Exact dyadic distribution contraction from authenticated saved controllers.
No new orbit data and no eventual ROOT oracle. Full receipt hashes retained.
"""
import argparse,hashlib,json
from fractions import Fraction
from itertools import product
from pathlib import Path
import collatz_interval_atlas_20261016e as A
from collatz_clock_tail_20261017b import tail
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def histogram_tail(hist,t):return sum(c for v,c in hist if v>=t)
def assess(X,old,new):
 hs=[1<<j for j in range(max(old['maximum'],new['maximum'],1).bit_length()+1)]
 rows=[];bad=[]
 for t in hs:
  u=histogram_tail(old['histogram'],t);v=histogram_tail(new['histogram'],2*t)
  rows.append([t,u,v,8*v<=u])
  if 8*v>u:bad.append(t)
 b=2*max(bad) if bad else 1
 ko=max(t*t*histogram_tail(old['histogram'],t) for t in hs)
 kn=max(t*t*histogram_tail(new['histogram'],t) for t in hs)
 need(kn<=max(4*b*b*X,Fraction(ko,2)),'good-lambda recurrence including integer threshold one')
 events=sorted(set(v for v,c in old['histogram'])|set(v for v,c in new['histogram']))
 oh=dict(old['histogram']);nh=dict(new['histogram']);u=v=X;debt=(0,0,0,0)
 for t in events:
  z=t*t*max(0,2*v-u)
  if z>debt[0]:debt=(z,t,u,v)
  u-=oh.get(t,0);v-=nh.get(t,0)
 need(u==v==0,'all histogram rows consumed without source renormalization')
 need(2*new['weak_numerator']<=old['weak_numerator']+debt[0],'same-height half-tail recurrence with exact unpaid debt')
 return dict(cutoff=b,cutoff_squared_over_X=str(Fraction(b*b,X)),old_dyadic_numerator=ko,new_dyadic_numerator=kn,thresholds=rows,
             half_tail_debt=dict(twice_numerator=debt[0],height=debt[1],old_rows=debt[2],new_rows=debt[3],normalized=str(Fraction(debt[0],2*X*X))))
def control():
 count=0
 for a in product((0,2,4),repeat=3):
  for b in product((0,2,4),repeat=3):
   if all(y<=x for x,y in zip(a,b)):
    assess(3,tail(a),tail(b));count+=1
 # One coherent constant coordinate contracts under a legal map (identity),
 # but has a macroscopic cutoff: local Lipschitz and l1 contraction are insufficient.
 bad=[]
 for X in (64,256,1024):
  q=tail([2*X]*X);r=assess(X,q,q)
  need(r['cutoff']==2*X,'constant-profile hostile forces a linear cutoff')
  bad.append(dict(X=X,cutoff=r['cutoff']))
 peaks=[]
 for m in (4,8,16):
  X=m**3;H=m*m;v=[max(H-4*j,0) for j in range(X)];q=tail(v)
  need(H**3<=32*q['weak_numerator'],'sharp-exponent peak hostile obeys universal peak bound')
  need(q['weak_numerator']<=6*X*X,'peak hostile has bounded critical tail')
  peaks.append(dict(m=m,X=X,maximum=H,energy=q['energy'],weak_numerator=q['weak_numerator']))
 return dict(exhaustive_monotone_pairs=count,constant_hostile=bad,peak_hostiles=peaks)
def run(path):
 raw=Path(path).read_bytes();d=json.loads(raw);out=[]
 for r in d['runs']:
  X=r['X'];previous=None;previous_profile=None;allanchors=[];stages=[]
  for st in r['stages']:
   sel=min(st['selections'],key=lambda q:q['bin_width']);allanchors+=sel['anchors']
   profile=A.cone_envelope(X,allanchors);q=tail(profile)
   need(q['energy']<=sel['cone_tail']['energy'],'retained older certificates cannot worsen the newest cone')
   stage=dict(stage=st['stage'],anchors=len(allanchors),energy=q['energy'],weak_numerator=q['weak_numerator'],maximum=q['maximum'])
   if previous is not None:
    need(all(v<=u for u,v in zip(previous_profile,profile)),'every literal row retains every earlier certificate')
    a=assess(X,previous,q);t=a['half_tail_debt']['height'];pieces=[];start=None;survivors=0
    for i,(u,v) in enumerate(zip(previous_profile+[0],profile+[0])):
     if i<X and u>=t and t>0:
      if start is None:start=i;survivors=0
      survivors+=v>=t
     elif start is not None:
      pieces.append([X+start,X+i,survivors,i-start-survivors]);start=None
    need(sum(s-e for l,h,s,e in pieces)*t*t==a['half_tail_debt']['twice_numerator'],'critical debt is signed survivors minus exits on the original row components')
    a['critical_components']=pieces
    a['component_positive_debt_twice_numerator']=t*t*sum(max(s-e,0) for l,h,s,e in pieces)
    stage['contraction']=a
   stages.append(stage);previous=q;previous_profile=profile
  out.append(dict(k=r['k'],X=X,stages=stages))
 return dict(source_sha256=hashlib.sha256(raw).hexdigest(),runs=out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);a=p.parse_args();ctl=control();out=[run(f) for f in a.input]
 print(json.dumps(dict(status='FINITE-EXACT dyadic good-lambda controls; uniform square-root cutoff OPEN',controls=ctl,receipts=out,checks=checks),sort_keys=True,separators=(',',':')))
