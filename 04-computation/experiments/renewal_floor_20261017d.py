#!/usr/bin/env python3
"""Join authenticated stages, retain old anchors, and locate unavoidable debt.
The lower bracket uses only already witnessed ROOT paths, never an oracle.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path
import collatz_interval_atlas_20261016e as A
from collatz_clock_tail_20261017b import tail
from collatz_good_lambda_20261017c import assess


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def controls():
    count=0
    for r,a,b in product(range(-2,3),repeat=3):
        lo,hi=r-abs(a)-abs(b),r+abs(a)+abs(b)
        for (d,e),(f,g) in product(product(range(3),(-1,1)),repeat=2):
            coeff=[r,0,0];coeff[d]+=e*a;coeff[f]+=g*b
            R=coeff[0];U=abs(coeff[1])+abs(coeff[2])
            need(lo<=R-U<=R+U<=hi,'every signed terminal merge preserves the interval')
            need(abs(R)+U<=abs(r)+abs(a)+abs(b),'all native pushforwards contract the upper bound')
            need(max(0,abs(R)-U)>=max(0,abs(r)-abs(a)-abs(b)),'lower certificate survives refinement')
            count+=1
    return count


def main(old_path,new_path,root_path):
    paths=[old_path,new_path,root_path];raw=[Path(p).read_bytes() for p in paths]
    old,new,root=[json.loads(r) for r in raw];o=old['runs'][0];n=new['runs'][0];r=root['runs'][0]
    need(o['k']==n['k']==r['k']==22,'same literal octave')
    need(o['stages'][-1]['audit']==n['stages'][0]['audit'],'all earlier stage19 clocks and endpoints reproduced')
    last=n['stages'][-1];rr=r['stages'][-1]
    for key in ('endpoint_sha256','clock_sha256'):
        need(last['audit'][key]==rr['audit'][key],'lower bracket uses the same fully audited stage27 sources')
    need(last['fine_tail']==rr['fine_tail'],'independent coefficient sweep reproduces the fine upper distribution')
    X=n['X'];anchors=[];rows=[];previous=None
    for st in o['stages'][:-1]+n['stages']:
        anchors+=min(st['selections'],key=lambda s:s['bin_width'])['anchors']
        q=tail(A.cone_envelope(X,anchors))
        row=dict(stage=st['stage'],anchors=len(anchors),energy=q['energy'],weak_numerator=q['weak_numerator'],
                 maximum=q['maximum'],critical_height=q['critical_height'])
        if previous is not None:
            row['half_debt']=assess(X,previous,q)['half_tail_debt']
        lower=rr['rooted_bracket']['lower_tail']
        row['unavoidable_future_half_debt']=assess(X,q,lower)['half_tail_debt']
        rows.append(row);previous=q
    return dict(status='FINITE-EXACT coverage and unavoidable future debt; all-scale bound OPEN',
                controls=controls(),inputs=[dict(path=p,sha256=hashlib.sha256(b).hexdigest()) for p,b in zip(paths,raw)],
                X=X,rows=rows,rooted_bracket=rr['rooted_bracket'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',required=True);p.add_argument('--new',required=True);p.add_argument('--root',required=True);a=p.parse_args()
    print(json.dumps(main(a.old,a.new,a.root),sort_keys=True,separators=(',',':')))
