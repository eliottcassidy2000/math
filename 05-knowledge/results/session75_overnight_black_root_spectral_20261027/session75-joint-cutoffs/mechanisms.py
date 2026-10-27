#!/usr/bin/env python3
"""Regenerate exact native mechanisms from complete source families."""
from collections import Counter
from pathlib import Path
import argparse,json
import literal_verify as a
ROOT=Path(__file__).resolve().parent

def save(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def need(ok,s):
    if not ok:raise RuntimeError(s)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--suffix',default='');args=ap.parse_args()
    rows=json.loads((ROOT/'results/native_sources_N4096.json').read_text())
    before=a.ledger(rows,4096,36,64);after=a.ledger(rows,4096,36,66);z=next(r for r in rows if r['s']==13759)
    additions=[r for r in rows if r['d']<=36 and 64<r['t']<=66]
    need(additions==[z],'unique boundary arrival')
    interactions=[]
    for r in rows:
        if r['u']==z['u'] and r['d']<=36 and r['t']<=64:
            gap=abs(r['m']-z['m']);q=gap&-gap;beta=a.cross(r['v'],z['v'])
            interactions.append(dict(q=q,beta=beta,source=r))
    nonzero=[x for x in interactions if x['beta']]
    need([(x['source']['s'],x['q'],x['beta']) for x in nonzero]==[(9,1,-1),(2831,8,1),(10065,1,-1)],'exact three interaction mechanism')
    difference=[b-c for c,b in zip(before['level_deltas'],after['level_deltas'])]
    need(difference==[-2,0,0,1]+[0]*9 and (before['D_plus'],after['D_plus'])==(495,494),'whole-level debt change')
    out=dict(N=4096,S_full=sum(a.norm(r['A']) for r in rows),D=36,before_T=64,after_T=66,added_source=z,before=before,after=after,
             same_terminal_pair_contributions_including_zero=interactions,same_terminal_nonzero_pair_contributions=nonzero)
    save(ROOT/'results'/f'N4096_clock_border_mechanism{args.suffix}.json',out)
    clock_only=a.ledger(rows,4096,67,64)
    excluded=[r for r in rows if r['t']<=64 and r['d']>36]
    norm_excluded=sum(a.norm(r['A']) for r in excluded)
    need(len(excluded)==16 and norm_excluded==22 and clock_only['D_plus']==439,'exact depth-clipping boundary')
    need(before['D_plus']-clock_only['D_plus']==56>norm_excluded,'native source-norm-only stability hostile')
    boundary=dict(N=4096,S_full=9949,T=64,before_D=67,after_D=36,clock_only=clock_only,joint=before,excluded_by_depth=excluded,
        excluded_norm=norm_excluded,debt_increase=56,level_delta_changes=[x-y for x,y in zip(before['level_deltas'],clock_only['level_deltas'])],
        scope='Post-hoc exact counterexample only to D(after restriction)<=D(before)+norm(excluded); not to R3 or joint S_full/2 target.')
    save(ROOT/'results'/f'N4096_depth_boundary_stability_hostile{args.suffix}.json',boundary)
    reports=[]
    for N in (4096,18275,65536):
        r=json.loads((ROOT/'results'/f'native_sources_N{N}.json').read_text());x=json.loads((ROOT/'results'/f'receipt_N{N}.json').read_text());D=x['max_joint']['D'];T=x['max_joint']['T'];counts=Counter();norms=Counter()
        for z in r:
            typ=('inside' if z['d']<=D and z['t']<=T else 'clock_only_excluded' if z['d']<=D else 'depth_only_excluded' if z['t']<=T else 'both_excluded')
            counts[typ]+=1;norms[typ]+=a.norm(z['A'])
        need(sum(counts.values())==x['sources'] and sum(norms.values())==x['S_full'],'complete exclusion partition')
        y=json.loads((ROOT/'results'/f'literal_receipt_N{N}.json').read_text())['max_joint_literal']
        reports.append(dict(N=N,D=D,T=T,counts=dict(counts),norms=dict(norms),same_port=y['same_port_level_deltas'],cross_port=y['cross_port_level_deltas'],per_terminal_positive_debt=y['per_terminal_positive_part_debt']))
    save(ROOT/'results'/f'resource_port_diagnostics{args.suffix}.json',dict(status='PASS',reports=reports))
    print('PASS: native clock-border and source-resolved resource/port diagnostics')
if __name__=='__main__':main()
