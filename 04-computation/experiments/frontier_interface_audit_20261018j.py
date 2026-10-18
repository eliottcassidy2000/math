#!/usr/bin/env python3
"""Independent original-source checks at the critical transient, plus formal
norm identities for the geometric links. Imports none of the research engines.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sympy as sp


def primitive(x):
    j=0
    while x%8==5:x=(x-1)//4;j+=1
    return x,j


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    here=Path(__file__).parent
    V={int(l.split()[0]) for l in (here/'first_arrival_20261018/V_closed_Y1000.txt').read_text().splitlines()
       if l and not l.startswith('#')}
    rows=[];clock_checks=0
    for N in (1204,5997):
        initial=defaultdict(lambda:[0,0,0]);land=defaultdict(lambda:[0,0,0]);events=defaultdict(lambda:[0,0,0])
        for m in range(1,4*N,2):
            weights=[([0,-3,3] if m<2*N else [-1,-1,2])[(m+r)%3] for r in range(3)]
            s,j0=primitive(m)
            for r in range(3):initial[s][r]+=weights[r]
            # Evolve the literal uncompressed odd source. The compressed
            # route is checked, not used as a substitute for this evolution.
            x=m;t=0;canonical=0;depth=0
            while s not in V:
                assert depth<10000
                sx=3*s+1;a=0
                while sx%2==0:sx//=2;a+=1
                ns,nj=primitive(sx);canonical+=a+2*nj
                x=3*x+1
                while x%2==0:x//=2;t+=1
                ss,jlast=primitive(x)
                assert (ss,jlast)==(ns,nj)
                assert t==canonical+2*j0-2*jlast
                clock_checks+=1;depth+=1;s=ss
            for r in range(3):
                value=(-1)**t*weights[r]
                land[s][r]+=value;events[depth,s][r]+=value
        src=sum(c*c for v in initial.values() for c in v);final=sum(c*c for v in land.values() for c in v)
        cur=defaultdict(lambda:[0,0,0]);peak=0;pd=0
        for d in range(max(k[0] for k in events)+1):
            for (dd,s),v in events.items():
                if dd==d:
                    for r in range(3):cur[s][r]+=v[r]
            en=sum(c*c for v in cur.values() for c in v)
            if en>peak:peak=en;pd=d
        assert en==final
        rows.append(dict(N=N,source_trace=src,final_trace=final,max_prefix_trace=peak,peak_depth=pd,
                         final_contracts=final<=src,prefix_contracts=peak<=src))
    assert rows[0]['final_contracts'] is False
    assert rows[1]['peak_depth']==66 and rows[1]['prefix_contracts']
    # Independent engines must agree on every saved exact Gram entry.
    receipts=here/'frontier_interfaces_20261018j'
    py=json.loads((receipts/'root_phase_prefixes.json').read_text())
    cc=json.loads((receipts/'root_phase_everyN.json').read_text())
    cp=json.loads((receipts/'root_phase_prefix_everyN.json').read_text())
    crows={v['N']:v for v in cc['dyadic_checks']};gram_checks=0
    for v in py['rows']:
        w=crows[v['N']]
        for a,b in [('source_gram','source'),('arrival_gram','arrival'),('clock_diagonal_gram','clock')]:
            flat=[x for r in v[a] for x in r]
            assert flat==w[b];gram_checks+=len(flat)
    for v in cp['dyadic_checks']:assert v==crows[v['N']]
    assert gram_checks==513 and len(cp['dyadic_checks'])==20
    assert cc['violations']==[0,0,0] and cp['prefix_violations']==[0,0]
    assert (cp['prefix_ratio_N'],cp['prefix_ratio_depth'])==(5997,66)
    # Conjugates are substituted symbolically using |q|=|s|=1.
    k,kb,q,s=sp.symbols('k kb q s',nonzero=True)
    qp=(q-s/k)/(1-q/(s*kb))
    qpb=(1/q-1/(s*kb))/(1-s/(q*k))
    delta=s+k*(qp-q);deltab=1/s+kb*(qpb-1/q)
    assert sp.factor(qp*qpb-1)==0
    assert sp.factor(delta*deltab-1)==0
    assert sp.factor((qp+s/k)/(1+qp/(s*kb))-q)==0
    out=dict(status='FINITE-EXACT and symbolic identities',critical_sources=rows,
             literal_integer_clock_checks=clock_checks,independent_gram_entry_checks=gram_checks,
             matching_C_mode_dyadics=len(cp['dyadic_checks']),reflection_formal_identities=3,
             nondegeneracy='k*kb != 1 and nonzero displayed denominators')
    args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
