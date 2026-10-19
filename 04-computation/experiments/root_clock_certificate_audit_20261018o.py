#!/usr/bin/env python3
"""Independent literal source/path replay of the bounded-clock certificate.

Checks every saved finite maximum, not a second exhaustive every-N sweep.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction
import json
from pathlib import Path

CHECKS=0
def check(x,why):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(why)
def norm(a,b):return a*a-a*b+b*b
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    root=Path(__file__).resolve().parents[2];res=root/'05-knowledge/results'
    tree=json.loads((res/'root_clock_tree_20261018o.json').read_text())
    receipt=json.loads((res/'root_clock_finite_certificate_20261018o.json').read_text())
    V=tree['ports'];vset=set(V);T=tree['clock'];routes=tree['routes'];source_set=set()
    for s,u,t,d,p in routes:
        check(s not in source_set,'source unique');source_set.add(s)
        x=s;clock=depth=last=0
        while x not in vset:
            last=x;x=3*x+1
            while x%2==0:x//=2;clock+=1
            while x%8==5:x=(x-1)//4;clock+=2
            depth+=1
            check(clock<=T,'literal replay clock budget')
        check((x,clock,depth,last)==(V[u],t,d,p),'all first arrivals and exterior ports')
    counts=Counter(u for s,u,t,d,p in routes);K2=sum(k*k for k in counts.values())
    bound=(8*K2+36)//5
    check(bound==tree['N_tail_cutoff']==receipt['N_end_exclusive'],'independent cutoff')
    check(K2==tree['terminal_fibre_square_sum']==receipt['tail_K2'],'fibre square sum')
    check(receipt['violations']==0 and receipt['target_factor']==1,'consumer accepted contractive target')
    check({tuple(p[:2]) for p in receipt['peaks']}=={(a,b) for a in range(T+1) for b in range(T+1)},'all clock/depth cap states')

    def source(s,N):
        row=[0,0,0]
        while s<4*N:
            for r in range(3):
                k=(s+r)%3
                row[r]+=(0,-3,3)[k] if s<2*N else (-1,-1,2)[k]
            s=4*s+1
        a=(row[0]-row[1])//3;b=-(row[0]+2*row[1])//3
        check([2*a-b,-a-b,-a+2*b]==row,'source inversion')
        return a,b

    cache={};peak_replays=[]
    for tc,dc,N,W,S in receipt['peaks']:
        if N not in cache:
            total=sum(norm(*source(s,N)) for s in range(1,4*N,2) if s%8!=5)
            cache[N]=(total,{s:source(s,N) for s in source_set})
        total,charges=cache[N];check(total==S,'complete literal source normalization')
        q=defaultdict(lambda:defaultdict(lambda:[0,0]))
        for s,u,t,d,p in routes:
            if t<=tc and d<=dc:
                a,b=charges[s];q[u][t][0]+=a;q[u][t][1]+=b
        D=0;lags=defaultdict(int)
        for poly in q.values():
            ts=sorted(poly)
            for i,t in enumerate(ts):
                a,b=poly[t];D+=norm(a,b)
                for v in ts[i+1:]:
                    c,e=poly[v];lags[v-t]+=2*a*c-a*e-b*c+2*b*e
        check(D+sum(abs(x) for x in lags.values())==W,'independent saved maximum replay')
        peak_replays.append([tc,dc,N,W,S])
    full=json.loads((res/'root_clock_channels_20261018o.json').read_text())
    sweep=json.loads((res/'root_wiener_scan_20261018o.json').read_text())
    samples={r['N']:r for r in sweep['samples']};shared=[]
    for r in full['rows']:
        if r['N'] not in samples:continue
        s=samples[r['N']]
        check((r['source_norm'],r['clock_diagonal'],r['frequency_norms']['pi'],r['wiener_budget'])==
              (s['source_norm'],s['clock_diagonal'],s['pi_norm'],s['wiener']), 'independent full-source engines')
        shared.append(r['N'])
    hostile=next(r for r in full['rows'] if r['N']==4096)
    z=(Fraction(-35,37),Fraction(12,37));zz=(Fraction(1),Fraction(0));value=Fraction(hostile['clock_diagonal'])
    lags={int(l):v for l,v in hostile['signed_lag_sums'].items()}
    for l in range(1,max(lags)+1):
        zz=(zz[0]*z[0]-zz[1]*z[1],zz[0]*z[1]+zz[1]*z[0]);value+=lags.get(l,0)*zz[0]
    h=hostile['exact_frequency_hostile']
    check(value==Fraction(int(h['spectral_numerator']),int(h['spectral_denominator'])),'rational frequency receipt')
    check(value>Fraction(1043,1000)*hostile['source_norm'],'strict actual-source hostile')
    args.out.write_text(json.dumps(dict(status='PASS',checks=CHECKS,literal_forward_routes=len(routes),
        all_saved_peak_replays=len(peak_replays),distinct_peak_scales=sorted(cache),
        independent_full_source_scales=shared,tail_cutoff=bound,
        scope='Independent paths and all saved maxima; not a second every-N scan'),indent=2)+'\n')
    print('PASS',CHECKS,'checks;',len(peak_replays),'maxima;',len(routes),'literal routes')
if __name__=='__main__':main()
