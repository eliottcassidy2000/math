#!/usr/bin/env python3
"""Exact boundary-pair invoice and sharp generic self-energy obstruction."""
from collections import defaultdict
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parent

def need(ok,s):
    if not ok:raise RuntimeError(s)
def Q(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def debt(es):return sum(max(x-y,0) for x,y in zip(es,es[1:]))
def scalar_energy(xs,q):return sum(sum(xs[r::q])**2 for r in range(q))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--suffix',default='');args=ap.parse_args()
    elementary=0
    for full in range(-40,41):
        for delta in range(-40,41):
            need(max(full-delta,0)-max(full,0)<=max(-delta,0),'one-level positive-part boundary repair')
            elementary+=1
    source_path=ROOT/'results/native_sources_N4096.json'
    source_hash=hashlib.sha256(source_path.read_bytes()).hexdigest()
    need(source_hash=='4ed6130bd0b368cab2af3252958b940f3d056749bf72f1d85784745914d044ba','native full source hash')
    rows=json.loads(source_path.read_text())
    native=json.loads((ROOT/'results/N4096_depth_boundary_stability_hostile.json').read_text())
    full=[r for r in rows if r['t']<=64];R={r['s'] for r in full if r['d']>36}
    groups=defaultdict(list)
    for r in full:groups[r['u']].append(r)
    boundary=defaultdict(int);by_ports=defaultdict(int);pair_records=[];pair_count=0;nonzero_count=0
    for u,group in groups.items():
        for j,r in enumerate(group):
            for z in group[:j]:
                if r['s'] not in R and z['s'] not in R:continue
                gap=abs(r['m']-z['m']);q=gap&-gap;b=beta(r['v'],z['v'])
                boundary[q]+=b;by_ports[(q,r['p'],z['p'])]+=b;pair_count+=1;nonzero_count+=bool(b)
                pair_records.append(dict(s1=z['s'],s2=r['s'],u=u,p1=z['p'],p2=r['p'],q=q,beta=b,both_removed=r['s'] in R and z['s'] in R))
    delta=[x-y for x,y in zip(native['clock_only']['level_deltas'],native['joint']['level_deltas'])]
    need(delta==[boundary[1<<k] for k in range(len(delta))],'terminal/LCA boundary ownership equals full-minus-retained level deltas')
    invoice=sum(max(-v,0) for v in delta);gain=native['debt_increase']
    need((invoice,gain,native['excluded_norm'])==(82,56,22),'authentic boundary invoice values')
    need(gain<=invoice,'native repaired inequality')
    generic=[]
    for h in range(2,13):
        n=1<<h;full=[-1]+[1]*(n-1);retained=[0]+[1]*(n-1)
        qs=[1<<k for k in range(h+1)];ef=[scalar_energy(full,q) for q in qs];er=[scalar_energy(retained,q) for q in qs]
        need(ef==[(n-4)*(n//q)+4 for q in qs],'generic full energy formula')
        need(er==[(n-2)*(n//q)+1 for q in qs],'generic retained energy formula')
        df=[x-y for x,y in zip(ef,ef[1:])];dr=[x-y for x,y in zip(er,er[1:])]
        dd=[x-y for x,y in zip(df,dr)]
        need(all(x>=0 for x in df+dr),'generic nonnegative whole-level drops')
        need(dd==[-n//q for q in qs[:-1]],'generic negative boundary delta formula')
        cost=sum(max(-x,0) for x in dd)
        need(debt(er)-debt(ef)==cost==2*(n-1),'sharp boundary repair and unbounded removed-self-norm ratio')
        generic.append(dict(N_leaves=n,removed_norm=1,full_energies=ef,retained_energies=er,full_debt=debt(ef),retained_debt=debt(er),boundary_invoice=cost,debt_increase=2*(n-1)))
    out=dict(status='PASS',source_ledger_sha256=source_hash,elementary_integer_pairs_checked=elementary,
             native=dict(N=4096,T=64,full_D=67,retained_D=36,S_full=9949,removed_sources=len(R),removed_norm=22,full_debt=439,retained_debt=495,debt_increase=56,
                         delta_full_minus_retained=delta,boundary_invoice=invoice,slack=invoice-gain,boundary_pairs=pair_count,nonzero_boundary_pairs=nonzero_count,pair_records=pair_records,
                         by_ordered_ports=[dict(q=q,p1=p1,p2=p2,beta=b) for (q,p2,p1),b in sorted(by_ports.items())]),
             generic=dict(scope='Generic scalar balanced low-bit tree only; not an authentic A_N source family.',formula='N=2^d>=4; one -1 and N-1 copies +1. Zero the -1 but retain its conceptual leaf. E_full(a)=(N-4)a+4; E_ret(a)=(N-2)a+1. Debt increase and exact boundary invoice both2(N-1), removed norm1.',controls=generic))
    (ROOT/'results'/f'boundary_invoice_receipt{args.suffix}.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps(dict(status='PASS',native_invoice=invoice,native_gain=gain,native_removed_norm=22,generic_sizes=[z['N_leaves'] for z in generic]),sort_keys=True))
if __name__=='__main__':main()
