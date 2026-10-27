#!/usr/bin/env python3
"""Exact source-preserving joint depth/integer-clock cutoff scan (stdlib only)."""
import argparse, hashlib, json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PIN='96e35642110104cc547516a63f37eebcf930f707'
LOWER=((1,2),(-2,-1),(1,-1)); UPPER=((0,1),(-1,-1),(1,0))

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,obj): path.write_text(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n')
def Q(v): a,b=v; return a*a-a*b+b*b
def edge(s):
    a=2 if s%8==1 else 1
    y=(3*s+1)//(1<<a); t=a
    while y%8==5: y=(y-1)//4; t+=2
    need((1<<t)*(3*y+1)==3*(3*s+1)+(1<<a),'exact affine clock')
    return y,t

def charge(s,N):
    a=b=0; x=s
    while x<4*N:
        c,d=(LOWER if x<2*N else UPPER)[x%3]
        a+=c; b+=d; x=4*x+1
    need(Q((a,b)) in (0,1,3),'native palette')
    return [a,b]

def frontier():
    pins=json.loads((ROOT/'SOURCE_PINS.json').read_text())
    need(pins['commit']==PIN,'source commit')
    for z in pins['sources']: need(sha(ROOT/z['local_path'])==z['sha256'],'source file '+z['local_path'])
    V={int(x.split()[0]) for x in (ROOT/'frontier.txt').read_text().splitlines() if x.strip() and not x.startswith('#')}
    need(len(V)==627,'frontier cardinality')
    need(all(edge(s)[0] in V for s in V),'frontier closure')
    return V

def records(N,V):
    memo={s:(s,0,0,0) for s in V}; out=[]; max_state=0
    for s in range(1,4*N,2):
        if s%8==5: continue
        x=s; path=[]; seen=set()
        while x not in memo:
            need(x not in seen and len(path)<10000,'unresolved route; no incomplete source may be omitted')
            seen.add(x); max_state=max(max_state,x)
            y,dt=edge(x); path.append((x,y,dt)); x=y
        u,p,t,d=memo[x]
        for x,y,dt in reversed(path):
            if d==0: p=x
            d+=1;t+=dt;memo[x]=(u,p,t,d)
        u,p,t,d=memo[s]; A=charge(s,N); sign=-1 if t%2 else 1
        out.append(dict(s=s,u=u,p=p,t=t,d=d,A=A,v=[sign*A[0],sign*A[1]],m=(s-1)//2))
    need(len(out)==sum(s%8!=5 for s in range(1,4*N,2)),'source family complete')
    return out,max_state

def literal(rows,N,D,T,merge_terminals=False,address='m',per_parent=False,absolute=False):
    eligible=[r for r in rows if r['d']<=D and r['t']<=T]
    es=[]; levels=[]; q=1
    while True:
        bins=defaultdict(lambda:[0,0])
        for r in eligible:
            v=r['A'] if absolute else r['v']; key=(0 if merge_terminals else r['u'],r[address]%q)
            bins[key][0]+=v[0]; bins[key][1]+=v[1]
        es.append(sum(Q(v) for v in bins.values()));levels.append(bins)
        if q>=2*N: break
        q*=2
    delta=[a-b for a,b in zip(es,es[1:])]
    pp=0
    if per_parent:
        for k,(coarse,fine) in enumerate(zip(levels,levels[1:])):
            q=1<<k
            for (u,r),v in coarse.items():
                pp+=max(Q(v)-Q(fine.get((u,r),(0,0)))-Q(fine.get((u,r+q),(0,0))),0)
    return dict(energies=es,level_deltas=delta,D_plus=sum(max(0,x) for x in delta),
                arrived=len(eligible),S_arrived=sum(Q(r['A']) for r in eligible),per_parent=pp)

def sweep(rows,N):
    depths=sorted({r['d'] for r in rows}); clocks=sorted({r['t'] for r in rows})
    clock_groups=defaultdict(list)
    for r in rows: clock_groups[r['t']].append(r)
    S=sum(Q(r['A']) for r in rows)
    L=(2*N-1).bit_length(); masks=[(1<<k)-1 for k in range(L+1)]
    best=None; all_max=[]; controls_depth=[]; all_rows=[]; clocks_control=[]
    for D in depths:
        bins=[{} for _ in masks]; energies=[0]*len(masks); count=norm=0
        for T in clocks:
            for r in clock_groups[T]:
                if r['d']>D: continue
                a,b=r['v']; qv=Q((a,b)); m=r['m'];u=r['u'];count+=1;norm+=qv
                if not qv: continue
                for j,mask in enumerate(masks):
                    key=(u,m&mask); x,y=bins[j].get(key,(0,0))
                    energies[j]+=qv+(2*x-y)*a+(2*y-x)*b
                    bins[j][key]=(x+a,y+b)
            deltas=[a-b for a,b in zip(energies,energies[1:])]; debt=sum(max(0,x) for x in deltas)
            need(energies[-1]==norm,'last level source norm')
            entry=[D,T,count,norm,energies[0],debt]
            all_rows.append(entry)
            if best is None or debt>best['D_plus']:
                best=dict(D=D,T=T,D_plus=debt,arrived=count,S_arrived=norm,energies=list(energies),level_deltas=deltas)
                all_max=[[D,T]]
            elif debt==best['D_plus']: all_max.append([D,T])
            if T==clocks[-1]: controls_depth.append(entry)
            if D==depths[-1]: clocks_control.append(entry)
    direct=literal(rows,N,best['D'],best['T'],per_parent=True)
    for k in ('energies','level_deltas','D_plus','arrived','S_arrived'): need(best[k]==direct[k],'literal maximum '+k)
    depth_max=max(controls_depth,key=lambda x:x[-1]);clock_max=max(clocks_control,key=lambda x:x[-1])
    return dict(N=N,S_full=S,sources=len(rows),zero_charge_sources=sum(Q(r['A'])==0 for r in rows),
                terminal_count=len({r['u'] for r in rows}),entry_port_count=len({(r['u'],r['p']) for r in rows}),
                distinct_depths=depths,distinct_clocks=clocks,rectangle_count=len(all_rows),
                empty_cutoff_D_plus=0,summary_columns=['D','T','arrived','S_arrived','terminal_energy','D_plus'],
                max_joint=best,maximizing_cutoff_pairs=all_max,max_joint_2D_over_S=str(Fraction(2*best['D_plus'],S)),
                full_denominator_margin=S-2*best['D_plus'],max_depth_only=depth_max,max_clock_only=clock_max,
                holds_all_rectangles=2*best['D_plus']<=S,rows=all_rows,
                hostile_controls=dict(correct=direct,
                    merged_terminals=literal(rows,N,best['D'],best['T'],merge_terminals=True),
                    lost_clock_sign=literal(rows,N,best['D'],best['T'],absolute=True),
                    wrong_low_bit_address=literal(rows,N,best['D'],best['T'],address='s')))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('N',type=int,nargs='+');ap.add_argument('--suffix',default='');args=ap.parse_args()
    prereg=json.loads((ROOT/'PREREGISTRATION.json').read_text()); allowed=set(prereg['primary_N']+prereg['structured_N']+prereg['optional_extension_N'])
    phase2_path=ROOT/'PREREGISTRATION.phase2.json'
    additional=set(json.loads(phase2_path.read_text())['additional_N']) if phase2_path.exists() else set()
    need(set(args.N)<=allowed|additional,'N must be preregistered')
    V=frontier()
    for N in args.N:
        rows,max_state=records(N,V)
        if N==4096:
            need(sha(ROOT/'input/native_sources_N4096.json')=='4ed6130bd0b368cab2af3252958b940f3d056749bf72f1d85784745914d044ba','authenticated source ledger hash')
            need(sha(ROOT/'input/producer_revision2_manifest.json')=='ad8d60b20884230f7bd2180e22115b16f970d4cd0aa284a1149d45cd51f321e8','revision2 manifest hash')
            need(rows==json.loads((ROOT/'input/native_sources_N4096.json').read_text()),'complete authenticated N4096 ledger equality')
        path=ROOT/'results'/f'native_sources_N{N}.json';save(path,rows)
        out=sweep(rows,N);out.update(status='PASS',source_commit=PIN,preregistration_sha256=sha(ROOT/'PREREGISTRATION.json'),source_ledger_sha256=sha(path),max_exterior_state=max_state,unresolved_sources=0,route_arithmetic='Python arbitrary precision')
        if N in additional: out['additional_preregistration_sha256']=sha(phase2_path)
        save(ROOT/'results'/f'receipt_N{N}{args.suffix}.json',out)
        print(json.dumps({k:out[k] for k in ('N','sources','S_full','rectangle_count','max_joint','max_joint_2D_over_S','max_depth_only','max_clock_only','holds_all_rectangles')},sort_keys=True),flush=True)
if __name__=='__main__':main()
