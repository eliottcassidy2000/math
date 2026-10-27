#!/usr/bin/env python3
"""Fresh literal reconstruction and pair-birth audit. Does not import scan.py.
This is a second implementation by the producer, not an independent human/agent audit.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parent

def require(ok,label):
    if not ok: raise ValueError(label)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,sort_keys=True,separators=(',',':'))+'\n')
def norm(v):return v[0]**2-v[0]*v[1]+v[1]**2
def plus(x,y):return (x[0]+y[0],x[1]+y[1])
def cross(x,y):return norm(plus(x,y))-norm(x)-norm(y)

def reconstruct(N):
    raw=(ROOT/'frontier.txt').read_bytes()
    require(hashlib.sha256(raw).hexdigest()=='df88f03b576b13c98c0c5378fc1cf61685eeb6817073681266d2a781dc329b26','frontier binding')
    V={int(line.split()[0]) for line in raw.decode().splitlines() if line and not line.startswith('#')}
    out=[]; h=hashlib.sha256(); max_state=0
    for s in range(1,4*N,2):
        if s%8==5:continue
        x=s; real=[0,0,0]; chain=[]
        while x<4*N:
            w=(0,-3,3) if x<2*N else (-1,-1,2)
            term=[w[(x+r)%3] for r in range(3)]; real=[a+b for a,b in zip(real,term)]; chain.append([x,term]);x=4*x+1
        require((real[0]-real[1])%3==0 and (real[0]+2*real[1])%3==0,'integral source reconstruction')
        A=[(real[0]-real[1])//3,-(real[0]+2*real[1])//3]
        require(real==[2*A[0]-A[1],-A[0]-A[1],-A[0]+2*A[1]],'all real rows')
        require(sum(c*c for c in real)==6*norm(A),'three-phase norm identity')
        x=s;t=0;route=[];seen=set()
        while x not in V:
            require(x not in seen and len(route)<10000,'unresolved source route, not dropped');seen.add(x);max_state=max(max_state,x)
            a=1+(x%8==1); y,rem=divmod(3*x+1,2**a);require(rem==0,'native initial divisibility');j=0
            while y%8==5:
                y,rem=divmod(y-1,4);require(rem==0,'inverse braid divisibility');j+=1
            dt=a+2*j;require(2**dt*(3*y+1)==9*x+3+2**a,'exact affine identity')
            route.append([x,y,a,j,dt]);x=y;t+=dt
        r=dict(s=s,u=x,p=route[-1][0] if route else 0,t=t,d=len(route),A=A,v=[((-1)**t)*c for c in A],m=(s-1)//2)
        out.append(r);h.update((json.dumps(dict(r,B_chain=chain,exterior_route=route),sort_keys=True,separators=(',',':'))+'\n').encode())
    return out,h.hexdigest(),max_state

def ledger(rows,N,D,T):
    arrived=[r for r in rows if r['d']<=D and r['t']<=T]
    energies=[];port_energies=[];terminal_deltas=[];previous=None;L=(2*N-1).bit_length()
    for k in range(L+1):
        by=defaultdict(lambda:(0,0)); ports=defaultdict(lambda:(0,0)); per_u=defaultdict(int)
        for r in arrived:
            residue=r['m']%(2**k);key=(r['u'],residue);pk=(r['u'],r['p'],residue)
            by[key]=plus(by[key],r['v']);ports[pk]=plus(ports[pk],r['v'])
        for (u,_),v in by.items():per_u[u]+=norm(v)
        energies.append(sum(per_u.values()));port_energies.append(sum(norm(v) for v in ports.values()))
        if previous is not None: terminal_deltas.append({u:previous.get(u,0)-per_u.get(u,0) for u in set(previous)|set(per_u)})
        previous=per_u
    delta=[a-b for a,b in zip(energies,energies[1:])];same=[a-b for a,b in zip(port_energies,port_energies[1:])]
    require(sum(delta)==energies[0]-energies[-1],'telescoping')
    require(energies[-1]==sum(norm(r['A']) for r in arrived),'full leaf endpoint')
    terminal_positive=sum(max(0,x) for d in terminal_deltas for x in d.values())
    return dict(energies=energies,level_deltas=delta,D_plus=sum(max(0,x) for x in delta),arrived=len(arrived),S_arrived=energies[-1],
                same_port_level_deltas=same,cross_port_level_deltas=[x-y for x,y in zip(delta,same)],
                per_terminal_positive_part_debt=terminal_positive,
                positive_debt_terminal_contributors=[sorted([[u,v] for u,v in d.items() if v],key=lambda x:(-x[1],x[0]))[:12] for d in terminal_deltas])

def pair_birth_verify(rows,N,receipt):
    """Each same-terminal source pair enters at(max d,max t); prefix its signed beta."""
    depths=receipt['distinct_depths'];clocks=receipt['distinct_clocks'];di={d:i for i,d in enumerate(depths)};ti={t:i for i,t in enumerate(clocks)}
    L=(2*N-1).bit_length();nd=len(depths);nt=len(clocks)
    birth=[[[0]*nt for _ in range(nd)] for _ in range(L)];groups=defaultdict(list)
    for r in rows:groups[r['u']].append(r)
    pairs=0;nonzero=0
    for group in groups.values():
        for i,r in enumerate(group):
            for z in group[:i]:
                gap=abs(r['m']-z['m']);require(gap>0,'distinct source labels');k=(gap&-gap).bit_length()-1
                b=cross(r['v'],z['v']);birth[k][di[max(r['d'],z['d'])]][ti[max(r['t'],z['t'])]]+=b
                pairs+=1;nonzero+=bool(b)
    for mat in birth:
        for i in range(nd):
            for j in range(nt):
                mat[i][j]+=(mat[i-1][j] if i else 0)+(mat[i][j-1] if j else 0)-(mat[i-1][j-1] if i and j else 0)
    for D,T,count,Sarr,E0,debt in receipt['rows']:
        ds=[mat[di[D]][ti[T]] for mat in birth]
        require(sum(max(0,x) for x in ds)==debt,'all-rectangle pair-birth debt')
        require(sum(ds)==E0-Sarr,'all-rectangle pair-birth telescope')
    return dict(rectangles_checked=len(receipt['rows']),same_terminal_pairs=pairs,nonzero_signed_pairs=nonzero)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('N',type=int,nargs='+');ap.add_argument('--suffix',default='');args=ap.parse_args()
    for N in args.N:
        claimed=json.loads((ROOT/'results'/f'receipt_N{N}.json').read_text());rows,h,maxstate=reconstruct(N)
        stored=json.loads((ROOT/'results'/f'native_sources_N{N}.json').read_text());require(rows==stored,'complete source reconstruction equality')
        require(sum(norm(r['A']) for r in rows)==claimed['S_full'],'full source norm')
        maxrow=claimed['max_joint'];actual=ledger(rows,N,maxrow['D'],maxrow['T'])
        for k in ('energies','level_deltas','D_plus','arrived','S_arrived'):require(actual[k]==maxrow[k],'literal maximum '+k)
        cuts={(maxrow['D'],maxrow['T']),(claimed['max_depth_only'][0],max(claimed['distinct_clocks'])),(max(claimed['distinct_depths']),claimed['max_clock_only'][1])}
        for d in [claimed['distinct_depths'][len(claimed['distinct_depths'])*j//3] for j in range(3)]:
            for t in [claimed['distinct_clocks'][len(claimed['distinct_clocks'])*j//3] for j in range(3)]:cuts.add((d,t))
        rows_by_cut={(r[0],r[1]):r for r in claimed['rows']};samples=[]
        for D,T in sorted(cuts):
            check=ledger(rows,N,D,T);r=rows_by_cut[D,T]
            require([check['arrived'],check['S_arrived'],check['energies'][0],check['D_plus']]==r[2:],'literal selected control cutoff')
            samples.append(dict(D=D,T=T,D_plus=check['D_plus']))
        pairs=pair_birth_verify(rows,N,claimed) if N==4096 else None
        out=dict(status='PASS',N=N,source_records_compared=len(rows),detailed_source_routes_and_Bchains_sha256=h,max_exterior_state=maxstate,
                 source_ledger_sha256=digest(ROOT/'results'/f'native_sources_N{N}.json'),max_joint_literal=actual,cutoff_controls=samples,all_rectangle_pair_birth_check=pairs,
                 independence_scope='Fresh producer implementation; not independently authored agent audit')
        dump(ROOT/'results'/f'literal_receipt_N{N}{args.suffix}.json',out)
        print(json.dumps({k:out[k] for k in ('N','status','source_records_compared','all_rectangle_pair_birth_check')},sort_keys=True),flush=True)
if __name__=='__main__':main()
