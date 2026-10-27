#!/usr/bin/env python3
"""Independent-agent audit; no imports from the producer packet, stdlib only.

Reads frozen inputs without changing them. All checks survive python -O.
Native scales are exactly those declared in the producer preregistrations.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import subprocess

MANIFEST = '09f5e8372a8133109cec03447284c49881440b36521eb22d2ef2ea7c67cf0ef7'
PIN = '96e35642110104cc547516a63f37eebcf930f707'

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def enc(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'))+'\n').encode()

def read(p):
    return json.loads(p.read_text())

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def Q(v):
    a,b = v
    return a*a+b*b-a*b

def mul(x,y):
    a,b=x; c,d=y
    return a*c-b*d, a*d+b*c-b*d

def B(v,w):
    a,b=v; c,d=w
    return (2*a-b)*c+(2*b-a)*d

def step(x):
    check(type(x) is int and x > 0 and x%2 and x%8 != 5, 'primitive route state')
    power = 2 if x%8 == 1 else 1
    y,r = divmod(3*x+1, 2**power)
    check(r == 0, 'initial division exact')
    braids=0
    while y%8 == 5:
        y,r=divmod(y-1,4)
        check(r == 0, 'braid division exact')
        braids+=1
    clock=power+2*braids
    check(2**clock*(3*y+1) == 9*x+3+2**power, 'affine clock identity')
    return y,clock

def reconstruct(N, frontier):
    """Uncached routes and Eisenstein polynomial multiplication, not source tables."""
    records=[]; trace=sha256(); peak=0; witness=None
    units=((1,0),(0,1),(-1,-1))
    for s in range(1,4*N,2):
        if s%8 == 5:
            continue
        a=b=0; z=s; chain=[]
        while z<4*N:
            factor=(1,2) if z<2*N else (0,1)
            c,d=mul(units[z%3],factor)
            a+=c; b+=d; chain.append(z); z=4*z+1
        A=[a,b]
        check(Q(A) in (0,1,3), 'native norm palette')
        z=s; seen=set(); t=0; route=[]
        while z not in frontier:
            check(z not in seen and len(route)<10000, 'unresolved route: audit incomplete, never erase mass')
            seen.add(z)
            if z>peak:
                peak=z; witness=s
            target,dt=step(z)
            route.append([z,target,dt]); t+=dt; z=target
        sign=1 if t%2==0 else -1
        rec={'s':s,'u':z,'p':route[-1][0] if route else 0,'t':t,'d':len(route),
             'A':A,'v':[sign*a,sign*b],'m':(s-1)//2}
        records.append(rec)
        trace.update(enc([rec,chain,route]))
    return records, {'uncached_route_Bchain_hash':trace.hexdigest(),'max_exterior_state':peak,'peak_source':witness}

def direct(rows,N,D,T,merge=False,address='m',unsigned=False):
    """Complete regrouping from leaves at every level; no energy updates."""
    selected=[r for r in rows if r['d']<=D and r['t']<=T]
    energies=[]; port_energies=[]; layers=[]
    q=1
    while True:
        bins={}; ports={}
        for r in selected:
            key=(0 if merge else r['u'],r[address]%q)
            a,b=bins.get(key,(0,0)); c,d=r['A'] if unsigned else r['v']
            bins[key]=a+c,b+d
            pk=(r['u'],r['p'],r[address]%q)
            a,b=ports.get(pk,(0,0)); ports[pk]=a+c,b+d
        energies.append(sum(map(Q,bins.values())))
        port_energies.append(sum(map(Q,ports.values())))
        layers.append(bins)
        if q>=2*N:
            break
        q*=2
    delta=[a-b for a,b in zip(energies,energies[1:])]
    same=[a-b for a,b in zip(port_energies,port_energies[1:])]
    parent_positive=0
    for k in range(len(layers)-1):
        q=2**k; fine=layers[k+1]
        for (u,r),v in layers[k].items():
            parent_positive+=max(0,Q(v)-Q(fine.get((u,r),(0,0)))-Q(fine.get((u,r+q),(0,0))))
    count=len(selected); norm=sum(Q(r['A']) for r in selected)
    if address == 'm':
        check(energies[-1]==norm, 'leaf endpoint')
    return {'energies':energies,'level_deltas':delta,'D_plus':sum(max(v,0) for v in delta),
            'arrived':count,'S_arrived':norm,'same_port_level_deltas':same,
            'cross_port_level_deltas':[x-y for x,y in zip(delta,same)],'per_parent':parent_positive}

def prefix(flat,nd,nt):
    # Two one-dimensional passes, not the producer's inclusion-exclusion recurrence.
    for i in range(nd):
        base=i*nt
        for j in range(1,nt):
            flat[base+j]+=flat[base+j-1]
    for i in range(1,nd):
        for j in range(nt):
            flat[i*nt+j]+=flat[(i-1)*nt+j]

def rectangle_audit(rows,N,receipt):
    ds=sorted(set(r['d'] for r in rows)); ts=sorted(set(r['t'] for r in rows))
    nd=len(ds); nt=len(ts); ix={d:i for i,d in enumerate(ds)}; iy={t:i for i,t in enumerate(ts)}
    L=(2*N-1).bit_length(); cells=nd*nt
    birth=[[0]*cells for _ in range(L)]; counts=[0]*cells; norms=[0]*cells
    groups=defaultdict(list)
    for r in rows:
        i=ix[r['d']]*nt+iy[r['t']]; counts[i]+=1; norms[i]+=Q(r['A'])
        groups[r['u']].append((r['m'],r['d'],r['t'],*r['v']))
    pairs=nonzero=0
    for group in groups.values():
        for x,y in combinations(group,2):
            m,d,t,a,b=x; n,e,u,c,f=y
            difference=m^n  # first differing bit equals v2(m-n), for distinct m,n.
            check(difference>0,'distinct named source addresses')
            level=(difference & -difference).bit_length()-1
            beta=(2*a-b)*c+(2*b-a)*f
            birth[level][ix[max(d,e)]*nt+iy[max(t,u)]]+=beta
            pairs+=1; nonzero+=bool(beta)
    prefix(counts,nd,nt); prefix(norms,nd,nt)
    for plane in birth:
        prefix(plane,nd,nt)
    surface=sha256(); peak=-1; maxcuts=[]
    for i,(D,T) in enumerate(product(ds,ts)):
        deltas=[plane[i] for plane in birth]
        debt=sum(max(x,0) for x in deltas)
        expected=[D,T,counts[i],norms[i],norms[i]+sum(deltas),debt]
        check(expected==receipt['rows'][i],f'all rectangle fields at N={N},D={D},T={T}')
        check(2*debt<=receipt['S_full'],'half-source target finite rectangle')
        if debt>peak:
            peak=debt; maxcuts=[[D,T]]
        elif debt==peak:
            maxcuts.append([D,T])
        surface.update(enc([expected,deltas]))
        if D==receipt['max_joint']['D'] and T==receipt['max_joint']['T']:
            check(deltas==receipt['max_joint']['level_deltas'],'full maximizer level-delta vector')
    check(maxcuts==receipt['maximizing_cutoff_pairs'],'all maximizers')
    return {'rectangles':cells,'same_terminal_pairs':pairs,'nonzero_pairs':nonzero,'surface_with_level_deltas_sha256':surface.hexdigest(),'maximum':peak}

def boundary_audit(rows,root):
    full=direct(rows,4096,67,64); kept=direct(rows,4096,36,64)
    selected=[r for r in rows if r['t']<=64]
    removed=[r for r in selected if r['d']>36]; removed_ids={r['s'] for r in removed}
    saved=read(root/'results/N4096_depth_boundary_stability_hostile.json')
    check(removed==saved['excluded_by_depth'],'every removed native label including all metadata')
    for key in ('energies','level_deltas','D_plus','arrived','S_arrived','same_port_level_deltas','cross_port_level_deltas'):
        check(full[key]==saved['clock_only'][key] and kept[key]==saved['joint'][key],'clipping literal paths '+key)
    pairs=[]; sums=Counter(); ports=Counter(); groups=defaultdict(list)
    for r in selected:
        groups[r['u']].append(r)
    for u,group in groups.items():
        for x,y in combinations(group,2):
            if x['s'] not in removed_ids and y['s'] not in removed_ids:
                continue
            v=x['m']^y['m']; q=v&-v; beta=B(x['v'],y['v'])
            pairs.append({'s1':x['s'],'s2':y['s'],'u':u,'p1':x['p'],'p2':y['p'],'q':q,'beta':beta,
                          'both_removed':x['s'] in removed_ids and y['s'] in removed_ids})
            sums[q]+=beta; ports[q,x['p'],y['p']]+=beta
    order=lambda p:(p['s1'],p['s2'])
    claimed=read(root/'results/boundary_invoice_receipt.json')['native']
    check(sorted(pairs,key=order)==sorted(claimed['pair_records'],key=order),'every boundary pair: owner, signs, terminals and both ports')
    check(len({(p['s1'],p['s2']) for p in pairs})==len(pairs),'no duplicate boundary pairs')
    check(Counter({(z['q'],z['p1'],z['p2']):z['beta'] for z in claimed['by_ordered_ports']})==ports,'all ordered port aggregates')
    delta=[x-y for x,y in zip(full['level_deltas'],kept['level_deltas'])]
    check(delta==[sums[2**k] for k in range(len(delta))]==claimed['delta_full_minus_retained'],'boundary pair low-bit ownership')
    invoice=sum(max(-x,0) for x in delta); gain=kept['D_plus']-full['D_plus']; norm=sum(Q(r['A']) for r in removed)
    check((len(removed),norm,full['D_plus'],kept['D_plus'],invoice,gain)==(16,22,439,495,82,56),'authentic clipping values')
    check(gain<=invoice and gain>norm and sum(Q(r['A']) for r in rows)==9949,'repair versus norm-only hostile, unchanged full denominator')
    return {'removed_labels':len(removed),'removed_norm':norm,'boundary_pairs':len(pairs),'nonzero_boundary_pairs':sum(p['beta']!=0 for p in pairs),
            'all_pair_metadata_checked':True,'delta':delta,'invoice':invoice,'gain':gain,'S_full':9949,'removed_s':[r['s'] for r in removed]}

def generic_audit(root):
    sizes=[]
    for saved in read(root/'results/boundary_invoice_receipt.json')['generic']['controls']:
        n=saved['N_leaves']; check(n>=4 and n&(n-1)==0,'generic dyadic size')
        full=[-1]+[1]*(n-1); kept=[0]+[1]*(n-1); ef=[]; ek=[]; q=1
        while q<=n:
            ef.append(sum(sum(full[r::q])**2 for r in range(q)))
            ek.append(sum(sum(kept[r::q])**2 for r in range(q)))
            check(ef[-1]==(n-4)*(n//q)+4 and ek[-1]==(n-2)*(n//q)+1,'generic energy formulas')
            q*=2
        df=[x-y for x,y in zip(ef,ef[1:])]; dk=[x-y for x,y in zip(ek,ek[1:])]
        invoice=sum(max(y-x,0) for x,y in zip(df,dk)); gain=sum(map(lambda x:max(x,0),dk))-sum(map(lambda x:max(x,0),df))
        check(ef==saved['full_energies'] and ek==saved['retained_energies'],'generic recorded energies')
        check(gain==invoice==2*(n-1)==saved['boundary_invoice']==saved['debt_increase'],'generic sharpness')
        sizes.append(n)
    check(sizes==[2**h for h in range(2,13)],'only declared generic controls')
    return {'scope':'NON-NATIVE scalar family; not a counterexample to R3','checked_leaf_sizes':sizes,'all_sharp':True}

def packet_audit(root,repo):
    check(digest(root/'MANIFEST.json')==MANIFEST,'frozen manifest binding')
    manifest=read(root/'MANIFEST.json'); actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    check(actual==set(manifest['files'])|{'MANIFEST.json'},'exact package membership')
    for p,h in manifest['files'].items():
        check(digest(root/p)==h,'manifest '+p)
    pins=read(root/'SOURCE_PINS.json'); check(pins['commit']==PIN,'source pin')
    pincheck=[]
    for s in pins['sources']:
        data=subprocess.check_output(['git','-C',str(repo),'show',PIN+':'+s['repository_path']])
        check(sha256(data).hexdigest()==s['sha256']==digest(root/s['local_path']),'repository pinned bytes '+s['local_path'])
        pincheck.append(s)
    declared=read(root/'PREREGISTRATION.json'); phase2=read(root/'PREREGISTRATION.phase2.json')
    Ns=sorted(declared['primary_N']+declared['structured_N']+declared['optional_extension_N']+phase2['additional_N'])
    check(len(Ns)==len(set(Ns))==19 and Ns==manifest['declared_N'],'all and only preregistered native scales')
    for p in manifest['files']:
        if p.endswith('.optimized.json'):
            check((root/p).read_bytes()==(root/p.replace('.optimized.json','.json')).read_bytes(),'frozen normal/-O equality '+p)
    return Ns,{'manifest_sha256':MANIFEST,'member_files':len(manifest['files']),'all_files_verified':True,'pinned_git_sources':pincheck}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--packet',type=Path,required=True);ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();root=args.packet
    Ns,binding=packet_audit(root,args.repo)
    ranks={int(a):int(b) for line in (root/'frontier.txt').read_text().splitlines() if line.strip() and not line.startswith('#') for a,b in [line.split()]}
    V=set(ranks); check(len(V)==627 and ranks[1]==0,'frontier grounded size')
    for s in V:
        z,_=step(s);check(z in V and (s==1 or ranks[z]<ranks[s]),'frontier closure and grounding')
    summary=read(root/'SUMMARY.json'); check(summary['declared_N']==Ns==summary['completed_N'] and not summary['pending_N'],'summary exact universe')
    report={'status':'PASS','binding':binding,'native_scales':Ns,'source_reconstruction':[],'all_rectangle_pair_audits':{},'direct_cutoff_audits':[]}
    total_rectangles=total_sources=0; ratios=[]
    for N in Ns:
        rows,route=reconstruct(N,V); saved=read(root/f'results/native_sources_N{N}.json'); rec=read(root/f'results/receipt_N{N}.json')
        check(rows==saved,'every independently reconstructed source field N'+str(N))
        check(sha256(enc(rows)).hexdigest()==rec['source_ledger_sha256'],'serialized source ledger')
        ds=sorted({r['d'] for r in rows});ts=sorted({r['t'] for r in rows}); S=sum(Q(r['A']) for r in rows)
        check(ds==rec['distinct_depths'] and ts==rec['distinct_clocks'],'all distinct global resource breaks')
        check([(r[0],r[1]) for r in rec['rows']]==list(product(ds,ts)),'complete Cartesian surface, ordering, no duplicate or missing rectangles')
        check(len(rec['rows'])==rec['rectangle_count']==len(ds)*len(ts),'rectangle count')
        check(S==rec['S_full'] and len(rows)==rec['sources'] and route['max_exterior_state']==rec['max_exterior_state'],'full denominator, labels, route maximum')
        check(rec['unresolved_sources']==0 and rec['empty_cutoff_D_plus']==0,'empty and unresolved policy')
        check(rec['preregistration_sha256']==digest(root/'PREREGISTRATION.json'),'preregistration receipt binding')
        if N in read(root/'PREREGISTRATION.phase2.json')['additional_N']:
            check(rec['additional_preregistration_sha256']==digest(root/'PREREGISTRATION.phase2.json'),'phase2 registration receipt binding')
        best=rec['max_joint']; check(best['D_plus']==max(r[-1] for r in rec['rows']),'saved global peak')
        check(rec['max_depth_only']==max((r for r in rec['rows'] if r[1]==ts[-1]),key=lambda x:x[-1]),'depth-only maximum')
        check(rec['max_clock_only']==max((r for r in rec['rows'] if r[0]==ds[-1]),key=lambda x:x[-1]),'clock-only maximum')
        check(2*best['D_plus']<=S and rec['full_denominator_margin']==S-2*best['D_plus'],'recorded finite bound')
        cuts={(best['D'],best['T'])}
        if N in (4096,18275,65536,262144):
            cuts|={(rec['max_depth_only'][0],ts[-1]),(ds[-1],rec['max_clock_only'][1]),(ds[-1],ts[-1]),(ds[len(ds)//2],ts[len(ts)//2]),(-1,ts[-1]),(ds[-1],-1)}
        cutdict={(r[0],r[1]):r for r in rec['rows']}
        for D,T in sorted(cuts):
            direct_result=direct(rows,N,D,T)
            if D<0 or T<0:
                check(direct_result['D_plus']==0 and direct_result['arrived']==0,'literal empty cutoffs')
            else:
                check([D,T,direct_result['arrived'],direct_result['S_arrived'],direct_result['energies'][0],direct_result['D_plus']]==cutdict[D,T],'direct selected rectangle')
            if (D,T)==(best['D'],best['T']):
                for key in ('energies','level_deltas','D_plus','arrived','S_arrived'):
                    check(direct_result[key]==best[key],'direct maximum '+key)
                lit=read(root/f'results/literal_receipt_N{N}.json')['max_joint_literal']
                check(direct_result['same_port_level_deltas']==lit['same_port_level_deltas'] and direct_result['cross_port_level_deltas']==lit['cross_port_level_deltas'],'port split maximum')
            report['direct_cutoff_audits'].append({'N':N,'D':D,'T':T,'D_plus':direct_result['D_plus'],'energies':direct_result['energies']})
        if N in (4096,18275,65536):
            report['all_rectangle_pair_audits'][str(N)]=rectangle_audit(rows,N,rec)
        if N==4096:
            check(rows==read(root/'input/native_sources_N4096.json'),'inherited authenticated N4096 ledger')
            report['native_boundary_invoice']=boundary_audit(rows,root)
            correct=direct(rows,N,36,64); merged=direct(rows,N,36,64,merge=True); unsigned=direct(rows,N,36,64,unsigned=True); address=direct(rows,N,36,64,address='s')
            check((correct['D_plus'],correct['per_parent'],merged['D_plus'],unsigned['D_plus'],address['D_plus'])==(495,11227,56715,209,495),'semantic hostile values')
            check(2*correct['per_parent']>S and 2*merged['D_plus']>S and address['energies']!=correct['energies'],'manufactured violations rejected; debt-only wrong-address blind spot')
            report['negative_controls']={'correct_debt':495,'positive_parent_false_debt':11227,'terminal_erasing_false_debt':56715,'lost_sign_false_debt':209,'wrong_address_debt':495,'wrong_address_energy_path_detected':True}
        total_sources+=len(rows);total_rectangles+=len(rec['rows']);ratios.append((Fraction(best['D_plus'],S),N))
        report['source_reconstruction'].append({'N':N,'records':len(rows),'S_full':S,'source_sha256':digest(root/f'results/native_sources_N{N}.json'),**route})
        print(json.dumps({'N':N,'sources':len(rows),'rectangles_structurally_checked':len(rec['rows']),'independent_pair_surface':N in (4096,18275,65536),'source_and_direct_checks':'PASS'}),flush=True)
    check((total_sources,total_rectangles)==(1000908,245826)==(summary['total_source_labels'],summary['total_rectangles']),'aggregate scope totals')
    check(max(ratios)==(Fraction(4132,44381),18275),'largest finite ratio')
    report.update(total_source_records=total_sources,structurally_checked_rectangles=total_rectangles,independently_recomputed_rectangles=sum(x['rectangles'] for x in report['all_rectangle_pair_audits'].values()),generic_sharpness=generic_audit(root),worst_ratio='4132/44381',worst_N=18275)
    check(report['source_reconstruction'][-1]['max_exterior_state']==28167930305>2**32-1,'exact large-field safeguard')
    args.out.write_bytes(enc(report))
    print('PASS independent audit; receipt '+str(args.out),flush=True)

if __name__=='__main__':
    main()
