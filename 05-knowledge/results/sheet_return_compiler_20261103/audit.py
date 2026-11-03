"""Replay the new certificates; inherited results are pinned, not re-proved."""
import hashlib,json
from bulk import *
from phase import run as phase_run,mpower
from singer import run as singer_run
from switch_search import run as switch_run,hpaths
from return_program import bin_values,field,induced,opposite
from sympy import Matrix


def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def run():
    from incoming import run as incoming_run
    incoming_run()
    literal=[]
    for kind in ('straight','revisit'):
        for n in range(2,7):
            T,ad,pp=paths(3,n);p=dict(pp)[kind]
            need(len(p)==len(T) and len(set(p))==len(T) and all(T[a]>>b&1 for a,b in zip(p,p[1:])),'actual path legality')
            U=native.native(T,p)
            need(native.native(opposite(T),p[::-1])==opposite(U),'native converse transport')
            V,ad=target(3,n,kind);h=tables(V)[0][-1]
            hh=[1,1,1]
            while len(hh)<=n+1:hh.append(sum(hh[-3:]))
            expected=97*hh[n+1] if kind=='straight' else 54*hh[n+1]+62*hh[n]+31*hh[n-1]
            need(h==expected and strong(V),'formula and literal strongness')
            P,E,S=bin_values(V,ad.index((0,2)),2)
            obs=(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2)
            if n in (2,3):
                hp=list(hpaths(V));need(len(hp)==h,'independent permutation-tree count')
                cycle=next((p for p in hp if V[p[-1]]>>p[0]&1),None)
                need(cycle is not None,'base strong cycle')
            else:cycle=None
            literal.append(dict(kind=kind,n=n,H=h,observer=obs,phase=field(obs),cycle=cycle))
            tables.cache_clear()
            if n>=3:
                W,wad=target(3,n+1,kind)
                def address(a):return (5,a[1]+1) if a[0]==5 and a[1]>=1 else a
                order=[wad.index(address(a)) for a in ad]
                pulled=induced(W,order);a=ad.index((5,0));b=ad.index((5,1))
                flipped=list(V);flipped[a]^=1<<b;flipped[b]^=1<<a
                need(pulled==tuple(flipped),'inductive cycle insertion changes exactly the old consecutive edge')
                v=wad.index((5,1));need(W[order[a]]>>v&1 and W[v]>>order[b]&1,'new two-edge cycle segment')
    # Exact literal identification with the previous canonical 161 escape.
    from branch_search import canon_order
    T,_,_=paths(3,3);_,order=canon_order(T)
    old=json.loads((HERE.parent/'defect_offset_dynamics_20261103/native_escape_161.json').read_text())
    for kind,h in (('straight',485),('revisit',487)):
        V,_=target(3,3,kind)
        need(induced(V,tuple(order)+(10,))==tuple(old[str(h)]['target']),'old source labels and root transported exactly')
    bulk_records=[audit(blueprint(3,k)) for k in ('straight','revisit')]
    (HERE/'bulk.json').write_text(json.dumps(bulk_records,indent=2)+'\n')
    phases=phase_run();singer=singer_run()
    for record in phases['records']:
        hseq=record['H_initial'][:]
        for _ in range(48):hseq.append(sum(hseq[-3:]))
        zeros=sorted({(i+4)%16 for i,h in enumerate(hseq[:48]) if h%7==0})
        need(zeros==([8] if record['kind']=='straight' else [1,4,5,10]),'complete H divisibility phases')
        C=[[0,1,0],[0,0,1],[1,1,1]]
        need(mpower(C,16,7)==[[4*int(i==j) for j in range(3)] for i in range(3)],'projective16 return')
        state=[x%7 for x in record['H_initial']];first=None;projective=[]
        for i in range(48):
            inv=pow(next(x for x in state if x),-1,7);normal=tuple(inv*x%7 for x in state)
            projective.append(normal);state=[state[1],state[2],sum(state)%7]
        need(len(set(projective))==16,'actual H projective period16')
        # Crosscheck the low-size joint response against the proven recurrence.
        from phase import recurrence
        sseq=record['S_initial'][:]
        s3=(sseq[2]-2*sseq[1]-sseq[0])//2
        s2=(sseq[1]-2*sseq[0]-s3)//2
        need([row['observer'][1]+row['observer'][2] for row in literal if row['kind']==record['kind']][:2]==[s2,s3],'joint bases n2,n3')
    switches=[switch_run(m) for m in (1,2,3)]
    need([r['shared_paths'] for r in switches]==[62,104,188] and all(not r['paired'] for r in switches),'complete common-program switch obstruction')
    (HERE/'switch_search.json').write_text(json.dumps(switches,indent=2)+'\n')
    # Same H/current F49 phase need not identify an admissible attachment source.
    bad=(250,40,51,52,2,16,62,126)
    good,_,_=paths(3,1)
    control=[]
    for V in (bad,good):
        P,E,S=bin_values(V,7,2);v=(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2)
        control.append(dict(observer=v,phase=field(v)));tables.cache_clear()
    need(control[0]['phase']==control[1]['phase']==(2,2,2) and control[0]['observer']!=control[1]['observer'],'inherited current-phase hostile')
    from next_escapes import profile
    continuations=json.loads((HERE/'next_escapes.json').read_text());next_checks=[]
    for row in continuations:
        T=tuple(row['source']);h=row['H']
        need(set(map(int,row['witnesses']))=={3*h+2,3*h+4},'both actual next offsets')
        for count,w in row['witnesses'].items():
            p=w['path'];U=T if p is None else native.native(T,tuple(p))
            V=birth(U,w['mask'])
            need(V==tuple(w['target']) and strong(V),'literal continuation and strong target')
            need(tables(V)[0][-1]==int(count)==sum(1 for _ in hpaths(V)),'independent next-target path count')
            need(profile(V)==dict(w['profile'],observer=tuple(w['profile']['observer']),phase=tuple(w['profile']['phase'])),'retained next-target root and joint defects')
            need(0 in w['profile']['root_defects'],'each child retains a zero guard')
            next_checks.append(dict(source=h,target=int(count),root_phase=w['profile']['phase']))
            tables.cache_clear()
    direct=continuations[1]['witnesses']['1465'];V=tuple(direct['target'])
    repeated=birth(V,direct['mask']|(1<<11))
    need(tables(repeated)[0][-1]==3745!=3*1465+4,'repeating the same boundary ear is not a return')
    tables.cache_clear()
    dependencies=[
      'defect_offset_dynamics_20261103/wing.py','defect_offset_dynamics_20261103/ESCAPES.md',
      'defect_offset_dynamics_20261103/native_escape_161.json',
      'overnight_rank_return_20261029/RUN_RECURRENCE.md','overnight_rank_return_20261029/defect_transfer.py',
      'overnight_rank_return_20261029/defect_transfer.json',
      'session35_dimensions_dualities_20261023/arithmetic/SINGER_LINE_FLIPS.md',
      'session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md',
      'session167_recurring_structures_20261103/norm/NOTE.md',
      'six_roots_hidden_bits_20261103/NATIVE_BRANCHES.md',
      'joint_port_geometry_20261103/JOINT_EFFECTS.md',
      'golden_inverse_defects_20261103/inverse.py','ready_return_pentagon_20261103/return_program.py']
    out=dict(literal=literal,phase_hostile=control,next_escapes=next_checks,
             inherited_81_27_scope='unchanged twelve-vertex two/three-run separator, not a member of this specialized carrier',
             dependencies={p:digest(HERE.parent/p) for p in dependencies},
             files={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.json') and p.name!='audit.json'})
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: actual sources, independent counts, exact recurrences, joint phases, all precision lifts, Singer map, common-switch universe',flush=True)


if __name__=='__main__':run()
