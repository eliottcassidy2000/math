"""Exact signed-affine compiler for paths sharing indexed square-clock ports."""
from itertools import product, combinations
from math import isqrt
from pathlib import Path
import argparse
import json

HERE=Path(__file__).resolve().parent
CHECKS=0

def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(message)

def constants(sums):
    c=[0]
    for s in sums: c.append(s-c[-1])
    return c

def values(c,t): return [((-1)**i)*t+x for i,x in enumerate(c)]

def local_domain(sums,N):
    c=constants(sums)
    L=max(1-x if i%2==0 else x-N for i,x in enumerate(c))
    U=min(N-x if i%2==0 else x-1 for i,x in enumerate(c))
    holes=set()
    empty=L>U
    for i,j in combinations(range(len(c)),2):
        a=(-1)**i-(-1)**j
        b=c[j]-c[i]
        if a==0:
            if b==0: empty=True
        elif b%a==0 and L<=b//a<=U: holes.add(b//a)
    empty |= len(holes)==max(0,U-L+1)
    return dict(c=c,L=L,U=U,holes=sorted(holes),empty=empty)

def compile_network(paths,ports=(),marks=(),distinct=()):
    """distinct entries are pairs of indexed ports required NOT to coincide.

    Return independent component root domains plus explicitly retained
    intercomponent forbidden-equality factors. No cross-path distinctness is
    inferred beyond this supplied list.
    """
    loc=[local_domain(p['sums'],p['N']) for p in paths]
    adj=[[] for _ in paths]
    equations=[]
    for (a,i),(b,j) in ports:
        eta=(-1)**(i+j)
        q=(-1)**j*(loc[a]['c'][i]-loc[b]['c'][j])
        equations.append((a,b,eta,q))
        adj[a].append((b,eta,q))
        adj[b].append((a,eta,-eta*q))
    component={};sign={};offset={};groups=[]
    for root in range(len(paths)):
        if root in component: continue
        cid=len(groups);group=[];todo=[root]
        component[root]=cid;sign[root]=1;offset[root]=0
        while todo:
            a=todo.pop();group.append(a)
            for b,eta,q in adj[a]:
                if b not in component:
                    component[b]=cid;sign[b]=eta*sign[a];offset[b]=eta*offset[a]+q;todo.append(b)
        groups.append(sorted(group))
    result=[]
    for cid,group in enumerate(groups):
        lo=[];hi=[];holes=set();pins=[];empty=any(loc[a]['empty'] for a in group)
        for a in group:
            if sign[a]==1:
                lo.append(loc[a]['L']-offset[a]);hi.append(loc[a]['U']-offset[a])
            else:
                lo.append(offset[a]-loc[a]['U']);hi.append(offset[a]-loc[a]['L'])
            holes.update(sign[a]*(x-offset[a]) for x in loc[a]['holes'])
        for a,b,eta,q in equations:
            if component[a]!=cid: continue
            coeff=sign[b]-eta*sign[a]
            rhs=eta*offset[a]+q-offset[b]
            if coeff==0:
                empty |= rhs!=0
            elif rhs%coeff:
                empty=True
            else: pins.append(rhs//coeff)
        for (a,i),value in marks:
            if component[a]!=cid: continue
            eps=(-1)**i
            pins.append((value-eps*offset[a]-loc[a]['c'][i])*(eps*sign[a]))
        for (a,i),(b,j) in distinct:
            if component[a]!=cid or component[b]!=cid: continue
            coeff=(-1)**i*sign[a]-(-1)**j*sign[b]
            rhs=(-1)**j*offset[b]+loc[b]['c'][j]-(-1)**i*offset[a]-loc[a]['c'][i]
            if coeff==0:
                empty |= rhs==0
            elif rhs%coeff==0: holes.add(rhs//coeff)
        L=max(lo);U=min(hi);holes=sorted(x for x in holes if L<=x<=U)
        empty |= L>U or len(set(pins))>1 or len(holes)==max(0,U-L+1)
        pin=pins[0] if pins else None
        if pin is not None:
            empty |= not(L<=pin<=U) or pin in holes
        count=0 if empty else (1 if pin is not None else U-L+1-len(holes))
        result.append(dict(paths=group,L=L,U=U,holes=holes,pin=pin,empty=bool(empty),count=count))
    cross=[]
    for (a,i),(b,j) in distinct:
        if component[a]==component[b]: continue
        cross.append(dict(left_component=component[a],left_sign=(-1)**i*sign[a],
                          right_component=component[b],right_sign=(-1)**j*sign[b],
                          rhs=(-1)**j*offset[b]+loc[b]['c'][j]-(-1)**i*offset[a]-loc[a]['c'][i]))
    return dict(local=loc,components=result,node_component=[component[a] for a in range(len(paths))],
                node_sign=[sign[a] for a in range(len(paths))],node_offset=[offset[a] for a in range(len(paths))],cross_forbidden=cross)

def root_values(c):
    if c['empty']: return []
    if c['pin'] is not None: return [c['pin']]
    holes=set(c['holes'])
    return [x for x in range(c['L'],c['U']+1) if x not in holes]

def compiled_solutions(record):
    answer=[]
    for roots in product(*(root_values(c) for c in record['components'])):
        if any(f['left_sign']*roots[f['left_component']]-f['right_sign']*roots[f['right_component']]==f['rhs']
               for f in record['cross_forbidden']): continue
        answer.append(tuple(s*roots[c]+b for s,c,b in zip(record['node_sign'],record['node_component'],record['node_offset'])))
    return sorted(answer)

def brute_solutions(paths,ports=(),marks=(),distinct=()):
    local=[]
    for p in paths:
        candidates=[]
        for start in range(1,p['N']+1):
            seq=[start]
            for s in p['sums']: seq.append(s-seq[-1])
            if all(1<=x<=p['N'] for x in seq) and len(seq)==len(set(seq)):
                candidates.append(seq)
        local.append(candidates)
    answer=[]
    for seqs in product(*local):
        if not all(seqs[a][i]==seqs[b][j] for (a,i),(b,j) in ports): continue
        if not all(seqs[a][i]==value for (a,i),value in marks): continue
        if not all(seqs[a][i]!=seqs[b][j] for (a,i),(b,j) in distinct): continue
        answer.append(tuple(s[0] for s in seqs))
    return sorted(answer)

def audit_case(paths,ports=(),marks=(),distinct=()):
    record=compile_network(paths,ports,marks,distinct)
    direct=brute_solutions(paths,ports,marks,distinct)
    check(compiled_solutions(record)==direct,('compiler mismatch',paths,ports,marks,distinct))
    if not record['cross_forbidden']:
        count=1
        for c in record['components']:count*=c['count']
        check(count==len(direct),'component product count')
    return record,direct

def run():
    cases=0
    # All clock words with0..3 edges, roots2..5, each host1..12.
    for k in range(4):
        for roots in product(range(2,6),repeat=k):
            for N in range(1,13):
                audit_case([dict(sums=[r*r for r in roots],N=N)]);cases+=1
    # Exhaust every subset of possible indexed attachments between two one-edge paths.
    port_pairs=[((0,i),(1,j)) for i in (0,1) for j in (0,1)]
    for s,t,N in product((4,9,16),repeat=3):
        paths=[dict(sums=[s],N=N),dict(sums=[t],N=N)]
        for mask in range(16):
            ports=[p for bit,p in enumerate(port_pairs) if mask>>bit&1]
            audit_case(paths,ports);cases+=1
    # All triangle clock assignments and independent host sizes1..8.
    cycle_ports=[((0,1),(1,0)),((1,1),(2,0)),((2,1),(0,0))]
    for sums in product((4,9,16),repeat=3):
        for Ns in product(range(1,9),repeat=3):
            audit_case([dict(sums=[s],N=N) for s,N in zip(sums,Ns)],cycle_ports);cases+=1
    hostiles={}
    specs={
      'marginal_ports':([dict(sums=[4],N=15),dict(sums=[16],N=15)],port_pairs[::3]),
      'odd_pin_nonintegral':([dict(sums=[s],N=15) for s in (4,9,16)],cycle_ports),
      'integer_pin_hits_collision_hole':([dict(sums=[s],N=8) for s in (4,9,9)],cycle_ports),
      'triangle_positive':([dict(sums=[s],N=30) for s in (25,49,36)],cycle_ports),
      'triangle_outside_host':([dict(sums=[s],N=29) for s in (25,49,36)],cycle_ports),
    }
    for name,(paths,ports) in specs.items():
        record,solutions=audit_case(paths,ports)
        check(all(not c['empty'] and c['count'] for c in compile_network(paths)['components']),name+' local nonempty')
        check(bool(solutions)==(name=='triangle_positive'),name+' global outcome')
        hostiles[name]=dict(paths=paths,ports=ports,compiled=record,solutions=solutions)
    # Each separate equality of the two-clock marginal hostile is feasible.
    for port in port_pairs[::3]:
        check(bool(audit_case(specs['marginal_ports'][0],[port])[1]),'separate port feasible')
    for M in range(1,33):
        x=[1,3];y=[1,(M+2)**2-1]
        check(all((a-b)%M==0 for a,b in zip(x,y)) and x!=y,'fixed congruence loses exact clock')
    # Optional same-component and cross-component distinctness guards.
    paths=[dict(sums=[9],N=8),dict(sums=[9],N=8)]
    inside,inside_solutions=audit_case(paths,[((0,0),(1,0))],distinct=[((0,1),(1,1))])
    check(not inside_solutions,'identical derived owners cannot be distinct')
    independent=[dict(sums=[],N=3),dict(sums=[],N=3)]
    cross,cross_solutions=audit_case(independent,distinct=[((0,0),(1,0))])
    check([c['count'] for c in cross['components']]==[3,3] and len(cross_solutions)==6,'cross-component factor')
    hostiles['cross_component_distinctness']=dict(compiled=cross,solutions=cross_solutions)
    # Filtered start domains do not inherit pairwise Helly from cyclic cosets.
    three_paths=[dict(sums=[],N=2),dict(sums=[4],N=3),dict(sums=[9],N=7)]
    pair_starts=[]
    for a,b in combinations(range(3),2):
        rec,sol=audit_case([three_paths[a],three_paths[b]],[((0,0),(1,0))])
        check(bool(sol),'every filtered-domain pair is feasible')
        pair_starts.append(dict(paths=[a,b],solutions=sol))
    three,three_sol=audit_case(three_paths,[((0,0),(1,0)),((1,0),(2,0))])
    check(not three_sol and three['components'][0]['empty'],'filtered domains lack pairwise Helly')
    hostiles['pairwise_not_joint_filtered_domains']=dict(paths=three_paths,pairs=pair_starts,compiled=three,solutions=three_sol)
    # Inherited full square-clock witness: local fiber2, indexed mark restores1.
    p=[3,33,16,9,27,22,14,11,25,24,1,35,29,20,5,31,18,7,2,23,26,10,6,19,30,34,15,21,4,12,13,36,28,8,17,32]
    sums=[a+b for a,b in zip(p,p[1:])]
    check(all(isqrt(s)**2==s for s in sums),'inherited witness square clocks')
    full,starts=audit_case([dict(sums=sums,N=36)])
    check(starts==[(3,),(9,)],'inherited even twin starts')
    for pos in range(36):
        marked,ms=audit_case([dict(sums=sums,N=36)],marks=[((0,pos),p[pos])])
        check(ms==[(3,)],'every indexed vertex fixes incoming twin')
    hostiles['incoming_full_clock_twin']=dict(roots=list(map(isqrt,sums)),compiled=full,solutions=starts)
    # A free component has a five-integer interval and one genuine collision hole.
    free_paths=[dict(sums=[9],N=8),dict(sums=[16],N=12)]
    free,fs=audit_case(free_paths,[((0,1),(1,0))])
    check(fs==[(2,7),(3,6),(4,5),(5,4)],'free-root positive control')
    hostiles['free_root_positive']=dict(compiled=free,solutions=fs)
    return dict(status='PASS',checks=CHECKS,exhaustive_network_cases=cases,
                universe='local words roots2..5 length0..3 hosts1..12; every two-one-edge port subset at sums/hosts4,9,16; triangle clocks sums4,9,16 and hosts1..8',
                inherited_clock_source='cf7c4fc327:05-knowledge/results/session42_shared_lattices_square_colors_20261023/twin_witness.json',
                hostiles=hostiles)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    r=json.loads(json.dumps(run()));target=HERE/'shared_clocks_receipt.json'
    if args.write:target.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    else:check(r==json.loads(target.read_text(encoding='utf-8')),'stored receipt')
    print(json.dumps(dict(status='PASS',checks=CHECKS,exhaustive_network_cases=r['exhaustive_network_cases']),sort_keys=True))
