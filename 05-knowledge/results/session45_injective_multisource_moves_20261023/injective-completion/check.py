#!/usr/bin/env python3
"""Exact partial completion, literal owner normalization, and resource controls."""
import argparse
from array import array
from collections import Counter
import csv
import gzip
import importlib.util
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(t, why):
    global CHECKS
    CHECKS += 1
    if not t:
        raise RuntimeError(why)

def matrix(c, n=8):
    a = [[0]*n for _ in range(n)]
    for bit, (i,j) in enumerate(combinations(range(n), 2)):
        a[i][j] = (c >> (n*(n-1)//2-1-bit)) & 1
        a[j][i] = 1-a[i][j]
    return tuple(map(tuple,a))

def paths(a):
    ans=[]
    def visit(p, left):
        if not left:
            ans.append(tuple(p)); return
        for v in sorted(left):
            if not p or a[p[-1]][v]: visit(p+[v], left-{v})
    visit([], set(range(len(a))))
    return ans

def valid(a,p):
    return sorted(p)==list(range(len(a))) and all(a[u][v] for u,v in zip(p,p[1:]))

def k(port):
    a,p=port; n=len(a)
    need(valid(a,p), 'K applied only to actual Hamiltonian port')
    adj={frozenset((u,v)) for u,v in zip(p,p[1:])}
    b=tuple(tuple(0 if u==v else a[u][v] if frozenset((u,v)) in adj else 1-a[u][v] for v in range(n)) for u in range(n))
    return b,p

def c(port):
    a,p=port
    return tuple(zip(*a)), tuple(reversed(p))

def orbit(port):
    kp=k(port); cp=c(port); kcp=k(cp)
    need(k(kp)==port, 'literal K involution')
    need(c(cp)==port, 'literal converse involution')
    need(c(kp)==kcp, 'literal K/converse square')
    result=frozenset((port,kp,cp,kcp))
    need(len(result)==4, 'free four-port orbit at n>=3')
    return result

def digest_port(port):
    return sha256(json.dumps(port,separators=(',',':')).encode()).hexdigest()

def match(neighbors):
    """Simple augmenting-path matching; returns assignments or None."""
    owner={}
    def aug(i,seen):
        for r in sorted(neighbors[i]):
            if r in seen: continue
            seen.add(r)
            if r not in owner or aug(owner[r],seen):
                owner[r]=i; return True
        return False
    for i in range(len(neighbors)):
        if not aug(i,set()): return None
    answer=[None]*len(neighbors)
    for r,i in owner.items(): answer[i]=r
    return answer

def brute(neighbors):
    return [a for a in product(*[sorted(x) for x in neighbors]) if len(set(a))==len(a)]

def laminar(neighbors):
    return all(not (a & b) or a<=b or b<=a for a,b in combinations(neighbors,2))

def laminar_check(neighbors):
    if not laminar(neighbors): raise ValueError('nonlaminar input')
    return all(sum(t<=s for t in neighbors)<=len(s) for s in set(map(frozenset,neighbors)))

def tree_join(tables, scopes, edges):
    """Exact tree semijoin; tables contain literal shared owner values."""
    a=[list(t) for t in tables]
    if len(edges)!=len(a)-1:raise ValueError('not a tree')
    graph=[[] for _ in a]
    for u,v in edges: graph[u].append(v); graph[v].append(u)
    parent={0:-1}; order=[0]
    for u in order:
        for v in graph[u]:
            if v!=parent[u]:
                if v in parent:raise ValueError('cycle in join tree')
                parent[v]=u; order.append(v)
    if len(order)!=len(a):raise ValueError('disconnected join tree')
    for var in set().union(*scopes):
        owners={i for i,s in enumerate(scopes) if var in s}
        reached={min(owners)};queue=list(reached)
        for i in queue:
            for j in graph[i]:
                if j in owners and j not in reached:reached.add(j);queue.append(j)
        if reached!=owners:raise ValueError('running intersection fails for '+var)
    for child in reversed(order[1:]):
        par=parent[child]; sep=set(scopes[par])&set(scopes[child])
        a[par]=[r for r in a[par] if any(all(r[v]==s[v] for v in sep) for s in a[child])]
    if not a[0]: return None
    chosen={0:a[0][0]}
    for child in order[1:]:
        par=parent[child]; sep=set(scopes[par])&set(scopes[child])
        chosen[child]=next(s for s in a[child] if all(chosen[par][v]==s[v] for v in sep))
    return [chosen[i] for i in range(len(a))]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text()); blobs={}
    for pin in pins['files']:
        b=(args.repo/pin['path']).read_bytes()
        need(len(b)==pin['bytes'] and sha256(b).hexdigest()==pin['sha256'], 'frozen original input '+pin['path'])
        blobs[Path(pin['path']).name]=b
    for pin in pins['inherited_results']:
        b=(args.repo/pin['path']).read_bytes()
        need(len(b)==pin['bytes'] and sha256(b).hexdigest()==pin['sha256'],'frozen inherited result '+pin['path'])
    inherited_path=args.repo/'05-knowledge/results/session44_guarded_native_moves_20261023/coupled-moves/check.py'
    spec=importlib.util.spec_from_file_location('session44_coupled_verified',inherited_path)
    inherited=importlib.util.module_from_spec(spec);spec.loader.exec_module(inherited)
    native_lifts=inherited.local_lifting()
    need(native_lifts['native_lifts']==27 and native_lifts['three_cell_native_lifts']==81,'authenticated actual C9 source and all whole-source lifts')
    rows=list(csv.DictReader(blobs['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'))
    mid={int(r[key]):int(r['merged_id']) for r in rows for key in ('class_a','class_b')}
    owner=array('H');owner.frombytes(gzip.decompress(blobs['owner_map.u16.gz']))
    if sys.byteorder!='little':owner.byteswap()
    x=json.loads(blobs['receipt.json'])['order8']['higher_layers'][-1]['lifts'][0]
    full=(1<<21)-1; chords=[(i,j) for i in range(8) for j in range(i+2,8)]
    def annotate(a,p):
        need(valid(a,p), 'actual directed Hamiltonian arcs')
        mask=sum(a[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords))
        s,t=mid[owner[mask]],mid[owner[full^mask]]
        need(k((a,p))[0]==matrix_from_path(a,p,full^mask), 'original complement word retains literal partner')
        return dict(path=list(p),source=s,target=t,mask=mask,phase=(x[s]-x[t])%8)
    def matrix_from_path(a,p,mask):
        b=[[0]*8 for _ in range(8)]
        for u,v in zip(p,p[1:]):b[u][v]=1
        for bit,(i,j) in enumerate(chords):b[p[i]][p[j]]=(mask>>bit)&1
        for i,j in combinations(range(8),2):
            u,v=p[i],p[j];b[v][u]=1-b[u][v]
        return tuple(map(tuple,b))
    image=[0]*len(x);ec=0
    for row in csv.DictReader(gzip.decompress(blobs['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t'):
        u,v,w=int(row['source']),int(row['target']),int(row['blue'])+int(row['black']);ec+=1
        image[u]+=w*(x[u]-x[v]);image[v]+=w*(x[v]-x[u])
    need(ec==438216 and all(v%256==0 for v in image),'entire original frozen harmonic character modulo256')

    a=matrix(int(rows[2242]['canon_a'])); bank=paths(a)
    brute_bank=[p for p in permutations(range(8)) if valid(a,p)]
    need(set(bank)==set(brute_bank) and len(bank)==407,'two independent complete native source2242 enumerations')
    annotated=[annotate(a,p) for p in bank]
    pool=[r for r in annotated if r['target']==3442]
    expected={(1,0,5,3,2,4,6,7),(1,5,4,3,2,0,6,7),(4,1,5,3,2,0,6,7)}
    need({tuple(r['path']) for r in pool}==expected,'complete original three-port target3442 fibre')
    pool=sorted(pool,key=lambda r:r['path'])
    eligible=[i for i,r in enumerate(pool) if r['path'][0]==1 and r['path'][3:5]==[3,2] and r['path'][6:8]==[6,7]]
    need(eligible==[0,1], 'all literal pins give exact two-port neighbor set')
    neighbors=[set(eligible)]*3
    need(all(match([neighbors[i],neighbors[j]]) is not None for i,j in combinations(range(3),2)), 'every pair separately injective')
    need(len(pool)>=3 and all(len(s)>=1 for s in neighbors), 'unrestricted target capacity and singleton capacity pass')
    need(match(neighbors) is None and not brute(neighbors), 'joint injection fails in two independent solvers')
    need(all(r['phase']==5 for r in pool),'original fixed mod8 guard retained')
    orb=[orbit((a,tuple(r['path']))) for r in pool]
    need(sum(len(s) for s in orb)==len(set().union(*orb))==12, 'three disjoint literal K/converse resources')
    # Remainder occupies p0's entire resource, while the remaining two slots use p1,p2.
    unavailable=orb[0]
    available=[i for i,s in enumerate(orb) if not s&unavailable]
    positive=match([set(available)]*2)
    need(positive is not None and match([set(available)]*3) is None, 'unavailable closed remainder is subtracted before matching')
    # Differently labelled local presentations, explicitly normalized to physical owners.
    relabelings=[tuple(range(8)),(1,2,3,4,5,6,7,0),(7,6,5,4,3,2,1,0)]
    normalized=[]
    for f in relabelings:
        inv=tuple(f.index(i) for i in range(8))
        local=tuple(tuple(a[f[u]][f[v]] for v in range(8)) for u in range(8))
        localports=[]
        for j in eligible:
            lp=tuple(inv[v] for v in pool[j]['path'])
            need(valid(local,lp), 'relabeled source retains actual path')
            physical_a=tuple(tuple(local[inv[u]][inv[v]] for v in range(8)) for u in range(8))
            physical_p=tuple(f[v] for v in lp)
            need((physical_a,physical_p)==(a,tuple(pool[j]['path'])), 'explicit physical owner normalization')
            localports.append(digest_port((physical_a,physical_p)))
        normalized.append(localports)
    need(normalized[0]==normalized[1]==normalized[2], 'labels cannot multiply physical resources')
    # Exhaust every ordered collection of at most four nonempty subsets of four resources.
    subsets=[set(j for j in range(4) if mask>>j&1) for mask in range(1,16)]
    tests=0; laminar_tests=0
    for count in range(1,5):
        for ns in product(subsets,repeat=count):
            tests+=1; m=match(ns)
            need((m is not None)==bool(brute(ns)), 'augmenting matcher agrees with exhaustive assignment')
            if laminar(ns):
                laminar_tests+=1
                need(laminar_check(ns)==(m is not None),'compressed laminar inequalities exact')
    # Join-tree tests retain all 27 native C9 local lifts, not one quotient representative.
    table=[{'a':u,'b':v} for u,v in product(range(9),repeat=2) if (u-v)%3==0]
    table2=[{'b':u,'c':v} for u,v in product(range(9),repeat=2) if (u-v)%3==0]
    joined=[(r,s) for r in table for s in table2 if r['b']==s['b']]
    need(len(table)==len(table2)==27 and len(joined)==81,'all local kernel lifts survive exactly when owners agree')
    need(tree_join([table,table2],[{'a','b'},{'b','c'}],[(0,1)]) is not None,'join-tree witness nonvacuous')
    bad1=[r for r in table if r['b']==0];bad2=[r for r in table2 if r['b']==3]
    need(len(bad1)==len(bad2)==3 and tree_join([bad1,bad2],[{'a','b'},{'b','c'}],[(0,1)]) is None,'equal quotient but wrong literal owner rejected')
    nonlaminar=[{0,1},{0,1},{0,2},{0,2}]
    need(all(sum(t<=s for t in nonlaminar)<=len(s) for s in nonlaminar) and match(nonlaminar) is None,'without laminarity compressed neighborhood tests fail')
    cycle_tables=[[{'a':v,'b':v} for v in (0,1)],[{'b':v,'c':v} for v in (0,1)],[{'a':v,'c':1-v} for v in (0,1)]]
    scopes=[{'a','b'},{'b','c'},{'a','c'}]
    need(all(any(all(r[v]==s[v] for v in scopes[i]&scopes[j]) for r in cycle_tables[i] for s in cycle_tables[j]) for i,j in combinations(range(3),2)),'cyclic shared-owner hostile pairwise feasible')
    need(not any(all(all(r[v]==s[v] for v in set(r)&set(s)) for r,s in combinations(rows,2)) for rows in product(*cycle_tables)),'cyclic shared-owner hostile has no global assignment')
    rejected=False
    try:tree_join(cycle_tables,scopes,[(0,1),(1,2)])
    except ValueError:rejected=True
    need(rejected,'join-tree solver rejects missing running intersection instead of false certificate')
    # Candidate independent audit. The second source is reversed to balance target incidence.
    candidate=json.loads((HERE/'candidate-input.json').read_text()); states=[[],[]]; origins=[]; targetcounts=[]
    for i,slot in enumerate(candidate['slots']):
        source=slot['source']; aa=matrix(int(rows[source]['canon_a']))
        for state in (0,1):
            key=('before','after')[state ^ (i==1)]
            rec=slot[key]; pp=tuple(rec['path']); rr=annotate(aa,pp)
            need(rr['source']==source and rr['target']==rec['target'] and rr['mask']==rec['mask'],'candidate authenticates original owner/mask/target')
            p0,p7,j,u,v=slot['pins']
            need(pp[0]==p0 and pp[-1]==p7 and pp[j:j+2]==(u,v),'candidate retains every listed actual pin')
            states[state].append(orbit((aa,pp)))
            origins.append(rr)
    for state in states:
        need(sum(map(len,state))==len(set().union(*state))==8,'each candidate endpoint uses eight distinct literal K/converse ports')
    need(len(set().union(*(states[0]+states[1])))==16,'all four candidate resource orbits distinct')
    debt=Counter()
    for i,slot in enumerate(candidate['slots']):
        s=slot['source']
        for state in (0,1):
            t=slot[('before','after')[state ^ (i==1)]]['target']; sign=(-1,1)[state]
            debt[s]+=sign*2*(x[s]-x[t]);debt[t]+=sign*2*(x[t]-x[s])
    need(all(v==0 for v in debt.values()),'candidate balanced converse-paired exact integral frozen current debt zero')
    def audit_physical_packet(filename):
        blob=(HERE/filename).read_bytes(); packet=json.loads(blob)
        states=[[],[]]; records=[]; maps={}; currents=[Counter(),Counter()]; physical_matrices=[]; shared_owner_sets=[]
        for slot in packet['slots']:
            if 'adjacency' in slot:
                slot=dict(slot)
                slot['canonical_matrix']=slot['adjacency']
                slot['local_to_physical']=slot['local_to_physical_owners']
                slot['physical_matrix']=slot['adjacency']
                slot['shared_physical_pins']=slot['preserved_owner_pins']
                for key,old in [('before','old'),('after','new')]:
                    slot[key]=dict(slot[old],physical_path=slot[old]['path'])
            source=slot['source']; aa=matrix(int(rows[source]['canon_a']))
            need(tuple(map(tuple,slot['canonical_matrix']))==aa,'packet source canonical matrix authenticates atlas')
            f=tuple(slot['local_to_physical']);need(sorted(f)==list(range(8)),'packet physical owner map bijective')
            if source in maps: need(maps[source]==f,'one physical owner map per repeated actual source')
            maps[source]=f
            inv=tuple(f.index(i) for i in range(8))
            physical=tuple(tuple(aa[inv[u]][inv[v]] for v in range(8)) for u in range(8))
            need(physical==tuple(map(tuple,slot['physical_matrix'])),'packet carries its full physical source matrix')
            physical_matrices.append(physical)
            shared_owner_sets.append({owner for clock,owner in slot['shared_physical_pins']})
            slotrec=[]
            for state,key in enumerate(('before','after')):
                record=slot[key];pp=tuple(record['path']);phys=tuple(record['physical_path'])
                need(phys==tuple(f[v] for v in pp),'packet physical path is exact owner transport')
                rr=annotate(physical,phys)
                need(rr['source']==source and rr['target']==record['target'] and rr['mask']==record['mask'],'physical packet retains original source/target/word')
                need(all(phys[clock]==owner for clock,owner in slot['shared_physical_pins']),'all shared physical indexed pins retained')
                need(all(phys[clock]==owner for clock,owner in slot.get('all_fixed_physical_clock_owners',[])),'every additional unchanged physical pin retained')
                resource=orbit((physical,phys));states[state].append(resource)
                need(all(c((physical,phys))[1][7-clock]==owner for clock,owner in slot['shared_physical_pins']),'converse transports every shared pin to reversed clock')
                target=record['target']; currents[state][source]+=2*(x[source]-x[target]);currents[state][target]+=2*(x[target]-x[source])
                slotrec.append(rr)
            need(slotrec[0]['phase']==slotrec[1]['phase'],'physical packet per-slot frozen mod8 guard preserved')
            records.append(slotrec)
        count=len(packet['slots'])
        for state in states:
            need(len(set().union(*state))==4*count,'physical packet full K/converse resource injection')
        union=set().union(*(states[0]+states[1]))
        need(len(union)==8*count,'before/after physical resource orbits all distinct')
        current_debt={node:currents[1][node]-currents[0][node] for node in currents[0].keys()|currents[1].keys()}
        need(all(v%256==0 for v in current_debt.values()),'physical packet full frozen current debt zero modulo256')
        if filename=='packet-four-source.json':need(all(v==0 for v in current_debt.values()),'four-source frozen integer current debt zero')
        shared_owners=sorted(set.intersection(*shared_owner_sets))
        conflicts=[dict(owner_pair=[u,v],source_values=[a[u][v] for a in physical_matrices]) for u,v in combinations(shared_owners,2) if len({a[u][v] for a in physical_matrices})>1]
        if filename=='packet-four-source.json':
            need(all(a[3][4]==1 for a in physical_matrices),'declared shared actual arc consistent in every source')
            need(any(r['owner_pair']==[0,7] and r['source_values']==[1,1,1,0] for r in conflicts),'shared pin charts are not a common-host induced tournament')
        return dict(input=filename,input_sha256=sha256(blob).hexdigest(),slots=count,source_owner_maps={str(s):list(f) for s,f in maps.items()},ports_per_endpoint=4*count,union_raw_ports=len(union),resource_orbits_per_endpoint=count,all_retained_records=records,current_debt=dict(sorted(current_debt.items())),integer_current_debt_zero=all(v==0 for v in current_debt.values()),current_debt_mod256_zero=True,resource_keys_per_endpoint=[[min(digest_port(p) for p in resource) for resource in state] for state in states],shared_owner_overlap_conflicts=conflicts,common_host_amalgamation_claimed=False,unavailable_remainder_scope='Any unchanged remainder disjoint from the proposed endpoint resources; a sufficient common condition is disjointness from the union of both endpoint resources.')
    physical_audits=[audit_physical_packet('packet-four-source.json')]
    if (HERE/'packet-three-source.json').is_file():
        physical_audits.append(audit_physical_packet('packet-three-source.json'))
    out=dict(status='FINITE-EXACT PASS',checks=CHECKS,source_commit=pins['commit'],native_hostile=dict(source=2242,target=3442,full_source_count=407,complete_target_fibre=pool,slots=3,eligible=eligible,pairwise_injective=True,global_injective=False,minimal_slot_count=3,original_phase_mod8=5,resource_orbit_sizes=[len(s) for s in orb],normalized_owner_maps=[list(f) for f in relabelings],physical_neighbor_hashes=normalized),unavailable_remainder=dict(resource=0,available_resources=available,two_slot_witness=positive,three_slot_completion=False),finite_matching_universe=dict(resources=4,max_slots=4,nonempty_neighborhood_tuples=tests,laminar_tuples=laminar_tests),full_kernel_join=dict(local_lifts=[27,27],full_joint_lifts=81,equal_quotient_hostile_counts=[3,3],equal_quotient_hostile_joint=0),candidate_audit=dict(input_sha256=sha256((HERE/'candidate-input.json').read_bytes()).hexdigest(),source_ids=[r['source'] for r in candidate['slots']],literal_ports_per_endpoint=8,union_literal_ports=16,resource_orbits_per_endpoint=2,all_retained_records=origins,integer_current_debt=dict(sorted(debt.items())),scope='Partial selected K/converse-closed resources. Extends only if proposed resources avoid unavailable remainder. No full-atlas incidence rewiring.'))
    out['physical_packet_audits']=physical_audits
    out['inherited_native_lift_authentication']=dict(source_path=str(inherited_path.relative_to(args.repo)),checks=inherited.CHECKS,local_lifts=native_lifts['native_lifts'],whole_source_three_cell_lifts=native_lifts['three_cell_native_lifts'],gluing_cases=native_lifts['gluing_cases'])
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(status=out['status'],checks=CHECKS,receipt_sha256=sha256(args.out.read_bytes()).hexdigest())))

if __name__=='__main__': main()
