#!/usr/bin/env python3
"""Exact selected-observer packets, shared-port lifting and native matching.

No claim that a selected packet is an elementary metagraph edge, or can be
inserted into an arbitrary already occupied global port inventory.
"""
import argparse
from array import array
from collections import Counter
import csv
import gzip
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def matword(A):
    return tuple(A[i][j] for i in range(len(A)) for j in range(len(A)))


def path_ok(A, p):
    return len(p) == len(A) and len(set(p)) == len(A) and all(A[a][b] for a, b in zip(p, p[1:]))


def partner(A, p):
    need(path_ok(A, p), 'actual Hamiltonian path before partner operation')
    keep = set(zip(p, p[1:]))
    return [[0 if i == j else (A[i][j] if (i,j) in keep or (j,i) in keep else 1-A[i][j])
             for j in range(len(A))] for i in range(len(A))]


def opposite(A):
    return [list(row) for row in zip(*A)]


def enumerate_paths(A):
    out = []
    def extend(p, unused):
        if not unused:
            out.append(tuple(p)); return
        for v in sorted(unused):
            if not p or A[p[-1]][v]:
                extend(p+[v], unused-{v})
    extend([], set(range(len(A))))
    return out


def witness_isomorphism(A, B):
    n=len(A)
    for f in permutations(range(n)):
        if all(A[i][j] == B[f[i]][f[j]] for i in range(n) for j in range(i+1,n)):
            return list(f)
    return None


def primitive_bank(A):
    kernel = [g for g in product((-1,0,1),repeat=4) if any(g) and all(sum(a*b for a,b in zip(row,g))==0 for row in A)]
    def below(h,g):
        return all(a*b>=0 and abs(a)<=abs(b) for a,b in zip(h,g))
    return [g for g in kernel if not any(h!=g and below(h,g) for h in kernel)]


def local_lifting():
    # The inherited actual two-nine-cell tournament, with 27 full native lifts.
    n=18
    A=[[0]*n for _ in range(n)]
    for u,v in combinations(range(n),2):
        a = int(1 <= (v-u)%9 <= 4) if u//9==v//9 else int((u-v)%3==0)
        A[u][v]=a;A[v][u]=1-a
    p=tuple(range(9))+tuple(9+(2+i)%9 for i in range(9))
    Q=partner(A,p)
    records=[]
    for a,b in product(range(9),repeat=2):
        f=tuple((i+a)%9 for i in range(9))+tuple(9+(i+b)%9 for i in range(9))
        actual=all(A[u][v]==A[f[u]][f[v]] for u in range(n) for v in range(n))
        need(actual==((a-b)%3==0), 'actual whole-block transporter congruence')
        if not actual:continue
        p2=tuple(f[u] for u in p);Q2=partner(A,p2)
        need(all(Q[u][v]==Q2[f[u]][f[v]] for u in range(n) for v in range(n)), 'partner transport retains actual owners')
        need(partner(opposite(A),tuple(reversed(p2)))==opposite(Q2),'C9 actual converse square')
        records.append((a,b))
    need(len(records)==27,'27 native lifts retained')
    need(len({(a,(1+b)%9) for a,b in records})==27,'27 distinct endpoint pairs')
    # Realize the shared-frame join on three named actual cells A,B,C.
    T3=[[0]*27 for _ in range(27)]
    for u,v in combinations(range(27),2):
        ci,cj=u//9,v//9
        forward=(1 <= (v-u)%9 <= 4) if ci==cj else ((u-v)%3==0 if cj-ci==1 else True)
        T3[u][v]=int(forward);T3[v][u]=int(not forward)
    p3=tuple(range(9))+tuple(9+(2+i)%9 for i in range(9))+tuple(18+(1+i)%9 for i in range(9))
    need(path_ok(T3,p3),'actual three-cell Hamiltonian carrier')
    full_lifts=0
    for phases in product(range(9),repeat=3):
        f=tuple(9*(u//9)+(u%9+phases[u//9])%9 for u in range(27))
        actual=all(T3[u][v]==T3[f[u]][f[v]] for u,v in combinations(range(27),2))
        need(actual==(phases[0]%3==phases[1]%3==phases[2]%3),'actual global source exactly shared local frame join')
        if actual:
            full_lifts+=1
            need(path_ok(T3,tuple(f[u] for u in p3)),'every shared global lift transports actual path')
    need(full_lifts==81,'81 global native frames, not729 independent pair choices')
    cases=0;counted=[]
    for q,r in product(range(3),repeat=2):
        for pin_left,pin_right in product([None]+list(range(9)),repeat=2):
            L=[x for x in records if x[0]%3==q and (pin_left is None or x[1]==pin_left)]
            R=[x for x in records if x[0]%3==r and (pin_right is None or x[0]==pin_right)]
            if not L or not R:continue
            ell,rr=L[0],R[0]
            KL=[tuple((x[i]-ell[i])%9 for i in range(2)) for x in L]
            KR=[tuple((x[i]-rr[i])%9 for i in range(2)) for x in R]
            image=Counter((u[1]-v[0])%9 for u in KL for v in KR)
            defect=(rr[0]-ell[1])%9
            direct=[(u,v) for u in L for v in R if u[1]==v[0]]
            need(bool(direct)==(defect in image),'exact obstruction class membership')
            need(len(direct)==image.get(defect,0),'exact compatible kernel-lift multiplicity')
            if image:
                need(len(set(image.values()))==1,'homomorphism fibres have equal kernel size')
            for u,v in direct:
                # Local physical B-owner is at clock16 in first pair and clock0 in second.
                need(u[1]==v[0], 'one actual B owner, not independent edge copies')
                need((u[0]-u[1])%3==0 and (v[0]-v[1])%3==0,'integer local quotient certificate')
            cases+=1
            if q==r==0 and (pin_left,pin_right) in [(None,None),(0,None),(0,0),(0,3)]:
                counted.append(dict(pins=[pin_left,pin_right],local_sizes=[len(L),len(R)],defect=defect,image=sorted(image),global_lifts=len(direct)))
    need(cases==144,'complete nonempty quotient-and-pin universe')
    need(next(x for x in counted if x['pins']==[None,None])['global_lifts']==27,'kernel product81 contracts to27')
    need(next(x for x in counted if x['pins']==[0,3])['global_lifts']==0,'same quotient and incompatible actual shared owner')
    return dict(native_lifts=27,quotient_states=3,local_kernel=9,gluing_cases=cases,controls=counted,
                quotient_return=[3,3],return_endpoint_owners=[3,13],three_cell_native_lifts=full_lifts,shared_owner_clocks=[16,0])


def load_atlas(repo,pins):
    blobs={}
    for item in pins['files']:
        data=(repo/item['path']).read_bytes()
        need(sha256(data).hexdigest()==item['sha256'],'frozen input hash '+item['path'])
        blobs[item['path'].split('/')[-1]]=data
    rows=list(csv.DictReader(blobs['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'))
    owner=array('H'); owner.frombytes(gzip.decompress(blobs['owner_map.u16.gz']))
    if sys.byteorder!='little':owner.byteswap()
    mid={}
    for row in rows:
        for key in ['class_a','class_b']:mid[int(row[key])]=int(row['merged_id'])
    original=json.loads(blobs['receipt.json'])
    x=original['order8']['higher_layers'][-1]['lifts'][0]
    image=[0]*len(x); degree=Counter();sep_ok=True;edge_count=0
    for row in csv.DictReader(gzip.decompress(blobs['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t'):
        i,j=int(row['source']),int(row['target']);w=int(row['blue'])+int(row['black']);edge_count+=1
        if i==j:continue
        image[i]+=w*(x[i]-x[j]);image[j]+=w*(x[j]-x[i])
        if i==2874 or j==2874:
            degree['incident']+=1;degree['weighted']+=w
            sep_ok &= w%2==0
    need(edge_count==438216,'full original edge export')
    need(all(v%256==0 for v in image),'original full character harmonic modulo256')
    need(sep_ok and degree==dict(incident=273,weighted=898),'sparse harmonic separator4e2874')
    chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
    def annotate(A,p):
        need(path_ok(A,p),'every listed actual arc')
        mask=sum(A[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(chords))
        source=mid[owner[mask]];ordinary=int(owner[full^mask]);target=mid[ordinary]
        Q=partner(A,p)
        qmask=sum(Q[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(chords))
        need(qmask==full^mask,'literal partner equals original forward-bit complement')
        rev=tuple(reversed(p));op=opposite(A)
        rmask=sum(op[rev[i]][rev[j]]<<bit for bit,(i,j) in enumerate(chords))
        need(mid[owner[rmask]]==source and mid[owner[full^rmask]]==target,'actual reversed-converse merged owners')
        need(partner(op,rev)==opposite(Q),'literal actual converse square')
        return dict(path=list(p),source=source,target=target,ordinary_target=ordinary,mask=mask,phase=(x[source]-x[target])%8)
    banks={};matrices={}
    for source,count in [(1238,409),(2242,407),(2246,523)]:
        canon=int(rows[source]['canon_a']);A=[[0]*8 for _ in range(8)]
        for bit,(i,j) in enumerate(combinations(range(8),2)):
            A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
        paths=enumerate_paths(A)
        direct=[p for p in permutations(range(8)) if all(A[a][b] for a,b in zip(p,p[1:]))]
        need(set(paths)==set(direct) and len(paths)==count,'two independent complete Hamilton enumerations')
        aut=sum(all(A[i][j]==A[f[i]][f[j]] for i in range(8) for j in range(i+1,8)) for f in permutations(range(8)))
        need(aut==1,'trivial source automorphism group')
        banks[source]=[annotate(A,p) for p in paths];matrices[source]=A
        need(all(r['source']==source for r in banks[source]),'owner map authenticates every source path')
    return banks,matrices,x,dict(edges=edge_count,modulus=256,separator=dict(node=2874,value=4,incident_edges=273,weighted_degree=898))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path);parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    repo=args.repo
    if repo is None:
        repo=next((p for p in HERE.parents if (p/'05-knowledge/results/metagraph_order8_20261023/owner_map.u16.gz').is_file()),None)
    if repo is None:raise RuntimeError('Supply --repo pointing to the research checkout containing the pinned original atlas files.')
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    native=local_lifting();banks,As,x,electrical=load_atlas(repo,pins)
    # Strongest packet: every slot's target stays fixed, so every character survives.
    source=1238;keys=[(3,2,2397),(2,3,2790)]
    pools=[[r for r in banks[source] if (r['path'][0],r['path'][-1],r['target'])==k] for k in keys]
    need([len(s) for s in pools]==[2,2],'complete target-pinned local fibres, not invented candidate lists')
    states=[(a,b) for a,b in product(*pools) if a['path'][1:7]==b['path'][1:7]]
    need(len(states)==2,'shared interior-owner join has exactly two states')
    need(all(a['path'][2:4]==b['path'][2:4]==[6,7] and a['path'][5:7]==b['path'][5:7]==[5,4] for a,b in states),'fixed literal shared arc pins')
    need(all([a['phase'],b['phase']]==[5,4] for a,b in states),'original per-slot phases')
    need(all(len({tuple(a['path']),tuple(b['path'])})==2 for a,b in states),'no recycled path inside selected packet')
    need(len({tuple(r['path']) for state in states for r in state})==4,'four distinct actual path ports')
    need(len({matword(partner(As[source],r['path'])) for state in states for r in state})==4,'four distinct literal complement partners')
    costs=[sum(a!=b for a,b in zip(states[0][i]['path'],states[1][i]['path'])) for i in range(2)]
    need(costs==[2,2],'each slot changes two actual owner clocks')
    need(all(states[0][i]['path']!=states[1][i]['path'] for i in range(2)),'minimum two-slot budget')
    need(sum(a!=b for a,b in zip(states[0][0]['path'][1:7],states[1][0]['path'][1:7]))==2,'two shared-frame owner variables change')
    need(len([(a,b) for a,b in states if a['path'][1]==0])==1,'pinning a changing shared owner freezes the move')
    isos=[]
    for i in range(2):
        Qa=partner(As[source],states[0][i]['path']);Qb=partner(As[source],states[1][i]['path'])
        iso=witness_isomorphism(Qa,Qb);branch='direct'
        if iso is None:iso=witness_isomorphism(Qa,opposite(Qb));branch='converse'
        need(iso is not None,'explicit target isomorphism, not merely same scalar phase')
        isos.append(dict(branch=branch,permutation=iso))
    Aeq=[[1,1,0,0],[0,0,1,1],[1,0,-1,0]];rhs=[1,1,0]
    solutions=[z for z in product((0,1),repeat=4) if all(sum(a*b for a,b in zip(row,z))==c for row,c in zip(Aeq,rhs))]
    atoms=primitive_bank(Aeq)
    need(set(solutions)=={(1,0,1,0),(0,1,0,1)},'exact one-hot integer lift')
    need(set(atoms)=={(1,-1,1,-1),(-1,1,-1,1)},'complete bounded primitive kernel bank')
    need(all(any(sum(a*b for a,b in zip(row,g)) for row in Aeq) for g in [(-1,1,0,0),(0,0,-1,1)]),'uncoupled moves fail actual shared-port row')
    target_packet=dict(source=source,keys=keys,states=states,source_matrix=As[source],fixed_arc_clocks=[2,5],shared_owner_clocks=list(range(1,7)),slot_budget=2,slot_clock_costs=costs,shared_variable_cost=2,target_isomorphisms=isos,onehot_matrix=Aeq,rhs=rhs,solutions=solutions,primitive_bank=atoms)
    # A weaker frozen-character packet changes each slot's target while preserving incidence.
    source=2242;keys2=[(6,4,2),(6,1,2)]
    pools2=[[r for r in banks[source] if (r['path'][0],r['path'][-1],r['phase'])==k] for k in keys2]
    need([len(s) for s in pools2]==[2,2],'complete frozen-phase local fibres')
    states2=[(a,b) for a,b in product(*pools2) if sorted([a['target'],b['target']])==[1957,2874]]
    need(len(states2)==2,'fixed target-incidence margin leaves two states')
    need(all(a['path'][:2]==b['path'][:2]==[6,7] for a,b in states2),'fixed literal shared initial arc')
    need(len({tuple(r['path']) for st in states2 for r in st})==4,'not an exchange of identical literal path inventory')
    need([r['target'] for r in states2[0]]==list(reversed([r['target'] for r in states2[1]])),'two slots swap target owners')
    sep_phases=[[(-4*int(r['target']==2874))%8 for r in st] for st in states2]
    need(sep_phases[0]!=sep_phases[1],'universal per-slot bank forbids frozen-character-only trade')
    # Selected Laplacians agree as INTEGER matrices, stronger than modulo-eight aggregate equality.
    def incidence(st):return Counter((min(r['source'],r['target']),max(r['source'],r['target'])) for r in st)
    need(incidence(states2[0])==incidence(states2[1]),'entire integer merged-edge multiset unchanged')
    weak_packet=dict(source=source,keys=keys2,states=states2,fixed_arc_clock=0,separator_phases=sep_phases,frozen_values={str(i):x[i] for i in [2242,1957,2874]})
    # A length-three matching cycle arises from actual clock pins, not a chosen edge list.
    source=2242;target=3442
    ps=[(1,0,5,3,2,4,6,7),(1,5,4,3,2,0,6,7),(4,1,5,3,2,0,6,7)]
    clocksets=[[0,3,4,6,7],[3,4,5,6,7],[2,3,4,6,7]];representatives=[ps[0],ps[1],ps[0]]
    pool=[r for r in banks[source] if r['target']==target]
    eligible=[]
    for clocks,p in zip(clocksets,representatives):
        matches=[r for r in pool if all(r['path'][j]==p[j] for j in clocks)]
        need(len(matches)==2,'complete target fibre yields exactly two pin-eligible native ports')
        eligible.append([ps.index(tuple(r['path'])) for r in matches])
    need([set(s) for s in eligible]==[{0,1},{1,2},{0,2}],'literal native cyclic eligibility')
    matchings=[s for s in product(*eligible) if len(set(s))==3]
    need(len(matchings)==2 and all(a!=b for a,b in zip(*matchings)),'sharp three-slot injective exchange')
    need(all(all(ps[i][j:j+2]==ps[0][j:j+2] for i in range(3)) for j in [3,6]),'two actual fixed shared arc flags')
    annotations=[next(r for r in pool if tuple(r['path'])==p) for p in ps]
    need(all(r['phase']==5 for r in annotations),'original character retained on all ternary matching ports')
    ternary=dict(source=source,target=target,complete_target_fibre_size=len(pool),paths=annotations,clocksets=clocksets,pin_values=[[p[j] for j in c] for p,c in zip(representatives,clocksets)],eligible=eligible,matchings=matchings,minimum_slot_budget=3,fixed_arc_clocks=[3,6],unlabelled_port_inventory_unchanged=True)
    out=dict(status='FINITE-EXACT PASS',source_commit=pins['commit'],checks=CHECKS,guard_scope='original frozen character and exact target-pinned universal preservation separately',native_lifting=native,original_electrical=electrical,target_pinned_packet=target_packet,frozen_character_packet=weak_packet,ternary_matching=ternary,global_extension='Requires all selected new literal ports unoccupied outside the packet, or an explicit global injective reassignment. No full-atlas rewiring asserted.')
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.out:args.out.write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()
