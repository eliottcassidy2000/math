#!/usr/bin/env python3
"""Complete exact repair of the six n12/14 failed source cells.

Pinned inherited GF(2) solver + independent local integer-row replay.
No assertion gates, random sampling, tournament quotient, or HP census.
"""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
from collections import Counter
import importlib.util
import json
import sys

HERE=Path(__file__).resolve().parent
PINS={'audit_core.py':'9e0f86d30773b46dad359c04ead55dc07e606cad54156fb4eb0025b42a7339a4',
      'independent_audit.py':'376782e9c9eab3f3dae563a22735271d064fe73f4ef0361de044beb56cdc85b0'}
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
for name,digest in PINS.items():
    need(sha256((HERE/'inherited'/name).read_bytes()).hexdigest()==digest,('dependency changed',name))
sys.path.insert(0,str(HERE/'inherited'))
import audit_core as gf
import independent_audit as bank

SIX=((12,(2,7,9),True,True,1),(12,(6,7,9),True,True,1),
     (14,(3,5,10),False,True,3),(14,(3,7,10),False,True,3),
     (14,(5,9,10),False,False,3),(14,(7,9,10),False,False,3))
def arc(A,u,v):return (A[u]>>v)&1
def edge_set(q):return {tuple(sorted((u,v))) for u,v in zip(q,q[1:])}
def directed(A,q):return all(arc(A,u,v) for u,v in zip(q,q[1:]))
def witness(A,q):
    r=dict(zip(q,q[::-1]));n=len(A)
    for u,v in combinations(range(n),2):
        if arc(A,u,v)!=arc(A,r[v],r[u]):
            return {'pair':[u,v],'reflected_reversed_pair':[r[v],r[u]],
                    'arc_bits':[arc(A,u,v),arc(A,r[v],r[u])]}
    return None
def move(A,q):
    n=len(A);out=[((1<<n)-1)^A[v]^(1<<v) for v in range(n)]
    for u,v in zip(q,q[1:]):out[u]^=1<<v;out[v]^=1<<u
    return tuple(out)
def neighbourhoods(A,v):
    return {'vertex':v,'in':[u for u in range(len(A)) if u!=v and arc(A,u,v)],
            'out':[u for u in range(len(A)) if u!=v and arc(A,v,u)]}
def clone_data(A):
    """Exact maximal transitive modules, retaining their literal orders."""
    n=len(A);twins=[];adj=[set() for _ in A]
    for u,v in combinations(range(n),2):
        if all(arc(A,u,w)==arc(A,v,w) for w in range(n) if w not in (u,v)):
            twins.append([u,v]);adj[u].add(v);adj[v].add(u)
    remaining=set(range(n));blocks=[]
    while remaining:
        stack=[min(remaining)];component=set()
        while stack:
            v=stack.pop()
            if v in component:continue
            component.add(v);stack.extend(adj[v]-component)
        remaining-=component
        ordered=sorted(component,key=lambda u:-sum(arc(A,u,v) for v in component-{u}))
        need(all(arc(A,u,v) for i,u in enumerate(ordered) for v in ordered[i+1:]),'clone component not transitive')
        need(all(len({arc(A,u,w) for u in component})==1 for w in set(range(n))-component),'clone component not module')
        blocks.append(ordered)
    # Every transitive module's consecutive vertices are global twins; hence
    # it lies in one component. Verified components therefore are maximal.
    return {'twin_pairs':twins,'ordered_maximal_transitive_modules':blocks}
def state_record(A):
    return {'rows':list(A),'ports_at_0':neighbourhoods(A,0),
            'ports_at_z':neighbourhoods(A,len(A)-1),'clones':clone_data(A)}
def source_equations(n,I,opened,chord):
    raw=gf.raw_matrix(n);eq=gf.constraints(n,raw)
    eq += [raw[0][v]^int(v not in I) for v in range(1,n-1)]
    eq += [raw[1][n-2]^int(opened),raw[2][n-2]^int(chord)]
    return raw,eq

def enumerate_models(raw,eq):
    basis=gf.eliminate(eq);need(basis is not None,'residual unexpectedly empty')
    d=len(raw)*(len(raw)-1)//2;free=sorted(set(range(1,d+1))-set(basis));out=[]
    for bits in range(1<<len(free)):
        assignment=sum(((bits>>j)&1)<<v for j,v in enumerate(free))
        for pivot,equation in sorted(basis.items()):
            if (equation&1)^((equation&~(1<<pivot)&assignment).bit_count()&1):assignment|=1<<pivot
        need(all(((e&1)^((e&assignment).bit_count()&1))==0 for e in eq),'model violates equations')
        A=tuple(sum(((raw[u][v]&1)^((raw[u][v]&assignment).bit_count()&1))<<v for v in range(len(raw)) if u!=v) for u in range(len(raw)))
        out.append(A)
    need(len(set(out))==len(out),'duplicate affine model')
    return free,out

def audit_source_cell(n,I,opened,chord):
    name,qs=bank.words(n,set(I),opened,chord)
    raw,eq=source_equations(n,I,opened,chord)
    basis=gf.eliminate(eq);need(basis is not None,'empty base cell')
    src=[[gf.reduce_form(f,basis) for f in row] for row in raw]
    at=src;bad=[];xor=set()
    for j,q in enumerate(qs):
        need(sorted(q)==list(range(n)),'old word not permutation')
        need(all(at[u][v]==1 for u,v in zip(q,q[1:])),'old join not forced')
        if gf.eliminate(gf.blue_equations(at,q)) is not None:bad.append(j+1)
        at=gf.k_move(at,q);xor^=edge_set(q)
    need(at==gf.k_move(src,list(range(n))) and xor==edge_set(list(range(n))),'old ledger target')
    return name,qs,bad

def compact_words(n,I,bit13=0,bit14=0):
    if n==12:
        if not bit13:W=(1,3,4,2,5,6,7,8,9,10,11)
        elif not bit14:W=(1,4,5,6,7,2,3,8,9,10,11)
        elif 2 in I:W=(1,6,7,8,2,3,4,5,9,10,11)
        else:W=(1,6,2,3,4,5,7,8,9,10,11)
    elif 3 in I:W=(1,3,4,5,2,6,7,8,9,10,11,12,13)
    else:W=(1,3,4,5,6,2,7,8,9,10,11,12,13)
    return [tuple(range(1,n))+(0,),W+(0,),(0,)+W]


def main():
    found=json.loads((HERE/'found.json').read_text());need(len(found)==20,'twenty certificates expected')
    indexed={tuple(r['source_rows']):r for r in found};need(len(indexed)==20,'duplicate certificate rows')
    failed=[];counts={}
    for n in (12,14):
        tested=deferred=0
        for I in combinations(range(2,n-1),3):
            for opened,chord in ((True,False),(True,True),(False,False),(False,True)):
                if not opened and {n-3,n-2}<=set(I):deferred+=1;continue
                name,qs,bad=audit_source_cell(n,I,opened,chord);tested+=1
                if bad:
                    need(len(bad)==1,'unexpected multiple old failures')
                    failed.append((n,I,opened,chord,bad[0]))
        counts[n]={'all_branch_cells':tested+deferred,'base_cells':tested,'residual_endpair_cells':deferred}
    need(tuple(sorted(failed))==tuple(sorted(SIX)),'wrong complete six-failure list')
    need(counts=={12:{'all_branch_cells':336,'base_cells':322,'residual_endpair_cells':14},14:{'all_branch_cells':660,'base_cells':642,'residual_endpair_cells':18}},'partition counts')
    cells=[];certificates=[];seen=set();compact_cells=[]
    for n,I,opened,chord,step in SIX:
        name,old,bad=audit_source_cell(n,I,opened,chord)
        raw,eq=source_equations(n,I,opened,chord);at=raw
        for q in old[:step-1]:at=gf.k_move(at,q)
        eq+=gf.blue_equations(at,old[step-1])
        free,models=enumerate_models(raw,eq)
        branches=([(0,None),(1,0),(1,1)] if n==12 else [(None,None)])
        for bit13,bit14 in branches:
            equations=eq[:]
            if bit13 is not None:equations.append(raw[1][3]^bit13)
            if bit14 is not None:equations.append(raw[1][4]^bit14)
            basis=gf.eliminate(equations);need(basis is not None,'empty compact repair cell')
            symbolic=[[gf.reduce_form(f,basis) for f in row] for row in raw]
            words=compact_words(n,I,bit13,bit14)
            gf.replay_family(symbolic,words)
            compact_cells.append({'n':n,'ports':list(I),'bit_T13':bit13,'bit_T14':bit14,
                                  'affine_dimension':n*(n-1)//2-len(basis),'words':words})
        need(len(free)==(3 if n==12 else 0),'unexpected residual dimension')
        need(all(A in indexed for A in models),'missing affine hostile completion')
        cells.append({'n':n,'ports':list(I),'middle_open':opened,'x_to_c':chord,
                      'old_branch':name,'old_words':old,'old_blue_step':step,
                      'residual_free_raw_variables':free,'residual_affine_dimension':len(free),
                      'residual_tournament_count':len(models)})
        for A in sorted(models):
            seen.add(A);rec=indexed[A];P=tuple(range(n));B=move(A,P)
            need(all(A[v]>>v&1==0 for v in range(n)) and all(arc(A,u,v)^arc(A,v,u)==1 for u,v in combinations(range(n),2)),'not tournament')
            need(directed(A,P) and witness(A,P) is None and arc(A,n-1,0),'not closed blue source')
            need(A[0].bit_count()==n-5 and A[n-1].bit_count()==4,'endpoint scores')
            need(set(neighbourhoods(A,0)['in'])==set(I)|{n-1},'actual incoming ports')
            need(set(neighbourhoods(A,n-1)['out'])=={0}|{n-1-v for v in I},'reflected outgoing ports')
            flags=[];at=A
            for q in old:
                need(directed(at,q),'old literal join');flags.append(witness(at,q) is None);at=move(at,q)
            need(flags==[j==step-1 for j in range(3)] and at==B,'old literal blue flags')
            qs=[tuple(q) for q in rec['words']];need(len(qs)==3,'three moves')
            need(qs==compact_words(n,I,arc(A,1,3),arc(A,1,4)),'compact recipe mismatch')
            need(qs[0]==tuple(range(1,n))+(0,),'cyclic-shift bridge first word')
            need(qs[1][-1]==0 and qs[2]==(0,)+qs[1][:-1] and qs[1][0]==1 and qs[1][-2]==n-1,'bridge shape')
            states=[state_record(A)];steps=[];xor=set();at=A;allpairs=set(combinations(range(n),2))
            for j,q in enumerate(qs):
                need(sorted(q)==list(range(n)) and directed(at,q),'repair not legal Hamilton word')
                black=witness(at,q);need(black is not None,'repair blue')
                nxt=move(at,q);xor^=edge_set(q)
                before=states[-1]['clones'];after=clone_data(nxt)
                oldtw={tuple(e) for e in before['twin_pairs']};newtw={tuple(e) for e in after['twin_pairs']}
                overlaps=[[sorted(set(u)&set(v)) for v in after['ordered_maximal_transitive_modules']] for u in before['ordered_maximal_transitive_modules']]
                steps.append({'word':list(q),'blackness_witness':black,
                              'word_endpoint_ports':[neighbourhoods(at,q[0]),neighbourhoods(at,q[-1])],
                              'flipped_unordered_pairs':[list(e) for e in sorted(allpairs-edge_set(q))],
                              'clone_transport':{'new_twin_pairs':[list(e) for e in sorted(newtw-oldtw)],
                                                 'lost_twin_pairs':[list(e) for e in sorted(oldtw-newtw)],
                                                 'literal_block_intersections':overlaps}})
                states.append(state_record(nxt));at=nxt
            need(xor==edge_set(P) and at==B,'repair literal target or edge XOR')
            certificates.append({'n':n,'ports':list(I),'source_rows':list(A),
                                 'old_blue_flags':flags,'steps':steps,'states':states,
                                 'literal_target_rows':list(B),'edge_XOR_equals_P':True})
    need(seen==set(indexed),'certificate outside complete residual universe')
    # Positive merger/split control: a triangle module is not a transitive clone.
    transitive=(6,4,0);cycle=move(transitive,(0,1,2))
    need(cycle==(2,4,1) and move(cycle,(0,1,2))==transitive,'clone merge/split literal K control')
    tc=clone_data(transitive);cc=clone_data(cycle)
    need(tc['ordered_maximal_transitive_modules']==[[0,1,2]],'transitive clone chain lost')
    need(cc['ordered_maximal_transitive_modules']==[[0],[1],[2]],'cyclic module wrongly contracted')
    need(all(len(b)==1 for c in certificates for s in c['states'] for b in s['clones']['ordered_maximal_transitive_modules']),'new clone unexpectedly appears in residual ledger')
    # A canonical real obstruction is not accidentally accepted as black-connected.
    H=(6,4,8,19,103,79,15);T=move(H,tuple(range(7)))
    need(T==(122,124,121,112,32,64,16),'order-seven hostile changed')
    def paths(A):
        def extend(q,remaining):
            if not remaining:yield tuple(q)
            for v in sorted(remaining):
                if not q or arc(A,q[-1],v):yield from extend(q+[v],remaining-{v})
        return list(extend([],set(range(len(A)))))
    hp=paths(T);need(len(hp)==9 and all(witness(T,q) is None for q in hp),'genuine hostile black obstruction lost')
    # Deliberately wrong orientation and word mutations must fail local checks.
    row=indexed[next(iter(indexed))];A=tuple(row['source_rows']);q=tuple(row['words'][0])
    need(not directed(A,q[::-1]),'reversed directed word hostile')
    need(move(A,q)!=move(A,tuple(range(len(A)))),'one shift wrongly equals target')
    receipt={'status':'PASS','method':'pinned inherited Gaussian elimination; local complete affine enumeration and literal rows',
             'dependency_sha256':PINS,'base_partition':counts,'six_failed_cells':cells,
             'complete_residual_sources':20,'direct_three_black_detours':20,'compact_affine_repair_cells':compact_cells,
             'scope':'six whole affine cells repaired by old ledger off the blue locus and complete 20-source certificates on it; inherited residual h families still required for all-order conclusion',
             'hostile_controls':{'order7_target_HP':9,'order7_target_black_HP':0,'reversed_word_rejected':True,'one_shift_target_rejected':True,
                                 'clone_merge_split_control':{'transitive_rows':transitive,'cycle_rows':cycle,'transitive_clones':tc,'cycle_clones':cc},
                                 'all_80_repair_states_clone_free':True},
             'certificates':certificates}
    print(json.dumps(receipt,sort_keys=True,indent=2))

if __name__=='__main__':main()
