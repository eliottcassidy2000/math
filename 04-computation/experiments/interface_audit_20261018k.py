#!/usr/bin/env python3
"""Independent numerical, support-set, and marked-minor hostile controls."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from pysat.solvers import Solver


def sums(A,B,q):return {(a+b)%q for a in A for b in B}


def compile_support(Ds,q):
    """Forward Horn support propagation, independently reconstructed."""
    top=0;clauses=[]
    def row():
        nonlocal top
        out=list(range(top+1,top+q+1));top+=q;return out
    vectors=[]
    for D in Ds:
        v=row();vectors.append(v)
        clauses.extend([[v[x]] for x in D])
    acc=vectors[0]
    for B in vectors[1:]:
        C=row()
        for x in range(q):
            for y in range(q):clauses.append([-acc[x],-B[y],C[(x+y)%q]])
        acc=C
    clauses.append([-acc[0]])
    with Solver(name='m22',bootstrap_with=clauses) as sol:return sol.solve()


def support_audit():
    controls=[]
    for Ds in ([{1}]*6,[{1}]*5+[{2}],[set(range(6)),{1},{2}],
               [{0,2,4},{0,2,4},{1,3,5}],[{0,2,4}]*3,[{0,1}]*5,[{0,1}]*4+[{1}],
               [{0,3},{0,3},{1,4}],[{0,3},{0,3},{0,3}]):
        literal={sum(xs)%6 for xs in itertools.product(*Ds)}
        got=compile_support(Ds,6);assert got==(0 not in literal)
        controls.append(dict(sets=[sorted(D) for D in Ds],extends=got))
    masks=[{x for x in range(6) if m>>x&1} for m in range(1,64)]
    # Exhaustive state reduction over ALL 63^6 support sequences, retaining
    # the largest number of rich factors for each exact reachable sumset.
    state={frozenset({0}):0};maxima=[]
    for arity in range(1,7):
        new={}
        for S,n in state.items():
            for D in masks:
                T=frozenset(sums(S,D,6));r=n+any((x+1)%6 in D for x in D)
                new[T]=max(new.get(T,-1),r)
        state=new;maxima.append(max(n for S,n in state.items() if 0 not in S))
    assert maxima==[1,2,3,4,4,4]
    # Independent literal enumeration in the smaller alphabet, including
    # six-way circuits, checks the Minkowski-sum recurrence.
    checked=0
    m3=[{x for x in range(3) if m>>x&1} for m in range(1,8)]
    for arity in range(2,7):
        for Ds in itertools.product(m3,repeat=arity):
            acc={0}
            for D in Ds:acc=sums(acc,D,3)
            assert acc=={sum(xs)%3 for xs in itertools.product(*Ds)};checked+=1
    return dict(sat_controls=controls,exhaustive_ternary_sequences=checked,
                exact_six_colour_mask_states=len(state),max_rich_factors_avoiding_zero=maxima)


def minor_control():
    vertices=list(itertools.combinations(range(5),2))
    E=[(i,j) for i in range(10) for j in range(i) if set(vertices[i]).isdisjoint(vertices[j])]
    assert len(E)==15
    subdiv=[edge for k,(i,j) in enumerate(E) for edge in ((i,10+k),(10+k,j))]
    contraction={i:i for i in range(10)}|{10+k:i for k,(i,j) in enumerate(E)}
    recovered={tuple(sorted((contraction[i],contraction[j]))) for i,j in subdiv if contraction[i]!=contraction[j]}
    assert recovered=={tuple(sorted(e)) for e in E}
    count=proper=0
    def visit(pat):
        nonlocal count,proper
        if len(pat)==10:
            full=pat+[next(c for c in range(6) if c not in (pat[i],pat[j])) for i,j in E]
            assert all(full[i]!=full[j] for i,j in subdiv)
            proper+=all(pat[i]!=pat[j] for i,j in E);count+=1;return
        for c in range(min(5,max(pat,default=-1)+1)+1):visit(pat+[c])
    visit([])
    assert all(([0]*10+[1]*15)[i]!=([0]*10+[1]*15)[j] for i,j in subdiv)
    return dict(ports=10,original_edges=15,subdivision_vertices=25,subdivision_edges=30,
                six_colour_canonical_assignments=count,original_extending=proper,subdivision_extending=count,
                hostile='All ten ports equal: excluded by Petersen; extends to its once-subdivided graph.',
                scope='Abstract marked graphs. Not a geometric unit-distance substitution.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,required=True);args=ap.parse_args();p=args.directory
    a=json.loads((p/'source_tensor.json').read_text());b=json.loads((p/'source_concentration_everyN.json').read_text())
    lookup={r['N']:r for r in a['rows']};checked=0
    for r in b['dyadic_checks']:
        if r['N'] not in lookup:continue
        x=lookup[r['N']];g=np.array(x['prefix_grams']);n=len(g)
        assert r['packet_source_energy']==x['packet_source_gram'][0][0]
        assert r['prefix_actual'][:n]==g[:,0,0].tolist()
        assert r['prefix_trace'][:n]==np.trace(g,axis1=1,axis2=2).tolist()
        assert r['prefix_actual'][n:]==[g[-1,0,0]]*(len(r['prefix_actual'])-n);checked+=2*n+1
    boundary=json.loads((p/'boundary_automaton.json').read_text());audit=json.loads((p/'boundary_everyN.json').read_text())
    assert boundary['maximum_D0']==audit['maximum']==1428 and boundary['minimum_D0']==audit['minimum']==-72
    assert boundary['first_maximum_N']==audit['maximum_N']==66475
    assert boundary['first_minimum_N']==audit['minimum_N']==171
    for N,D in audit['dyadics']:
        if N in lookup:
            g=np.array(lookup[N]['prefix_grams'][0]);assert D==3*g[0,0]-np.trace(g)
    old=json.loads((p/'mixed_quads.json').read_text());new=json.loads((p/'support_compiler_control.json').read_text())
    assert old['coordinate_sizes']==new['coordinate_sizes'] and old['blocked']==new['blocked']==[]
    assert sorted(old['selected_families'].values())==sorted(t['cartesian_constraints'] for t in new['families'].values())
    # Recheck OLD, separately generated literal-edge witnesses using the NEW
    # support compiler's exact group indices. No solver involved here.
    crosschecks=0
    for w in old['witnesses']:
        for group in new['support_groups']:
            acc={0}
            for axis,ps in group['factors']:
                D={(w['colours'][axis][j]-w['colours'][axis][i])%6 for i,j in ps};acc=sums(acc,D,6)
            assert 0 not in acc;crosschecks+=1
    increases=[]
    for N in (4096,5997,1048576):
        g=np.array(lookup[N]['prefix_grams']);D=3*g[:,0,0]-np.trace(g,axis1=1,axis2=2);up=np.diff(D)
        increases.append(dict(N=N,max_local_increase=int(up.max()),at_depth=int(up.argmax()+1)))
    out=dict(status='FINITE-EXACT independent audit',source_engine_comparisons=checked,
             old_witness_new_support_checks=crosschecks,source_local_monotonicity_hostiles=increases,
             support_compiler=support_audit(),marked_minor=minor_control())
    (p/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
