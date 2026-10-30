#!/usr/bin/env python3
"""Complete R344(73) root-case list; bounded native shared-edge SAT probes."""
from pathlib import Path
from itertools import combinations
from threading import Timer
import json,time
from ramsey_checks import decode,need
P=Path(__file__).parent
def case_list():
    cases=[]
    for r,patterns in ((15,((28,29),)),(16,((27,29),(28,28))),(17,((26,29),(27,28)))):
        lines=(P/'data'/f'r44_{r}.g6').read_text().splitlines()
        for index,line in enumerate(lines):
            for b,g in patterns:
                cases.append(dict(id=len(cases),degrees=[r,b,g],red_port_graph_index=index,
                                  red_port_graph6=line,global_red_cap=r))
    need(len(cases)==646,'case count')
    return cases
def cliques4(adj):
    n=len(adj)
    for a in range(n):
        bs=adj[a]&~((1<<(a+1))-1)
        while bs:
            bbit=bs&-bs;bs-=bbit;b=bbit.bit_length()-1;cs=bs&adj[b]
            while cs:
                cbit=cs&-cs;cs-=cbit;c=cbit.bit_length()-1;ds=cs&adj[c]
                while ds:
                    dbit=ds&-ds;ds-=dbit
                    yield a,b,c,dbit.bit_length()-1
def probe(case,seconds=25):
    from pysat.solvers import Minicard
    n=73;edges=list(combinations(range(n),2));ids={e:3*i+1 for i,e in enumerate(edges)}
    def var(i,j,c):return ids[tuple(sorted((i,j)))]+c
    S=Minicard();start=time.monotonic();initial=0
    for i,j in edges:
        vs=[var(i,j,c) for c in range(3)]
        S.add_clause(vs);S.add_atmost(vs,1)
    for i in range(n):
        for c,cap in enumerate((case['global_red_cap'],29,29)):
            S.add_atmost([var(i,j,c) for j in range(n) if i!=j],cap)
    # Pick a root of maximum red degree, and keep all native edge identities.
    r,b,g=case['degrees']
    for c,range_ in enumerate((range(1,r+1),range(r+1,r+b+1),range(r+b+1,n))):
        for j in range_:S.add_clause([var(0,j,c)])
    port=decode(case['red_port_graph6'].encode())
    for i,j in combinations(range(r),2):S.add_clause([var(i+1,j+1,1 if port[i]>>j&1 else 2)])
    for tri in combinations(range(n),3):S.add_clause([-var(i,j,0) for i,j in combinations(tri,2)]);initial+=1
    added=0;rounds=0;remaining=None;verdict='UNKNOWN';timer=Timer(seconds,S.interrupt);timer.start()
    try:
        while time.monotonic()-start<seconds:
            S.conf_budget(5000);sat=S.solve_limited(expect_interrupt=True)
            if sat is None:break
            if not sat:
                verdict='SOLVER_UNSAT_UNCERTIFIED';break
            rounds+=1;positive={x for x in S.get_model() if x>0};adj=[[0]*n for _ in range(3)]
            for i,j in edges:
                chosen=[c for c in range(3) if var(i,j,c) in positive];need(len(chosen)==1,'edge color')
                c=chosen[0];adj[c][i]|=1<<j;adj[c][j]|=1<<i
            bad=[]
            for c in (1,2):
                for clique in cliques4(adj[c]):bad.append([-var(i,j,c) for i,j in combinations(clique,2)])
            remaining=len(bad)
            if not bad:
                verdict='SAT_CANDIDATE_REQUIRES_INDEPENDENT_AUDIT'
                (P/f'ramsey_candidate_{case["id"]}.json').write_text(json.dumps(adj)+'\n');break
            for clause in bad:S.add_clause(clause)
            added+=len(bad)
        stats=S.accum_stats()
    finally:
        timer.cancel();S.delete()
    return dict(case=case['id'],status=verdict,degree_caps=[case['global_red_cap'],29,29],
                model_rounds=rounds,added_native_K4_clauses=added,last_model_K4_count=remaining,
                seconds_budget=seconds,conflicts_per_solve=5000,solver_stats=stats)
def run():
    cases=case_list();(P/'ramsey_root_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    selected=[cases[0],cases[1],*cases[640:]];out=[]
    for c in selected:
        z=probe(c);out.append(z);print(z,flush=True)
        (P/'ramsey_root_probes.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
