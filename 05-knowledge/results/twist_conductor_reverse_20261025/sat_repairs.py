#!/usr/bin/env python3
"""Witness search only; accepted words are replayed by standard-library code."""
import json
import sys
from pathlib import Path
import z3
from search_repairs import arc, blue, directed, edge, move


def solve(T,timeout=60000,first=None):
    n=len(T);solver=z3.Solver();solver.set(timeout=timeout)
    q=[[z3.Int(f'q{s}_{i}') for i in range(n)] for s in range(3)]
    pos=[[z3.Int(f'p{s}_{i}') for i in range(n)] for s in range(3)]
    Q=[z3.Function(f'Q{s}',z3.IntSort(),z3.IntSort()) for s in range(3)]
    F=[z3.Function(f'F{s}',z3.IntSort(),z3.IntSort(),z3.BoolSort()) for s in range(3)]
    E={}
    for s in range(3):
        solver.add(z3.Distinct(q[s]),z3.Distinct(pos[s]))
        for i in range(n):
            solver.add(q[s][i]>=0,q[s][i]<n,pos[s][i]>=0,pos[s][i]<n,Q[s](i)==q[s][i])
            for v in range(n):solver.add(z3.Implies(q[s][i]==v,pos[s][v]==i))
        for u in range(n):
            solver.add(z3.Not(F[s](u,u)))
            for v in range(u+1,n):
                E[s,u,v]=z3.Or(pos[s][u]-pos[s][v]==1,pos[s][v]-pos[s][u]==1)
                value=z3.BoolVal(bool(arc(T,u,v))) if s==0 else z3.Xor(F[s-1](u,v),z3.Not(E[s-1,u,v]))
                solver.add(F[s](u,v)==value,F[s](v,u)==z3.Not(value))
        solver.add(*[F[s](q[s][i],q[s][i+1]) for i in range(n-1)])
        solver.add(z3.Or(*[F[s](u,v)!=F[s](Q[s](n-1-pos[s][v]),Q[s](n-1-pos[s][u]))
                            for u in range(n) for v in range(u+1,n)]))
    for u in range(n):
        for v in range(u+1,n):solver.add(z3.Xor(z3.Xor(E[0,u,v],E[1,u,v]),E[2,u,v])==z3.BoolVal(v==u+1))
    if first is not None:solver.add(*[q[0][i]==v for i,v in enumerate(first)])
    result=solver.check()
    if result!=z3.sat:return {'solver_result':str(result),'reason':solver.reason_unknown() if result==z3.unknown else None}
    model=solver.model();words=[[model.eval(v).as_long() for v in row] for row in q]
    current=T;target=move(T,tuple(range(n)));xor=frozenset()
    for word in words:
        if sorted(word)!=list(range(n)) or not directed(current,word) or blue(current,word):raise RuntimeError('invalid solver witness')
        xor^=edge(word);current=move(current,word)
    if current!=target or xor!=edge(tuple(range(n))):raise RuntimeError('wrong solver target')
    return {'solver_result':'sat','words':words,'literal_replay':'PASS'}


def main():
    base=Path(__file__).resolve().parent
    rows=json.loads((base.parent/'port_residual_integral_20261025/low_orders.normal.json').read_text())['failed_ledger_hostiles']
    index=int(sys.argv[1]) if len(sys.argv)>1 else 0
    row=[r for r in rows if r['n']==14][index]
    first=None
    if len(sys.argv)>3:
        first=json.loads((base/'sat_repair_0.json').read_text())['words'][0]
    result=solve(tuple(row['source_rows']),int(sys.argv[2]) if len(sys.argv)>2 else 60000,first)
    result.update({'n':14,'ports':row['ports'],'source_rows':row['source_rows'],'z3_version':z3.get_version_string(),
                   'scope':'literal three-black detour search; unknown is not an obstruction'})
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
