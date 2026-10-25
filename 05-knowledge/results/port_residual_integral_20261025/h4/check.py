"""Exact h4 residual cells plus independent literal K-move replay."""
from pathlib import Path
from itertools import combinations
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'ports_regular_actions_hn_20261025/ports/check.py'
spec = importlib.util.spec_from_file_location('h4_affine_source', SOURCE)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def words(n):
    V = tuple(range(5,n))
    return [(0,)+V+(1,2,3,4), (2,0)+V+(1,3,4), (2,0,1,3,4)+V]


def edges(q):
    return {frozenset((a,b)) for a,b in zip(q,q[1:])}


def move(arcs, q):
    keep = edges(q)
    return {(a,b) if frozenset((a,b)) in keep else (b,a) for a,b in arcs}


def literal(arcs, qs):
    n = len(qs[0])
    P = tuple(range(n))
    target = move(arcs,P)
    for j,q in enumerate(qs):
        need(sorted(q) == list(range(n)), ('Hamilton word',n,j))
        need(all((a,b) in arcs for a,b in zip(q,q[1:])), ('literal join',n,j))
        reverse = dict(zip(q,q[::-1]))
        need(any((reverse[b],reverse[a]) not in arcs for a,b in arcs), ('literal blue',n,j))
        arcs = move(arcs,q)
    need(arcs == target, ('literal target',n))


def geometry(n):
    P = tuple(range(n))
    q1,q2,q3 = words(n)
    X = {frozenset(e) for e in ((0,1),(4,5),(0,5),(1,n-1))}
    Y = {frozenset(e) for e in ((1,2),(2,3),(0,2),(1,3))}
    need(edges(q1) == edges(P)^X, ('first rectangle side',n))
    need(edges(q2) == edges(P)^X^Y, ('second rectangle side',n))
    need(edges(q3) == edges(P)^Y, ('third rectangle side',n))
    reverse = [dict(zip(q,q[::-1])) for q in (q1,q2,q3)]
    need((reverse[0][3],reverse[0][0]) == (5,4), ('Q1 witness image',n))
    need((reverse[1][n-1],reverse[1][3]) == (6,0), ('Q2 witness image',n))
    need((reverse[2][3],reverse[2][0]) == (n-4,n-2), ('Q3 witness image',n))


def main():
    cells = models = 0
    for n in range(12,121):
        qs = words(n)
        src,report = base.symbolic(n,{4,n-3,n-2},qs,[(3,1)],True)
        need(src.forced(0,3,1) and src.forced(5,4,0), ('Q1 source values',n))
        need(src.forced(3,n-1,1) and src.forced(6,0,0), ('Q2 source values',n))
        need(frozenset((3,n-1)) not in edges(qs[0]) and
             frozenset((6,0)) not in edges(qs[0]), ('Q2 omitted pairs',n))
        need(src.forced(0,3,1) and src.forced(n-4,n-2,0), ('Q3 source values',n))
        for seed in (17,233,30030):
            literal(src.model(seed),qs)
            models += 1
        cells += 1
    for n in range(12,1001):
        geometry(n)
    old = base.source(12,{4,9,10},[(3,1)])
    need(old.ok and old.forced(0,4,0), 'inherited false first join')
    extra, report = base.symbolic(12,{4,5,9,10},words(12),[(3,1)],False)
    need(extra.ok and extra.forced(0,5,0) and [0,0,5] in report['illegal'],
         'outside-scope added-port hostile')
    print(json.dumps({'status':'PASS','complete_affine_cells':cells,
        'affine_orders':[12,120],'literal_models':models,'seeds':[17,233,30030],
        'edge_rectangle_and_reversal_orders':[12,1000],
        'source':'closed blue P, I={4,n-3,n-2},3->1; all source completions',
        'hostiles':{'old_ledger':'0->4 is false',
                    'extra_port_5':'first new join 0->5 becomes false'},
        'scope':'specified three-word certificate only; no Hamiltonian path census'},
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
