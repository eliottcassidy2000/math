"""Exact h3-forward certificates and independent permutation/degree controls.

Uses the previously audited signed-equation engine, then literal arc sets.
All gates are explicit and remain enabled under python -O.
"""
import importlib.util
import json
from pathlib import Path


PRIOR = Path(__file__).resolve().parents[2] / 'ports_regular_actions_hn_20261025/ports/check.py'
spec = importlib.util.spec_from_file_location('signed_source', PRIOR)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
need = base.need


def words(n):
    u, r = tuple(range(3,n-3)), tuple(range(n-3,n))
    return [(0,1)+u+r+(2,), (1,)+u+(0,)+r+(2,), (1,2)+u+(0,)+r]


def action(q,n):
    rev = dict(zip(q,reversed(q)))
    return {v:rev[n-1-v] for v in range(n)}


def cycles(g):
    unseen, result = set(g), []
    while unseen:
        at, cycle = min(unseen), []
        while at in unseen:
            unseen.remove(at)
            cycle.append(at)
            at = g[at]
        result.append(cycle)
    return result


def geometry(n):
    q1,q2,q3 = words(n)
    g3 = action(q3,n)
    expected3 = {0:1,1:2,2:3,n-5:0,n-4:n-3,n-3:n-2,n-2:n-1,n-1:4}
    expected3.update({v:v+2 for v in range(3,n-5)})
    need(g3 == expected3, ('last permutation',n))
    final_cycles = cycles(g3)
    outgoing = {1,3,n-3,n-2,n-1}
    if n%2:
        need(len(final_cycles)==1, ('last odd transitivity',n))
    else:
        orbit = set(final_cycles[0])
        need(len(orbit)==n//2 and orbit&outgoing=={1,3}, ('last even induced score',n))

    g = action(q2,n)
    expected = {0:3,1:4,2:5,n-7:0,n-6:n-3,n-5:n-2,n-4:n-1,
                n-3:1,n-2:2,n-1:6}
    expected.update({v:v+4 for v in range(3,n-7)})
    need(g == expected, ('middle permutation',n))
    delta = {2:-1,n-1:1}
    increment = {v:-delta.get(n-1-v,0)-delta.get(g[v],0) for v in range(n)}
    need(increment == {v:-1 if v in {0,n-4} else 1 if v in {n-3,n-2} else 0
                       for v in range(n)}, ('middle degree increment',n))
    cs = cycles(g)
    if n%2:
        need(len(cs)==1, ('middle odd transitivity',n))
        offset,total = 0,0
        for v in cs[0]:
            total += offset
            offset += increment[v]
        need(offset==0 and total==(0 if n%4==1 else -n), ('odd degree total',n,total))
    elif n%4==2:
        expected_cycle = (0,)+tuple(range(3,n-6,4))
        need(tuple(cs[0])==expected_cycle, ('nonclosing degree cycle',n))
        need(sum(increment[v] for v in cs[0])==-1, ('degree loses one on return',n))
    else:
        t = n//4
        expected_cycle = (0,)+tuple(range(3,n-4,4))+(n-2,2)+tuple(range(5,n-6,4))
        need(tuple(cs[0])==expected_cycle and len(cs[0])==2*t, ('induced cycle',n))
        orbit = set(cs[0])
        reflected = {n-1-v for v in orbit}
        need(reflected==set(range(n))-orbit, ('reflected complementary orbit',n))
        need({g[v] for v in reflected}==reflected, ('complement invariant',n))
        f = {(1,2),(2,3),(1,3),(2,n-1)}
        need({e for e in f if set(e)<=orbit}=={(2,3)}, ('one internal defect',n))
        need(not any(set(e)<=reflected for e in f), ('no complement defect',n))
        need(orbit&outgoing=={3,n-2}, ('induced actual outdegree',n))
        offset,total = 0,0
        for v in cs[0]:
            total += offset
            offset += -1 if g[v]==3 else 1 if g[v]==2 else 0
        need(offset==0 and total==-t, ('induced offset total',n))
        need(4*t+total != t*(2*t-1), ('induced contradiction requires n>=12',n))


def move(arcs,q):
    path = set(zip(q,q[1:]))
    need(path<=arcs, 'literal joins')
    return path | {(v,u) for u,v in arcs-path}


def is_blue(arcs,q):
    rev = dict(zip(q,reversed(q)))
    return all((rev[v],rev[u]) in arcs for u,v in arcs)


def certify(n):
    qs = words(n)
    src,_ = base.symbolic(n,{3,n-3,n-2},qs,[(1,3)],True)
    for seed in (0,1917):
        original = src.model(seed)
        b = move(original,tuple(range(n)))
        s = move(original,qs[0])
        f = {base.pair(u,v) for u,v in b^s}
        need(f=={(1,2),(2,3),(1,3),(2,n-1)}, ('literal four-pair defect',n))
        for v in range(n):
            change=sum(a==v for a,_ in s)-sum(a==v for a,_ in b)
            need(change==(-1 if v==2 else 1 if v==n-1 else 0), ('literal endpoint degrees',n,v))
        at=original
        for q in qs:
            need(not is_blue(at,q), ('literal black',n,q))
            at=move(at,q)
        need(at==b, ('literal labelled target',n))


def blue_hostile(n,step):
    qs=words(n)
    src,report=base.symbolic(n,{3,n-3,n-2},qs,[(1,3)],False)
    need(not report['illegal'] and report['target'], ('hostile must have exact legal ledger',n))
    need(report['blue_feasible'][step], ('hostile blue step feasible',n))
    flips=set()
    for q in qs[:step]:
        flips ^= set(src.edges)^base.pathset(q)
    src.anti(qs[step],flips)
    need(src.ok, ('hostile source exists',n))
    arcs=src.model(0)
    at=arcs
    for q in qs[:step]:
        at=move(at,q)
    need(is_blue(at,qs[step]), ('literal hostile is blue',n))
    return {'n':n,'blue_move':step+1,'source_rows':[
        sum(1<<b for a,b in arcs if a==v) for v in range(n)]}


def main():
    orders=list(range(12,101))+[120,160]
    for n in orders:
        certify(n)
    for n in range(12,501):
        geometry(n)
    original=[words(12)[0],(1,0)+tuple(range(3,12))+(2,),
              (1,2,0)+tuple(range(3,12))]
    old,report=base.symbolic(12,{3,9,10},original,[(1,3)],False)
    need([1,1,0] in report['illegal'] and [2,2,0] in report['illegal']
         and not report['target'], 'old extracted-predecessor hostile')
    print(json.dumps({'status':'PASS','complete_affine_orders':orders,
        'complete_affine_cells':len(orders),'literal_models_per_cell':2,
        'permutation_and_degree_orders':[12,500],
        'lower_order_blue_hostiles':[blue_hostile(8,1),blue_hostile(10,2)],
        'old_template_hostile':report,
        'scope':'all h3,1->3 endpoint-four sources n>=12; exact labelled target'},indent=2))


if __name__=='__main__':
    main()
