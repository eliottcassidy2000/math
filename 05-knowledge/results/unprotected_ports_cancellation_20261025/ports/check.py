"""Exact new-ledger controls; every gate remains active under python -O."""
import importlib.util
import json
from pathlib import Path


PRIOR = Path(__file__).resolve().parents[2] / 'ports_regular_actions_hn_20261025/ports'
spec = importlib.util.spec_from_file_location('prior_ports', PRIOR / 'check.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
spec2 = importlib.util.spec_from_file_location('prior_boundary', PRIOR / 'boundary_check.py')
old = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(old)
need = base.need


def words(n, h, forward):
    if forward:
        u, r = tuple(range(3, h)), tuple(range(h, n))
        return [(0, 1) + u + r + (2,), (1,) + u + (0,) + r + (2,),
                (1, 2) + u + (0,) + r]
    v = tuple(range(4, n))
    return [(0,) + v + (1, 2, 3), (2, 0) + v + (1, 3), (2, 0, 1, 3) + v]


def one_cell(n, interior, forward):
    k, h = len(interior) + 1, min(interior)
    if n == 2*k+3 and interior == set(range(k+3, n-1)):
        # The existing theorem covers this entire source, split by 2 versus 4.
        old.cell(n, interior, False)
        old.cell(n, interior, True)
        return 1
    qs = words(n, h, forward)
    guard = [(1, 3)] if forward else [(3, 1)]
    src, _ = base.symbolic(n, interior, qs, guard, True)
    for seed in (0, 2261):
        base.literal(src.model(seed), qs)
    return 0


def geometry(n):
    q = words(n, 5, False)[1]
    r = dict(zip(q, reversed(q)))
    g = {v: r[n-1-v] for v in range(n)}
    un, cycles = set(range(n)), []
    while un:
        at, cycle = min(un), []
        while at in un:
            un.remove(at)
            cycle.append(at)
            at = g[at]
        cycles.append(cycle)
    defect = {3: -1, n-1: 1}
    inc = {v: -defect.get(n-1-v, 0)-defect.get(g[v], 0) for v in range(n)}
    need(inc == {v: (-1 if v in {0,n-5} else 1 if v in {n-4,n-3} else 0)
                 for v in range(n)}, ('degree increments', n))
    if n % 2:
        need(len(cycles) == 1, ('odd transitivity', n))
        offset, total = 0, 0
        for v in cycles[0]:
            total += offset
            offset += inc[v]
        need(offset == 0, 'degree cycle closes')
        need(total == (-n if n % 4 == 1 else 0), ('offset sum', n, total))
    else:
        if n % 4 == 0:
            need(sorted(map(len, cycles)) == [n//2,n//2], ('two cycles',n))
        else:
            t = (n-2)//4
            need(sorted(map(len, cycles)) == [t,t,t+1,t+1], ('four cycles',n))
        chosen = next(cycle for cycle in cycles if len(cycle) % 2 == 0 and len(cycle)>=4)
        u, v = chosen[0], chosen[len(chosen)//2]
        start = (u, v)
        f = {(0,1),(3,4),(0,4),(1,n-1)}
        for _ in range(len(chosen)//2):
            need(base.pair(n-1-u,n-1-v) not in f, ('reflected antipodal defect',n))
            u,v = g[u],g[v]
            need(base.pair(u,v) not in f, ('antipodal defect',n))
        need((u,v) == tuple(reversed(start)), ('antipodal endpoint reversal',n))


def hostile(n, interior):
    qs = words(n, min(interior), False)
    src, report = base.symbolic(n, interior, qs, [(3,1)], False)
    need(not report['illegal'] and report['target'], 'hostile ledger must be legal and exact')
    need(report['blue_feasible'] == [False,True,False], 'hostile blue step')
    flips = set(src.edges) ^ base.pathset(qs[0])
    src.anti(qs[1], flips)
    need(src.ok, 'critical source exists')
    arcs = src.model(0)
    intermediate = {edge if base.pair(*edge) not in flips else tuple(reversed(edge)) for edge in arcs}
    rev = dict(zip(qs[1], reversed(qs[1])))
    need(all((rev[b], rev[a]) in intermediate for a,b in intermediate), 'literal blue hostile')
    return {'n':n,'interior_ports':sorted(interior),
            'rows':[sum(1<<b for a,b in arcs if a==v) for v in range(n)]}


def main():
    endpoint = general = inherited = h4_forward = 0
    for n in range(12,41):
        for h in range(5,n-3):
            interior = {h,n-3,n-2}
            for forward in (False,True):
                inherited += one_cell(n, interior, forward)
                endpoint += 1
    for n in range(12,18):
        fixed = {n-3,n-2}
        free = sorted(set(range(5,n-1))-fixed)
        for mask in range(1<<len(free)):
            interior = fixed | {v for i,v in enumerate(free) if (mask>>i)&1}
            k = len(interior)+1
            if k<4 or (k%2==0 and n in {2*k+1,2*k+3}):
                continue
            for forward in (False,True):
                inherited += one_cell(n, interior, forward)
                general += 1
        fixed4 = {4,n-3,n-2}
        free4 = sorted(set(range(5,n-1))-fixed4)
        for mask in range(1<<len(free4)):
            interior = fixed4 | {v for i,v in enumerate(free4) if (mask>>i)&1}
            qs = words(n,4,True)
            src,_ = base.symbolic(n,interior,qs,[(1,3)],True)
            flips = set(src.edges)^base.pathset(qs[0])
            need(src.forced(3,4,1^(base.pair(3,4) in flips)), 'h4 retained path witness')
            need(src.forced(n-3,n-1,0^(base.pair(n-3,n-1) in flips)), 'h4 reflected witness')
            for seed in (0,2261):
                base.literal(src.model(seed),qs)
            h4_forward += 1
    for n in range(12,301):
        geometry(n)
    critical = [hostile(13,{5,7,9,10,11}),hostile(15,{6,8,10,12,13})]
    qs = words(12,4,False)
    bad, report = base.symbolic(12,{4,9,10},qs,[(3,1)],False)
    need(bad.forced(0,4,0) and [0,0,4] in report['illegal'], 'h4 false first join')
    h3, report3 = base.symbolic(12,{3,9,10},words(12,3,False),[(3,1)],False)
    need(not report3['illegal'] and report3['target'] and
         report3['blue_feasible']==[False,False,True], 'h3 possible blue last word')
    print(json.dumps({'status':'PASS','endpoint_four_guard_cells':endpoint,
        'endpoint_four_order_range':[12,40],'general_guard_cells':general,
        'general_order_range':[12,17],'inherited_tail_routes':inherited,
        'h4_forward_affine_cells':h4_forward,
        'permutation_degree_antipodal_orders':[12,300],
        'literal_models_per_new_cell':2,'critical_blue_sources':critical,
        'h4_hostile':'n12,I={4,9,10},3->1 forces a->4 false',
        'h3_hostile':'n12,I={3,9,10},3->1 permits blue last word',
        'scope':'endpoint-four h>=5 and h4 with1->3; exact general degree exclusions'},indent=2))


if __name__ == '__main__':
    main()
