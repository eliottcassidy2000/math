"""Native regular-host block moves with pinned actual owners.

All guard controls here are explicitly defined integer-weight filters, not
surrogates for the repository's electrical mod-8 predicate. Standard library.
"""
from array import array
from collections import Counter
import csv
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

CHECKS = 0
HERE = Path(__file__).resolve().parent


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def phases(word, r, q, alpha=0):
    a = [alpha]
    for i in range(r):
        a.append((a[-1] + (word >> i & 1)) % q)
    return tuple(a)


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


def host(q, t):
    need(q >= 3 and q % 2 == 1 and t >= 3 and t % 2 == 1, 'regular-host parameters')
    def arc(x, y):
        i, u = divmod(x, q)
        j, v = divmod(y, q)
        if i == j:
            return (v-u) % q in range(1, (q+1)//2)
        if (j-i) % t in range(1, (t+1)//2):
            return (v-u) % q in (1, 2)
        return (u-v) % q not in (1, 2)
    n = q*t
    need(all(arc(x,y) != arc(y,x) for x,y in combinations(range(n),2)), 'actual tournament')
    need(all(sum(arc(x,y) for y in range(n) if x != y) == (n-1)//2 for x in range(n)),
         'actual regular scores')
    return arc


def actual_path(a, q):
    return tuple(i*q + (a[i]+u) % q for i in range(len(a)) for u in range(q))


def verify_path(a, q, arc):
    p = actual_path(a, q)
    need(len(set(p)) == len(p), 'literal owners occur once')
    need(all(arc(x,y) for x,y in zip(p,p[1:])), 'every actual Hamiltonian junction')
    return p


def connected(ids, bound, aa):
    todo = [ids[0]]
    seen = set(todo)
    while todo:
        x = todo.pop()
        for y in ids:
            if y not in seen and distance(aa[x], aa[y]) <= bound:
                seen.add(y)
                todo.append(y)
    return len(seen) == len(ids)


def minimum_budget(weights, r, q):
    if len(weights) >= 2:
        return max((q-1)*(b-a)//q for a,b in zip(weights, weights[1:]))
    k = weights[0]
    return int(0 < k < r)


def one_family(q, r, special=False):
    t = r+1
    arc = host(q,t)
    words = [w for w in range(1 << r) if not special or w.bit_count() % q == 0]
    aa = {w: phases(w,r,q) for w in words}
    actual = {}
    for w in words:
        p = verify_path(aa[w],q,arc)
        actual[w] = p
        # Global phase is retained as an actual lift, not silently quotiented.
        for alpha in range(1,q):
            pa = verify_path(phases(w,r,q,alpha),q,arc)
            need(all(pa[j]//q == p[j]//q and (pa[j]-p[j]) % q == alpha for j in range(len(p))),
                 'actual common-phase automorphism lift')
        # Each primitive is checked on all eligible positions.
        for i in range(r-1):
            if (w >> i & 1) != (w >> (i+1) & 1):
                z = w ^ (3 << i)
                need(distance(aa[w], phases(z,r,q)) == 1, 'one-cell adjacent swap')
                need(z.bit_count() == w.bit_count(), 'swap preserves integer weight')
        for i in range(r-q+1):
            bits = (w >> i) & ((1 << q)-1)
            if bits in (0, (1 << q)-1):
                z = w ^ (((1 << q)-1) << i)
                need(distance(aa[w], phases(z,r,q)) == q-1, 'q-run block support')
                need(phases(z,r,q)[-1] == aa[w][-1], 'block keeps actual endpoint phase')

        # Deterministic legal sorting certificate inside the integer-weight layer.
        z = w
        sorting_steps = 0
        while z != (1 << w.bit_count())-1:
            i = next(i for i in range(r-1) if (z >> i & 3) == 2)
            z ^= 3 << i
            verify_path(phases(z,r,q),q,arc)
            sorting_steps += 1
        need(sorting_steps <= w.bit_count()*(r-w.bit_count()), 'sorting bound')

    filters = 0
    endpoint_reports = []
    for delta in range(q):
        ids = [w for w in words if w.bit_count() % q == delta]
        if not ids:
            continue
        weights = sorted({w.bit_count() for w in ids})
        for k,l in combinations(weights,2):
            lower = min(distance(aa[x], aa[y]) for x in ids if x.bit_count() == k
                        for y in ids if y.bit_count() == l)
            need(lower == (q-1)*(l-k)//q, 'sharp distance between weight layers')
        for subset in range(1,1 << len(weights)):
            accepted = [k for i,k in enumerate(weights) if subset >> i & 1]
            kept = [w for w in ids if w.bit_count() in accepted]
            budget = minimum_budget(accepted,r,q)
            need(connected(kept,budget,aa), 'guarded budget sufficient')
            if budget:
                need(not connected(kept,budget-1,aa), 'guarded budget necessary')
            filters += 1
        endpoint_reports.append({'delta': delta, 'weights': weights, 'states': len(ids),
                                 'unfiltered_minimum_budget': minimum_budget(weights,r,q)})
    return {'q': q, 'r': r, 'regular_host_order': q*t, 'words': len(words),
            'weight_filters_checked': filters, 'endpoint_fibres': endpoint_reports}


def two_interval_hostile():
    q,r = 3,8
    arc = host(q,r+1)
    pins = (0,3,5,8)
    accepted = []
    for w in range(1 << r):
        a = phases(w,r,q)
        if any(a[i] != 0 for i in pins):
            continue
        weights = tuple(sum(w >> j & 1 for j in range(l,h)) for l,h in zip(pins,pins[1:]))
        if weights[0] == weights[2]:
            verify_path(a,q,arc)
            accepted.append((w,a,weights))
    need(len(accepted) == 2, 'correlated weight guard keeps exactly two states')
    need(distance(accepted[0][1], accepted[1][1]) == 4, 'two intervals require coupled support four')
    p0,p1 = [actual_path(x[1],q) for x in accepted]
    for i in pins:
        # One indexed owner in a C3 fixes its entire phase.
        need(p0[i*q] == p1[i*q], 'actual shared indexed owner pin retained')
    need(all(not arc(y,x) for x,y in zip(p0,p0[1:])), 'converse reversal control')
    reversed0 = tuple(reversed(p0))
    need(all(arc(y,x) for x,y in zip(reversed0,reversed0[1:])), 'literal reversed path in converse host')
    def backward_word(path):
        return sum(1 << k for k,(i,j) in enumerate(combinations(range(len(path)),2))
                   if not arc(path[i],path[j]))
    need(backward_word(p0) != backward_word(p1), 'native backward word genuinely changes')
    return {'regular_host_order': 27, 'pinned_cell_indices': pins,
            'guard': 'the integer weights in the first and last pinned intervals agree; not an electrical guard',
            'accepted_difference_words': [x[0] for x in accepted],
            'accepted_weight_vectors': [x[2] for x in accepted], 'sharp_coupled_phase_budget': 4}


def original_electrical_pins():
    """Independent realization and lookup in the ORIGINAL order-eight atlas."""
    results = HERE.parents[1]
    frozen = json.loads((HERE.parent/'guards'/'receipt.json').read_text(encoding='utf-8'))['original_eight']
    folder = results/'metagraph_order8_20261023'
    potential_bytes = (results/'metagraph_order8_electrical_20261023'/'receipt.json').read_bytes()
    x = json.loads(potential_bytes)['order8']['higher_layers'][-1]['lifts'][0]
    need(sha256(potential_bytes).hexdigest() == frozen['original_potential_sha256'],
         'original frozen electrical potential bytes')
    nodes = list(csv.DictReader((folder/'merged_nodes.tsv').open(), delimiter='\t'))
    merged = [None]*6880
    for row in nodes:
        for key in ('class_a','class_b'):
            merged[int(row[key])] = int(row['merged_id'])
    raw = gzip.decompress((folder/'owner_map.u16.gz').read_bytes())
    need(sha256(raw).hexdigest() == frozen['owner_bytes_sha256'], 'original native owner map bytes')
    owners = array('H')
    owners.frombytes(raw)
    if sys.byteorder != 'little':
        owners.byteswap()
    need(len(owners) == 1 << 21 and len(x) == 3528, 'original atlas dimensions')
    def arc(u,v):
        if u == v:
            return False
        if u == 0 or v == 7:
            return True
        if v == 0 or u == 7:
            return False
        i,a = divmod(u-1,3)
        j,b = divmod(v-1,3)
        if i == j:
            return (b-a) % 3 == 1
        return (b-a) % 3 in (1,2) if i < j else (a-b) % 3 not in (1,2)
    need(all(arc(u,v) != arc(v,u) for u,v in combinations(range(8),2)), 'padded actual tournament')
    chords = [(i,j) for i in range(8) for j in range(i+2,8)]
    full = (1 << 21)-1
    def lookup(path, converse=False):
        # This atlas stores FORWARD bits; the n<=6 atlas uses another gauge.
        mask = sum((arc(path[j],path[i]) if converse else arc(path[i],path[j])) << bit
                   for bit,(i,j) in enumerate(chords))
        source = merged[owners[mask]]
        target = merged[owners[full ^ mask]]
        return mask,source,target,(x[source]-x[target]) % 8,int(owners[mask])
    aa, rows = {}, []
    for a in range(3):
        for b in range(3):
            path = (0,)+tuple(1+(a+u) % 3 for u in range(3))+tuple(4+(b+u) % 3 for u in range(3))+(7,)
            valid = all(arc(u,v) for u,v in zip(path,path[1:]))
            need(valid == ((b-a) % 3 in (0,1)), 'all nine actual phase cases')
            if not valid:
                continue
            mask,source,target,value,ordinary = lookup(path)
            need(path[0] == 0 and path[7] == 7, 'actual indexed endpoint pins')
            need(source == 3517, 'actual original merged source')
            delta = (b-a) % 3
            need((mask,target,value) == ((1765183,532,1) if delta == 0 else (1799615,1290,5)),
                 'actual original electrical phase sectors')
            rev = tuple(reversed(path))
            need(all(arc(v,u) for u,v in zip(rev,rev[1:])), 'actual converse Hamiltonian path')
            need(lookup(rev,True)[1:4] == (source,target,value), 'actual converse transports original phase')
            aa[(a,b)] = (a,b)
            rows.append({'a':a,'b':b,'path':path,'ordinary_source':ordinary,
                         'normalized_forward_mask':mask,'source_merged':source,'target_merged':target,
                         'phase_mod8':value})
    need(connected(list(aa),1,aa), 'unfiltered native two-cell move graph connected')
    need(len({row['ordinary_source'] for row in rows}) == 1, 'fixed exact ordinary source')
    for value in (1,5):
        ids = [(row['a'],row['b']) for row in rows if row['phase_mod8'] == value]
        need(len(ids) == 3 and not connected(ids,1,aa) and connected(ids,2,aa),
             'original guard has sharp simultaneous budget two')
    for delta in (1,2):
        perm = [0]+[1+((u-1)//3)*3+((u-1)%3+delta) % 3 for u in range(1,7)]+[7]
        need(all(arc(u,v) == arc(perm[u],perm[v]) for u,v in combinations(range(8),2)),
             'shared common-phase move is actual fixed-port automorphism')
        need(perm[0] == 0 and perm[7] == 7, 'shared actual endpoint owners fixed')
    return {'n':8,'regular_source':False,'actual_pins':[[0,0],[7,7]],
            'original_potential_sha256':frozen['original_potential_sha256'],
            'owner_bytes_sha256':frozen['owner_bytes_sha256'],
            'phase_values':[1,5],'sharp_guarded_phase_budget':2,'rows':rows,
            'scope':'original frozen native electrical character; every common-phase kernel lift retained'}


if __name__ == '__main__':
    reports = [one_family(q,r) for q in (3,5,7) for r in (2,4,6,8)]
    reports.append(one_family(9,10,special=True))
    hostile = two_interval_hostile()
    electrical = original_electrical_pins()
    result = {'status': 'FINITE-EXACT', 'checks': CHECKS, 'families': reports,
              'correlated_pin_guard_hostile': hostile,
              'original_electrical_pinned_two_cells':electrical,
              'scope': 'actual regular cyclic-cell families and explicitly chosen weight filters; separate actual pinned order-eight original electrical character control'}
    receipt = Path(__file__).with_name('receipt.json')
    if '--write' in sys.argv:
        receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    else:
        need(json.loads(json.dumps(result)) == json.loads(receipt.read_text(encoding='utf-8')), 'saved receipt')
    print(json.dumps({'status':'PASS','checks':CHECKS,'families':len(reports),
                      'weight_filters':sum(x['weight_filters_checked'] for x in reports),
                      'correlated_pin_guard_hostile':hostile},indent=2))
