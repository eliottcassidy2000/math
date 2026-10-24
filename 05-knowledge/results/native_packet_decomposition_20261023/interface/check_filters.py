"""Exact endpoint-filter controls, including full native charge tables."""
from itertools import combinations, product
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def closure(n, edges):
    reach = [(1 << i) for i in range(n)]
    for i, j in edges:
        reach[i] |= 1 << j
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return reach


def upper(mask, reach):
    return all(not (mask >> i & 1) or not (r & ~mask)
               for i, r in enumerate(reach))


def zero_charge(mask, charges, moduli):
    for j, modulus in enumerate(moduli):
        v = sum(q[j] for i, q in enumerate(charges) if mask >> i & 1)
        if (v % modulus if modulus else v) != 0:
            return False
    return True


def bottleneck(accepted, n):
    """Exact monotone subset recurrence; slot cost, not arbitrary-route cost."""
    full = (1 << n)-1
    distance = {0: 0}
    for v in sorted(accepted, key=lambda a: (a.bit_count(), a)):
        if v == 0:
            continue
        distance[v] = min(max(distance[u], (v ^ u).bit_count())
                          for u in distance if u & ~v == 0)
    return distance[full]


def exhaustive_four():
    pairs = list(combinations(range(4), 2))
    models = 0
    histogram = {}
    for types in product(range(4), repeat=6):
        edges = []
        for (i, j), kind in zip(pairs, types):
            if kind & 1:
                edges.append((i, j))  # forbid 10
            if kind & 2:
                edges.append((j, i))  # forbid 01
        reach = closure(4, edges)
        direct = [mask for mask in range(16) if all(
            not (kind & 1 and (mask >> i & 1) and not (mask >> j & 1))
            and not (kind & 2 and not (mask >> i & 1) and (mask >> j & 1))
            for (i, j), kind in zip(pairs, types))]
        need(direct == [m for m in range(16) if upper(m, reach)], 'literal/closure')
        for q in range(16):
            if q.bit_count() % 2:
                continue
            accepted = [m for m in direct if (q & m).bit_count() % 2 == 0]
            charges = [(q >> i & 1,) for i in range(4)]
            need(accepted == [m for m in range(16) if upper(m, reach)
                              and zero_charge(m, charges, [2])], 'joint charge')
            # Independently ask the reachability consequence for each support limit.
            for bound in range(1, 5):
                seen = {0}
                for v in sorted(accepted, key=lambda a: (a.bit_count(), a)):
                    if any(u & ~v == 0 and (v ^ u).bit_count() <= bound for u in seen):
                        seen.add(v)
                if 15 in seen:
                    break
            need(bound == bottleneck(accepted, 4), 'exact monotone bottleneck')
            histogram[str(bound)] = histogram.get(str(bound), 0)+1
            models += 1
    return {'pair_systems': 4096, 'charged_models': models,
            'monotone_bottleneck_histogram': histogram}


def compatible(a, b):
    return not (a['ports'] & b['ports']) and all(
        b['pins'].get(k, v) == v for k, v in a['pins'].items())


def full_record_controls():
    rng = random.Random(20261023)
    count = 500
    for trial in range(count):
        n = rng.randrange(2, 9)
        records = [[{'ports': {('own', i, b)}, 'pins': {},
                     'kernel': (i, b, 'retained raw lift')} for b in range(2)]
                   for i in range(n)]
        for i, j in combinations(range(n), 2):
            mode = rng.randrange(5)
            if mode in (1, 2, 3):
                for a, b in ((1, 0), (0, 1)):
                    if mode == 3 or (mode == 1 and a == 1) or (mode == 2 and a == 0):
                        key = ('shared normalized resource', i, j, a)
                        records[i][a]['ports'].add(key)
                        records[j][b]['ports'].add(key)
            elif mode == 4:
                for b in range(2):
                    records[i][b]['pins'][(i, j)] = b
                    records[j][b]['pins'][(i, j)] = b
        for endpoint in (0, 1):
            need(all(compatible(records[i][endpoint], records[j][endpoint])
                     for i, j in combinations(range(n), 2)), 'feasible endpoints')
        edges = []
        for i, j in combinations(range(n), 2):
            if not compatible(records[i][1], records[j][0]):
                edges.append((i, j))
            if not compatible(records[i][0], records[j][1]):
                edges.append((j, i))
        charges = [(rng.randrange(-4, 5), rng.randrange(4), rng.randrange(2))
                   for _ in range(n-1)]
        charges.append(tuple(-sum(q[j] for q in charges) for j in range(3)))
        reach = closure(n, edges)
        direct, compiled = [], []
        for mask in range(1 << n):
            charge_ok = zero_charge(mask, charges, [0, 4, 2])
            if charge_ok and all(compatible(records[i][mask >> i & 1],
                                            records[j][mask >> j & 1])
                                 for i, j in combinations(range(n), 2)):
                direct.append(mask)
            if charge_ok and upper(mask, reach):
                compiled.append(mask)
        need(direct == compiled, 'actual resources/shared owners/integer charges')
        need(bottleneck(direct, n) <= n, 'finite monotone decomposition')
    # Two raw lifts of one normalized key must still collide.
    a = {'ports': {('normalized', 344062)}, 'pins': {}, 'kernel': 'lift A'}
    b = {'ports': {('normalized', 344062)}, 'pins': {}, 'kernel': 'different lift B'}
    need(not compatible(a, b), 'kernel decoration is not normalized capacity')
    return {'seed': 20261023, 'full_record_models': count,
            'charge_group': 'Z x Z/4 x Z/2', 'normalized_kernel_hostile': True}


def native_controls():
    old = RESULTS/'injective_native_blocks_20261023'
    doc = json.loads((old/'search/selected_blocks.json').read_text(encoding='utf-8'))
    electrical = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text(encoding='utf-8'))['order8']
    x = electrical['higher_layers'][-1]['lifts'][0]
    ys = electrical['radical_lifts_mod4'][1:]
    summaries = []
    larger = json.loads((old/'search/larger16_selected.json').read_text(encoding='utf-8'))['circuits'][0]
    for slots, shared, expected in ((doc['slots'], doc['shared_owner_equalities'], [0, 17, 238, 255]),
                                    (larger['slots'], [], [0, 65535])):
        n = len(slots)
        edges = []
        for (i, ci), (j, cj) in shared:
            for ai, bj in ((1, 0), (0, 1)):
                ki, kj = ('before', 'after')[ai], ('before', 'after')[bj]
                if slots[i][ki]['path'][ci] != slots[j][kj]['path'][cj]:
                    edges.append((i, j) if ai else (j, i))
        charges = [(s['sign'], s['sign']*x[s['source']],
                    *(y[s['source']] % 2 for y in ys)) for s in slots]
        reach = closure(n, edges)
        accepted = [m for m in range(1 << n) if upper(m, reach)
                    and zero_charge(m, charges, [0, 128, 2, 2, 2, 2])]
        need(accepted == expected, 'inherited native endpoint interval')
        exact = [m for m in accepted if zero_charge(m, charges, [0, 0, 2, 2, 2, 2])]
        need(exact == expected, 'integer X retained in native interval')
        summaries.append({'slots': n, 'implications': edges, 'accepted_masks': accepted,
                          'exact_integer_accepted': exact,
                          'sharp_monotone_slot_budget': bottleneck(accepted, n)})
    return summaries


def main():
    out = {'status': 'FINITE-EXACT controls of proved endpoint-filter theorem',
           'exhaustive': exhaustive_four(), 'records': full_record_controls(),
           'native': native_controls()}
    out['checks'] = CHECKS
    frozen = json.loads(json.dumps(out))
    path = HERE/'receipt.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(frozen, indent=2)+'\n', encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == frozen, 'frozen receipt')
    print(json.dumps(frozen, sort_keys=True))


if __name__ == '__main__':
    main()
