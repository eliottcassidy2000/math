"""Exact path counts and all arc changes on iterated marked multiplicative constructions."""
import argparse
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, permutations, product
from math import comb, prod
from pathlib import Path
import sys

from rail_arcs import fillings, two_slot_kernel


ARMS = {
    1353: (8, 5, 73, 50, 39, 71, 27),
    3963: (6, 12, 72, 129, 15, 31, 59, 119),
    5551: (6, 20, 136, 67, 45, 79, 151, 59),
    4681: (10, 12, 65, 164, 15, 23, 59, 119),
    8275: (136, 37, 73, 82, 39, 141, 51, 94),
    2**61+1: 61,
}


def endpoint_counts(rows):
    n = len(rows)
    counts = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        counts[1 << i][i] = 1
    for mask, row in enumerate(counts):
        for i, ways in enumerate(row):
            if not ways:
                continue
            targets = rows[i] & ~mask
            while targets:
                bit = targets & -targets
                counts[mask | bit][bit.bit_length()-1] += ways
                targets ^= bit
    return counts


def subset_parts(mask):
    subset = mask
    while subset:
        yield subset
        subset = (subset-1) & mask


def tables(rows):
    """Number of ordered nonempty runs with their optional end/start pinned in distinct runs."""
    n = rows if isinstance(rows, int) else len(rows)
    indices = (None, *range(n))
    if isinstance(rows, int):
        return {(a,b): tuple(fillings(n,r,a,b) for r in (1,2,3)) for a in indices for b in indices}
    full = (1 << n)-1
    assert all(0 <= row <= full and not row & (1 << i) for i,row in enumerate(rows))
    assert all(((rows[i] >> j) ^ (rows[j] >> i)) & 1 for i,j in combinations(range(n),2))
    ends = endpoint_counts(rows)
    starts = endpoint_counts(tuple(full ^ (1<<i) ^ row for i,row in enumerate(rows)))
    h = list(map(sum,ends))

    @cache
    def count(mask, free, a, b):
        if free < 0:
            return 0
        if a is not None and b is not None:
            if a == b:
                return 0
            return sum(ends[s][a]*count(mask ^ s, free, None, b) for s in subset_parts(mask) if s >> a & 1 and not s >> b & 1)
        if a is not None or b is not None:
            data, vertex = (ends,a) if a is not None else (starts,b)
            return sum(data[s][vertex]*count(mask ^ s, free, None, None) for s in subset_parts(mask) if s >> vertex & 1)
        if free == 0:
            return int(mask == 0)
        if free == 1:
            return h[mask]
        return sum(h[s]*count(mask ^ s, free-1, None, None) for s in subset_parts(mask))

    return {(a,b): tuple(count(full,r-int(a is not None)-int(b is not None),a,b) for r in (1,2,3))
            for a in indices for b in indices}


def quotient(length):
    rows, ports, root = (0,), [], 0
    for _ in range(length):
        n = len(rows)
        rows = tuple(row | (int(i == root) << (n+1)) for i,row in enumerate(rows)) + (
            (1 << n)-1, ((1 << (n+1))-1) ^ (1 << root))
        ports.append(n)
        root = n+1
    return rows, tuple(ports)


def contract(kernel, ports, perm, table):
    """Integrate out every port other than the one or two that contain the joining endpoints."""
    port_num = {v:i for i,v in enumerate(ports)}
    out = {}
    for (i,j), terms in kernel.items():
        a,b = port_num.get(i),port_num.get(j)
        outside = [k for k in range(len(ports)) if k not in (a,b)]
        # Preaggregate at the two pinned slots. All labels of a run-multiplicity tuple are <4.
        c = defaultdict(int)
        for runs, ways in terms.items():
            for k in outside:
                ways *= table[perm[k]][None,None][runs[k]-1]
            x,y = (runs[a]-1 if a is not None else 0),(runs[b]-1 if b is not None else 0)
            c[x,y] += ways
        out[i,j] = c
    return out


def changes(rows, ports, perm, table, kernel, arms=ARMS):
    coeff = contract(kernel, ports, perm, table)
    factors = dict(zip(ports,perm))
    profiles = {i:table[t] for i,t in factors.items()}
    sizes = {i:((arms[t] if isinstance(arms[t],int) else len(arms[t]))) for i,t in factors.items()}
    names = [(i,a) for i in range(len(rows)) for a in range(sizes.get(i,1))]
    def joined(i,a,j,b):
        c=coeff.get((i,j))
        if not c:return 0
        if i==j:
            return sum(val*profiles[i][a,b][t] for (t,_),val in c.items()) if i in factors else 0
        left=profiles[i][a,None] if i in factors else (1,)
        right=profiles[j][None,b] if j in factors else (1,)
        return sum(val*left[t]*right[s] for (t,s),val in c.items())
    for (i,a),(j,b) in combinations(names,2):
        if i==j:
            mod=arms[factors[i]]
            forward=isinstance(mod,int) or bool(mod[a]>>b & 1)
        else:
            forward=bool(rows[i]>>j & 1)
        delta=joined(j,b,i,a)-joined(i,a,j,b)
        yield ((i,a),(j,b)), delta if forward else -delta


def class_count():
    """Count uniquely decodable words from the full small class catalogue."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_arc_fiber_grammar_20261122"))
    from fibers import catalog
    goal = 2**60 - 1
    bank = Counter()
    for n in range(2, 9):
        for rows in (catalog(n) if n >= 3 else ((2,0),)):
            paths = list(map(sum, endpoint_counts(rows)))
            paths[0] = 1
            factor = paths[-1] + sum(a*b for a,b in zip(paths, reversed(paths)))
            if goal % factor == 0:
                bank[n+1, factor] += 1

    @cache
    def words(target):
        if target == 1:
            return {(0,0): 1}
        answer = Counter()
        for (cost, factor), multiplicity in bank.items():
            if target % factor == 0:
                for (n,k), count in words(target // factor).items():
                    answer[n+cost,k+1] += multiplicity * count
        return answer

    totals = Counter()
    for (n,k), count in words(goal).items():
        # Choose a position for the one large arm. The initial seed has order one.
        totals[n+63] += (k+1)*count
    assert (len(bank), sum(bank.values())) == (66,243)
    assert totals == {107: 645120, 108: 737280}
    print('small arm (cost,factor) entries, types:', len(bank), sum(bank.values()))
    print('strong classes by order:', sorted(totals.items()))


def check_small():
    """Compare every arc in a nontransitive two-arm instance with literal subset DP."""
    arms = {15:(2,4,1), 5:(2,0)}
    perm = tuple(arms)
    rows, ports = quotient(len(perm))
    table = {f:tables(mod) for f,mod in arms.items()}
    kernel = two_slot_kernel(rows,ports)
    blocks = dict(zip(ports, arms.values()))
    names = [(i,a) for i in range(len(rows)) for a in range(len(blocks.get(i, (0,))))]
    graph = [sum((1 << k) for k,(j,b) in enumerate(names)
                 if (rows[i] >> j & 1) or (i == j and blocks.get(i,(0,))[a] >> b & 1)) for i,a in names]
    base = sum(endpoint_counts(graph)[-1])
    assert base == prod(arms)
    for (i,j), delta in changes(rows, ports, perm, table, kernel, arms):
        a,b = names.index(i),names.index(j)
        changed = graph.copy()
        changed[a] ^= 1 << b
        changed[b] ^= 1 << a
        assert sum(endpoint_counts(changed)[-1]) == base + delta


def dyadic_spines():
    """All normalized one-apex words with nontrivial count dividing 32767."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_e31_physical_20261122"))
    from check import p, two_slots

    divisors = Counter()
    hits = []
    for v in range(2, 1 << 16, 2):
        h = 1 + p(v)
        if h == 1 or 32767 % h:
            continue
        divisors[h] += 1
        q = two_slots(v)
        if q % h == 0 and ((q//h) & (q//h - 1)) == 0:
            hits.append((v,h,q))
    assert divisors == {31:1, 151:2, 217:2, 4681:2, 32767:4}
    assert not hits
    print('one-apex normalized words by height', sorted(divisors.items()), 'dyadic hits', hits)


def four_layer_check():
    """Universal signed run coefficients of the one exceptional four-layer crossing."""
    rows, ports = quotient(4)
    kernel = two_slot_kernel(rows, ports)
    delta = Counter(kernel[2,6])
    delta.subtract(kernel[6,2])
    expected = (60,22,33,11,33,12,18,6,22,8,12,4,11,4,6,2)
    assert {k:v for k,v in delta.items() if v} == dict(zip(product((1,2), repeat=4), expected))
    root_cases = {
        2: {(0,2): (-6,-2,-3,-1), (2,4): (-6,-3,-2,-1), (4,0): (4,3,3,2)},
        3: {(4,0): (16,6,9,3,11,4,6,2), (6,2): (16,11,9,6,6,4,3,2),
            (2,4): (-18,-6,-9,-3,-6,-2,-3,-1)}
    }
    for depth, cases in root_cases.items():
        rows, ports = quotient(depth)
        kernel = two_slot_kernel(rows, ports)
        for (a,b), coefficients in cases.items():
            delta = Counter(kernel[b,a])
            delta.subtract(kernel[a,b])
            assert {k:v for k,v in delta.items() if v} == dict(zip(product((1,2), repeat=depth), coefficients))
    rows, ports = quotient(3)
    updates = changes(rows, ports, (3,3,3), {3:tables(1)}, two_slot_kernel(rows, ports), {3:1})
    assert sum(d == 2 for _,d in updates) == 2


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-orders',action='store_true',help='use all 720 orders of the six distinct arms')
    parser.add_argument('--class-count', action='store_true', help='count words of all catalogued small arms')
    parser.add_argument('--dyadic-spines', action='store_true', help='exclude every normalized one-apex dyadic arm')
    args=parser.parse_args()
    if args.dyadic_spines:
        dyadic_spines()
        return
    if args.class_count:
        class_count()
        return
    check_small()
    four_layer_check()
    table={k:tables(mod) for k,mod in ARMS.items()}
    for factor,t in table.items():
        h,two,_=t[None,None]
        assert 3*h+two==factor
    lower=2**60*(2**61-1)-1
    assert prod(ARMS)==lower
    rows,ports=quotient(len(ARMS))
    kernel=two_slot_kernel(rows,ports)
    assert all(max(r)<=3 for k in kernel.values() for r in k)
    totals=Counter(); hits=[]; minimum=None; total=0
    for perm in (permutations(ARMS) if args.all_orders else (tuple(ARMS),)):
        ctr=Counter()
        for edge,delta in changes(rows,ports,perm,table,kernel):
            ctr[(delta>0)-(delta<0)]+=1
            if delta==2:
                hits.append((perm,edge))
            if delta>0 and (minimum is None or delta<minimum):
                minimum=delta
        assert sum(ctr.values())==comb(107,2)
        totals.update(ctr);total+=1
    print('strong lower isomorphism classes tested',total,'directed arc reversals',sum(totals.values()))
    print('signed totals',sorted(totals.items()),'least positive',minimum,'+2 hits',hits)


if __name__=='__main__':
    main()
