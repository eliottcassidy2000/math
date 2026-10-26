"""Exact integration checks; does not rerun the incoming nauty census."""
import csv
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[4]


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def mul(a, b, n):
    return [sum(a[j] * b[i-j] for j in range(i+1)) for i in range(n+1)]


def main():
    incoming = ROOT / '05-knowledge/results/a000568_intransitivity_families_20261026'
    rows = list(csv.DictReader((incoming / 'data/prime_cycle_types_n9.tsv').read_text().splitlines(), delimiter='\t'))
    a, s = [0] * 11, [0] * 11
    a[0] = a[1] = s[1] = 1
    prime = {}
    for row in rows:
        k = int(row['k'])
        prime[k] = prime.get(k, 0) + int(row['number_of_prime_classes'])
    for n in range(2, 11):
        total = Q(7906161 if n == 10 else 0)  # Incoming n10 prime census input.
        for row in rows:
            if int(row['k']) > n:
                continue
            weighted, mass = Q(0), 0
            for term in row['Aut_cycle_types(count)'].split():
                cycles, copies = term.rsplit('x', 1)
                copies = int(copies)
                mass += copies
                polynomial = [1] + [0] * n
                for factor in cycles.split(','):
                    length, count = map(int, factor.split('^'))
                    base = [a[i // length] if i and i % length == 0 else 0 for i in range(n+1)]
                    for _ in range(count):
                        polynomial = mul(polynomial, base, n)
                weighted += copies * polynomial[n]
            need(mass == int(row['|Aut|']), 'cycle-index mass')
            total += Q(int(row['number_of_prime_classes']), mass) * weighted
        need(total.denominator == 1, 'integral strong count')
        s[n] = int(total)
        a[n] = sum(s[j] * a[n-j] for j in range(1, n+1))
    need(a == [1,1,1,2,4,12,56,456,6880,191536,9733056], 'prime-species reconstruction')
    own = json.loads((ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/order9-normal.json').read_text())
    need(all(a[r['order']] == r['classes'] for r in own['census']), 'incoming prime families versus saved clone families')
    need(prime[9] == 136655 and own['census'][-1]['cores'] == 147095, 'prime and pair-free are distinct')

    z = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    pairs = [(0,3),(1,5),(2,4)]
    reflect = [0,5,4,3,2,1]
    norm = lambda u,v: (u[0]-v[0])**2+(u[0]-v[0])*(u[1]-v[1])+(u[1]-v[1])**2
    need(all(norm(z[i],z[j]) == norm(z[reflect[i]],z[reflect[j]]) for i,j in combinations(range(6),2)), 'whole metric reflection')
    need([norm(z[i],z[j]) for i,j in pairs] == [4,3,3], 'non-antipodal pairing gate')
    for x in range(8):
        chosen = [pairs[j][x>>j & 1] for j in range(3)]
        target = [pairs[j][(x ^ 6)>>j & 1] for j in range(3)]
        need([reflect[v] for v in chosen] == target and x != x ^ 6, 'fixed-free even translation')
    support = total = 0
    for beta in product(range(6), repeat=6):
        weight = (6-len(set(beta))) if all(beta[i] != beta[(i+1)%6] for i in range(6)) else 0
        support += weight > 0
        total += weight
    need((support,total) == (14910,24600), 'full wheel boundary remains colourable')

    edges = list(combinations(range(4),2))
    def adj(t,i,j):
        if i>j: return not adj(t,j,i)
        return bool(t>>edges.index((i,j)) & 1)
    h = [sum(all(adj(t,u,v) for u,v in zip(p,p[1:])) for p in permutations(range(4))) for t in range(64)]
    t = sum(1<<edges.index(e) for e in [(0,1),(1,2),(0,3),(1,3),(2,3)])
    op = t ^ 63
    def c3(t):
        return tuple(tr for tr in combinations(range(4),3) if all(sum(adj(t,v,w) for w in tr if w!=v)==1 for v in tr))
    def scores(t):
        return sorted(sum(adj(t,v,w) for w in range(4) if w!=v) for v in range(4))
    need(c3(t) == c3(op) and all(h[t^f] == h[op^f] for f in range(64)), 'unsigned cycles plus response still loses converse')
    need(scores(t) != scores(op), 'even simultaneous unsigned records fail up to isomorphism')
    print(json.dumps({'status':'PASS','scope':'independent aggregation of incoming census inputs; no nauty replay',
                      'prime_cycle_index_rows':len(rows),'A000568_through10':a,'strong_through10':s,
                      'order9_prime_cores':prime[9],'order9_pair_free_cores':147095,
                      'mirror_pair_squared_lengths':[4,3,3],'mirror_gain':6,
                      'wheel_boundary_support':support,'wheel_six_colorings':total,
                      'unsigned_C3_response_hostile':{'masks':[t,op],'scores':[scores(t),scores(op)]}},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
