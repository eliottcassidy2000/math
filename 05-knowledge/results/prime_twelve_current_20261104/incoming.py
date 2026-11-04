"""Reconcile S187's complete modular cubic and gate clocks with retained integers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def cutoff(a):
    n = len(a)
    full = (1 << n)-1
    g = [[0]*n for _ in range(full+1)]
    starts = [0]*(full+1)
    for m in range(1, full+1):
        for u in bits(m):
            rest = m ^ (1 << u)
            ends = 1 << u if not rest else 0
            for v in bits(rest & a[u]):
                ends ^= g[rest][v]
            g[m][u] = ends
            starts[m] |= (ends.bit_count() % 2) << u
    # Store endpoint parities together. Keep the two starts ordered until
    # extracting a named coefficient; no endpoint-support marginal product.
    k = [[0]*n for _ in range(n)]
    partitions = 0
    for left in range(1, full+1):
        rem = full ^ left
        right = rem
        while right:
            size = (left | right).bit_count()
            if size % 2 and size <= n-3:
                partitions += 1
                for u in bits(starts[left]):
                    for v in bits(right):
                        k[u][v] ^= g[right][v]
            right = (right-1) & rem
    coeff = {}
    for u in range(n):
        for v in range(u+1, n):
            ends = k[u][v] ^ k[v][u]
            for w in range(n):
                pair, single = (1 << u) | (1 << v), 1 << w
                coeff[pair, single] = coeff[single, pair] = (ends >> w) & 1
    return coeff, starts, partitions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'incoming.json')
    args = parser.parse_args()
    inputs = ['response_tables.json', 'response.json',
              '../session187_named_response_20261104/prime12/cubic.json',
              '../homogeneous_pair_response_20261103/response_tables.json']
    pins = {name: sha256((HERE/name).read_bytes()).hexdigest() for name in inputs+['incoming.py']}
    tables = json.loads((HERE/inputs[0]).read_text())['cases']
    cards = {x['name']: x['cards'] for x in json.loads((HERE/inputs[1]).read_text())['cases']}
    cases = []
    for t in tables:
        masks = t['masks']
        c = {(x,y): t['contrast_coefficients'][i*len(masks)+j]
             for i,x in enumerate(masks) for j,y in enumerate(masks)}
        cases.append((t['name'], t['rows'], c))
    incoming = json.loads((HERE/inputs[2]).read_text())
    cases.append(('incoming-S187-P', incoming['rows'], {(x,y): c for x,y,c in incoming['coefficients']}))
    records = []
    shared_checks = 0
    for name, rows, exact in cases:
        modular, starts, partitions = cutoff(rows)
        need(len(modular) == 1584, 'complete named modular cubic')
        for key, value in modular.items():
            need(exact[key] % 2 == 0 and (exact[key]//2) % 2 == value, 'incoming cutoff versus exact integer cubic')
        if name in cards:
            for card in cards[name]:
                v = card['deleted']
                short = 0
                for support in range(1, 4096):
                    if not (support >> v) & 1 and support.bit_count() in (2,4,6,8):
                        short ^= starts[support]
                for u, current in zip(card['names'], card['one_port_current']):
                    need(current % 2 == 0 and (current//2) % 2 == (short >> u) & 1, 'short-even-support deletion current')
                    shared_checks += 1
        records.append(dict(name=name, cubic_coordinates=len(modular),
                            nonzero_modular_cubics=sum(modular.values()),
                            ordered_disjoint_support_pairs=partitions))
    gates = json.loads((HERE/inputs[3]).read_text())['cases']
    rank_checks = 0
    for t in gates:
        j, r, masks = t['state'], t['r'], t['masks']
        pos = {m:i for i,m in enumerate(masks)}
        h = 4*3**r+9*2**r-3+2*j
        b,c,d = [(44,-5,-4),(42,-1,-3),(42,3,-3)][j]
        for k in range(1,r+1):
            score = 18*3**(r-k)*4**(k-1)+b*2**(r-k)*3**(k-1)+c*2**(k-1)+d*(k==r)
            actual = t['values'][pos[1 << (k+4)]*len(masks)]
            need(actual == h+2*score, 'incoming rank clock versus retained actual query')
            rank_checks += 1
        key = pos[1]*len(masks)
        delta = t['values'][key]-t['converse_values'][key]
        expected = 2*[4*3**r+6*2**r-6,10*3**r+12*2**r-9,10*3**r-3*2**r-9][j]
        need(delta == expected > 0, 'incoming fixed exterior clock')
    need(pins == {name: sha256((HERE/name).read_bytes()).hexdigest() for name in pins}, 'inputs changed')
    result = dict(status='PASS scoped integration; no prime12 or spectrum closure', source_sha256=pins,
                  cases=records, named_modular_cubic_checks=1584*len(cases),
                  short_even_current_checks=shared_checks, gate_rank_queries=rank_checks,
                  gate_whole_converse_queries=len(gates))
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','cases')},indent=2))


if __name__ == '__main__':
    main()
