"""Independent Burnside fixed-point audit; does not import parent solver."""
from itertools import combinations, product
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent


def regular(p, d):
    rows = [tuple(int(j in s) for j in range(p)) for s in combinations(range(p), d)]
    return [B for B in product(rows, repeat=p)
            if all(sum(row[j] for row in B) == d for j in range(p))]


def fixed_counts(p, d):
    blocks = regular(p, d)
    return {(a, b): sum(all(B[i][j] == B[(i+a) % p][(j+b) % p]
                                for i in range(p) for j in range(p)) for B in blocks)
            for a, b in product(range(p), repeat=2)}


def count(p, n, edges, degrees, local):
    burnside_sum = 0
    nonidentity_pairs = 0
    for shifts in product(range(p), repeat=n):
        fixed = 1
        for (i, j), d in zip(edges, degrees):
            fixed *= local[p, d][shifts[i], shifts[j]]
        burnside_sum += fixed
        if any(shifts):
            nonidentity_pairs += fixed
    if burnside_sum % (p**n):
        raise ValueError('nonintegral independent Burnside sum')
    symmetric_labelled = nonidentity_pairs//(p-1)
    if symmetric_labelled*(p-1) != nonidentity_pairs:
        raise ValueError('nonintegral stabilizer multiplicity')
    return {'total': burnside_sum//p**n,
            'symmetric': symmetric_labelled//p**(n-1),
            'burnside_sum': burnside_sum,
            'nonidentity_fixed_incidence': nonidentity_pairs}


def main():
    local = {(3, 1): fixed_counts(3, 1), (5, 1): fixed_counts(5, 1),
             (5, 2): fixed_counts(5, 2)}
    source = json.loads((HERE.parent/'networks.json').read_text())
    results = []
    for row in source['complete_network_censuses']:
        observed = count(row['p'], row['vertices'], row['edges'], [1]*len(row['edges']), local)
        if (observed['total'], observed['symmetric']) != (row['total'], row['symmetric']):
            raise ValueError('independent fixed-point sum differs from parent orbit enumeration')
        results.append({'p': row['p'], 'vertices': row['vertices'], 'edges': row['edges'], **observed})
    mixed = count(5, 3, [(0, 1), (1, 2), (0, 2)], [1, 2, 1], local)
    target = source['mixed_degree_triangle_formula']
    if (mixed['total'], mixed['symmetric']) != (target['total'], target['symmetric']):
        raise ValueError('mixed-degree formula failed independent sum')
    out = {'status': 'FINITE-EXACT independent audit',
           'method': 'enumerate actual fixed binary matrices for each local translation, then sum over every shared vertex translation',
           'networks': results, 'mixed_degree_triangle': mixed,
           'local_fixed_counts': {str(k): {str(ab): v for ab, v in D.items()} for k, D in local.items()}}
    text = json.dumps(out, indent=2)+'\n'
    dest = HERE/'audit_networks.json'
    if '--write' in sys.argv:
        dest.write_text(text, encoding='utf-8')
    elif json.loads(dest.read_text()) != out:
        raise ValueError('independent audit receipt changed')
    print(json.dumps({'total_orbits': [r['total'] for r in results], 'mixed': mixed}, sort_keys=True))


if __name__ == '__main__':
    main()
