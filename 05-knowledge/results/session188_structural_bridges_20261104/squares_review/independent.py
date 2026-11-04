#!/usr/bin/env python3
from math import isqrt
from itertools import combinations
import json


def need(x, msg):
    if not x:
        raise RuntimeError(msg)


def main():
    oo = ee = cross = 0
    rows = []
    # Add edges by their later endpoint; oracle tests square values directly.
    for n in range(1, 10001):
        born = []
        k = 1
        while k*k < 2*n:
            other = k*k-n
            if 1 <= other < n:
                born.append((other, n))
            k += 1
        need(len(born) == isqrt(2*n-1)-isqrt(n), ('birth', n))
        for a, b in born:
            if a % 2 and b % 2:
                oo += 1
            elif a % 2 == 0 and b % 2 == 0:
                ee += 1
            else:
                cross += 1
        formula = isqrt((n-1)//2) if n % 2 else isqrt((n+1)//4)
        need(oo-ee == formula, ('inventory', n))
        if n in (1, 2, 15, 16, 17, 25, 26, 100, 10000):
            rows.append(dict(n=n, odd_odd=oo, even_even=ee, cross=cross, births=len(born)))
    # Endpoint identities tested against all possible parity arrangements.
    paths = cycles = 0
    for n in range(2, 19):
        odds = (n+1)//2
        for positions in combinations(range(n), odds):
            labels = [0]*n
            for p in positions:
                labels[p] = 1
            po = sum(a == b == 1 for a, b in zip(labels, labels[1:]))
            pe = sum(a == b == 0 for a, b in zip(labels, labels[1:]))
            eo = labels[0]+labels[-1]
            ee = 2-eo
            need(2*(po-pe) == 2*(2*odds-n)-(eo-ee), ('path incidence', n, labels))
            need(po-pe in ((0, 1, 2) if n % 2 else (-1, 0, 1)), ('path interval', n))
            paths += 1
            if n >= 3:
                po += labels[-1] == labels[0] == 1
                pe += labels[-1] == labels[0] == 0
                need(po-pe == 2*odds-n, ('cycle incidence', n, labels))
                cycles += 1
    # Explicit genuine square-sum Hamilton paths, with odd/even graph orders.
    base = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    witnesses = [base, base+[16], [17]+base+[16]]
    witness_rows = []
    for path in witnesses:
        n = len(path)
        need(sorted(path) == list(range(1,n+1)), 'not spanning')
        need(all(isqrt(a+b)**2 == a+b for a,b in zip(path,path[1:])), 'non-square edge')
        all_edges = {(i,j) for j in range(2,n+1) for i in range(1,j) if isqrt(i+j)**2 == i+j}
        used = {tuple(sorted((a,b))) for a,b in zip(path,path[1:])}
        omitted = all_edges-used
        missed_oo = sum(i%2 == j%2 == 1 for i,j in omitted)
        missed_ee = sum(i%2 == j%2 == 0 for i,j in omitted)
        inventory = sum(i%2 == j%2 == 1 for i,j in all_edges)-sum(i%2 == j%2 == 0 for i,j in all_edges)
        eo = path[0]%2+path[-1]%2
        need(missed_oo-missed_ee == inventory-(n%2)+(eo-1), 'endpoint-specific omitted inventory')
        need(missed_oo >= max(0, inventory-(2 if n%2 else 1)), 'omission bound')
        witness_rows.append(dict(n=n, path=path, omitted_odd_odd=missed_oo, omitted_even_even=missed_ee))
    print(json.dumps(dict(status='PASS', inventory_and_births_through=10000,
                          parity_path_cases=paths, parity_cycle_cases=cycles,
                          samples=rows, actual_hamilton_paths=witness_rows), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
