"""Complete small-base search for the requested 99*2^k+23 tail ray."""
from pathlib import Path
import sys
import json
from itertools import combinations

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'modular_centers_collatz_pell_20261111'
sys.path.insert(0, str(PRIOR))
import content_checks as cc
import tournament_checks as tc


def oddpart(x):
    if not x:
        return 0, 0
    j = 0
    while x % 2 == 0:
        x //= 2
        j += 1
    return x, j


def main():
    bank, counts, hits = [[]], [], []
    for n in range(1, 8):
        codes = {cc.canonical(cc.child(rows, mask)) for rows in bank for mask in range(1 << (n - 1))}
        bank = [cc.decode(n, code) for code in sorted(codes)]
        counts.append(len(bank))
        for code, rows in zip(sorted(codes), bank):
            h, u, w, _, _, _ = tc.path_tables(rows)
            if h < 23:
                continue
            for selected in range(1, 1 << n):
                vertices = [v for v in range(n) if selected >> v & 1]
                offset = h - sum(w[a][b] for a, b in combinations(vertices, 2))
                if offset != 23:
                    continue
                slope = sum(u[v] for v in vertices)
                c, j = oddpart(slope)
                if c == 99:
                    hits.append(dict(n=n, code=code, rows=rows, selected=selected,
                                     H=h, offset=offset, U=slope, first_exponent=j+1))
        if hits:
            break
    result = {'ordinary_base_counts': counts, 'hits': hits}
    (HERE / 'ray_probe.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
