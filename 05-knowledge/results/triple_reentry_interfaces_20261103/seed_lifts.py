#!/usr/bin/env python3
"""Rooted first-ear slack on the six scalar seeds, via subset path counts."""
import argparse
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    seeds = json.loads((HERE / 'scalar_seeds.json').read_text())
    out = []
    for h, seed in seeds.items():
        rows = seed['rows']; n = len(rows)
        dp = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                todo = rows[v] & ~mask
                while todo:
                    bit = todo & -todo; todo -= bit
                    dp[mask | bit][bit.bit_length()-1] += dp[mask][v]
        hs = [sum(a) for a in dp]; hs[0] = 1
        full = (1 << n) - 1
        if hs[-1] != int(h):
            raise RuntimeError('seed raw count')
        K = [sum(dp[mask][q] * hs[full ^ mask] for mask in range(1, 1 << n)
                 if mask >> q & 1) for q in range(n)]
        J = [x - int(h) for x in K]
        if min(J) < 0:
            raise RuntimeError('negative rooted cover slack')
        out.append({'H': int(h), 'rows': rows, 'J': J, 'minimum_J': min(J),
                    'minimum_lift_H': 9*int(h)+6*min(J),
                    'best_ports': [q for q,j in enumerate(J) if j == min(J)]})
    text = json.dumps(out, indent=2, sort_keys=True)+'\n'
    if args.write:
        (HERE / 'seed_lifts.json').write_text(text)
    elif (HERE / 'seed_lifts.json').read_text() != text:
        raise RuntimeError('seed lift receipt mismatch')
    print('PASS: six rooted slack profiles and exact best-port lift counts')
if __name__ == '__main__':
    main()
