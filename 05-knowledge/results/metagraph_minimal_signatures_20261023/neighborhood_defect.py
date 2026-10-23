"""Rooted five-vertex overlap defect for doubly regular neighborhoods.

Independent of the frozen adaptive-signature API. Default checks receipt.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def paley(p):
    squares = {x*x % p for x in range(1, p)}
    return [[int((j-i) % p in squares) for j in range(p)] for i in range(p)]


def double(a):
    n = len(a); out = [[0]*(2*n+1) for _ in range(2*n+1)]
    for i in range(n):
        out[n][i] = out[n+1+i][n] = 1
        for j in range(n):
            out[i][j] = a[i][j]
            out[n+1+i][n+1+j] = a[j][i]
            out[i][n+1+j] = int(i == j or a[i][j])
            out[n+1+j][i] = 1-out[i][n+1+j]
    return out


def examine(a):
    n = len(a); need(n % 8 == 7, 'order must be 8h+7'); h = (n-7)//8
    need(all(a[i][i] == 0 for i in range(n)), 'zero diagonal')
    need(all(a[i][j]+a[j][i] == 1 for i,j in combinations(range(n),2)), 'tournament')
    need(all(sum(row) == (n-1)//2 for row in a), 'regular source')
    need(all(sum(a[i][k]*a[j][k] for k in range(n)) == (n-3)//4
             for i,j in combinations(range(n),2)), 'doubly regular source')
    roots = []
    for v in range(n):
        vertices = [u for u in range(n) if a[v][u]]
        need(len(vertices) == 4*h+3, 'neighborhood size')
        need(all(sum(a[u][w] for w in vertices) == 2*h+1 for u in vertices), 'regular neighborhood')
        counts = [sum(a[u][w]*a[z][w] for w in vertices) for u,z in combinations(vertices,2)]
        pair_count = len(vertices)*(len(vertices)-1)//2
        need(sum(counts) == h*pair_count, 'mean common-out count')
        defect = sum((k-h)**2 for k in counts)
        blocks = sum(k*(k-1)//2 for k in counts)
        # Independent literal four-set partition count inside the neighborhood.
        literal = 0; transitive_fours = 0
        for four in combinations(vertices,4):
            transitive_fours += sorted(sum(a[u][w] for w in four) for u in four) == [0,1,2,3]
            for sources in combinations(four,2):
                targets = [w for w in four if w not in sources]
                literal += all(a[u][w] for u in sources for w in targets)
        need(literal == blocks, 'rooted five-vertex overlap count')
        need(transitive_fours == blocks, 'rooted transitive-five count')
        need(defect == 2*blocks-h*(h-1)*pair_count, 'variance/overlap identity')
        need((defect == 0) == all(k == h for k in counts), 'zero iff doubly regular neighborhood')
        roots.append({'root':v, 'J':defect, 'B':blocks, 'common_out_histogram':sorted(Counter(counts).items())})
    total_transitive_five = sum(sorted(sum(a[u][w] for w in five) for u in five) == [0,1,2,3,4]
                                for five in combinations(range(n),5))
    need(total_transitive_five == sum(r['B'] for r in roots), 'global transitive-five count')
    four_counts = Counter(tuple(sorted(sum(a[u][w] for w in four) for u in four))
                          for four in combinations(range(n),4))
    t = (n-3)//4
    source_cycle = n*(comb(2*t+1,3)-(2*t+1)*comb(t,2))
    strong_four = comb(n,2)*comb(t+1,2)
    predicted = Counter({(0,1,2,3):comb(n,4)-strong_four-2*source_cycle,
                         (1,1,1,3):source_cycle, (0,2,2,2):source_cycle, (1,1,2,2):strong_four})
    need(four_counts == predicted, 'all four induced-motif classes')
    return {'n':n, 'h':h, 'J_histogram':sorted(Counter(r['J'] for r in roots).items()),
            'B_histogram':sorted(Counter(r['B'] for r in roots).items()),
            'transitive_five_total':total_transitive_five,
            'induced_four_motifs_by_score_sequence':sorted(four_counts.items()),
            'zero_roots':sum(r['J'] == 0 for r in roots), 'roots':roots}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    cases = {name:examine(a) for name,a in [('Paley7_control',paley(7)),('Paley31',paley(31)),
                                          ('double_double_Paley7',double(double(paley(7))))]}
    need(cases['Paley7_control']['J_histogram'] == [(0,7)], 'positive zero control')
    need(cases['Paley31']['zero_roots'] == 0 and cases['double_double_Paley7']['zero_roots'] == 9,
         'authenticated pair zero-root counts')
    out = {'status':'PROVED variance/overlap identity for DRTs of order 8h+7; FINITE-EXACT cases',
           'cases':cases, 'checks':CHECKS,
           'source_lf_sha256':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
    raw = (json.dumps(out,sort_keys=True,indent=2)+'\n').encode(); path = HERE/'neighborhood_defect.json'
    if args.write:
        path.write_bytes(raw)
    else:
        need(path.read_bytes().replace(b'\r\n',b'\n') == raw, 'frozen receipt')
    print(json.dumps({'checks':CHECKS, 'receipt_lf_sha256':sha256(raw).hexdigest(),
                      'J_histograms':{name:case['J_histogram'] for name,case in cases.items()}}))


if __name__ == '__main__':
    main()
