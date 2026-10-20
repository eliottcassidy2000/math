"""Exact center-deficit controls; finite census is not an all-order proof."""
import csv
import gzip
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def load(name, path):
    spec = spec_from_file_location(name, ROOT / path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load('center_blue', '05-knowledge/results/defect_transitions_20261020/blue_check.py')
C = load('center_cycles', '05-knowledge/results/black_descent_20261020c/cycles.py')
E = load('center_enum', '05-knowledge/results/metagraph_growth_20261020/enumeration.py')


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def invariants(adj):
    paths = B.hamilton_paths(adj)
    aut = E.canonical(adj, aut=True)[1]
    blue = sum(B.blue(adj, p) for p in paths)
    need(len(paths) % aut == blue % aut == 0, 'free path action')
    return len(paths) // aut, blue // aut


def main():
    known = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}
    word = [0] * 30
    word[0] = 1
    for j in range(1, len(word)):
        word[j] = word[j-1] + (word[j-3] if j >= 3 else 0)
    for n in range(3, 10):
        if n < 9:
            path = ROOT / f'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv'
            stream = path.open(encoding='utf-8')
        else:
            stream = gzip.open(HERE.parent / 'blue_quotient_20261020b/n9_classes.tsv.gz', 'rt')
        count = sc_count = strong_sc = 0
        low = Counter()
        deficit_two = []
        strong_two = []
        with stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                count += 1
                code = int(row['canon' if n < 9 else 'code'])
                sc = int(row['SC']) if n < 9 else code == int(row['converse_code'])
                if not sc:
                    continue
                sc_count += 1
                h = int(row['tc']) if n < 9 else int(row['SC_H']) // int(row['Aut'])
                b = int(row['bm' if n < 9 else 'blue_mass'])
                adj = tuple(C.rows(code, n))
                strong = C.strongly_connected(adj)
                strong_sc += strong
                if strong and h - b <= 30:
                    low[h-b] += 1
                if h-b == 2:
                    need(invariants(adj) == (h, b), 'independent deficit replay')
                    deficit_two.append(code)
                    if strong:
                        strong_two.append(code)
        need(count == known[n], 'complete class count')
        expected = word[(n-4)//2] if n >= 4 and n % 2 == 0 else 0
        need(len(deficit_two) == expected, 'conditional center-series finite control')
        need(len(strong_two) == int(n == 4), 'unique center finite control')
        print(f'n={n} classes={count} self_converse={sc_count} strong_self_converse={strong_sc} deficit2={len(deficit_two)} strong_deficit2={len(strong_two)} low_center_deficits={dict(sorted(low.items()))}')

    r4 = (2, 4, 9, 3)
    need(B.strong(r4) and invariants(r4) == (5, 3), 'R4 control')
    tests = 0
    for mass in range(7):
        for w in B.words(mass):
            pieces = [(0,) if a == 1 else B.cyclic(3) for a in w]
            left = B.ordered_join(pieces)
            for center in ((), (0,), B.cyclic(3), r4, B.cyclic(5)):
                full = B.ordered_join([left, center, E.converse(left)])
                if not full:
                    continue
                hc, bc = invariants(center) if center else (1, 1)
                need(invariants(full) == (hc, bc), 'rigid padding preserves center deficit')
                tests += 1
    nonrigid = B.cyclic(5)
    pair = B.ordered_join([nonrigid, E.converse(nonrigid)])
    need(invariants(pair) == (9, 3), 'nonrigid half is not harmless padding')
    print(f'rigid_padding_controls={tests}; R5->R5op has (h,b)=(9,3), deficit=6')
    print('PASS; local controls reach center order nine; all-order D(z)=z^4 is inherited separately.')


if __name__ == '__main__':
    main()
