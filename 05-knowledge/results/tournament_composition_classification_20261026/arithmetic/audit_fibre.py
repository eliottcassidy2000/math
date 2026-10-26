"""Independent arithmetic-agent audit of the four order-twelve fibre classes."""
import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

from check_arithmetic import automorphisms, cycle, need, strong, substitute, transitive


def cyclic_triples(A):
    return [list(vertices) for vertices in combinations(range(len(A)), 3)
            if all(sum(A[v][w] for w in vertices if w != v) == 1 for v in vertices)]


def hamilton_paths(A):
    n = len(A)
    counts = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        counts[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for j in range(n):
            if mask & (1 << j):
                smaller = mask ^ (1 << j)
                if smaller:
                    counts[mask][j] = sum(counts[smaller][i] for i in range(n)
                                          if smaller & (1 << i) and A[i][j])
    return sum(counts[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    c3, one = cycle(3), transitive(1)
    A = substitute(transitive(2), [c3, one])
    # Converse on exactly the same labels preserves the triangle support.
    B = [[A[j][i] for j in range(4)] for i in range(4)]
    report = {'universe': 'four stated representatives on the same twelve labels', 'representatives': {}}
    common = None
    for b_count, name in enumerate(('AAA', 'AAB', 'ABB', 'BBB')):
        T = substitute(c3, [A]*(3-b_count)+[B]*b_count)
        triples = cyclic_triples(T)
        if common is None:
            common = triples
        need(triples == common, 'literal equality of the labelled cyclic triples')
        group = automorphisms(T)
        expected = 81 if b_count in (0, 3) else 27
        need(len(group) == expected, 'full automorphism order from independent arc code')
        hp = hamilton_paths(T)
        need(hp == 462915 and hp % len(group) == 0, 'Hamilton count and free quotient')
        need(strong(T), 'all representatives are strong')
        report['representatives'][name] = {'automorphism_order': len(group),
                                         'hamilton_paths': hp, 'path_orbits': hp//len(group),
                                         'strong': True}
    report['common_triples'] = common
    need(len(common) == 67, 'triple count')
    codegrees = Counter()
    for triple in common:
        codegrees.update(combinations(triple, 2))
    report['pair_codegree_histogram'] = dict(sorted(Counter(codegrees[pair]
                                                    for pair in combinations(range(12), 2)).items()))
    need(report['pair_codegree_histogram'] == {0: 9, 1: 9, 4: 48}, 'intrinsic block recovery')
    proof = Path(__file__).resolve().parent.parent/'fibre'/'RESULT.md'
    report['audited_proof_sha256'] = hashlib.sha256(proof.read_bytes()).hexdigest()
    with args.out.open('w', encoding='utf-8', newline='\n') as output:
        output.write(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report['representatives'], sort_keys=True))


if __name__ == '__main__':
    main()
