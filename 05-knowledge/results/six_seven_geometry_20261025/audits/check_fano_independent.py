"""Independent finite Fano audit; no producer imports or label-permutation generation."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def main():
    triples = list(combinations(range(7), 3))
    systems = []
    def extend(candidates, chosen):
        systems.append(tuple(chosen))
        for pos, index in enumerate(candidates):
            t = set(triples[index])
            next_candidates = [j for j in candidates[pos+1:] if len(t.intersection(triples[j])) < 2]
            extend(next_candidates, chosen+[index])
    extend(list(range(35)), [])
    fano = [set(triples[j] for j in indices) for indices in systems if len(indices) == 7]
    need(len(fano) == 30, 'independent exact-cover Fano generation')
    need(len(systems) == 5596, 'partial Steiner system count')
    min_free = min(sum(all(triples[j] not in f for j in indices) for f in fano) for indices in systems)
    need(min_free == 8, 'universal disjoint Fano bound')
    star = {(0, 1, 2), (0, 3, 4), (0, 5, 6)}
    star_hist = Counter(len(star.intersection(f)) for f in fano)
    need(star_hist == {0: 16, 1: 12, 3: 2}, 'hourglass star distribution')
    allwords = set(range(128))
    reconstructed = []
    for f in fano:
        # A Fano incidence structure supplies its 7 complement words, without
        # choosing vector labels, an incidence kernel, or a simplex generator.
        code = {0} | {127 ^ sum(1 << j for j in t) for t in f}
        need(len(code) == 8 and all(a ^ b in code for a in code for b in code), 'simplex closure')
        even_canonical = [127 ^ (1 << j) for j in range(7)]
        need(all(any((a & char).bit_count() % 2 for a in code) for char in even_canonical), 'no invariant canonical form')
        for zeros in [set()] + [set(t) for size in (1, 2) for t in combinations(range(7), size)] + [set(t) for t in star]:
            nonzero = set(range(7))-zeros
            stabilizer = [a for a in code if len({(a >> j) & 1 for j in nonzero}) == 1]
            need(len(stabilizer) == (2 if tuple(sorted(zeros)) in f else 1), 'fixed locus criterion')
        reconstructed.append(tuple(sorted(code)))
    producer = json.loads((HERE.parent/'fano_codes/receipt.json').read_text(encoding='utf-8'))
    need(set(reconstructed) == {tuple(row['simplex_words']) for row in producer['all_hourglass_markings']}, 'producer subgroup set')
    need(producer['universal_partial_Steiner_systems']['minimum_free_subgroups'] == min_free, 'producer universal lower bound')
    receipt = {'status': 'PASS', 'method': 'all pair-compatible triple systems; Fano structures selected as exact7 covers',
               'partial_system_count': len(systems), 'Fano_structures': len(fano), 'minimum_disjoint_Fano': min_free,
               'hourglass_bad_triple_histogram': dict(sorted(star_hist.items())),
               'all30_group_sets_match_producer': True, 'canonical_invariants': 0,
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'producer_source_sha256': hashlib.sha256((HERE.parent/'fano_codes/check.py').read_bytes()).hexdigest()}
    out = json.dumps(receipt, indent=2, sort_keys=True)+'\n'
    (HERE/'fano_independent_receipt.json').write_text(out, encoding='utf-8')
    print(out, end='')


if __name__ == '__main__':
    main()
