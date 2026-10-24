"""Independent Boolean group-algebra census of every receiver pin set.

No producer imports, count subtraction, or supplied short charge labels.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
FULL = (1 << 2048) - 1
LO = [sum(1 << q for q in range(2048) if not q & (1 << j)) for j in range(4)]
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def translate(bits, x, parity):
    for j in range(4):
        if parity & (1 << j):
            bits = ((bits & LO[j]) << (1 << j)) | ((bits & (FULL ^ LO[j])) >> (1 << j))
    shift = 16 * (x % 128)
    return ((bits << shift) | (bits >> (2048-shift))) & FULL


def grades(labels):
    d = [1] + [0] * len(labels)
    for i, label in enumerate(labels):
        for k in range(i+1, 0, -1):
            d[k] |= translate(d[k-1], *label)
    return d


def main():
    receipt = json.loads((HERE/'receipt.json').read_text())
    e = json.loads((RESULTS/'metagraph_order8_electrical_20261023/receipt.json').read_text())['order8']
    x, y = e['higher_layers'][-1]['lifts'][0], e['radical_lifts_mod4']
    ids = receipt['sources']
    labels = {s: (x[s]%128, sum((y[j][s]%2) << (j-1) for j in range(1,5))) for s in ids+[1841]}
    need([16*labels[s][0]+labels[s][1] for s in ids] == receipt['labels'],
         'original electrical potentials recover producer charge coordinates')
    # All 2048 one-point translations for two generators test bit order and wrap.
    for q in range(2048):
        for a,b in ((127,15),(1,1)):
            expected = ((((q>>4)+a)%128)<<4) | ((q&15)^b)
            need(translate(1<<q,a,b) == 1<<expected, 'literal point translation')
    for erased_count in range(4):
        profiles = Counter()
        center_bad = []
        for erased in combinations(ids, erased_count):
            remaining = [labels[s] for s in ids if s not in erased]
            d = grades(remaining)
            profiles[tuple(v.bit_count() for v in d)] += 1
            if erased_count == 3:
                repaired = grades(remaining+[labels[1841]])
                need(all(repaired[k] == FULL for k in (8,9,10)),
                     'independent Boolean DP: every triple repaired by actual1841')
                if d[8] != FULL or d[9] != FULL:
                    center_bad.append(list(erased))
        expected = {tuple(r['coverage']):r['pin_sets'] for r in receipt['levels'][erased_count]['profiles']}
        need(profiles == expected, 'all erased-set grade profiles agree independently')
    need(center_bad == [r['erased'] for r in receipt['adaptive_repair']['exceptional_original_centers']],
         'exact six exceptional central sections retained')
    print(f'PASS {CHECKS}; independent Boolean census, all1351 pin sets, all1140 repairs')


if __name__ == '__main__':
    main()
