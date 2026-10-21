#!/usr/bin/env python3
"""For each port pattern of the (n-5,4) layer that needed the labelled search (n = 9, 10), test whether a single word
triple (one of the found certificates, read as position words) is a valid three-black detour for EVERY presentation of
that pattern; report the best template per pattern and its coverage."""
import sys, json
from collections import defaultdict
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK)); sys.path.insert(0, '.')
import probe as P
from endpoint_four_bank import endpoint_k
for n in (9, 10):
    data = json.load(open(f'endpoint_four_n{n}.json'))
    certs = data['search_certificates']
    bypat = defaultdict(list)
    for c in certs: bypat[(c['kind'], tuple(c['incoming']))].append(tuple(tuple(q) for q in c['words']))
    pres = defaultdict(list)
    for t, u, v in endpoint_k(n, 4):
        pres[(('open' if P.arc(t, 1, n - 2) else 'closed'), u)].append(t)
    for key, words_list in sorted(bypat.items()):
        universe = pres[key]; best = (0, None)
        for w in sorted(set(words_list)):
            cov = sum(1 for t in universe if P.evaluate(t, w)['success'])
            if cov > best[0]: best = (cov, w)
        print(f"n = {n} pattern {key}: {len(universe)} presentations, {len(set(words_list))} distinct certificate templates; best template covers {best[0]}/{len(universe)}: {best[1]}", flush=True)
