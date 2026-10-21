#!/usr/bin/env python3
"""Consecutive-port patterns of the (n-5,4) layer at n = 9, 10: split each pattern into the half covered by the inherited bank
and the residual half; find (greedily) the fewest certificate templates covering the residual; report the bank template that
covers the bank half, the separating arcs between the two halves (orientation constant on each half and opposite), and whether
one bank template + one certificate template suffice -- the shape of a two-ledger structural statement."""
import sys, json, itertools
from collections import defaultdict, Counter
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK)); sys.path.insert(0, '.')
import probe as P, probe_refined as R
from endpoint_four_bank import endpoint_k
for n in (9, 10):
    data = json.load(open(f'endpoint_four_n{n}.json'))
    bypat = defaultdict(set)
    for c in data['search_certificates']: bypat[(c['kind'], tuple(c['incoming']))].add(tuple(tuple(q) for q in c['words']))
    bank = list(P.bank(n)) + list(R.new_bank(n))
    pres = defaultdict(list)
    for t, u, v in endpoint_k(n, 4): pres[('open' if P.arc(t, 1, n - 2) else 'closed', u)].append(t)
    for key, temps in sorted(bypat.items()):
        if len(temps) < 2 and len(pres[key]) < 256: continue
        universe = pres[key]
        bankhit = {}
        for i, t in enumerate(universe):
            for name, w in bank:
                if P.evaluate(t, w)['success']: bankhit[i] = (name, w); break
        covered = set(bankhit); residual = [i for i in range(len(universe)) if i not in covered]
        bank_names = Counter(v[0] for v in bankhit.values())
        # greedy cover of the residual by certificate templates
        temps = sorted(temps); cov = {w: {i for i in residual if P.evaluate(universe[i], w)['success']} for w in temps}
        left = set(residual); chosen = []
        while left:
            w = max(temps, key=lambda w: len(cov[w] & left))
            if not (cov[w] & left): break
            chosen.append((w, len(cov[w] & left))); left -= cov[w]
        # separating arcs between bank half and residual half
        A = [universe[i] for i in covered]; B = [universe[i] for i in residual]; sep = []
        for i in range(n):
            for j in range(i + 1, n):
                a_ij = {P.arc(t, i, j) for t in A}; b_ij = {P.arc(t, i, j) for t in B}
                if len(a_ij) == 1 and len(b_ij) == 1 and a_ij != b_ij: sep.append((i, j, 'bank half has i->j' if True in a_ij else 'bank half has j->i'))
        print(f"n = {n} {key}: {len(universe)} presentations; bank covers {len(covered)} via {dict(bank_names)}; residual {len(residual)}; greedy certificate cover: {[(c[1]) for c in chosen]} templates {len(chosen)}, uncovered {len(left)}; separating arcs: {sep}", flush=True)
        for w, k in chosen[:3]: print(f"      template covering {k}: {w}", flush=True)
