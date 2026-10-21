#!/usr/bin/env python3
"""For each consecutive-port pattern of the (n-5,4) layer that the inherited bank misses (n = 9, 10), search for a PAIR of word
templates (from the found certificates, read as position words) that together cover the whole pattern; report the pair, the
split sizes, and the arc that separates the two halves (an arc (i,j) whose orientation agrees with the template choice on
every presentation of the pattern) -- the candidate guard for a two-ledger structural theorem."""
import sys, json, itertools
from collections import defaultdict
from pathlib import Path
PK = Path('/Users/e/Documents/GitHub/math-wt-mersenne-561/05-knowledge/results/endpoint_three_repairs_20261020h')
sys.path.insert(0, str(PK)); sys.path.insert(0, '.')
import probe as P
from endpoint_four_bank import endpoint_k
for n in (9, 10):
    data = json.load(open(f'endpoint_four_n{n}.json'))
    bypat = defaultdict(set)
    for c in data['search_certificates']: bypat[(c['kind'], tuple(c['incoming']))].add(tuple(tuple(q) for q in c['words']))
    pres = defaultdict(list)
    for t, u, v in endpoint_k(n, 4): pres[('open' if P.arc(t, 1, n - 2) else 'closed', u)].append(t)
    for key, temps in sorted(bypat.items()):
        if len(temps) < 2: continue
        universe = pres[key]; temps = sorted(temps)
        cover = {w: frozenset(i for i, t in enumerate(universe) if P.evaluate(t, w)['success']) for w in temps}
        full = frozenset(range(len(universe))); best = None
        for w1, w2 in itertools.combinations(temps, 2):
            if cover[w1] | cover[w2] == full:
                if best is None or (len(cover[w1] & cover[w2]) < best[0]): best = (len(cover[w1] & cover[w2]), w1, w2)
        if best is None:
            print(f"n = {n} pattern {key}: {len(universe)} presentations, {len(temps)} templates, NO pair covers; best single {max(len(c) for c in cover.values())}", flush=True); continue
        ov, w1, w2 = best; A = cover[w1]; B = cover[w2]
        # separating arcs: pairs (i,j) with arc i->j in every presentation of A\B and j->i in every presentation of B\A (or vice versa)
        sep = []
        onlyA = [universe[i] for i in A - B]; onlyB = [universe[i] for i in B - A]
        for i in range(n):
            for j in range(i + 1, n):
                a_ij = {P.arc(t, i, j) for t in onlyA}; b_ij = {P.arc(t, i, j) for t in onlyB}
                if len(a_ij) == 1 and len(b_ij) == 1 and a_ij != b_ij: sep.append((i, j, 'A: i->j' if True in a_ij else 'A: j->i'))
        print(f"n = {n} pattern {key}: {len(universe)} presentations, {len(temps)} templates; covering pair with overlap {ov}: |A| = {len(A)}, |B| = {len(B)}\n   w1 = {w1}\n   w2 = {w2}\n   separating arcs (A\\B vs B\\A): {sep}", flush=True)
