"""Bounded exact scout: length-four excursions and depth-seven mutation balls."""
from collections import Counter, deque
from itertools import permutations
from pathlib import Path
import json
import sys
import verify as v


def representatives(n):
    seen = set()
    reps = []
    for mask in range(1 << (n*(n-1)//2)):
        if mask in seen:
            continue
        B = v.matrix(n, mask)
        reps.append(mask)
        seen.update(v.maskof(tuple(tuple(B[i][j] for j in p) for i in p))
                    for p in permutations(range(n)))
    return reps


def four_excursions():
    census = []
    witnesses = []
    for n in (3, 4, 5):
        counts = Counter()
        for mask in range(1 << (n*(n-1)//2)):
            B = v.matrix(n, mask)
            h = v.hamilton(B)
            for k in range(n):
                C = v.mutate(B, k)
                if v.tournament(C):
                    continue
                for l in range(n):
                    if l == k:
                        continue
                    D = v.mutate(C, l)
                    for m in range(n):
                        if m == l:
                            continue
                        E = v.mutate(D, m)
                        if v.tournament(E):
                            continue
                        for q in range(n):
                            if q == m:
                                continue
                            F = v.mutate(E, q)
                            if not v.tournament(F):
                                continue
                            word = (k, l, m, q)
                            names = tuple(dict.fromkeys(word))
                            shape = ''.join(str(names.index(x)) for x in word)
                            dh = v.hamilton(F)-h
                            counts[(shape, dh)] += 1
                            if dh:
                                witnesses.append(dict(n=n, mask=mask, word=word,
                                                      old_H=h, new_H=h+dh,
                                                      initial=B, final=F))
        expected = {3: {'0120': 24}, 4: {'0101': 96, '0120': 96},
                    5: {'0101': 1920, '0120': 480}}[n]
        v.need(dict(counts) == {(shape, 0): count for shape, count in expected.items()},
               'length-four exact census')
        census.append(dict(n=n, excursion_shapes=expected, H_changes=0))
    return census, witnesses


def balls(depth=7):
    census = []
    witnesses = []
    for n in (3, 4, 5):
        reps = representatives(n)
        v.need(len(reps) == {3: 2, 4: 4, 5: 12}[n], 'all ordinary seed types')
        for seed in reps:
            B = v.matrix(n, seed)
            h = v.hamilton(B)
            todo = deque([(B, ())])
            seen = {B}
            layer = Counter({0: 1})
            returned = set()
            max_weight = 1
            while todo:
                C, word = todo.popleft()
                if len(word) == depth:
                    continue
                for k in range(n):
                    if word and k == word[-1]:
                        continue
                    D = v.mutate(C, k)
                    if D in seen:
                        continue
                    seen.add(D)
                    w = word+(k,)
                    layer[len(w)] += 1
                    max_weight = max(max_weight, max(abs(x) for row in D for x in row))
                    if v.tournament(D):
                        returned.add(v.maskof(D))
                        dh = v.hamilton(D)-h
                        if dh:
                            witnesses.append(dict(n=n, seed=seed, old_H=h, new_H=h+dh,
                                                  word=w, initial=B, final=D))
                    todo.append((D, w))
            census.append(dict(n=n, seed_mask=seed, H=h, depth=depth,
                               distinct_matrices=len(seen), layers=dict(sorted(layer.items())),
                               returned_tournaments=len(returned), max_abs_entry=max_weight))
    return census, witnesses


def main():
    four, fw = four_excursions()
    bfs, bw = balls()
    result = dict(status='FINITE-EXACT; bounded search, not full mutation-class closure',
                  four_step=four, depth_seven=bfs, H_changing_witnesses=fw+bw,
                  weight_cutoff=None, mutation_words_cutoff=7)
    name = sys.argv[1] if len(sys.argv) > 1 else 'scout_receipt.json'
    output = (v.HERE/name).resolve()
    v.need(output.parent == v.HERE, 'owned output directory')
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', seeds=len(bfs),
                          matrices=sum(r['distinct_matrices'] for r in bfs),
                          H_changing_witnesses=len(fw+bw), four_step=four), sort_keys=True))


if __name__ == '__main__':
    main()
