#!/usr/bin/env python3
"""Exact native prefix-bank control; not a Mal'tsev CSP implementation."""
import argparse
import json
from itertools import permutations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


G = list(product(range(3), repeat=4))
ID = (0, 0, 0, 0)


def mul(g, h):
    return ((g[0]+h[0]) % 3,) + tuple((h[a+1]+g[(a+h[0]) % 3+1]) % 3 for a in range(3))


def inv(g):
    return (-g[0] % 3,) + tuple(-g[(a-g[0]) % 3+1] % 3 for a in range(3))


def perm(g):
    return tuple(3*((a+g[0]) % 3)+(x+g[a+1]) % 3 for a, x in product(range(3), repeat=2))


def arc(u, v):
    a, x = divmod(u, 3)
    b, y = divmod(v, 3)
    return (y-x) % 3 == 1 if a == b else (b-a) % 3 == 1


def run():
    maps = {g: perm(g) for g in G}
    found = {p for p in permutations(range(9))
             if all(arc(u, v) == arc(p[u], p[v]) for u in range(9) for v in range(u+1, 9))}
    need(found == set(maps.values()), 'complete internal automorphism group')
    need(all(tuple(maps[g][maps[h][u]] for u in range(9)) == maps[mul(g, h)]
             for g, h in product(G, repeat=2)), 'permutation multiplication')
    need(all(mul(g, inv(g)) == ID and mul(inv(g), g) == ID for g in G), 'group inverses')
    need(any(mul(g, h) != mul(h, g) for g, h in product(G, repeat=2)), 'nonabelian positive control')

    relations = []
    for c in (0, 1):
        # Source block v_outer-u_outer=0, target block =c.
        rel = {(g, h) for g, h in product(G, repeat=2)
               if all(((v//3-u//3) % 3 == 0) ==
                      ((maps[h][v]//3-maps[g][u]//3) % 3 == c)
                      for u in range(9) for v in range(9))}
        need(rel == {(g, h) for g, h in product(G, repeat=2) if (h[0]-g[0]) % 3 == c},
             'actual block quotient relation')
        need(len(rel) == 2187, 'one edge relation size')
        relations.append(rel)

    banks = []
    for i in range(3):
        bank = {}
        for g in G:
            if i and g[0] != 0:
                continue
            lift = [(g[0], 0, 0, 0) if i == 0 else ID for _ in range(3)]
            lift[i] = g
            bank[g] = tuple(lift)
        banks.append(bank)
    need([len(b) for b in banks] == [81, 27, 27], 'prefix image sizes')
    count = twisted = 0
    decompositions = set()
    for state in product(G, repeat=3):
        flat = all((state[i], state[j]) in relations[0] for i, j in ((0, 1), (1, 2), (0, 2)))
        twisted += ((state[0], state[1]) in relations[0] and
                    (state[1], state[2]) in relations[0] and
                    (state[0], state[2]) in relations[1])
        if not flat:
            continue
        count += 1
        remainder = state
        rebuilt = (ID,)*3
        word = []
        for i, bank in enumerate(banks):
            key = remainder[i]
            need(key in bank, 'complete prefix image')
            word.append(key)
            lift = bank[key]
            remainder = tuple(mul(inv(a), b) for a, b in zip(lift, remainder))
            rebuilt = tuple(mul(a, b) for a, b in zip(rebuilt, lift))
            need(all(g == ID for g in remainder[:i+1]), 'prefix retained')
        need(remainder == (ID,)*3 and rebuilt == state, 'exact factorization')
        decompositions.add(tuple(word))
    need(count == 59049 and len(decompositions) == count, 'unique complete bank')
    need(twisted == 0, 'local feasibility does not imply joint compatibility')

    # Restore an arbitrary nontrivial target frame: every relation is a coset.
    r = (G[17], G[32], G[68])
    shifts = tuple((r[j][0]-r[i][0]) % 3 for i, j in ((0, 1), (1, 2), (0, 2)))
    need(shifts[2] == (shifts[0]+shifts[1]) % 3, 'target phase cocycle')
    for word in decompositions:
        h = (ID,)*3
        for i, key in enumerate(word):
            h = tuple(mul(a, b) for a, b in zip(h, banks[i][key]))
        target = tuple(mul(a, b) for a, b in zip(r, h))
        need(all((target[j][0]-target[i][0]) % 3 == c
                 for (i, j), c in zip(((0, 1), (1, 2), (0, 2)), shifts)), 'source-to-target coset')

    return dict(status='FINITE-EXACT native control; cited Mal\'tsev algorithm not implemented',
                cell='C3[C3,C3,C3]', full_internal_group_order=len(found),
                candidate_tuples=len(G)**3, flat_solutions=count, twisted_solutions=twisted,
                local_relation_sizes=[len(r) for r in relations],
                prefix_image_sizes=[len(b) for b in banks],
                stored_nonidentity_tuples=sum(len(b)-1 for b in banks),
                unique_factorizations=len(decompositions), checks=CHECKS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = run()
    path = HERE/'nonabelian.json'
    if args.write:
        path.write_text(json.dumps(output, indent=2)+'\n')
    else:
        need(json.loads(path.read_text()) == output, 'nonabelian receipt mismatch')
    print(json.dumps(output, indent=2))
