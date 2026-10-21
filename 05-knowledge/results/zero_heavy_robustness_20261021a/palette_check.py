#!/usr/bin/env python3
"""Exact physical ports, root gluing, colour witnesses and repair-radius controls."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from math import comb
import hashlib
import json
import os


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def repo_root():
    candidates = [Path(os.environ['MATH_REPO'])] if os.environ.get('MATH_REPO') else []
    for p in Path(__file__).resolve().parents:
        candidates.extend((p, p/'math'))
    for p in candidates:
        if (p/'05-knowledge/results/first_exit_complement_games_20261020/palette-glue.json').is_file():
            return p
    raise RuntimeError('Set MATH_REPO to the pinned repository checkout')


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def neg(a):
    return tuple(-x for x in a)


def mul(a, b):
    out = [F(0)]*4
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i ^ j] += x*y*(-3 if i & j & 1 else 1)*(-11 if i & j & 2 else 1)
    return tuple(out)


def norm(a):
    return mul(a, (a[0], -a[1], -a[2], a[3]))


def main():
    root = repo_root()
    manifest = json.loads((Path(__file__).parent/'manifest.json').read_text())
    for p, expected in manifest['dependencies'].items():
        need(hashlib.sha256((root/p).read_bytes()).hexdigest() == expected, 'pinned dependency '+p)
    src = json.loads((root/'05-knowledge/results/first_exit_complement_games_20261020/palette-glue.json').read_text())
    points = [tuple(map(F, p)) for p in src['points'][:16]]
    zero = (F(0),)*4; one = (F(1),F(0),F(0),F(0))
    need(points[0] == zero, 'actual fixed source root')
    es = [(a,b) for a,b in combinations(range(16),2) if norm(add(points[a],neg(points[b]))) == one]
    need(es == [tuple(e) for e in src['edges'] if max(e) < 16], 'reconstructed entire induced source graph')
    need(len(es) == 31, 'source contact count')
    adj = [set() for _ in range(16)]
    for a,b in es:
        adj[a].add(b);adj[b].add(a)
    base = {v:(v-7)//3 for v in range(7,16)}
    repaired = dict(base); repaired[7] = 1
    def enumerate_extensions(boundary):
        need(all(boundary[a] != boundary[b] for a,b in es if a in boundary and b in boundary), 'proper actual boundary')
        output = []; witnesses = []
        for target in range(6):
            cs = dict(boundary); cs[0] = target; witness = None
            def go():
                nonlocal witness
                if len(cs) == 16:
                    if witness is None:
                        witness = tuple(cs[i] for i in range(16))
                    return 1
                allowed,i = min([(set(range(6))-{cs[j] for j in adj[i] if j in cs},i)
                                for i in range(16) if i not in cs],
                                key=lambda p:(len(p[0]),-len(adj[p[1]]),p[1]))
                total = 0
                for c in sorted(allowed):
                    cs[i] = c; total += go()
                cs.pop(i, None)
                return total
            output.append(go()); witnesses.append(witness)
        return output,witnesses
    normal_counts,normal_witnesses = enumerate_extensions(base)
    repair_counts,repair_witnesses = enumerate_extensions(repaired)
    need(normal_counts == [24,24,24,0,0,0], 'source palette relation')
    need(repair_counts == [60,84,84,8,8,8], 'one literal port repair relation')
    for witness in [w for w in normal_witnesses+repair_witnesses if w]:
        need(all(witness[a] != witness[b] for a,b in es), 'literal extension witness')

    # Reconstruct one exact six-copy geometry without importing a drawing.
    union = list(points); embeddings = [tuple(range(16))]; rotations = [None]
    expected_edges = set(es)
    for j in range(1,6):
        found = None
        for t in range(1,900*j+2):
            rho = (F(1-3*t*t,1+3*t*t),F(2*t,1+3*t*t),F(0),F(0))
            need(norm(rho) == one, 'rational rotation norm')
            copied = [mul(rho,p) for p in points]
            if any(norm(add(a,neg(b))) in (zero,one) for a in union[1:] for b in copied[1:]):
                continue
            found = t; break
        need(found is not None, 'proved finite avoidance bound')
        emb = (0,)+tuple(range(len(union),len(union)+15))
        union.extend(copied[1:]); embeddings.append(emb); rotations.append(found)
        expected_edges.update(tuple(sorted((emb[a],emb[b]))) for a,b in es)
    actual_edges = {(a,b) for a,b in combinations(range(len(union)),2)
                    if norm(add(union[a],neg(union[b]))) == one}
    need(len(set(union)) == 91, 'all actual vertices retained separately')
    need(actual_edges == expected_edges and len(actual_edges) == 186, 'no missing or extra unit contacts')
    c4 = src['four_colouring'][:16]; global4 = [-1]*len(union)
    for emb in embeddings:
        for i,v in enumerate(emb):
            need(global4[v] in (-1,c4[i]), 'root four-colour alignment')
            global4[v] = c4[i]
    need(all(global4[a] != global4[b] for a,b in actual_edges), 'actual four-colour hostile')
    near = []
    for n in range(1,7):
        palettes = [(0,1,2)]*(n//2)+[(3,4,5)]*(n-n//2)
        loads = [sum(c in a for a in palettes) for c in range(6)]
        target = max(range(6),key=lambda c:loads[c])
        colours = [-1]*(15*n+1); boundary_changes = 0
        for i,a in enumerate(palettes):
            b = tuple(c for c in range(6) if c not in a)
            perm = a+b
            changed = target not in a
            source_target = perm.index(target)
            witness = (repair_witnesses if changed else normal_witnesses)[source_target]
            need(witness is not None, 'nearest-repair witness exists')
            for u,v in enumerate(embeddings[i]):
                c = perm[witness[u]]
                need(colours[v] in (-1,c), 'shared actual output agrees')
                colours[v] = c
            boundary_changes += sum(colours[embeddings[i][v]] != a[(v-7)//3] for v in range(7,16))
        used_edges = {(a,b) for a,b in actual_edges if max(a,b) < len(colours)}
        need(all(colours[a] != colours[b] for a,b in used_edges), 'nearest repair on full induced union')
        radius = n-max(loads)
        need(radius == n//2 == boundary_changes, 'sharp physical-port radius')
        near.append(dict(copies=n,palettes=palettes,loads=loads,root=target,radius=radius,
                         vertices=len(colours),target_extensions=24**loads[target]*8**radius,
                         colouring=colours))
    palette_bank = list(combinations(range(6),3)); universe = []
    for n in range(1,5):
        checked = 0; hist = Counter()
        for palettes in combinations_with_replacement(palette_bank,n):
            formula = n-max(sum(c in a for a in palettes) for c in range(6))
            subset_cost = n
            for mask in range(1 << n):
                support = set(range(6))
                for i,a in enumerate(palettes):
                    if not mask >> i & 1:
                        support.intersection_update(a)
                if support:
                    subset_cost = min(subset_cost,mask.bit_count())
            need(subset_cost == formula, 'independent repaired-copy subset optimization')
            checked += 1; hist[formula] += 1
        universe.append(dict(copies=n,palette_multisets=checked,radius_histogram=dict(hist)))
    print(json.dumps(dict(status='FINITE-EXACT actual-port and geometry controls',
                         base_edges=es,source_counts=normal_counts,repair_counts=repair_counts,
                         repair_witnesses=repair_witnesses,rotation_parameters=rotations,
                         points=[[str(x) for x in p] for p in union],
                         actual_edges=sorted(actual_edges),four_colouring=global4,
                         unordered_distances_checked=comb(len(union),2),
                         nearest_repairs=near,palette_multiset_controls=universe),
                     indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
