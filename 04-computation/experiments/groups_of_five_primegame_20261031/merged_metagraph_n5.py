#!/usr/bin/env python3
"""mac-mini-2026-10-31-fivegroups: the merged tiling metagraph at n = 5, in the owner's own gauge.

Replicates the owner's explorer `inbox/processed/2026-03-05/new/tournament_n5_v5.html` exactly:
  * vertices labelled VERTS = [n, n-1, ..., 1] (index 0 = label n);
  * fixed base path n -> n-1 -> ... -> 1;
  * tiles TILES = [(x, y) for y = 1..n-2 for x = n down to y+2]; bit k of a mask is tile k;
  * bit 0: arc x -> y (the transitive direction); bit 1: arc y -> x;
  * class index #ci = rank of the lexicographically least row-major adjacency string among sorted classes;
  * grid transpose (x, y) -> (n-y+1, n-x+1) on tiles (downstairs: the converse); grid-symmetric = blue tiling;
  * flip = complement mask; a line is a flip pair {t, flip t}; blue iff both ends are grid-symmetric;
  * merged mode: a class and its transpose class become one node (12 classes -> 10 nodes at n = 5).

Prints the class table, the 32 lines, the merged multigraph and the three points of symmetry, and checks
the census row of THM-643 (01-canon/theorems/THM-643-blue-black-line-parity-structure.md) for n = 3..6.
Writes merged_metagraph_n5.json for the correspondence tests.

usage: python3 merged_metagraph_n5.py
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def build(n):
    verts = [n - i for i in range(n)]
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1)]
    tidx = {t: i for i, t in enumerate(tiles)}
    trans = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
    m = len(tiles)
    perms = list(itertools.permutations(range(n)))

    def adj_of(mask):
        A = [[0] * n for _ in range(n)]
        for k in range(n - 1):
            A[k][k + 1] = 1
        for i, (x, y) in enumerate(tiles):
            xi, yi = verts.index(x), verts.index(y)
            if (mask >> i) & 1:
                A[yi][xi] = 1
            else:
                A[xi][yi] = 1
        return A

    def canon(A):
        return min(''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n)) for p in perms)

    def gridsym(mask):
        return all(((mask >> i) & 1) == ((mask >> trans[i]) & 1) for i in range(m))

    def transmask(mask):
        out = 0
        for i in range(m):
            if (mask >> i) & 1:
                out |= 1 << trans[i]
        return out

    full = (1 << m) - 1
    tilings = []
    for mask in range(1 << m):
        A = adj_of(mask)
        tilings.append(dict(mask=mask, adj=A, canon=canon(A), gs=gridsym(mask),
                            flip=full ^ mask, trans=transmask(mask)))
    sigs = sorted(set(t['canon'] for t in tilings))
    ci_of = {s: i for i, s in enumerate(sigs)}
    for t in tilings:
        t['ci'] = ci_of[t['canon']]
    by_mask = {t['mask']: t for t in tilings}
    for t in tilings:
        t['tci'] = by_mask[t['trans']]['ci']
    return dict(n=n, verts=verts, tiles=tiles, trans=trans, m=m, tilings=tilings, sigs=sigs)


def sig_to_adj(sig, n):
    return [[int(sig[i * n + j]) for j in range(n)] for i in range(n)]


def ham_paths(A):
    n = len(A)
    return sum(1 for p in itertools.permutations(range(n)) if all(A[p[i]][p[i + 1]] for i in range(n - 1)))


def aut_order(A):
    n = len(A)
    return sum(1 for p in itertools.permutations(range(n))
               if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))


def cycles(A, k):
    """number of directed k-cycles (as vertex sets with a cyclic order)"""
    n = len(A)
    cnt = 0
    for S in itertools.combinations(range(n), k):
        s0 = S[0]
        for rest in itertools.permutations(S[1:]):
            cyc = (s0,) + rest
            if all(A[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                cnt += 1
    return cnt


def scc_shape(A):
    """strong components in dominance order, as sizes (a tournament's condensation is a transitive order)"""
    n = len(A)
    reach = [[A[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]:
                    reach[i][j] = True
    comps = []
    seen = set()
    for i in range(n):
        if i in seen:
            continue
        c = sorted(j for j in range(n) if reach[i][j] and reach[j][i])
        seen.update(c)
        comps.append(c)
    # order: component X before Y if X beats Y
    comps.sort(key=lambda c: -sum(A[c[0]][j] for j in range(n)))
    return [len(c) for c in comps]


def analyse(n, verbose=True):
    D = build(n)
    T, sigs, m = D['tilings'], D['sigs'], D['m']
    classes = []
    for ci, s in enumerate(sigs):
        mem = [t for t in T if t['ci'] == ci]
        A = sig_to_adj(s, n)
        H, aut = ham_paths(A), aut_order(A)
        tci = {t['tci'] for t in mem}
        assert len(tci) == 1
        tci = tci.pop()
        nblue = sum(t['gs'] for t in mem)
        typ = 'PB' if nblue == len(mem) else ('PK' if nblue == 0 else 'MX')
        scores = sorted(sum(r) for r in A)
        classes.append(dict(ci=ci, sig=s, scores=scores, fibre=len(mem), H=H, aut=aut, sc=(tci == ci), tci=tci,
                            blue=nblue, type=typ, c3=cycles(A, 3), c4=cycles(A, 4) if n >= 4 else 0,
                            c5=cycles(A, 5) if n >= 5 else 0, scc=scc_shape(A),
                            masks=[t['mask'] for t in mem], blue_masks=[t['mask'] for t in mem if t['gs']]))
        assert len(mem) * aut == H, 'LEM-003 fibre law'
    # lines
    lines = []
    for t in T:
        if t['mask'] < t['flip']:
            u = T[t['flip']]
            blue = t['gs'] and u['gs']
            assert t['gs'] == u['gs']  # flip preserves Fix(sigma)
            lines.append(dict(a=t['mask'], b=u['mask'], ca=t['ci'], cb=u['ci'], blue=blue))
    # merged nodes in the explorer's display order
    disp, ci2d = [], {}
    for c in classes:
        if c['ci'] in ci2d:
            continue
        if c['sc']:
            ci2d[c['ci']] = len(disp)
            disp.append(dict(members=[c['ci']]))
        else:
            a, b = sorted((c['ci'], c['tci']))
            ci2d[a] = ci2d[b] = len(disp)
            disp.append(dict(members=[a, b]))
    medges = {}
    for L in lines:
        x, y = sorted((ci2d[L['ca']], ci2d[L['cb']]))
        key = (x, y, 'blue' if L['blue'] else 'black')
        medges[key] = medges.get(key, 0) + 1
    census = dict(classes=len(classes), SC=sum(c['sc'] for c in classes),
                  PB=sum(c['type'] == 'PB' for c in classes), MX=sum(c['type'] == 'MX' for c in classes),
                  PK=sum(c['type'] == 'PK' for c in classes), gridsym=sum(t['gs'] for t in T),
                  blue_lines=sum(L['blue'] for L in lines), black_lines=sum(not L['blue'] for L in lines),
                  blue_self_loops=sum(1 for L in lines if L['blue'] and L['ca'] == L['cb']),
                  merged_nodes=len(disp))
    return D, classes, lines, disp, ci2d, medges, census


THM643 = {3: (2, 2, 2, 0, 0, 2, 1, 0, 0), 4: (4, 2, 1, 1, 2, 4, 2, 2, 1), 5: (12, 8, 3, 5, 4, 16, 8, 24, 0),
          6: (56, 12, 2, 10, 44, 64, 32, 480, 2)}


def main():
    print('CENSUS CHECK against the THM-643 table (classes, SC, PB, MX, PK, gridsym, blue lines, black lines, '
          'blue self-loops)')
    for n in (3, 4, 5, 6):
        _, _, _, _, _, _, cen = analyse(n, verbose=False)
        row = (cen['classes'], cen['SC'], cen['PB'], cen['MX'], cen['PK'], cen['gridsym'], cen['blue_lines'],
               cen['black_lines'], cen['blue_self_loops'])
        ok = row == THM643[n]
        print(f'  n={n}: {row}  merged nodes {cen["merged_nodes"]}  THM-643 row {THM643[n]}  match={ok}')
        assert ok
    n = 5
    D, classes, lines, disp, ci2d, medges, cen = analyse(n)
    print()
    print('n = 5, owner gauge: tiles (bit order) =', D['tiles'], ' transpose map =', D['trans'])
    print()
    print('CLASS TABLE (#ci = the owner explorer label)')
    print('  ci  scores       fibre  H   |Aut| SC  T(ci) blue type c3 c4 c5 scc      masks')
    for c in classes:
        print(f"  #{c['ci']:<2} {str(c['scores']):12} {c['fibre']:>4} {c['H']:>4} {c['aut']:>4}  "
              f"{'Y' if c['sc'] else 'n'}  #{c['tci']:<3} {c['blue']:>3}  {c['type']}  {c['c3']:>2} {c['c4']:>2} "
              f"{c['c5']:>2} {str(c['scc']):8} {c['masks']}")
    print()
    print('LINES (flip pairs {t, 63-t}); colour blue iff t is grid-symmetric')
    for L in lines:
        print(f"  {L['a']:>2} -- {L['b']:>2}   #{L['ca']} -- #{L['cb']}   {'BLUE' if L['blue'] else 'black'}")
    print()
    print('MERGED NODES (explorer display order)')
    for i, d in enumerate(disp):
        mem = d['members']
        typ = '/'.join(classes[c]['type'] for c in mem)
        print(f'  node {i}: classes {["#%d" % c for c in mem]}  type {typ}  '
              f'scores {[classes[c]["scores"] for c in mem]}')
    print()
    print('MERGED MULTIGRAPH (node, node, colour): multiplicity')
    for k in sorted(medges):
        print(f'  {k}: {medges[k]}')
    pts = [c for c in classes if c['type'] == 'PB']
    print()
    print('POINTS OF SYMMETRY (pure-blue classes):')
    for c in pts:
        print(f"  #{c['ci']} scores {c['scores']} scc {c['scc']} H={c['H']} |Aut|={c['aut']} fibre={c['fibre']} "
              f"c3={c['c3']} blue masks {c['blue_masks']}")
    print()
    print("THE OWNER'S {7, 17, 55} READ AS TILING MASKS IN THIS GAUGE:")
    by_mask = {t['mask']: t for t in D['tilings']}
    for mk in (7, 17, 55):
        t = by_mask[mk]
        c = classes[t['ci']]
        bits = [(mk >> k) & 1 for k in range(D['m'])]
        print(f"  mask {mk:>2} bits {bits} tiles on {[D['tiles'][k] for k in range(D['m']) if bits[k]]}: "
              f"class #{t['ci']} ({c['type']}, scores {c['scores']}), grid-symmetric={t['gs']}, "
              f"flip partner {t['flip']} in #{by_mask[t['flip']]['ci']}")
    out = dict(n=n, tiles=D['tiles'], trans=D['trans'], census=cen,
               classes=[{k: v for k, v in c.items()} for c in classes],
               lines=lines, merged_nodes=disp, ci_to_node={str(k): v for k, v in ci2d.items()},
               merged_edges=[dict(u=k[0], v=k[1], colour=k[2], mult=v) for k, v in sorted(medges.items())])
    with open(os.path.join(HERE, 'merged_metagraph_n5.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print()
    print('wrote merged_metagraph_n5.json')
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    sys.exit(main())
