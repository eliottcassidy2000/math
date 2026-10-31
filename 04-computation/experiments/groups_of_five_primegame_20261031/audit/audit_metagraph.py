#!/usr/bin/env python3
"""Independent audit of the n = 5 merged metagraph replication in groups_of_five_primegame_merged_metagraph_20261031.md.

Own code (nothing imported from the session).  Semantics read from the owner's explorer
inbox/processed/2026-03-05/new/tournament_n5_v5.html (buildData, bitsToAdj, canonicalize, isGridSym, flip lines,
buildOverview merged mode):
  VERTS = [n..1]; TILES = (x,y) for y = 1..n-2, x = n down to y+2; bit value 0 -> arc x -> y, 1 -> y -> x;
  base path n -> n-1 -> ... -> 1; class label #ci = rank of the lexicographically least row-major adjacency string;
  TRANS_MAP (x,y) -> (n-y+1, n-x+1); blue tiling = transpose-fixed; line = {t, complement t}, blue iff both ends blue;
  merged display order = first encounter in #ci order.
Classes are computed with a different complete invariant (least string over score-sorted orderings); the explorer's
own least string is then computed per class only to recover the owner's labels.  H by subset DP.

Run from the repository root:
    python3 04-computation/experiments/groups_of_five_primegame_20261031/audit/audit_metagraph.py
"""
import itertools
import sys

FAILS = []


def check(name, ok, info=''):
    print(f"  [{'OK' if ok else 'FAIL'}] {name}" + (f"  -- {info}" if info else ''))
    if not ok:
        FAILS.append(name)


def tiles_of(n):
    return [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1)]


def adj_owner(n, mask):
    """explorer gauge: index i <-> label n-i."""
    T = tiles_of(n)
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1
    for k, (x, y) in enumerate(T):
        xi, yi = n - x, n - y
        if (mask >> k) & 1:
            A[yi][xi] = 1
        else:
            A[xi][yi] = 1
    return A


def adj_s126(n, mask):
    """session126 gauge: vertices 1..n (index v-1), path 1 -> 2 -> ... -> n, mask 0 = ascending transitive,
    a set bit reverses the nonconsecutive arc of cell (x, y)."""
    T = tiles_of(n)
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = 1
    for k, (x, y) in enumerate(T):
        lo, hi = y - 1, x - 1
        if (mask >> k) & 1:
            A[hi][lo] = 1
        else:
            A[lo][hi] = 1
    return A


def transpose_perm(n):
    T = tiles_of(n)
    idx = {t: i for i, t in enumerate(T)}
    return [idx[(n - y + 1, n - x + 1)] for (x, y) in T]


def apply_bits(mask, perm):
    out = 0
    for i, j in enumerate(perm):
        if (mask >> i) & 1:
            out |= 1 << j
    return out


def string_of(A, order):
    return ''.join(str(A[a][b]) for a in order for b in order)


def canon_scores(A):
    """complete invariant: least string over orderings with nondecreasing scores."""
    n = len(A)
    sc = [sum(r) for r in A]
    groups = {}
    for v in range(n):
        groups.setdefault(sc[v], []).append(v)
    keys = sorted(groups)
    best = None
    for parts in itertools.product(*[itertools.permutations(groups[k]) for k in keys]):
        order = [v for part in parts for v in part]
        s = string_of(A, order)
        if best is None or s < best:
            best = s
    return best


def explorer_string(A):
    n = len(A)
    return min(string_of(A, p) for p in itertools.permutations(range(n)))


def ham_paths_dp(A):
    n = len(A)
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        f[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            if f[S][v]:
                for w in range(n):
                    if not (S >> w) & 1 and A[v][w]:
                        f[S | (1 << w)][w] += f[S][v]
    return sum(f[(1 << n) - 1])


def aut_count(A):
    n = len(A)
    return sum(1 for p in itertools.permutations(range(n))
               if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)))


def palindromic_paths(A):
    n = len(A)
    cnt = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[i + 1]] for i in range(n - 1)):
            rev = {p[i]: p[n - 1 - i] for i in range(n)}
            if all(A[u][w] == A[rev[w]][rev[u]] for u in range(n) for w in range(n) if u != w):
                cnt += 1
    return cnt


def converse(A):
    n = len(A)
    return [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]


def analyse(n, labels=False):
    m = (n - 1) * (n - 2) // 2
    full = (1 << m) - 1
    tp = transpose_perm(n)
    cls_of, members = {}, {}
    for t in range(1 << m):
        key = canon_scores(adj_owner(n, t))
        cls_of[t] = key
        members.setdefault(key, []).append(t)
    keys = list(members)
    # owner labels: rank of the explorer's least string
    if labels:
        exs = {k: explorer_string(adj_owner(n, members[k][0])) for k in keys}
        order = sorted(keys, key=lambda k: exs[k])
    else:
        order = sorted(keys)
    ci = {k: i for i, k in enumerate(order)}
    gs = {t: apply_bits(t, tp) == t for t in range(1 << m)}
    # converse class via transpose, and directly via converse tournament
    conv_ok = True
    tci = {}
    for k in keys:
        s = {ci[cls_of[apply_bits(t, tp)]] for t in members[k]}
        direct = canon_scores(converse(adj_owner(n, members[k][0])))
        conv_ok &= (len(s) == 1 and ci[direct] == next(iter(s)))
        tci[ci[k]] = next(iter(s))
    classes = []
    for k in order:
        mem = members[k]
        nb = sum(gs[t] for t in mem)
        typ = 'PB' if nb == len(mem) else ('PK' if nb == 0 else 'MX')
        A = adj_owner(n, mem[0])
        classes.append(dict(ci=ci[k], masks=mem, nb=nb, typ=typ, sc=(tci[ci[k]] == ci[k]), tci=tci[ci[k]],
                            scores=sorted(sum(r) for r in A), A=A))
    lines = []
    for t in range(1 << m):
        u = t ^ full
        if t < u:
            assert gs[t] == gs[u]
            lines.append((t, u, ci[cls_of[t]], ci[cls_of[u]], gs[t]))
    # merged display nodes
    disp, node = [], {}
    for c in classes:
        if c['ci'] in node:
            continue
        mem = [c['ci']] if c['sc'] else sorted([c['ci'], c['tci']])
        for x in mem:
            node[x] = len(disp)
        disp.append(mem)
    return dict(n=n, m=m, tp=tp, cls_of=cls_of, ci=ci, classes=classes, lines=lines, disp=disp, node=node, gs=gs,
                conv_ok=conv_ok, full=full)


def black_pairing(R):
    """sigma acts freely on black lines and twins join the same merged node pair; merged black multiplicities even."""
    n, tp, full = R['n'], R['tp'], R['full']
    line_of = {}
    for L in R['lines']:
        line_of[L[0]] = L
        line_of[L[1]] = L
    free = True
    same = True
    mult = {}
    for (a, b, ca, cb, blue) in R['lines']:
        key = tuple(sorted((R['node'][ca], R['node'][cb])))
        if not blue:
            mult[key] = mult.get(key, 0) + 1
            sa = apply_bits(a, tp)
            twin = line_of[sa]
            if twin[0] == a:
                free = False
            tkey = tuple(sorted((R['node'][twin[2]], R['node'][twin[3]])))
            same &= (tkey == key)
    fixed_tiles = sum(1 for i, j in enumerate(tp) if i == j)
    return free, same, all(v % 2 == 0 for v in mult.values()), fixed_tiles, sum(mult.values())


def main():
    print('=' * 100)
    print('B1  census n = 3..6 against THM-643 (classes, SC, PB, MX, PK, gridsym, blue lines, black lines, '
          'blue self-loops) and merged node counts; black-line pairing under the transpose')
    print('=' * 100)
    THM643 = {3: (2, 2, 2, 0, 0, 2, 1, 0, 0), 4: (4, 2, 1, 1, 2, 4, 2, 2, 1), 5: (12, 8, 3, 5, 4, 16, 8, 24, 0),
              6: (56, 12, 2, 10, 44, 64, 32, 480, 2)}
    MERGED = {3: 2, 4: 3, 5: 10, 6: 34}
    for n in (3, 4, 5, 6):
        R = analyse(n)
        C = R['classes']
        row = (len(C), sum(c['sc'] for c in C), sum(c['typ'] == 'PB' for c in C), sum(c['typ'] == 'MX' for c in C),
               sum(c['typ'] == 'PK' for c in C), sum(R['gs'].values()), sum(L[4] for L in R['lines']),
               sum(not L[4] for L in R['lines']), sum(1 for L in R['lines'] if L[4] and L[2] == L[3]))
        check(f'n={n}: census {row} = THM-643 row; merged nodes {len(R["disp"])} = {MERGED[n]}; transpose = converse '
              f'class', row == THM643[n] and len(R['disp']) == MERGED[n] and R['conv_ok'])
        free, same, even, ft, nblack = black_pairing(R)
        check(f'n={n}: transpose has {ft} = floor((n-1)/2) fixed tiles, acts freely on the {nblack} black lines, '
              f'twins join the same merged pair, all merged black multiplicities even',
              free and same and even and ft == (n - 1) // 2)

    print()
    print('=' * 100)
    print('B2  n = 5 class table in the owner gauge (labels = explorer least-string ranks)')
    print('=' * 100)
    R = analyse(5, labels=True)
    C = R['classes']
    for c in C:
        A = c['A']
        c['H'] = ham_paths_dp(A)
        c['aut'] = aut_count(A)
        c['pal'] = palindromic_paths(A)
        c['node'] = R['node'][c['ci']]
    s126 = {}
    least = sorted((min(m for ci_ in mem for m in C[ci_]['masks']), i) for i, mem in enumerate(R['disp']))
    for rank, (lm, i) in enumerate(least):
        s126[i] = (rank, lm)
    print('  #ci scores        H  Aut tilings blue type SC  conv node  s126(least mask)  palindromic paths')
    for c in C:
        print(f"  #{c['ci']:<2} {str(c['scores']):15} {c['H']:>2} {c['aut']:>3} {len(c['masks']):>6} {c['nb']:>4}  "
              f"{c['typ']}   {'Y' if c['sc'] else 'n'}  #{c['tci']:<3} {c['node']:>3}   {s126[c['node']]}   {c['pal']}")
    check('fibre law H = tilings * |Aut| and #palindromic Hamiltonian paths = #blue tilings * |Aut| for all 12 classes',
          all(c['H'] == len(c['masks']) * c['aut'] and c['pal'] == c['nb'] * c['aut'] for c in C))
    # the note's table (owner #ci -> scores, H, Aut, tilings, blue, type, display node, THM-4729/S70 node, s126 node)
    NOTE = {0: ([0, 1, 2, 3, 4], 1, 1, 1, 1, 'PB', 0, 0, (0, 0)),
            1: ([0, 1, 3, 3, 3], 3, 3, 1, 0, 'PK', 1, 1, (4, 4)), 7: ([1, 1, 1, 3, 4], 3, 3, 1, 0, 'PK', 1, 1, (4, 4)),
            2: ([0, 2, 2, 3, 3], 5, 1, 5, 0, 'PK', 2, 2, (2, 2)), 5: ([1, 1, 2, 2, 4], 5, 1, 5, 0, 'PK', 2, 2, (2, 2)),
            3: ([0, 2, 2, 2, 4], 3, 3, 1, 1, 'PB', 3, 3, (8, 16)),
            4: ([1, 1, 2, 3, 3], 9, 1, 9, 1, 'MX', 4, 6, (1, 1)),
            6: ([1, 1, 2, 3, 3], 9, 1, 9, 3, 'MX', 5, 7, (5, 5)),
            8: ([1, 2, 2, 2, 3], 11, 1, 11, 3, 'MX', 6, 5, (3, 3)),
            9: ([1, 2, 2, 2, 3], 15, 3, 5, 1, 'MX', 7, 4, (9, 17)),
            10: ([1, 2, 2, 2, 3], 13, 1, 13, 3, 'MX', 8, 8, (6, 10)),
            11: ([2, 2, 2, 2, 2], 15, 5, 3, 3, 'PB', 9, 9, (7, 11))}
    ok = True
    for c in C:
        e = NOTE[c['ci']]
        mine = (c['scores'], c['H'], c['aut'], len(c['masks']), c['nb'], c['typ'], c['node'])
        ok &= (mine == e[:7] and s126[c['node']] == e[8])
    check("note's section 5 table: scores, H, |Aut|, tilings, blue tilings, type, display node, session126 node "
          "(least mask) -- all 12 classes", ok)
    # lines
    blue = sorted(tuple(sorted((R['node'][L[2]], R['node'][L[3]]))) for L in R['lines'] if L[4])
    black = {}
    for L in R['lines']:
        if not L[4]:
            k = tuple(sorted((R['node'][L[2]], R['node'][L[3]])))
            black[k] = black.get(k, 0) + 1
    NOTE_BLUE = sorted([(0, 5), (3, 8), (4, 6), (5, 6), (5, 8), (6, 9), (7, 9), (8, 9)])
    NOTE_BLACK = {(1, 4), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (4, 5), (4, 8), (5, 7), (6, 6), (6, 8), (8, 8)}
    check("note's blue lines (display nodes), each simple, no loop", blue == NOTE_BLUE, str(blue))
    check("note's black lines, each x2", set(black) == NOTE_BLACK and all(v == 2 for v in black.values()))
    # blue graph cycle
    adjb = {}
    for a, b in blue:
        adjb.setdefault(a, set()).add(b)
        adjb.setdefault(b, set()).add(a)
    V = sorted(adjb)
    comps = 0
    seen = set()
    for v in V:
        if v not in seen:
            comps += 1
            st = [v]
            while st:
                x = st.pop()
                if x in seen:
                    continue
                seen.add(x)
                st += list(adjb[x])
    cyc_rank = len(blue) - len(V) + comps
    # find the cycle: strip leaves
    deg = {v: len(adjb[v]) for v in V}
    alive = set(V)
    changed = True
    while changed:
        changed = False
        for v in list(alive):
            if sum(1 for w in adjb[v] if w in alive) <= 1:
                alive.remove(v)
                changed = True
    check('blue graph: one component on 8 nodes, cycle rank 1, the cycle is 5-6-9-8', comps == 1 and cyc_rank == 1
          and alive == {5, 6, 8, 9})
    # support degree sequence
    sup = set(blue) | {k for k in black if k[0] != k[1]}
    loops = [k for k in black if k[0] == k[1]]
    dg = [0] * 10
    for a, b in sup:
        dg[a] += 1
        dg[b] += 1
    check('support: 18 non-loop edges, 2 loops, degree sequence (1,1,1,3,3,5,5,5,6,6)',
          len(sup) == 18 and len(loops) == 2 and sorted(dg) == [1, 1, 1, 3, 3, 5, 5, 5, 6, 6])

    print()
    print('=' * 100)
    print("B3  dictionary to S70 section 4 / THM-4729 n = 5 table and to session126's metagraph_ca README")
    print('=' * 100)
    S70_BLUE = [(0, 7), (3, 8), (4, 9), (5, 6), (5, 7), (5, 9), (7, 8), (8, 9)]
    S70_BLACK = [(1, 6), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (4, 7), (5, 5), (5, 8), (6, 7), (6, 8), (8, 8)]
    THM4729 = {0: (1, 1, 1, 1), 1: (3, 3, 2, 0), 2: (5, 1, 10, 0), 3: (3, 3, 1, 1), 4: (15, 3, 5, 1),
               5: (11, 1, 11, 3), 6: (9, 1, 9, 1), 7: (9, 1, 9, 3), 8: (13, 1, 13, 3), 9: (15, 5, 3, 3)}
    dic = [0, 1, 2, 3, 6, 7, 5, 4, 8, 9]       # note: display node -> S70 node
    tb = sorted(tuple(sorted((dic[a], dic[b]))) for a, b in blue)
    tk = sorted(tuple(sorted((dic[a], dic[b]))) for a, b in black)
    check("under the note's dictionary my blue/black lists equal S70's section 4 lists exactly",
          tb == sorted(S70_BLUE) and tk == sorted(S70_BLACK))
    # invariants per merged node: (H, Aut, tilings of the node, symmetric tilings)
    inv = {}
    for i, mem in enumerate(R['disp']):
        c0 = C[mem[0]]
        inv[i] = (c0['H'], c0['aut'], sum(len(C[x]['masks']) for x in mem), sum(C[x]['nb'] for x in mem))
    check("under the dictionary (H, |Aut|, t, s) agree with THM-4729's n = 5 table node by node",
          all(inv[i] == THM4729[dic[i]] for i in range(10)))
    # how forced is the dictionary? all colour- and multiplicity-preserving bijections onto S70's multigraph
    mine = {}
    for a, b in blue:
        mine[(a, b, 'B')] = 1
    for k, v in black.items():
        mine[(k[0], k[1], 'K')] = v
    s70 = {}
    for a, b in S70_BLUE:
        s70[(min(a, b), max(a, b), 'B')] = 1
    for a, b in S70_BLACK:
        s70[(min(a, b), max(a, b), 'K')] = 2
    def sig(G, v):
        return tuple(sorted((col, mult, a == b) for (a, b, col), mult in G.items() if v in (a, b)))
    gm, gs_ = {}, {}
    for v in range(10):
        gm.setdefault(sig(mine, v), []).append(v)
        gs_.setdefault(sig(s70, v), []).append(v)
    isos = []
    if sorted(gm) == sorted(gs_) and all(len(gm[k]) == len(gs_[k]) for k in gm):
        ks = sorted(gm)
        for parts in itertools.product(*[itertools.permutations(gs_[k]) for k in ks]):
            p = [None] * 10
            for k, part in zip(ks, parts):
                for a, b in zip(gm[k], part):
                    p[a] = b
            img = {}
            for (a, b, col), v in mine.items():
                x, y = sorted((p[a], p[b]))
                img[(x, y, col)] = v
            if img == s70:
                isos.append(tuple(p))
    inv_ok = [p for p in isos if all(inv[i] == THM4729[p[i]] for i in range(10))]
    check('the coloured multigraph alone determines the dictionary up to its automorphisms; invariants pin it',
          tuple(dic) in isos and inv_ok == [tuple(dic)],
          f'{len(isos)} colour-preserving isomorphisms, {len(inv_ok)} consistent with (H,|Aut|,t,s)')
    # session126 gauge: same merged fibres
    n = 5
    m = 6
    same = all(canon_scores(adj_s126(n, t)) == canon_scores(converse(adj_owner(n, t))) for t in range(1 << m))
    check('session126 gauge (path 1->...->5, mask 0 ascending) gives the converse tournament of the owner gauge for '
          'every mask, hence identical converse-merged fibres', same)
    node_of_mask = {}
    for c in C:
        for t in c['masks']:
            node_of_mask[t] = s126[R['node'][c['ci']]][0]
    pb = {i: sorted(t for t in node_of_mask if node_of_mask[t] == i) for i in (0, 7, 8)}
    check('session126: pure-blue fibres node 0 {0}, node 7 {11,46,53}, node 8 {16}; masks 7,17,55 in nodes 1,9,6; '
          'rule-150 fixed points 26,63,37 in nodes 5,5,6',
          pb == {0: [0], 7: [11, 46, 53], 8: [16]} and [node_of_mask[x] for x in (7, 17, 55)] == [1, 9, 6]
          and [node_of_mask[x] for x in (26, 63, 37)] == [5, 5, 6])
    # owner {7,17,55} as masks
    cl = {t: R['ci'][R['cls_of'][t]] for t in range(64)}
    info = [(x, cl[x], R['gs'][x], 63 - x, cl[63 - x]) for x in (7, 17, 55)]
    check("masks 7, 17, 55 in classes #4, #9, #10 (all mixed), grid-symmetric F/T/F; 17's line goes to 46 in #11",
          [(a, b, c_) for a, b, c_, _, _ in info] == [(7, 4, False), (17, 9, True), (55, 10, False)]
          and info[1][3:] == (46, 11), str(info))
    t5 = [L for L in R['lines'] if 11 in (L[2], L[3])]
    check('T5 (#11) masks {11,46,53}; its three lines go to three distinct nodes',
          sorted(C[11]['masks']) == [11, 46, 53] and len({L[2] if L[3] == 11 else L[3] for L in t5}) == 3)
    print()
    if FAILS:
        print('FAILED:', FAILS)
        return 1
    print('ALL AUDIT CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
