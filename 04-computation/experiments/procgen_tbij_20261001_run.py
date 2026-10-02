#!/usr/bin/env python3
"""procgen_tbij_20261001_run.py -- re-verifies every computational claim of
05-knowledge/results/procgen_tbij_20261001_natural_bijections.md   (lane tbij, OPEN-Q-060 bijective form).

Run:   nice python3 -u 04-computation/experiments/procgen_tbij_20261001_run.py [--full]
       (--full adds the exhaustive n = 10 analysis CLASS -> EVENE and n = 9 TOUR -> EVENG, about 10 minutes;
        each runs in a child process so memory is released.)
Prints to stdout only; ends with ALL CHECKS PASSED.

Sections
  A  closed forms (A049313, A002854, THM-479 branch split)
  B  P1: odd n, switching classes <-> mod-4-Eulerian tournaments
  C  Theorem G: twisted Mallows-Sloane for 8 submodules U (two independent counts); per-permutation identity
  D  Corollary G2: self-converse tournaments / switching classes = signed even-graph counts
  E  exhaustive natural-map (matching) analysis: CLASS->EVENE, TOUR->EVENG, TWOG<->EULER
  F  data behind the hand proofs (n = 5 cases; Mallows-Sloane families for even n <= 30; cycle index)
  G  block-lemma certificates (twist-rigid blocks; E1 for 5 <= n <= 100; DFGPR families)
  H  P2: bipartite Euler graphs are even
  I  P3: Higashitani-Ueyama s_{l,n} = t_{l,n}
  J  branch split: sum of rho(Aut C) = N_odd; n = 2 mod 4 matching-reversal involution
  K  automorphism-order multisets
"""
import itertools
import math
import os
import subprocess
import sys
import time
from collections import Counter
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_tbij_20261001_lib as L  # noqa: E402
import procgen_tbij_20261001_nat as N  # noqa: E402
import procgen_tbij_20261001_twist as TW  # noqa: E402
import procgen_tbij_20261001_hu as HU  # noqa: E402
import procgen_tbij_20261001_blocks as B  # noqa: E402

NCHECK = [0]
T0 = time.time()


def check(cond, msg):
    NCHECK[0] += 1
    if not cond:
        print('CHECK FAILED:', msg)
        sys.stdout.flush()
        raise SystemExit(1)


def say(*a):
    print(*a)
    sys.stdout.flush()


def el():
    return '%.1fs' % (time.time() - T0)


# ----------------------------------------------------------------------------------------------
# A
# ----------------------------------------------------------------------------------------------


def branch_split(n):
    odd = Fraction(0)
    lev = Fraction(0)
    for mu in L.partitions(n):
        if not L.is_level_type(mu):
            continue
        k = len(mu)
        o = L.pair_orbits(mu)
        if mu[0] % 2 == 1:
            odd += Fraction(2 ** (o - k + 1), L.zee(mu))
        else:
            lev += Fraction(2 ** (o - k), L.zee(mu))
    return odd, lev


def section_A():
    say('== A. closed forms', el())
    for n in range(1, 17):
        check(L.a049313_closed(n) == L.A049313[n], 'A049313 n=%d' % n)
        check(L.a002854_closed(n) == L.A002854[n], 'A002854 n=%d' % n)
    vals = {}
    for n in range(2, 13):
        o, l = branch_split(n)
        check(o + l == L.A049313[n], 'branch sum n=%d' % n)
        if n >= 3:
            check(o.denominator == 1 and l.denominator == 1, 'branch integrality n=%d' % n)
        vals[n] = (o, l)
    say('   A049313, A002854 closed forms = OEIS for n <= 16; branch split N_odd/N_lev (n=2..12):',
        [(n, str(vals[n][0]), str(vals[n][1])) for n in range(2, 13)])
    return vals


# ----------------------------------------------------------------------------------------------
# B
# ----------------------------------------------------------------------------------------------


def section_B():
    say('== B. P1: odd n, classes <-> mod-4-Eulerian tournaments', el())
    for n in (3, 5, 7, 9):
        Ts = L.gentourng(n)
        R = [T for T in Ts if L.mod4_euler(T)]
        check(len(R) == L.A049313[n], 'mod4 count n=%d' % n)
        lab = 0
        gd = N.dreadnaut_generators([(n, N.tour_nbrs(T)) for T in R], digraph=True)
        for T, (gens, gs) in zip(R, gd):
            lab += math.factorial(n) // N.group_order([g for g in gens if g != tuple(range(n))], n)
        check(lab == 2 ** math.comb(n - 1, 2), 'mod4 labelled n=%d' % n)
        # Aut(class) = Aut(member): class_orbits uses the mod-4 members as representatives
        C = N.class_orbits(n)
        byc = {o.canon: o.order for o in C}
        orders = sorted(o.order for o in C)
        orders2 = sorted(N.group_order([g for g in gens if g != tuple(range(n))], n) for (gens, gs) in gd)
        check(orders == orders2, 'Aut(class) = Aut(mod-4 member) n=%d' % n)
        check(all(o % 2 == 1 for o in orders), 'odd orders n=%d' % n)
        say('   n=%d: %d mod-4-Eulerian tournaments up to iso (= A049313), labelled 2^C(n-1,2), '
            'Aut(class) = Aut(member) (all odd order)' % (n, len(R)))
    # labelled uniqueness: every class contains exactly one mod-4-Eulerian member (n = 3, 5, 7)
    for n in (3, 5, 7):
        Tstar = None
        for T in L.gentourng(n):
            if L.mod4_euler(T):
                Tstar = T
                break
        # coset Tstar + Z (reversals along Euler graphs) = all mod-4-Eulerian labelled tournaments
        mods = TW.submodules(n)
        Z = mods['Z']
        x0 = L.flip_vector(Tstar)
        coset = set()
        for bits in range(1 << len(Z)):
            x = x0
            b = bits
            i = 0
            while b:
                if b & 1:
                    x ^= Z[i]
                b >>= 1
                i += 1
            coset.add(x)
        check(len(coset) == 2 ** math.comb(n - 1, 2), 'coset size n=%d' % n)
        bad = 0
        for x in coset:
            T = L.from_flip(x, n)
            check(L.mod4_euler(T), 'coset member mod4 n=%d' % n)
            cnt = 0
            for U in range(1 << (n - 1)):
                if L.mod4_euler(L.switch(T, U)):
                    cnt += 1
            if cnt != 1:
                bad += 1
        check(bad == 0, 'unique mod-4 member per class n=%d' % n)
        say('   n=%d: all %d labelled mod-4-Eulerian tournaments form one coset T* + Z, and each has '
            'exactly one mod-4-Eulerian member among its 2^(n-1) switchings' % (n, len(coset)))
    # 3-cycle count of the canonical member is constant mod 4
    for n in (5, 7, 9):
        R = [T for T in L.gentourng(n) if L.mod4_euler(T)]
        want = (math.comb(n, 3) - n * math.comb((n - 1) // 2, 2)) % 4
        check(all(L.cyclic_triangles(T) % 4 == want for T in R), 'c3 mod 4 n=%d' % n)
    say('   c3(mod-4-Eulerian member) = C(n,3) - n C((n-1)/2, 2) (mod 4) for every class, n = 5, 7, 9')


# ----------------------------------------------------------------------------------------------
# C
# ----------------------------------------------------------------------------------------------


def section_C(graph_cache):
    say('== C. Theorem G (twisted Mallows-Sloane) for 8 submodules U', el())
    names = ['0', 'J', 'C', 'C+J', 'Z', 'Z+J', 'CcapZ', 'C+Z']
    table = {}
    for n in range(2, 9):
        G = graph_cache[n]
        mods = TW.submodules(n)
        m = len(L.edge_list(n))
        row = []
        for name in names:
            Ub = mods[name]
            Up = TW.perp_basis(Ub, m)
            a = TW.torsor_orbits(n, Ub, Up)
            b = TW.even_graph_orbits_in(n, Up, G)
            check(a == b, 'Theorem G n=%d U=%s (%d vs %d)' % (n, name, a, b))
            row.append(a)
        table[n] = row
    for k, name in enumerate(names):
        say('   U=%-6s n=2..8: %s' % (name, [table[n][k] for n in range(2, 9)]))
    check([table[n][0] for n in range(2, 9)] == [1, 2, 4, 12, 56, 456, 6880], 'U=0 is A000568')
    check([table[n][2] for n in range(2, 9)] == [L.A049313[n] for n in range(2, 9)], 'U=C is A049313')
    # per-permutation identity, every permutation, n <= 5
    tot = 0
    for n in range(2, 6):
        mods = TW.submodules(n)
        m = len(L.edge_list(n))
        allx = range(1 << m)
        for name in ('0', 'J', 'C', 'C+J', 'Z'):
            Ub = mods[name]
            Up = TW.perp_basis(Ub, m)
            piv = TW.piv_dict(Up)
            members = [x for x in allx if TW.in_span(piv, x)]
            for g in itertools.permutations(range(n)):
                lhs = TW.fix_count_torsor(g, n, Ub, Up)
                emap = TW.perm_edge_map(g, n)
                rhs = 0
                for x in members:
                    if TW.apply_edges(emap, x) == x:
                        rhs += L.eps_orient(L.adj_from_mask(x, n), g)
                check(lhs == rhs, 'per-element identity n=%d U=%s g=%s' % (n, name, g))
                tot += 1
    say('   per-permutation identity #Fix_{T/U}(g) = sum_{X in (U^perp)^g} sgn_X(g): %d (n, U, g) triples, '
        'n <= 5, U in {0, J, C, C+J, Z}' % tot)
    return table


# ----------------------------------------------------------------------------------------------
# D
# ----------------------------------------------------------------------------------------------


def converse(T):
    n = len(T)
    full = (1 << n) - 1
    return tuple(full & ~T[i] & ~(1 << i) for i in range(n))


def section_D(graph_cache, euler_cache, class_cache):
    say('== D. Corollary G2 (converse refinements)', el())
    sc_t = []
    for n in range(2, 9):
        sc = TW.self_converse_tournaments(n)
        ev = [o for o in graph_cache[n] if o.info['even']]
        signed = sum((-1) ** o.info['edges'] for o in ev)
        check(sc == signed, 'self-converse tournaments n=%d' % n)
        sc_t.append(sc)
    say('   self-converse tournaments n=2..8 =', sc_t, '= sum over even graphs of (-1)^edges')
    sc_c = []
    for n in range(2, 10):
        C = class_cache[n]
        conv = N.class_canon_batch([converse(o.rep) for o in C])
        sc = sum(1 for o, c in zip(C, conv) if c == o.canon)
        ev = [o for o in euler_cache[n] if o.info['even']]
        signed = sum((-1) ** o.info['edges'] for o in ev)
        check(sc == signed, 'self-converse classes n=%d' % n)
        sc_c.append(sc)
    say('   self-converse switching classes n=2..9 =', sc_c, '= sum over even Euler graphs of (-1)^edges')
    ev9 = [o for o in euler_cache[9] if o.info['even']]
    e0 = sum(1 for o in ev9 if o.info['edges'] % 2 == 0)
    check((e0, len(ev9) - e0) == (484, 308), 'n=9 edge parity split')
    say('   n=9: even Euler graphs with even/odd edge count: %d / %d' % (e0, len(ev9) - e0))
    # Proposition S: odd n, self-converse classes on n = self-converse tournaments on n-1 (= A002785(n-1))
    a002785 = [1, 1, 2, 2, 8, 12, 88, 176, 2752]
    for n in (3, 5, 7, 9):
        check(sc_c[n - 2] == sc_t[n - 3] == a002785[n - 2], 'Proposition S n=%d' % n)
    check(sc_t == a002785[1:8], 'self-converse tournaments = A002785')
    say('   Proposition S: self-converse classes on n = self-converse tournaments on n-1 for n = 3,5,7,9: %s'
        % [sc_c[n - 2] for n in (3, 5, 7, 9)])


# ----------------------------------------------------------------------------------------------
# E
# ----------------------------------------------------------------------------------------------


def section_E(class_cache, euler_cache, twog_cache, tour_cache, graph_cache):
    say('== E. exhaustive natural-map analysis (Lemma N1 matching)', el())
    res = {}
    for n in range(2, 10):
        C = class_cache[n]
        Ev = [o for o in euler_cache[n] if o.info['even']]
        check(len(C) == len(Ev) == L.A049313[n], 'E1 counts n=%d' % n)
        r = N.analyse(C, Ev, n, 'euler')
        res[('CE', n)] = r
        check(r['perfect'] == (n <= 4), 'CLASS->EVENE perfect iff n<=4 (n=%d)' % n)
        say('   CLASS->EVENE n=%d: %s, symmetric classes %d, matched %d' %
            (n, 'perfect' if r['perfect'] else 'NO PERFECT MATCHING', r['nsym'], r['matched_sym']))
    for n in range(2, 9):
        T = tour_cache[n]
        Gv = [o for o in graph_cache[n] if o.info['even']]
        check(len(T) == len(Gv), 'DFGPR counts n=%d' % n)
        r = N.analyse(T, Gv, n, 'graph')
        res[('TG', n)] = r
        check(r['perfect'] == (n <= 4), 'TOUR->EVENG perfect iff n<=4 (n=%d)' % n)
        say('   TOUR->EVENG n=%d: %s, symmetric tournaments %d, matched %d' %
            (n, 'perfect' if r['perfect'] else 'NO PERFECT MATCHING', r['nsym'], r['matched_sym']))
    for n in range(2, 9):
        W = twog_cache[n]
        E = euler_cache[n]
        r3 = N.analyse(W, E, n, 'euler')
        r4 = N.analyse(E, W, n, 'twog')
        expect = (n % 2 == 1) or n == 2
        check(r3['perfect'] == expect and r4['perfect'] == expect, 'TWOG<->EULER n=%d' % n)
        say('   TWOG->EULER n=%d: %s ; EULER->TWOG: %s' % (n, 'perfect' if r3['perfect'] else 'NO',
                                                         'perfect' if r4['perfect'] else 'NO'))
    return res


# ----------------------------------------------------------------------------------------------
# F
# ----------------------------------------------------------------------------------------------


def section_F(euler_cache, tour_cache, graph_cache):
    say('== F. data behind the hand proofs', el())
    n = 5
    E = euler_cache[5]
    with3 = []
    with5 = []
    for o in E:
        els = N.group_elements(o.gens, n)
        if any(N.perm_order(e) == 3 for e in els):
            with3.append(o)
        if any(N.perm_order(e) == 5 for e in els):
            with5.append(o)
    d3 = sorted(L.nedges(o.rep) for o in with3)
    d5 = sorted(L.nedges(o.rep) for o in with5)
    check(d3 == [0, 3, 7, 10] and d5 == [0, 5, 10], 'n=5 Euler graphs with order 3/5 automorphisms')
    check([L.nedges(o.rep) for o in with3 if o.info['even']] == [0] and
          [L.nedges(o.rep) for o in with5 if o.info['even']] == [0], 'only the empty one is even')
    say('   n=5: Euler graphs with an automorphism of order 3: edge counts %s (empty, K3, K5-K3, K5); '
        'of order 5: %s (empty, C5, K5); even among them: only the empty graph' % (d3, d5))
    # DFGPR n = 5
    h = (1, 2, 0, 3, 4)
    inv = list(N.invariant_graph_masks([h], 5, False))
    check(len(inv) == 16, '16 invariant graphs')
    adjs = [L.adj_from_mask(m, 5) for m in inv]
    gd = N.dreadnaut_generators([(5, N.graph_nbrs(a)) for a in adjs], digraph=False)
    ev = [a for a, (g, s) in zip(adjs, gd) if all(N.sgn_graph(a, x) == 1 for x in g)]
    types = sorted(set(N.canon_graphs(ev)))
    degs = sorted(tuple(sorted(L.degrees(a))) for a in {c: a for c, a in zip(N.canon_graphs(ev), ev)}.values())
    check(len(types) == 4, 'four even (3,1,1)-invariant graph types')
    check(degs == [(0, 0, 0, 0, 0), (0, 1, 1, 1, 3), (1, 1, 1, 1, 4), (2, 2, 2, 3, 3)], 'types')
    T5 = tour_cache[5]
    o3 = [o for o in T5 if o.order == 3]
    o5 = [o for o in T5 if o.order == 5]
    check(len(o3) == 4 and len(o5) == 1 and all(o.order in (1, 3, 5) for o in T5), 'symmetric tournaments n=5')
    check(sorted(tuple(sorted(L.scores(o.rep), reverse=True)) for o in o3) ==
          sorted([(4, 3, 1, 1, 1), (4, 2, 2, 2, 0), (3, 2, 2, 2, 1), (3, 3, 3, 1, 0)]), 'C3 tournaments')
    say('   n=5 DFGPR: (abc)-invariant even graph types: empty, K13+K1, K14, K23; tournaments with Aut C3: '
        '4, with Aut C5: 1 (5 > 4)')
    # Mallows-Sloane families, even n <= 30
    for n in range(4, 31, 2):
        sn = [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n)) + [0])]
        sn1 = [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, n - 1)) + [0, n - 1])]
        h = n // 2
        # S_{n/2} wr S_2 on parts {0..h-1}, {h..n-1}
        sw = [tuple([1, 0] + list(range(2, n))), tuple(list(range(1, h)) + [0] + list(range(h, n))),
              tuple(list(range(h, n)) + list(range(0, h)))]
        for gens in (sn, sn1, sw):
            cls = list(N.invariant_class_masks(gens, n))
            check(len(cls) == 2, 'two invariant two-graphs n=%d' % n)
            # they are the empty class and the class of K_n
            full = L.mask_from_adj(tuple(((1 << n) - 1) & ~(1 << i) for i in range(n)))
            kn_norm = L.mask_from_adj(N.graph_switch(L.adj_from_mask(full, n), ((1 << n) - 1) & ~1))
            check(set(cls) == {0, kn_norm}, 'invariant two-graphs are E and K n=%d' % n)
        Bn = (tuple(((1 << n) - 1) & ~((1 << h) - 1) if i < h else (1 << h) - 1 for i in range(n))
              if n % 4 == 0 else
              tuple((((1 << h) - 1) & ~(1 << i)) if i < h else ((((1 << n) - 1) & ~((1 << h) - 1)) & ~(1 << i))
                    for i in range(n)))
        check(L.is_euler(Bn), 'B_n Euler n=%d' % n)
        Kn1 = tuple((((1 << (n - 1)) - 1) & ~(1 << i)) if i < n - 1 else 0 for i in range(n))
        check(L.is_euler(Kn1), 'K_{n-1}+K_1 Euler n=%d' % n)
        check(not L.is_euler(tuple(((1 << n) - 1) & ~(1 << i) for i in range(n))), 'K_n not Euler n=%d' % n)
        for G, gens in ((Bn, sw), (Kn1, sn1)):
            for g in gens:
                check(L.apply_perm_g(G, g) == G, 'group preserves the graph n=%d' % n)
    say('   Mallows-Sloane, even n = 4..30: S_n, S_{n-1}, S_{n/2} wr S_2 each fix exactly the two trivial '
        'two-graphs; the Euler graphs empty, K_{n-1}+K_1, B_n are invariant under them; K_n is not Euler')
    # cycle index equality (Brauer) per cycle type, n <= 8
    cnt = 0
    for n in range(2, 9):
        for mu in L.partitions(n):
            g = N.standard_perm(mu, n)
            a = sum(1 for _ in N.invariant_class_masks([g], n))
            b = sum(1 for _ in N.invariant_graph_masks([g], n, True))
            check(a == b, 'cycle index equality n=%d mu=%s' % (n, mu))
            check(b == L.euler_fixed_count(mu), 'fixed Euler formula')
            cnt += 1
    say('   |Fix_TWOG(g)| = |Fix_EULER(g)| for all %d cycle types with n <= 8 (equal cycle index)' % cnt)


# ----------------------------------------------------------------------------------------------
# G
# ----------------------------------------------------------------------------------------------


def section_G(e_res, class_cache, euler_cache):
    say('== G. block-lemma certificates', el())
    sizes = [b for b in B.RIGID_SIZES if b <= 100]
    for b in sizes:
        ok, cnt = B.check_twist_rigid(B.block_kind_for(b), b)
        check(ok, 'twist-rigid b=%d' % b)
    say('   twist-rigid block sizes <= 100 (every nonempty invariant graph has an odd automorphism):', sizes)
    counts = {}
    forced_sets = {}
    for n in range(5, 101):
        good = []
        for sz in B.constructions(n):
            T, gens = B.block_tournament(list(sz))
            ok, k = B.forced_to_empty(gens, n)
            check(ok, 'forced n=%d %s' % (n, sz))
            good.append((sz, T))
        cans = N.class_canon_batch([T for sz, T in good]) if good else []
        d = {}
        for (sz, T), c in zip(good, cans):
            d.setdefault(c, []).append(sz)
        check(len(d) >= 2, 'two distinct forced classes n=%d' % n)
        counts[n] = len(d)
        forced_sets[n] = set(d)
    say('   E1: every 5 <= n <= 100 has >= 2 non-isomorphic block classes forced to the empty graph; '
        'min count %d, per n (5..30): %s' % (min(counts.values()), [counts[n] for n in range(5, 31)]))
    # cross-check with the exhaustive forced lists (n <= 9)
    for n in range(5, 10):
        r = e_res[('CE', n)]
        Ev = [o for o in euler_cache[n] if o.info['even']]
        empty = [j for j, o in enumerate(Ev) if L.nedges(o.rep) == 0][0]
        forced = {class_cache[n][i].canon for i in r['sym'] if r['adj'][i] == {empty}}
        check(forced == forced_sets[n], 'exhaustive forced classes = block classes n=%d' % n)
    say('   n = 5..9: the classes forced to the empty graph (exhaustive) are exactly the block classes')
    cov = []
    for n in range(5, 101):
        rr = B.dfgpr_certificates_general(n)
        if any(x[-1] for x in rr):
            cov.append(n)
    check(set([5, 6, 7, 9, 10, 11, 13, 14]) <= set(cov), 'DFGPR certificates small n')
    say('   DFGPR block certificates (forced T0 + pattern family), n <= 100:', cov)
    # block-type Hall search: all patterns (b', 1^k), k <= 4, plus forced tournaments, with neighbourhood bounds
    cov2 = []
    defic = {}
    for n in range(5, 101):
        nl, nm, viol = B.dfgpr_block_hall(n)
        if viol is not None:
            check(viol[0] > viol[1], 'block Hall violator n=%d' % n)
            cov2.append(n)
            defic[n] = nl - nm
    check(set(range(5, 16)) <= set(cov2) and set(cov) <= set(cov2), 'block Hall coverage')
    for n in range(5, 10):
        r = e_res[('TG', n)] if ('TG', n) in e_res else None
        if r is not None:
            check(r['nsym'] - r['matched_sym'] >= defic[n], 'exhaustive deficit >= block deficit n=%d' % n)
    say('   DFGPR block-type Hall search (k <= 4), n <= 100: violation for', cov2)
    say('   block-restricted deficits (lower bounds for the true deficit):', [(n, defic[n]) for n in cov2 if n <= 15])
    return forced_sets, cov2


# ----------------------------------------------------------------------------------------------
# H
# ----------------------------------------------------------------------------------------------


def is_bipartite(adj):
    n = len(adj)
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        stack = [s]
        while stack:
            v = stack.pop()
            for u in range(n):
                if (adj[v] >> u) & 1:
                    if col[u] < 0:
                        col[u] = 1 - col[v]
                        stack.append(u)
                    elif col[u] == col[v]:
                        return False
    return True


def section_H(euler_cache):
    say('== H. P2: bipartite Euler graphs are even', el())
    rows = []
    for n in range(1, 10):
        E = euler_cache[n]
        bip = [o for o in E if is_bipartite(o.rep)]
        check(all(o.info['even'] for o in bip), 'bipartite => even n=%d' % n)
        rows.append((n, len(bip), sum(1 for o in E if o.info['even'])))
    say('   (n, #bipartite Euler graphs, #even Euler graphs):', rows)


# ----------------------------------------------------------------------------------------------
# I
# ----------------------------------------------------------------------------------------------


def section_I():
    say('== I. P3: Higashitani-Ueyama s_{l,n} = t_{l,n}', el())
    for l in (2, 3, 4, 6, 8, 9, 12):
        nmax = 7 if l <= 4 else 6
        row = []
        for n in range(1, nmax + 1):
            for mu in L.partitions(n):
                g = N.standard_perm(mu, n)
                check(HU.coker_order_A(g, n, l) == HU.fixed_eulerian(g, n, l), 'per-sigma l=%d n=%d %s' % (l, n, mu))
            s, t = HU.s_and_t(n, l)
            check(s == t, 's=t l=%d n=%d' % (l, n))
            row.append(s)
        if l == 2:
            check(row == [L.A002854[n] for n in range(1, 8)], 'l=2 is A002854')
        if l == 3:
            check(row == [1, 1, 2, 4, 14, 120, 3222], 'l=3 (A240973, Cheng-Wells)')
        if l == 4:
            check(row[:6] == [1, 1, 3, 8, 62, 1760], 'l=4 matches H-U table')
        say('   l=%2d: s = t = %s' % (l, row))
    check([HU.brute_s(n, 4) for n in range(1, 5)] == [1, 1, 3, 8], 'brute s l=4')
    check([HU.brute_t(n, 4) for n in range(1, 5)] == [1, 1, 3, 8], 'brute t l=4')
    check([HU.brute_s(n, 6) for n in range(1, 4)] == [1, 1, 4] == [HU.brute_t(n, 6) for n in range(1, 4)], 'l=6')
    say('   brute force agrees (l=4, n<=4; l=6, n<=3); per-permutation |A^s| = |K^s| for every cycle type')


# ----------------------------------------------------------------------------------------------
# J
# ----------------------------------------------------------------------------------------------


def reverse_pairs(T, pairs):
    T = list(T)
    for (x, y) in pairs:
        if (T[x] >> y) & 1:
            T[x] &= ~(1 << y)
            T[y] |= 1 << x
        else:
            T[y] &= ~(1 << x)
            T[x] |= 1 << y
    return tuple(T)


def tau_check(C, n):
    """n = 2 mod 4: classes with even-order Aut; Sylow 2-subgroup has order 2; tau: C -> C + M_g (g an
    involution of Aut C) is a well-defined fixed-point-free involution on their iso classes"""
    idx = {o.canon: i for i, o in enumerate(C)}
    ev = [i for i, o in enumerate(C) if o.order % 2 == 0]
    img = {}
    for i in ev:
        o = C[i]
        check(o.order % 4 != 0, 'Sylow 2 of order 2 (n=%d)' % n)
        els = N.group_elements(o.gens, n)
        invs = [g for g in els if N.perm_order(g) == 2]
        targets = set()
        for g in invs:
            check(all(g[x] != x for x in range(n)), 'involution fixed-point-free')
            pairs = [(x, g[x]) for x in range(n) if x < g[x]]
            targets.add(N.class_canon_batch([reverse_pairs(o.rep, pairs)])[0])
        check(len(targets) == 1, 'tau well defined')
        j = idx[targets.pop()]
        img[i] = j
    for i, j in img.items():
        check(j in img and img[j] == i, 'tau involution')
        check(j != i, 'tau fixed-point-free')
    return len(ev)


def section_J(class_cache, branch_vals):
    say('== J. branch split', el())
    for n in range(3, 10):
        C = class_cache[n]
        s = Fraction(0)
        for o in C:
            els = N.group_elements(o.gens, n)
            s += Fraction(sum(1 for e in els if N.perm_order(e) % 2 == 1), len(els))
        check(s == branch_vals[n][0], 'sum rho = N_odd n=%d' % n)
    say('   sum over classes of (odd-order fraction of Aut) = N_odd(n) for n = 3..9')
    k6 = tau_check(class_cache[6], 6)
    say('   n=6: %d classes with even-order Aut; tau (reverse the arcs of the perfect matching of an '
        'involution) pairs them without fixed points' % k6)


# ----------------------------------------------------------------------------------------------
# K
# ----------------------------------------------------------------------------------------------


def section_K(twog_cache, euler_cache, class_cache, tour_cache, graph_cache):
    say('== K. automorphism-order multisets', el())
    for n in range(3, 10):
        a = sorted(o.order for o in twog_cache[n])
        b = sorted(o.order for o in euler_cache[n])
        check((a == b) == (n % 2 == 1), 'TWOG vs EULER Aut multisets n=%d' % n)
    say('   |Aut| multisets of two-graphs and Euler graphs: equal for n = 3, 5, 7, 9, different for n = 4, 6, 8')
    for n in range(3, 10):
        a = sum(Fraction(1, o.order) for o in class_cache[n])
        b = sum(Fraction(1, o.order) for o in euler_cache[n] if o.info['even'])
        check(a == Fraction(2 ** math.comb(n - 1, 2), math.factorial(n)) and b < a, 'sum 1/|Aut| n=%d' % n)
    for n in range(3, 9):
        a = sum(Fraction(1, o.order) for o in tour_cache[n])
        b = sum(Fraction(1, o.order) for o in graph_cache[n] if o.info['even'])
        check(b < a, 'DFGPR sum 1/|Aut| n=%d' % n)
    say('   sum 1/|Aut| over CLASS exceeds that over EVENE (n = 3..9), same for TOUR vs EVENG (n = 3..8): '
        'no |Aut|-preserving bijection')


# ----------------------------------------------------------------------------------------------
# L  independent re-implementations of the automorphism groups and of the sign
# ----------------------------------------------------------------------------------------------


def section_L(class_cache, euler_cache, tour_cache, graph_cache, twog_cache):
    say('== L. independent re-implementations (backtracking automorphisms, cycle-form sign)', el())
    cnt = 0
    for n in range(3, 9):
        for o in class_cache[n]:
            a = N.group_elements(o.gens, n)
            b = set(L.class_stabilizer(o.rep))
            check(a == b, 'class Aut (S-digraph) = class Aut (backtracking) n=%d' % n)
            cnt += 1
    for n in range(1, 9):
        for o in euler_cache[n]:
            full = L.automorphisms_graph(o.rep)
            check(set(full) == N.group_elements(o.gens, n), 'Euler Aut n=%d' % n)
            ev = all(L.eps_cycle(o.rep, g) == 1 for g in full)
            ev2 = all(L.eps_arcs(o.rep, g) == 1 for g in full) if n <= 6 else ev
            check(ev == o.info['even'] == ev2, 'evenness by cycle form / arc sign n=%d' % n)
            cnt += 1
    for n in range(2, 8):
        for o in tour_cache[n]:
            check(set(L.automorphisms_tournament(o.rep)) == N.group_elements(o.gens, n), 'tournament Aut n=%d' % n)
            cnt += 1
    for n in range(2, 7):
        for o in graph_cache[n]:
            full = L.automorphisms_graph(o.rep)
            check(set(full) == N.group_elements(o.gens, n), 'graph Aut n=%d' % n)
            check(all(L.eps_cycle(o.rep, g) == 1 for g in full) == o.info['even'], 'graph evenness n=%d' % n)
            cnt += 1
    for n in range(3, 8):
        for o in twog_cache[n]:
            # two-graph automorphisms by brute force over S_n: g fixes the class iff G + gG is a cut
            x = L.mask_from_adj(o.rep)
            brute = set()
            for g in itertools.permutations(range(n)):
                if L.is_cut(x ^ L.mask_from_adj(L.apply_perm_g(o.rep, g)), n):
                    brute.add(tuple(g))
            check(brute == N.group_elements(o.gens, n), 'two-graph Aut n=%d' % n)
            cnt += 1
    say('   %d orbit representatives: nauty generators (S-digraph / incidence graph / dreadnaut) generate exactly the '
        'groups found by independent backtracking or brute force; evenness agrees with the cycle form and the arc-sign '
        'form' % cnt)


# ----------------------------------------------------------------------------------------------
# heavy parts (child processes)
# ----------------------------------------------------------------------------------------------


def part_e1n10():
    n = 10
    t = time.time()
    C = N.class_orbits(n)
    check(len(C) == L.A049313[n], 'classes n=10')
    check(sum(math.factorial(n) // o.order for o in C) == 2 ** math.comb(n - 1, 2), 'labelled n=10')
    E = N.euler_graph_orbits(n)
    check(len(E) == L.A002854[n], 'Euler n=10')
    Ev = [o for o in E if o.info['even']]
    check(len(Ev) == L.A049313[n], 'E1 n=10')
    r = N.analyse(C, Ev, n, 'euler')
    check(not r['perfect'], 'n=10 not perfect')
    empty = [j for j, o in enumerate(Ev) if L.nedges(o.rep) == 0][0]
    forced = {C[i].canon for i in r['sym'] if r['adj'][i] == {empty}}
    good = []
    for sz in B.constructions(n):
        T, gens = B.block_tournament(list(sz))
        good.append(T)
    bset = set(N.class_canon_batch(good))
    check(forced == bset, 'n=10 forced classes = block classes')
    say('   CLASS->EVENE n=10: NO PERFECT MATCHING, symmetric classes %d, matched %d; forced classes %d '
        '(= block classes (1,9), (3,7), (5,5)) [%.0fs]' % (r['nsym'], r['matched_sym'], len(forced), time.time() - t))
    k = tau_check(C, 10)
    say('   n=10: %d classes with even-order Aut, Sylow 2-subgroups of order 2; tau is a fixed-point-free '
        'involution on them' % k)
    o, lev = branch_split(10)
    check(Fraction(k, 2) == lev, 'N_lev(10) = k/2')
    say('   N_lev(10) = %s = (number of even-Aut classes)/2' % lev)


def part_dfgpr9():
    def log(s):
        pass
    r = N.tour_eveng_exhaustive(9, log)
    check(not r['perfect'] and r['ntour'] == r['neven'] == 191536, 'n=9 DFGPR')
    say('   TOUR->EVENG n=9: NO PERFECT MATCHING, symmetric tournaments %d, matched %d' % (r['nsym'], r['matched']))


def main():
    if '--part' in sys.argv:
        part = sys.argv[sys.argv.index('--part') + 1]
        {'e1n10': part_e1n10, 'dfgpr9': part_dfgpr9}[part]()
        say('PART OK %d checks' % NCHECK[0])
        return
    full = '--full' in sys.argv
    say('procgen_tbij_20261001_run.py  (full=%s)' % full)
    branch_vals = section_A()
    section_B()
    say('-- building orbit data', el())
    class_cache = {n: N.class_orbits(n) for n in range(2, 10)}
    euler_cache = {n: N.euler_graph_orbits(n) for n in range(1, 10)}
    twog_cache = {n: N.twograph_orbits(n) for n in range(2, 10)}
    tour_cache = {n: N.tournament_orbits(n) for n in range(2, 9)}
    graph_cache = {n: N.graph_orbits(n) for n in range(2, 9)}
    for n in range(2, 10):
        check(len(class_cache[n]) == L.A049313[n] and len(twog_cache[n]) == L.A002854[n] and
              len(euler_cache[n]) == L.A002854[n], 'orbit counts n=%d' % n)
        lab = 2 ** math.comb(n - 1, 2)
        check(sum(math.factorial(n) // o.order for o in class_cache[n]) == lab, 'lab classes n=%d' % n)
        check(sum(math.factorial(n) // o.order for o in twog_cache[n]) == lab, 'lab twog n=%d' % n)
        check(sum(math.factorial(n) // o.order for o in euler_cache[n]) == lab, 'lab euler n=%d' % n)
        check(sum(1 for o in euler_cache[n] if o.info['even']) == L.A049313[n], 'E1 n=%d' % n)
    for n in range(2, 9):
        check(sum(math.factorial(n) // o.order for o in tour_cache[n]) == 2 ** math.comb(n, 2), 'lab tour')
        check(sum(math.factorial(n) // o.order for o in graph_cache[n]) == 2 ** math.comb(n, 2), 'lab graph')
        check(sum(1 for o in graph_cache[n] if o.info['even']) == len(tour_cache[n]), 'DFGPR n=%d' % n)
    say('   orbit data consistent: classes = A049313, Euler = two-graphs = A002854, even Euler = A049313 '
        '(n <= 9), even graphs = tournaments (n <= 8); all labelled counts by orbit-stabiliser')
    section_C(graph_cache)
    section_D(graph_cache, euler_cache, class_cache)
    e_res = section_E(class_cache, euler_cache, twog_cache, tour_cache, graph_cache)
    section_F(euler_cache, tour_cache, graph_cache)
    section_G(e_res, class_cache, euler_cache)
    section_H(euler_cache)
    section_I()
    section_J(class_cache, branch_vals)
    section_K(twog_cache, euler_cache, class_cache, tour_cache, graph_cache)
    section_L(class_cache, euler_cache, tour_cache, graph_cache, twog_cache)
    if full:
        for part in ('e1n10', 'dfgpr9'):
            say('== heavy part %s (child process)' % part, el())
            out = subprocess.run([sys.executable, '-u', os.path.abspath(__file__), '--part', part],
                                 capture_output=True, text=True)
            sys.stdout.write(out.stdout)
            check(out.returncode == 0 and 'PART OK' in out.stdout, 'heavy part %s' % part)
            NCHECK[0] += int(out.stdout.split('PART OK ')[1].split()[0])
    say('total checks: %d, time %s' % (NCHECK[0], el()))
    say('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
