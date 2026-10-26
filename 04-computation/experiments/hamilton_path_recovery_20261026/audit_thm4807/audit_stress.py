"""Auditor's stress test of THM-4807 beyond exhaustive orders.

Generates strong tournaments with planted untraced arcs (Busch's Cor. 4.14 structure,
recursively nested: components may themselves be planted), embeds them as blocks of
ordered joins, relabels at random, and checks with subset DP:
  (a) both halves, (b) |C_k| >= 3 / HP u->v / no HP v->u / recipe cycle in T_D,
  (c) on every arc of every strong block, the two-step rule (true >= 2, false <= 1),
  the three decoders (endpoint in DP form, G1, G4), and for n <= 11 the endpoint decoder
  and cut positions computed from the explicit list of Hamilton words.
No code from the lane is imported.

usage: python3 audit_stress.py [cases_per_n=150] [seed=20261026]
"""

import random
import sys
from collections import Counter
from itertools import combinations

from audit_lib import (
    all_hps, bits_iter, comp_index, decode_G1, decode_G4, decode_endpoint_from_words,
    hp_tables, in_masks, is_ham_word, is_legal, ordered_join, path_from_tables, popcount,
    random_size, random_strong_on, relabel, reverse_arcs, scc_ordered, traced_from_tables,
)

FAIL = Counter()


def fail(tag, detail=None):
    FAIL[tag] += 1
    if FAIL[tag] <= 3:
        print("  FAILURE", tag, detail)


def planted(rng, size_budget, depth):
    """Strong tournament with a planted untraced arc (Busch Cor 4.14 structure).
    Components may recursively be planted strong tournaments."""
    k = rng.randint(4, 5)
    sizes = []
    for t in range(k):
        sizes.append(random_size(rng, 1, 4))
    # optionally replace some component by a recursively planted one
    blocks = []
    for s in sizes:
        if depth > 0 and rng.random() < 0.3:
            blocks.append(planted(rng, size_budget // 2, depth - 1))
        else:
            blocks.append(None)
    lens = [len(b) if b is not None else s for b, s in zip(blocks, sizes)]
    n0 = sum(lens)
    out = [0] * (n0 + 1)
    comps = []
    start = 0
    for L in lens:
        comps.append(list(range(start, start + L)))
        start += L
    for a in range(k):
        for b in range(a + 1, k):
            for x in comps[a]:
                for y in comps[b]:
                    out[x] |= 1 << y
    for t, C in enumerate(comps):
        if blocks[t] is None:
            random_strong_on(C, out, rng)
        else:
            sub = blocks[t]
            for i in range(len(sub)):
                for j in bits_iter(sub[i]):
                    out[C[i]] |= 1 << C[j]
    z = n0
    x = rng.choice(comps[0])
    y = rng.choice(comps[-1])
    for t, C in enumerate(comps):
        for a in C:
            if t == 0:
                to_z = (a != x)
            elif t == k - 1:
                to_z = (a == y)
            elif t == 1:
                to_z = True
            elif t == k - 2:
                to_z = False
            else:
                to_z = rng.random() < 0.5
            if to_z:
                out[a] |= 1 << z
            else:
                out[z] |= 1 << a
    return out


def embed(rng, blocks_list):
    """Ordered join of the given strong blocks (lists of out-masks)."""
    n = sum(len(b) for b in blocks_list)
    out = [0] * n
    offs = []
    s = 0
    for b in blocks_list:
        offs.append(s)
        s += len(b)
    for a in range(len(blocks_list)):
        for i in range(len(blocks_list[a])):
            for j in bits_iter(blocks_list[a][i]):
                out[offs[a] + i] |= 1 << (offs[a] + j)
        for b in range(a + 1, len(blocks_list)):
            for i in range(len(blocks_list[a])):
                for j in range(len(blocks_list[b])):
                    out[offs[a] + i] |= 1 << (offs[b] + j)
    return out


def ends_from(n, out, inn, s, A):
    G = {1 << s: 1 << s}
    m = A & ~(1 << s)
    subs = []
    sub = 0
    while True:
        sub = (sub - m) & m
        if sub == 0:
            break
        subs.append(sub)
    for sub in subs:
        S = sub | (1 << s)
        g = 0
        for v in bits_iter(sub):
            if G.get(S ^ (1 << v), 0) & inn[v]:
                g |= 1 << v
        G[S] = g
    return G.get(A, 0)


def check(n, out, stats, explicit):
    full = (1 << n) - 1
    inn = in_masks(n, out)
    F, B = hp_tables(n, out, inn)
    tr = traced_from_tables(n, out, F, B)
    comps, chain = scc_ordered(n, out)
    if not chain:
        fail("chain")
    # decoders from trace support
    d1, i1 = decode_G1(n, tr)
    d4, i4 = decode_G4(n, tr)
    stats["G1"] += (d1 == out and i1["bad_pairs"] == 0)
    stats["G4"] += (d4 == out and i4["bad_pairs"] == 0)
    for (x, y, cxy, cyx) in i1["stats"]:
        tc, fc = (cxy, cyx) if (out[x] >> y) & 1 else (cyx, cxy)
        stats["twostep"][(tc, fc)] += 1
        if not (tc >= 2 and fc <= 1):
            fail("two-step rule", (n, out, x, y))
    # endpoint decoder, DP form
    dE = tr[:]
    idx = comp_index(n, comps)
    for x in range(n):
        for y in range(n):
            if x != y and idx[x] < idx[y]:
                dE[x] |= 1 << y
    badE = 0
    cache = {}
    for C in comps:
        cl = list(bits_iter(C))
        for a in range(len(cl)):
            for b in range(a + 1, len(cl)):
                x, y = cl[a], cl[b]
                if (tr[x] >> y) & 1 or (tr[y] >> x) & 1:
                    continue
                for s in (x, y):
                    if (s, C) not in cache:
                        cache[(s, C)] = ends_from(n, out, inn, s, C)
                exy = (cache[(x, C)] >> y) & 1
                eyx = (cache[(y, C)] >> x) & 1
                if exy == eyx:
                    badE += 1
                elif exy:
                    dE[x] |= 1 << y
                else:
                    dE[y] |= 1 << x
    stats["E"] += (dE == out and badE == 0)
    stats["pairs"] += len(i1["stats"])
    if explicit:
        words = all_hps(n, out)
        dW, iW = decode_endpoint_from_words(n, words)
        stats["E_words"] += (dW == out and iW["bad_pairs"] == 0)
        sizes = [popcount(c) for c in comps]
        bounds = [0]
        for s_ in sizes:
            bounds.append(bounds[-1] + s_)
        if iW["cuts"] != bounds:
            fail("cut positions", (n, out))
        stats["explicit"] += 1
    # (a), (b), (c) inside every strong block
    for C in comps:
        if popcount(C) < 3:
            continue
        untr = [(u, v) for u in bits_iter(C) for v in bits_iter(out[u] & C) if not (tr[u] >> v) & 1]
        stats["untraced"] += len(untr)
        for (u, v) in untr:
            cu, _ = scc_ordered(n, out, C & ~(1 << u))
            cv, _ = scc_ordered(n, out, C & ~(1 << v))
            if len(cu) < 2 or not (cu[-1] >> v) & 1:
                fail("(a) tail half", (n, out, u, v))
                continue
            if len(cv) < 2 or not (cv[0] >> u) & 1:
                fail("(a) dual half", (n, out, u, v))
            if popcount(cu[-1]) < 3:
                fail("(b) |C_k|>=3")
            if not (ends_from(n, out, inn, u, C) >> v) & 1:
                fail("(b) no HP u->v")
            if (ends_from(n, out, inn, v, C) >> u) & 1:
                fail("(b) HP v->u exists")
            c = next(bits_iter(out[u] & cu[0]))
            P = []
            for t, Cc in enumerate(cu):
                if t == 0:
                    P += path_from_tables(out, inn, F, B, Cc, start=c)
                elif t == len(cu) - 1:
                    P += path_from_tables(out, inn, F, B, Cc, end=v)
                else:
                    P += path_from_tables(out, inn, F, B, Cc)
            if not is_legal(out, [u] + P) or P[-1] != v or len(P) + 1 != popcount(C):
                fail("(b) recipe path")
            others = [a for a in untr if a != (u, v)]
            for r in range(min(len(others), 3) + 1):
                for extra in combinations(others, r):
                    TD = reverse_arcs(n, out, [(u, v)] + list(extra))
                    if not is_legal(TD, P + [u]) or not (TD[u] >> P[0]) & 1:
                        fail("(b) cycle in T_D")
                    stats["flipsets"] += 1
            if popcount(out[u]) < popcount(out[v]):
                fail("Cor 4.8")
        # (c) on every arc of the block
        for u in bits_iter(C):
            cu, _ = scc_ordered(n, out, C & ~(1 << u))
            for v in bits_iter(out[u] & C):
                pred = False
                if len(cu) >= 2 and (cu[-1] >> v) & 1:
                    Ck = cu[-1]
                    cond2 = (cu[-2] & ~out[u]) == 0
                    cond3 = True
                    rest = Ck & ~(1 << v)
                    A = rest
                    while A:
                        if (F[A] & inn[u]) and (B[Ck & ~A] >> v) & 1:
                            cond3 = False
                            break
                        A = (A - 1) & rest
                    pred = cond2 and cond3
                actual = not (tr[u] >> v) & 1
                stats["c_arcs"] += 1
                if pred != actual:
                    fail("(c)", (n, out, u, v))


def main():
    cases = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261026
    rng = random.Random(seed)
    for n in range(9, 15):
        stats = Counter()
        stats["twostep"] = Counter()
        made = 0
        tries = 0
        while made < cases and tries < 200 * cases:
            tries += 1
            blocks = []
            budget = n
            while budget > 0:
                r = rng.random()
                if budget >= 6 and r < 0.7:
                    blk = planted(rng, budget, depth=2)
                    if len(blk) > budget:
                        continue
                elif budget >= 3 and r < 0.85:
                    s = rng.randint(3, min(budget, 6))
                    blk = [0] * s
                    random_strong_on(list(range(s)), blk, rng)
                else:
                    blk = [0]
                blocks.append(blk)
                budget -= len(blk)
            out = embed(rng, blocks)
            if len(out) != n:
                continue
            perm = list(range(n))
            rng.shuffle(perm)
            out = relabel(n, out, perm)
            check(n, out, stats, explicit=(n <= 11 and made % 3 == 0))
            made += 1
            stats["cases"] += 1
        print(f"n = {n}: cases {stats['cases']}, untraced arcs {stats['untraced']}, untraced pairs "
              f"{stats['pairs']}, (c) arcs {stats['c_arcs']}, flip sets {stats['flipsets']}; decoders exact: "
              f"endpoint(DP) {stats['E']}, G1 {stats['G1']}, G4 {stats['G4']}, endpoint(words) "
              f"{stats['E_words']}/{stats['explicit']}")
        print(f"    two-step (true,false): {dict(sorted(stats['twostep'].items()))}")
        sys.stdout.flush()
    print("FAILURES:", dict(FAIL) if FAIL else "none")


if __name__ == "__main__":
    main()
