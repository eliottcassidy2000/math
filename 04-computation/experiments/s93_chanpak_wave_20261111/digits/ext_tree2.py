#!/usr/bin/env python3
"""Lane "digits", part (4), two-vertex steps: the binary tree H -> 2H + 1, 2H - 1 when a step may add ONE or TWO
new vertices with arbitrary arcs (to the old tournament and between them).  This contains every one-vertex rule,
every two-vertex gadget rule with any marks, and "pad by a source or a sink, then add one vertex".

K_0 = all;  K_(j+1) = { T : some step with H = 2H(T) + 1 lands in K_j and some step with H = 2H(T) - 1 lands in K_j }.
Roots: every isomorphism class of order <= ROOTMAX.  Children are kept up to isomorphism (exact canonical forms).

Usage: ext_tree2.py ROOTMAX MAXDEPTH [NPROC]
"""
import sys
import time
from collections import Counter
from multiprocessing import Pool

import numpy as np

from tour import adj_from_key, canon, ext_values, extend, iso_classes


def children(key):
    """(h, sorted list of child keys with 2h+1, with 2h-1) for one- and two-vertex steps."""
    n = key[0]
    adj = adj_from_key(key)
    ev = ext_values(n, adj)
    h = int(ev[0])
    cp = set()
    cm = set()
    for O in np.nonzero(ev == 2 * h + 1)[0]:
        cp.add(canon(n + 1, extend(n, adj, int(O))))
    for O in np.nonzero(ev == 2 * h - 1)[0]:
        cm.add(canon(n + 1, extend(n, adj, int(O))))
    for Ov in range(1 << n):
        a1 = extend(n, adj, Ov)
        ev2 = ext_values(n + 1, a1)
        for O in np.nonzero(ev2 == 2 * h + 1)[0]:
            cp.add(canon(n + 2, extend(n + 1, a1, int(O))))
        for O in np.nonzero(ev2 == 2 * h - 1)[0]:
            cm.add(canon(n + 2, extend(n + 1, a1, int(O))))
    return key, h, sorted(cp), sorted(cm)


def main():
    rootmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    maxdepth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    t0 = time.time()
    classes = iso_classes(rootmax)
    roots = [k for n in range(1, rootmax + 1) for k in classes[n]]
    print(f"roots: {len(roots)} isomorphism classes of orders 1..{rootmax}")
    info = {}
    frontier = list(roots)
    # level j: nodes whose children are known = everything needed to decide K_j for the roots
    known = set()
    with Pool(nproc) as pool:
        for level in range(1, maxdepth + 1):
            todo = [k for k in frontier if k not in info]
            for key, h, cp, cm in pool.imap_unordered(children, todo, chunksize=4):
                info[key] = (h, cp, cm)
            print(f"level {level}: children computed for {len(todo)} tournaments (orders {Counter(k[0] for k in todo)}), total known {len(info)}  [{time.time() - t0:.1f} s]", flush=True)

            # which nodes are in K_1 .. K_level (relative to what is known)
            memo = {}

            def inK(key, j):
                if j == 0:
                    return True
                if (key, j) in memo:
                    return memo[(key, j)]
                if key not in info:
                    return None          # unknown
                h, cp, cm = info[key]
                okp = any(inK(c, j - 1) for c in cp)
                okm = any(inK(c, j - 1) for c in cm)
                res = bool(okp and okm)
                memo[(key, j)] = res
                return res

            for j in range(1, level + 1):
                cnt = Counter()
                for k in roots:
                    if inK(k, j):
                        cnt[k[0]] += 1
                print(f"   roots in K_{j} by order: {dict(sorted(cnt.items()))}")
            # next frontier: children (both signs) of the nodes that are in K_(level) candidates, i.e. of every node
            # reachable from a root in K_level through `level` steps; simply: children of roots in K_1, then their
            # children when the parent is in K_1, ...
            nxt = set()
            cur = [k for k in roots if inK(k, 1)]
            seen = set(cur)
            for step in range(level):
                new = []
                for k in cur:
                    h, cp, cm = info[k]
                    for c in cp + cm:
                        if c not in seen:
                            seen.add(c)
                            new.append(c)
                cur = [c for c in new if (c not in info) or (info[c][1] and info[c][2])]
                for c in new:
                    if c not in info:
                        nxt.add(c)
            frontier = sorted(nxt)
            if not frontier:
                print("   no further tournaments to examine")
                break
            # report the K_level roots
            deep = [k for k in roots if inK(k, level)]
            hv = sorted(set(info[k][0] for k in deep))
            print(f"   roots in K_{level}: {len(deep)}; their H-values: {hv[:60]}{' ...' if len(hv) > 60 else ''}")
            print(f"   next: {len(frontier)} new tournaments (orders {dict(Counter(k[0] for k in frontier))})", flush=True)
    print(f"total time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
