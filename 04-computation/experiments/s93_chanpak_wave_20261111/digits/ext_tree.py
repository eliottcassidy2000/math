#!/usr/bin/env python3
"""Lane "digits", part (4): the binary tree H -> 2H + 1, 2H - 1 under UNRESTRICTED one-vertex extensions.

For a tournament T let Ext(T) = { H(T + v) : all 2^|T| out-sets of the new vertex }.
K_0 = all tournaments;  K_(j+1) = { T : some extension with H = 2H(T) + 1 lies in K_j and some extension with
H = 2H(T) - 1 lies in K_j }.  depth(T) = max { j : T in K_j }.  The largest hereditary class is the intersection.

Usage: ext_tree.py ROOTMAX NMAX
  roots: every isomorphism class of order <= ROOTMAX; children are followed up to order NMAX (classes up to
  isomorphism, exact canonical forms).  The depth computed for a root of order n is min(depth(T), NMAX - n).
"""
import sys
import time
from collections import Counter, defaultdict

import numpy as np

from tour import canon, adj_from_key, ext_values, extend, iso_classes


def main():
    rootmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 11
    t0 = time.time()
    classes = iso_classes(rootmax)
    print(f"roots: isomorphism classes of orders 1..{rootmax}: {[len(classes[n]) for n in range(1, rootmax + 1)]}")
    # level[n] = dict key -> adj for the classes of order n that are roots or reachable
    level = {n: dict(classes.get(n, {})) for n in range(1, nmax + 1)}
    isroot = {n: set(classes.get(n, {}).keys()) for n in range(1, nmax + 1)}
    hval = {}
    child = {}          # key -> (list of keys with 2h+1, list of keys with 2h-1)
    nplus = {}
    nminus = {}
    for n in range(1, nmax):
        cur = level[n]
        for key, adj in cur.items():
            ev = ext_values(n, adj)
            h = int(ev[0])
            hval[key] = h
            plus = np.nonzero(ev == 2 * h + 1)[0]
            minus = np.nonzero(ev == 2 * h - 1)[0]
            nplus[key] = len(plus)
            nminus[key] = len(minus)
            cp = set()
            cm = set()
            for O in plus:
                a = extend(n, adj, int(O))
                k = canon(n + 1, a)
                cp.add(k)
                if k not in level[n + 1]:
                    level[n + 1][k] = adj_from_key(k)
            for O in minus:
                a = extend(n, adj, int(O))
                k = canon(n + 1, a)
                cm.add(k)
                if k not in level[n + 1]:
                    level[n + 1][k] = adj_from_key(k)
            child[key] = (sorted(cp), sorted(cm))
        print(f"order {n}: {len(cur)} classes examined ({len(isroot[n])} roots); order {n + 1} now holds {len(level[n + 1])} classes  [{time.time() - t0:.1f} s]", flush=True)
    for key, adj in level[nmax].items():
        hval[key] = int(ext_values(nmax, adj)[0]) if nmax <= 16 else -1
    # depth by backward induction
    depth = {}
    for key in level[nmax]:
        depth[key] = 0
    for n in range(nmax - 1, 0, -1):
        for key in level[n]:
            cp, cm = child[key]
            dp = max((depth[k] for k in cp), default=-1)
            dm = max((depth[k] for k in cm), default=-1)
            depth[key] = 1 + min(dp, dm)
    # reports
    print()
    print("ROOT STATISTICS (all isomorphism classes of the order):")
    print(" n  classes  has 2H+1  has 2H-1  both  | depth histogram (depth: classes), cap = NMAX - n")
    for n in range(1, rootmax + 1):
        keys = isroot[n]
        a = sum(1 for k in keys if nplus[k] > 0)
        b = sum(1 for k in keys if nminus[k] > 0)
        c = sum(1 for k in keys if nplus[k] > 0 and nminus[k] > 0)
        hist = Counter(depth[k] for k in keys)
        print(f"{n:2d} {len(keys):8d} {a:9d} {b:9d} {c:5d}  | {dict(sorted(hist.items()))}  cap {nmax - n}")
    print()
    print("Same statistics with the degenerate value H = 1 removed (2H - 1 = H there):")
    for n in range(1, rootmax + 1):
        keys = [k for k in isroot[n] if hval[k] > 1]
        c = sum(1 for k in keys if nplus[k] > 0 and nminus[k] > 0)
        hist = Counter(depth[k] for k in keys)
        print(f"{n:2d} {len(keys):8d} both {c:5d}  | {dict(sorted(hist.items()))}")
    print()
    # deepest roots
    for n in range(4, rootmax + 1):
        keys = sorted(isroot[n], key=lambda k: (-depth[k], hval[k]))
        top = keys[:12]
        print(f"order {n}: deepest roots (depth, H, #ext 2H+1, #ext 2H-1, key): ")
        for k in top:
            print(f"    depth {depth[k]}  H = {hval[k]:6d}  +:{nplus[k]:3d} -:{nminus[k]:3d}  key {k[1]}")
    # by H value: max depth over roots with that H
    print()
    best = defaultdict(int)
    for n in range(1, rootmax + 1):
        for k in isroot[n]:
            d = depth[k] - 0
            cap = nmax - n
            # normalise: a root that reaches its cap is 'alive'
            best[hval[k]] = max(best[hval[k]], d if d < cap else 100 + d)
    alive = sorted(h for h, d in best.items() if d >= 100)
    print(f"H-values of roots (order <= {rootmax}) that survive to the horizon (depth = NMAX - order): {alive}")
    dead = sorted((h, d) for h, d in best.items() if d < 100)
    print(f"H-values whose every root dies before the horizon, with the largest depth: {dead[:80]}")
    # level sizes and survivors per level
    print()
    print("order  classes reached  alive (depth = NMAX - order)")
    for n in range(1, nmax + 1):
        al = sum(1 for k in level[n] if depth[k] >= nmax - n)
        print(f"{n:5d} {len(level[n]):9d} {al:9d}")
    print(f"total time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
