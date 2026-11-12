#!/usr/bin/env python3
"""Lane "digits", part (4), detail: why does the binary tree H -> 2H +- 1 die at depth 2?

For every isomorphism class T of order <= ROOTMAX that has both a (2H+1)-extension and a (2H-1)-extension,
list the classes of its children and what each child admits.  Aggregate: for children reached by the map
2H+1 (type P) and by 2H-1 (type M), how many admit a further 2H'+1, a further 2H'-1, both.
"""
import sys
from collections import Counter

import numpy as np

from tour import canon, adj_from_key, ext_values, extend, iso_classes


def info(n, adj):
    ev = ext_values(n, adj)
    h = int(ev[0])
    plus = [int(O) for O in np.nonzero(ev == 2 * h + 1)[0]]
    minus = [int(O) for O in np.nonzero(ev == 2 * h - 1)[0]]
    return h, plus, minus


def scores(n, adj):
    return sorted(bin(int(a)).count("1") for a in adj)


def main():
    rootmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    classes = iso_classes(rootmax)
    agg = {"P": Counter(), "M": Counter()}
    shown = 0
    for n in range(2, rootmax + 1):
        for key, adj in classes[n].items():
            h, plus, minus = info(n, adj)
            if h == 1:
                continue
            for typ, lst in (("P", plus), ("M", minus)):
                seen = set()
                for O in lst:
                    a = extend(n, adj, O)
                    k = canon(n + 1, a)
                    if k in seen:
                        continue
                    seen.add(k)
                    h2, p2, m2 = info(n + 1, adj_from_key(k))
                    agg[typ][(n + 1, len(p2) > 0, len(m2) > 0)] += 1
            if plus and minus and shown < 25 and n >= 5:
                shown += 1
                print(f"root order {n}, H = {h}, scores {scores(n, adj)}:")
                for typ, lst in (("2H+1", plus), ("2H-1", minus)):
                    seen = set()
                    for O in lst:
                        a = extend(n, adj, O)
                        k = canon(n + 1, a)
                        if k in seen:
                            continue
                        seen.add(k)
                        h2, p2, m2 = info(n + 1, adj_from_key(k))
                        print(f"    child {typ}: out-set {O:0{n}b} (|O| = {bin(O).count('1')}), H = {h2}, scores {scores(n + 1, a)}: admits 2H'+1: {len(p2)}, 2H'-1: {len(m2)}")
    print()
    print("AGGREGATE over all roots with H > 1 (children counted per parent, up to isomorphism):")
    print("type  order  (has 2H'+1, has 2H'-1): count")
    for typ in ("P", "M"):
        for (n, a, b), c in sorted(agg[typ].items()):
            print(f"  {typ}   {n:3d}   ({a!s:5}, {b!s:5}): {c}")


if __name__ == "__main__":
    main()
