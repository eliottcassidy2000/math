#!/usr/bin/env python3
"""Lane "digits", part (1): rule sets that mix "below" and "above" vertices with a moving mark have NO finite state.

One mark z.  D: new v beaten by everything except z (v -> z), mark := v.
             U: new v beats everything except z (z -> v), mark := v.
(These are the two rules of the lane text, "v is beaten by everything except z" / "v beats everything except z",
with the mark moving to the new vertex.)  In (D U)^n the path  T-seg, v1, v2, T-seg, v3, v4, ...  uses n + 1
segments of the base T, so the count involves k-path covers of T for unbounded k (PROOFS.txt, P2).
Numerical confirmation: the rank of the Hankel matrix [H(word(T, z))] over bases of order <= n and all words of
length <= LEN grows with n.  For comparison the same computation with the mark kept fixed (rules D0, U0), and with
the two below-type rules (a = D0-like stay, b = D).

Usage: mixed_rules_rank.py NB LEN
"""
import itertools
import sys
import time

import numpy as np

from tour import hcount, iso_classes
from twomark import rank_mod


def apply(n, adj, z, word, rules):
    N = n + len(word)
    a = np.zeros(N, np.int64)
    for i in range(n):
        a[i] = int(adj[i])
    cur = n
    for o in word:
        kind, move = rules[o]
        v = cur
        if kind == "D":           # v -> z only
            a[v] = 1 << z
            for i in range(cur):
                if i != z:
                    a[i] |= 1 << v
        else:                     # v beats everything except z
            a[v] = ((1 << cur) - 1) & ~(1 << z)
            a[z] |= 1 << v
        if move:
            z = v
        cur += 1
    return N, a


def hankel_rank(bases, words, rules):
    Hk = [[int(hcount(*apply(n, adj, z, w, rules))) for w in words] for (n, adj, z) in bases]
    r, _ = rank_mod(Hk)
    return r


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ln = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    t0 = time.time()
    classes = iso_classes(nb)
    systems = {
        "D, U with moving mark": [("D", True), ("U", True)],
        "D, U with fixed mark": [("D", False), ("U", False)],
        "below only: D fixed (ray step a), D moving (tower step b)": [("D", False), ("D", True)],
    }
    words = [()]
    for k in range(1, ln + 1):
        words.extend(itertools.product(range(2), repeat=k))
    print(f"words: all {len(words)} words of length <= {ln} in two rules; bases: all marked tournaments of order <= n")
    for name, rules in systems.items():
        out = []
        for n in range(1, nb + 1):
            bases = [(m, adj, z) for m in range(1, n + 1) for key, adj in classes[m].items() for z in range(m)]
            r = hankel_rank(bases, words, rules)
            out.append((n, len(bases), r))
        print(f"{name}: (n, bases, Hankel rank) = {out}   [{time.time() - t0:.1f} s]", flush=True)


if __name__ == "__main__":
    main()
