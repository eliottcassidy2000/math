"""b5_depth_word_level.py - depth <= 2 and depth <= 3 hit sets of the pure chains, computed at WORD level (every
binary word is a block; no normal form, no q-parametrisation), level by level, below 2^K (K = 14).

Compared with the author's depth_cover() (scripts/chain_semigroup.py) at the same bound.  Together with
b4_depth_and_multiplicity.py (stripped normal form, K = 18) this checks C2, C3 and the sentence
'depth at most 3 misses 1225 odd numbers below 2^14, the full closure 1214'.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from b3_independent_closure import word_arrays          # noqa: E402


def expand(states, X, g, P1, Q1, pow_n, hit):
    """one more block: returns the set of live successor states; marks hits"""
    out = set()
    for (A, W) in states:
        t = 2 * W
        while A + t <= X:
            m = int(np.searchsorted(Q1, (X - A) // t, side="right"))
            vals = g[:m] * A + t * Q1[:m]
            ok = vals <= X
            A2 = vals[ok]
            hit[A2] = True
            W2 = P1[:m][ok] * A + t * pow_n[:m][ok]
            live = (A2 + 2 * W2) <= X
            out.update(zip(A2[live].tolist(), W2[live].tolist()))
            t *= 2
    return out


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    X = 1 << K
    g, P1, Q1, nn = word_arrays(K + 1)
    order = np.argsort(Q1, kind="stable")
    g, P1, Q1, nn = g[order], P1[order], Q1[order], nn[order]
    pow_n = np.int64(1) << nn
    hit1 = np.zeros(X + 1, dtype=bool)
    sel = g <= X
    hit1[g[sel]] = True
    level1 = set((a, w) for a, w in zip(g[sel].tolist(), P1[sel].tolist()) if a + 2 * w <= X)
    hit2 = np.zeros(X + 1, dtype=bool)
    level2 = expand(level1, X, g, P1, Q1, pow_n, hit2)
    hit3 = np.zeros(X + 1, dtype=bool)
    level3 = expand(level2, X, g, P1, Q1, pow_n, hit3)
    odd = np.arange(1, X + 1, 2)
    m2 = odd[~(hit1 | hit2)[odd]]
    m3 = odd[~(hit1 | hit2 | hit3)[odd]]
    print("word-level, below 2^%d: %d level-1 states, %d level-2 states, %d level-3 states" % (K, len(level1), len(level2), len(level3)))
    print("depth <= 2: %d odd numbers missed (%d in the top dyadic block); depth <= 3: %d missed (%d in the top block)"
          % (len(m2), int(np.sum(m2 >= X // 2)), len(m3), int(np.sum(m3 >= X // 2))))

    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "scripts")))
    import chain_semigroup as cs
    G = np.array(cs.g_list(X), dtype=np.int64)
    Lb = np.array([cs.L_of(p) for p in range(X + 1)], dtype=np.int64)
    Rb = np.array([cs.R_of(p) for p in range(X + 1)], dtype=np.int64)
    a1 = np.zeros(X + 1, dtype=np.uint8)
    a1[G[G <= X]] = 1
    a2 = np.zeros(X + 1, dtype=np.uint8)
    a3 = np.zeros(X + 1, dtype=np.uint8)
    cs.depth_cover(X, G, Lb, Rb, a2, a3)
    am2 = odd[(a1 | a2)[odd] == 0]
    am3 = odd[(a1 | a2 | a3)[odd] == 0]
    print("author's depth_cover at the same bound: depth <= 2: %d missed; depth <= 3: %d missed; lists identical to the word-level ones: %s, %s"
          % (len(am2), len(am3), am2.tolist() == m2.tolist(), am3.tolist() == m3.tolist()))


if __name__ == "__main__":
    main()
