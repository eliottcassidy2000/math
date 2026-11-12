#!/usr/bin/env python3
"""Lane twoapex - added by the FINISHING agent: independent re-computation of the finite facts used in
FINISHER_PROOFS.txt P5 (uniform limit of H modulo m over all 4^n words), for m = 3, 4, 8.

For each m: the orbit O of x(empty) modulo m under the four letters (states kept as byte strings, plain
breadth-first search - a different implementation from c1/c2), closure, strong connectivity of the transition
graph, the self-loop at the start state, and invariance of O under the translation by t * x_0, x_0 = (1, ..., 1)
(t = 1 for odd m, t = 2 for even m, where H is odd and the translation permutes the odd classes cyclically).
If all hold, the chain is irreducible and aperiodic and its stationary law is invariant under the translation,
so the admissible classes of H modulo m have equal limit shares.
"""
import json
import sys

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components


def main():
    d = json.load(open("minrep17.json"))
    r = d["dimension"]
    A = [np.array(M, dtype=np.int64) for M in d["A"]]
    x0 = np.array(d["x0"], dtype=np.int64)
    one = np.ones(r, np.int64)
    print("x_0 A_l = x_0 over Z for the four letters:", all((one @ A[l] == one).all() for l in range(4)),
          "; x(empty) A_0 = x(empty) over Z:", bool((x0 @ A[0] == x0).all()))
    for m in (3, 4, 8):
        start = (x0 % m).astype(np.uint8)
        index = {start.tobytes(): 0}
        states = [start]
        frontier = [0]
        depth = 0
        edges = [[], [], [], []]
        while frontier:
            X = np.array([states[i] for i in frontier], dtype=np.int64)
            nxt = []
            for l in range(4):
                Y = ((X @ A[l]) % m).astype(np.uint8)
                for row, src in zip(Y, frontier):
                    key = row.tobytes()
                    j = index.get(key)
                    if j is None:
                        j = len(states)
                        index[key] = j
                        states.append(row)
                        nxt.append(j)
                    edges[l].append((src, j))
            frontier = nxt
            depth += 1
        n = len(states)
        # every state was expanded exactly once, so the orbit is closed by construction
        src = np.array([e[0] for l in range(4) for e in edges[l]])
        dst = np.array([e[1] for l in range(4) for e in edges[l]])
        assert len(src) == 4 * n
        G = sp.csr_matrix((np.ones(len(src)), (src, dst)), shape=(n, n))
        ncomp, _ = connected_components(G, directed=True, connection="strong")
        loop = any(s == t == 0 for s, t in edges[0][:1])
        t = 1 if m % 2 else 2
        S = np.array(states, dtype=np.int64)
        inv = all(((row + t) % m).astype(np.uint8).tobytes() in index for row in S)
        hval = S[:, 0]
        cls = sorted(set(int(v) for v in hval))
        print(f"m = {m}: orbit {n} states (closed, depth {depth}); strongly connected components: {ncomp}; "
              f"self-loop at the start state under the letter 0: {loop}; O + {t} x_0 = O: {inv}; "
              f"classes of H on O: {cls}")


if __name__ == "__main__":
    main()
