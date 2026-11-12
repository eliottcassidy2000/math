#!/usr/bin/env python3
"""Lane thinorbit (S93).  Exhaustive check of Lemma D (PROOFS.txt): in every tournament, every connected component of
the odd-cycle conflict graph Omega(T) contains a cycle meeting all other cycles of the component (a dominating
vertex), realised by an odd cycle on at least |W| - 1 of the |W| vertices of the component; hence every component
has diameter <= 2 and no "chain of cliques" C_1 - C_2 - ... - C_k with k >= 4 occurs.

All labelled tournaments on n = 3..7 vertices (2^(n(n-1)/2) of them).  Vertex sets of odd cycles are found by a
Hamiltonian-cycle DP on every odd subset; the multiplicity of directed cycles on a vertex set does not matter for
adjacency.  Also tabulated: the number of tournaments whose Omega has a component that is a chain of exactly three
classes C_1 - C_2 - C_3 (C_1, C_3 nonempty and mutually disjoint, C_2 = the cycles meeting everything), and the
maximum diameter met.
"""
import sys

import numpy as np
from numba import njit, prange


@njit(cache=True)
def analyse(n, code, res):
    adj = np.zeros((n, n), np.int64)
    t = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> t) & 1:
                adj[i, j] = 1
            else:
                adj[j, i] = 1
            t += 1
    nsub = 1 << n
    sets = np.zeros(64, np.int64)
    ns = 0
    dp = np.zeros((nsub, n), np.int64)
    for S in range(1, nsub):
        pc = 0
        for v in range(n):
            if (S >> v) & 1:
                pc += 1
        if pc < 3 or pc % 2 == 0:
            continue
        s0 = 0
        while ((S >> s0) & 1) == 0:
            s0 += 1
        # paths from s0 inside S
        for m in range(nsub):
            if (m & ~S) == 0:
                for v in range(n):
                    dp[m, v] = 0
        dp[1 << s0, s0] = 1
        found = False
        for m in range(1, nsub):
            if (m & ~S) != 0 or ((m >> s0) & 1) == 0:
                continue
            for v in range(n):
                c = dp[m, v]
                if c == 0:
                    continue
                if m == S:
                    if adj[v, s0] == 1:
                        found = True
                else:
                    for w in range(n):
                        if ((S >> w) & 1) and ((m >> w) & 1) == 0 and adj[v, w] == 1:
                            dp[m | (1 << w), w] += c
        if found:
            sets[ns] = S
            ns += 1
    # components
    comp = np.arange(ns)
    changed = True
    while changed:
        changed = False
        for i in range(ns):
            for j in range(i + 1, ns):
                if (sets[i] & sets[j]) != 0 and comp[i] != comp[j]:
                    a = min(comp[i], comp[j])
                    b = max(comp[i], comp[j])
                    for k in range(ns):
                        if comp[k] == b:
                            comp[k] = a
                    changed = True
    maxdiam = 0
    ok_dom = 1
    ok_long = 1
    three = 0
    for c in range(ns):
        if comp[c] != c:
            continue
        W = 0
        for i in range(ns):
            if comp[i] == c:
                W |= sets[i]
        wsize = 0
        for v in range(n):
            if (W >> v) & 1:
                wsize += 1
        dom = False
        domlong = False
        ndom = 0
        nmem = 0
        for i in range(ns):
            if comp[i] != c:
                continue
            nmem += 1
            alladj = True
            for j in range(ns):
                if comp[j] == c and (sets[i] & sets[j]) == 0:
                    alladj = False
                    break
            if alladj:
                dom = True
                ndom += 1
                sz = 0
                for v in range(n):
                    if (sets[i] >> v) & 1:
                        sz += 1
                if sz >= wsize - 1:
                    domlong = True
        if not dom:
            ok_dom = 0
        if not domlong:
            ok_long = 0
        # diameter
        diam = 0
        for i in range(ns):
            if comp[i] != c:
                continue
            for j in range(i + 1, ns):
                if comp[j] != c:
                    continue
                if (sets[i] & sets[j]) != 0:
                    d = 1
                else:
                    d = 3
                    for k in range(ns):
                        if comp[k] == c and (sets[i] & sets[k]) != 0 and (sets[j] & sets[k]) != 0:
                            d = 2
                            break
                if d > diam:
                    diam = d
        if diam > maxdiam:
            maxdiam = diam
        # chain of three classes: non-dominating cycles split into exactly two cliques with no edges between them
        if ndom < nmem:
            # classes of the non-dominating members under "adjacent"
            first = -1
            okchain = True
            cls = np.zeros(ns, np.int64)
            for i in range(ns):
                if comp[i] != c:
                    continue
                isdom = True
                for j in range(ns):
                    if comp[j] == c and (sets[i] & sets[j]) == 0:
                        isdom = False
                        break
                if isdom:
                    cls[i] = 0
                    continue
                if first < 0:
                    first = i
                    cls[i] = 1
                elif (sets[i] & sets[first]) != 0:
                    cls[i] = 1
                else:
                    cls[i] = 2
            for i in range(ns):
                if comp[i] != c or cls[i] == 0:
                    continue
                for j in range(i + 1, ns):
                    if comp[j] != c or cls[j] == 0:
                        continue
                    adjij = (sets[i] & sets[j]) != 0
                    if (cls[i] == cls[j]) != adjij:
                        okchain = False
            if okchain:
                three = 1
    res[0] = ok_dom
    res[1] = ok_long
    res[2] = maxdiam
    res[3] = three
    res[4] = ns


@njit(cache=True, parallel=True)
def run(n):
    npairs = n * (n - 1) // 2
    N = 1 << npairs
    bad_dom = 0
    bad_long = 0
    d3 = 0
    chain3 = 0
    hist = np.zeros(4, np.int64)
    for code in prange(N):
        res = np.zeros(5, np.int64)
        analyse(n, code, res)
        if res[0] == 0:
            bad_dom += 1
        if res[1] == 0:
            bad_long += 1
        if res[2] >= 3:
            d3 += 1
        if res[3] == 1:
            chain3 += 1
    return bad_dom, bad_long, d3, chain3


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    for n in range(3, nmax + 1):
        bad_dom, bad_long, d3, chain3 = run(n)
        print(f"n = {n}: {1 << (n * (n - 1) // 2)} labelled tournaments; with a component lacking a dominating cycle: {bad_dom}; "
              f"lacking a dominating cycle on >= |W| - 1 vertices: {bad_long}; with a component of diameter >= 3: {d3}; "
              f"with a component that is a chain of exactly three classes: {chain3}")


if __name__ == "__main__":
    main()
