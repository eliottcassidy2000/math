#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The spectral horizon table  (monad-explorer-2026-06-14).

GOAL: sharpen claim 3 of `the-spectral-resolution-ladder-of-the-ocf.md` from an
ASSERTION into a verified onset table. The reflection claims alpha_1 is the UNIQUE
OCF-derived invariant with a "delayed" spectral break (spectral past n=5). The sharp
test: is c6 (the directed 6-cycle count) non-spectral ALREADY at n=6, or does it --
like alpha_1 -- survive cospectral splitting until n=7?

Two tournaments are COSPECTRAL iff they have the same characteristic polynomial iff
the same trace vector (tr A^k)_{k}. For a tournament tr A^1 = tr A^2 = 0, so the
spectral signature is sig = (tr A^3, tr A^4, ..., tr A^n).

An invariant X is "spectral at n" iff X is constant on every cospectral class of
n-tournaments. We DETECT non-spectrality without any isomorphism testing: if two
tournaments share sig but differ in X, then X is provably non-spectral at n (and the
two are automatically non-isomorphic).

We tabulate, for each n, which of these split within some cospectral class:
  c6  (directed 6-cycles)      p33 (intersecting directed-triangle pairs, the
  c7  (directed 7-cycles)            bowtie/figure-8 carrier of tr A^6)
  c8                            TQ  (overlapping (triangle,4-cycle) pairs, carrier of tr A^7)
  alpha_2 (vertex-disjoint odd-cycle pairs)    H = I(Omega,2) = #Hamiltonian paths
"""

import sys, itertools
import numpy as np
from collections import defaultdict


# ----------------------------------------------------------------------------
def gen_labeled_tournaments(n):
    """Yield adjacency matrices of all 2^C(n,2) labeled tournaments on n vertices."""
    pairs = list(itertools.combinations(range(n), 2))
    m = len(pairs)
    for mask in range(1 << m):
        A = np.zeros((n, n), dtype=np.int64)
        for b, (i, j) in enumerate(pairs):
            if (mask >> b) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A


def directed_cycles_by_len(A, n):
    """Return dict length->list of frozenset(vertices) of simple directed cycles.
    Counts each directed cycle once (canonical: smallest vertex first, follow arcs)."""
    adj = [[j for j in range(n) if A[i, j]] for i in range(n)]
    cycles = defaultdict(list)
    for start in range(n):
        # only paths whose minimum is `start`, to count each cycle once
        stack = [(start, 1 << start, [start])]
        while stack:
            v, vis, path = stack.pop()
            for w in adj[v]:
                if w == start and len(path) >= 3:
                    cycles[len(path)].append(frozenset(path))
                elif w > start and not (vis >> w) & 1:
                    stack.append((w, vis | (1 << w), path + [w]))
    return cycles


def all_4cycles_sets(A, n):
    out = []
    seen = set()
    for verts in itertools.combinations(range(n), 4):
        a = verts[0]
        for perm in itertools.permutations(verts[1:]):
            seq = (a,) + perm
            if all(A[seq[t], seq[(t + 1) % 4]] for t in range(4)):
                edges = frozenset((seq[t], seq[(t + 1) % 4]) for t in range(4))
                if edges not in seen:
                    seen.add(edges)
                    out.append(frozenset(seq))
    return out


def H_count(A, n):
    """Number of directed Hamiltonian paths (Redei's H)."""
    # DP over subsets: dp[mask][v] = # paths covering mask ending at v
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and A[v, w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full][v] for v in range(n))


def invariants(A, n):
    cyc = directed_cycles_by_len(A, n)
    tri = cyc.get(3, [])
    c3 = len(tri)
    c5 = len(cyc.get(5, []))
    c6 = len(cyc.get(6, []))
    c7 = len(cyc.get(7, []))
    c8 = len(cyc.get(8, []))
    # p33: intersecting directed-triangle pairs
    p33 = 0
    for a in range(len(tri)):
        ta = tri[a]
        for b in range(a + 1, len(tri)):
            if ta & tri[b]:
                p33 += 1
    # TQ: overlapping (triangle, 4-cycle) pairs
    quads = all_4cycles_sets(A, n)
    TQ = sum(1 for ts in tri for qs in quads if ts & qs)
    # alpha_2: vertex-disjoint odd-cycle (3 or 5) pairs
    odd = list(tri) + list(cyc.get(5, []))
    a2 = 0
    for a in range(len(odd)):
        oa = odd[a]
        for b in range(a + 1, len(odd)):
            if not (oa & odd[b]):
                a2 += 1
    # spectral signature
    Apow = A.astype(object)
    traces = []
    P = np.eye(n, dtype=object)
    for k in range(1, n + 1):
        P = P @ Apow
        traces.append(int(np.trace(P)))
    sig = tuple(traces[2:])  # tr A^3 .. tr A^n  (tr1=tr2=0)
    H = H_count(A, n)
    return sig, dict(c3=c3, c5=c5, c6=c6, c7=c7, c8=c8,
                     p33=p33, TQ=TQ, alpha2=a2, H=H)


def run(n, sample=None, seed=0):
    import random
    random.seed(seed)
    groups = defaultdict(lambda: defaultdict(set))  # sig -> invname -> set of values
    counts = defaultdict(int)
    total = 0
    if sample is None:
        gen = gen_labeled_tournaments(n)
    else:
        def gen_sample():
            for _ in range(sample):
                A = np.zeros((n, n), dtype=np.int64)
                for i in range(n):
                    for j in range(i + 1, n):
                        if random.getrandbits(1):
                            A[i, j] = 1
                        else:
                            A[j, i] = 1
                yield A
        gen = gen_sample()
    invnames = ['c3', 'c5', 'c6', 'c7', 'c8', 'p33', 'TQ', 'alpha2', 'H']
    for A in gen:
        sig, inv = invariants(A, n)
        counts[sig] += 1
        for name in invnames:
            groups[sig][name].add(inv[name])
        total += 1

    # which invariants split within some cospectral (multi-member) class?
    split = {name: False for name in invnames}
    split_examples = {}
    n_cospectral_classes = 0   # classes proven to contain >=2 non-isomorphic members
    for sig, gd in groups.items():
        multi = any(len(gd[name]) > 1 for name in invnames)
        if multi:
            n_cospectral_classes += 1
        for name in invnames:
            if len(gd[name]) > 1:
                if not split[name]:
                    split_examples[name] = (sig, sorted(gd[name]))
                split[name] = True

    print(f"==== n={n}  ({'EXHAUSTIVE '+str(total)+' labeled' if sample is None else 'SAMPLED '+str(total)+' random'} tournaments) ====", flush=True)
    print(f"  distinct spectral signatures: {len(groups)}", flush=True)
    print(f"  signatures proven to host non-isomorphic mates (some inv splits): {n_cospectral_classes}", flush=True)
    print(f"  SPECTRAL HORIZON (does invariant split within a cospectral class?):", flush=True)
    for name in invnames:
        tag = "NON-SPECTRAL (splits)" if split[name] else "spectral (constant on every cospectral class)"
        ex = ""
        if split[name]:
            sig, vals = split_examples[name]
            ex = f"   e.g. values {vals}"
        print(f"    {name:7s}: {tag}{ex}", flush=True)
    print(flush=True)
    return split


def main():
    ns = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [4, 5, 6]
    results = {}
    for n in ns:
        if n <= 6:
            results[n] = run(n)
        else:
            # n>=7 exhaustive is 2^21+, sample heavily
            results[n] = run(n, sample=200000, seed=7)
    # summary onset table
    print("==== ONSET SUMMARY: smallest tested n at which each invariant first splits ====", flush=True)
    invnames = ['c3', 'c5', 'c6', 'c7', 'c8', 'p33', 'TQ', 'alpha2', 'H']
    for name in invnames:
        onset = next((n for n in ns if results[n].get(name)), None)
        print(f"    {name:7s}: first non-spectral at n = {onset if onset else '> '+str(max(ns))}", flush=True)


if __name__ == "__main__":
    main()
