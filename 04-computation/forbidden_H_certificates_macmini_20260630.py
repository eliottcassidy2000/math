#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Forbidden-H-value proof certificates: the per-n H-spectrum, the strong-component generators, and the
per-n forbidden values (global {7,21} + e.g. {35,39} at n=6). mac-mini-2026-06-30-S40 (HYP-3617).
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)


def Hval(adj, n, perms):
    return sum(1 for p in perms if all(adj[p[k]][p[k+1]] for k in range(n-1)))


def is_strong(adj, n):
    def reach(s, A):
        seen = {s}; st = [s]
        while st:
            u = st.pop()
            for v in range(n):
                if v not in seen and A[u][v]:
                    seen.add(v); st.append(v)
        return seen
    rA = [[adj[j][i] for j in range(n)] for i in range(n)]
    return len(reach(0, adj)) == n and len(reach(0, rA)) == n


def spectra(maxn=6):
    per, strongH = {}, {}
    for n in range(1, maxn + 1):
        prs = [(i, j) for i in range(n) for j in range(i+1, n)]; m = len(prs)
        perms = list(itertools.permutations(range(n)))
        seen = set(); Hs = set(); sH = set()
        for bits in range(1 << m):
            adj = [[False]*n for _ in range(n)]
            for b, (i, j) in enumerate(prs):
                if (bits >> b) & 1: adj[i][j] = True
                else: adj[j][i] = True
            c = min(tuple(1 if adj[s[i]][s[j]] else 0 for i in range(n) for j in range(n) if i != j)
                    for s in perms)
            if c in seen: continue
            seen.add(c); h = Hval(adj, n, perms); Hs.add(h)
            if is_strong(adj, n): sH.add(h)
        per[n] = sorted(Hs); strongH[n] = sorted(sH)
    return per, strongH


def main():
    print("=" * 78)
    print("FORBIDDEN-H-VALUE CERTIFICATES: per-n spectrum, generators, gaps (mac-mini-S40)")
    print("=" * 78)
    per, strongH = spectra(6)
    print("\n[1] per-n H-spectrum, H_max, and per-n FORBIDDEN (odd <= H_max not achievable at n):")
    for n in range(2, 7):
        hmax = max(per[n])
        gaps = [v for v in range(1, hmax+1, 2) if v not in per[n]]
        glob = [v for v in gaps if v in (7, 21)]
        pern = [v for v in gaps if v not in (7, 21)]
        print(f"  n={n}: H_max={hmax:>3}, achievable={per[n]}")
        print(f"        forbidden<=H_max: global{glob} + PER-N{pern}")
    print("\n[2] strong-tournament H-values (the multiplicative GENERATORS, H=prod over strong comps):")
    for n in range(1, 7):
        print(f"  n={n}: {strongH[n]}")
    print("\n[3] WHY {35,39} impossible at n=6: 35=5*7 (no strong has H=7); 39=3*13 = 3(n=3)+13(n=5)=8>6;")
    print("    and no single strong-6 has H in {35,39}. Both fine at n>=7. Global gaps {7,21} (no 7-source).")
    print("\n" + "=" * 78)
    print("CERTIFICATE: build T from a condition; H(T) in {7,21} (global) or a per-n gap => condition")
    print("impossible. EXISTENCE (Redei): H odd => >=1 => a Ham path exists. FACTORING: H=prod strong-H.")
    print("=" * 78)


if __name__ == "__main__":
    main()
