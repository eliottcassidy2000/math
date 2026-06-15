#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THE OCF NON-SPECTRAL DEFECT  (monad-explorer-2026-06-15)

Builds on my THM-500/502 and kind-pasteur THM-499.  Thesis:

  The OCF Hamiltonian-path count H = I(Omega, 2) splits CANONICALLY into
     H = (a SPECTRAL skeleton, a function of the traces tr A^k)
         + (an integer-linear combination of NON-SPECTRAL cycle/overlap counts).

  The non-spectral carriers are exactly the Witt/census DEFECTS
     delta_k = W_k - (spectral cycle counts), W_k = (1/k) sum_{d|k} mu(d) tr A^{k/d},
  i.e. the simple-vs-overlap split that the spectrum (Euler product
     zeta(u) = 1/det(I-uA) = prod_k (1-u^k)^{-W_k}) cannot see.

CLEAN n=7 CLOSED FORM (proved from THM-500 + THM-502, here VERIFIED):
     H = [ 1 + 2 c3 + 2 c5 + 4 C(c3,2) - 4 W_6 ]   +   4 c6 + 2 c7
         \________________ SPECTRAL ________________/      \__ NON-SPECTRAL __/
  where  c3 = trA^3/3,  c5 = trA^5/5,  W_6 = (trA^6 - trA^3)/6  are spectral,
  and  c6 (#6-cycles), c7 (#Ham 7-cycles) are the two non-spectral carriers.

  Consequence: within any COSPECTRAL class at n=7,  Delta H = 4 Delta c6 + 2 Delta c7.
  NEW prediction (tested): c6 co-varies with c7 in the THM-500 split classes.

n=8 GENERALIZATION (tested):
     H = [ 1 + 2 c3 + 2 c5 + 4 C(c3,2) + 4 c3 c5 - 4 W_6 - 4 W_8 ]
         + 2 c7 + 4 c6 + 4 c8 + 4 Q44
  W_8 = (trA^8 - trA^4)/8 spectral; carriers c6,c7,c8,Q44 non-spectral.
"""

import sys, random
from collections import defaultdict
from math import comb
sys.setrecursionlimit(10000)


def random_tournament(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1):
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def matmul(X, Y, n):
    Z = [[0] * n for _ in range(n)]
    for i in range(n):
        Xi = X[i]
        Zi = Z[i]
        for k in range(n):
            a = Xi[k]
            if a:
                Yk = Y[k]
                for j in range(n):
                    Zi[j] += a * Yk[j]
    return Z


def trace(A, k):
    n = len(A)
    P = [row[:] for row in A]
    for _ in range(k - 1):
        P = matmul(P, A, n)
    return sum(P[i][i] for i in range(n))


def simple_cycles_by_len(A, n, maxL):
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    out = defaultdict(list)
    for start in range(n):
        stack = [(start, 1 << start, (start,))]
        while stack:
            v, vis, path = stack.pop()
            for w in adj[v]:
                if w == start and len(path) >= 3:
                    out[len(path)].append(path)
                elif w > start and not (vis >> w) & 1 and len(path) < maxL:
                    stack.append((w, vis | (1 << w), path + (w,)))
    return out


def opairs(la, lb, same):
    sa = [frozenset(x) for x in la]
    cnt = 0
    if same:
        for i in range(len(sa)):
            for j in range(i + 1, len(sa)):
                if sa[i] & sa[j]:
                    cnt += 1
    else:
        sb = [frozenset(x) for x in lb]
        for x in sa:
            for y in sb:
                if x & y:
                    cnt += 1
    return cnt


def disjoint_pairs(la, lb, same):
    """count vertex-DISJOINT pairs."""
    sa = [frozenset(x) for x in la]
    cnt = 0
    if same:
        for i in range(len(sa)):
            for j in range(i + 1, len(sa)):
                if not (sa[i] & sa[j]):
                    cnt += 1
    else:
        sb = [frozenset(x) for x in lb]
        for x in sa:
            for y in sb:
                if not (x & y):
                    cnt += 1
    return cnt


def H_dp(A, n):
    """Ground-truth Hamiltonian path count via bitmask DP."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            cnt = row[v]
            if cnt == 0 or not (mask & (1 << v)):
                continue
            Av = A[v]
            for u in range(n):
                if mask & (1 << u):
                    continue
                if Av[u]:
                    dp[mask | (1 << u)][u] += cnt
    return sum(dp[full][v] for v in range(n))


def analyze7(A):
    n = 7
    cyc = simple_cycles_by_len(A, n, 7)
    c3, c4, c5 = len(cyc[3]), len(cyc[4]), len(cyc[5])
    c6, c7 = len(cyc[6]), len(cyc[7])
    p33 = opairs(cyc[3], None, True)
    alpha2 = disjoint_pairs(cyc[3], None, True)        # disjoint triangle pairs (only odd disjoint pairs at n=7)
    tr3 = trace(A, 3); tr5 = trace(A, 5); tr6 = trace(A, 6)
    W6 = (tr6 - tr3) // 6
    H = H_dp(A, n)
    # identities
    id_ocf = (H == 1 + 2 * (c3 + c5 + c7) + 4 * alpha2)            # THM-500 OCF
    id_a2 = (alpha2 == comb(c3, 2) - p33)                          # disjoint = total - intersecting
    id_p33 = (p33 == W6 - c6)                                      # census defect
    skeleton = 1 + 2 * c3 + 2 * c5 + 4 * comb(c3, 2) - 4 * W6
    id_closed = (H == skeleton + 4 * c6 + 2 * c7)                  # THE CLEAN FORM
    sig = (tr3, trace(A, 4), tr5, tr6, trace(A, 7))                # spectral signature
    return dict(c3=c3, c5=c5, c6=c6, c7=c7, p33=p33, W6=W6, alpha2=alpha2,
                H=H, skeleton=skeleton, sig=sig,
                ok=(id_ocf, id_a2, id_p33, id_closed))


def analyze8(A):
    n = 8
    cyc = simple_cycles_by_len(A, n, 8)
    c3, c4, c5 = len(cyc[3]), len(cyc[4]), len(cyc[5])
    c6, c7, c8 = len(cyc[6]), len(cyc[7]), len(cyc[8])
    Q44 = opairs(cyc[4], None, True)
    D33 = disjoint_pairs(cyc[3], None, True)
    D35 = disjoint_pairs(cyc[3], cyc[5], False)
    alpha2 = D33 + D35                                             # all disjoint odd pairs at n=8
    tr3 = trace(A, 3); tr4 = trace(A, 4); tr5 = trace(A, 5)
    tr6 = trace(A, 6); tr8 = trace(A, 8)
    W6 = (tr6 - tr3) // 6
    W8 = (tr8 - tr4) // 8
    H = H_dp(A, n)
    id_ocf = (H == 1 + 2 * (c3 + c5 + c7) + 4 * alpha2)
    skeleton = (1 + 2 * c3 + 2 * c5 + 4 * comb(c3, 2) + 4 * c3 * c5 - 4 * W6 - 4 * W8)
    id_closed = (H == skeleton + 2 * c7 + 4 * c6 + 4 * c8 + 4 * Q44)
    sig = tuple(trace(A, k) for k in range(3, 9))
    return dict(c6=c6, c7=c7, c8=c8, Q44=Q44, H=H, skeleton=skeleton, sig=sig,
                ok=(id_ocf, id_closed))


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "7"
    NS = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
    rng = random.Random(20260615)

    if which == "7":
        n = 7
        print(f"=== OCF NON-SPECTRAL DEFECT  n={n}, {NS} random tournaments ===", flush=True)
        bad = [0, 0, 0, 0]
        classes = defaultdict(lambda: defaultdict(set))   # sig -> H -> {(c6,c7,alpha2,skeleton)}
        Hspread = defaultdict(set)
        for _ in range(NS):
            A = random_tournament(n, rng)
            r = analyze7(A)
            for i, b in enumerate(r["ok"]):
                if not b:
                    bad[i] += 1
            classes[r["sig"]][r["H"]].add((r["c6"], r["c7"], r["alpha2"], r["skeleton"]))
            Hspread[r["sig"]].add(r["H"])
        names = ["OCF H=1+2a1+4a2", "a2=C(c3,2)-p33", "p33=W6-c6", "CLOSED H=skel+4c6+2c7"]
        for nm, b in zip(names, bad):
            print(f"  [{ 'OK ' if b==0 else 'FAIL' }] {nm:28s}: {NS-b}/{NS} hold", flush=True)

        # cospectral split analysis
        split = {s: hs for s, hs in Hspread.items() if len(hs) >= 2}
        print(f"\n  cospectral classes sampled: {len(classes)};  with split H: {len(split)}", flush=True)
        print(f"  --- within-class check:  Delta H = 4 Delta c6 + 2 Delta c7 ---", flush=True)
        shown = 0
        all_ok = True
        c6_varies = 0
        for sig, hmap in classes.items():
            if len(hmap) < 2:
                continue
            # skeleton must be constant across the class
            skels = set(t[3] for hs in hmap.values() for t in hs)
            reps = []
            for Hval, tset in sorted(hmap.items()):
                c6v = sorted(set(t[0] for t in tset))
                c7v = sorted(set(t[1] for t in tset))
                reps.append((Hval, c6v, c7v))
            # verify the linear law on the (min) representatives
            base = reps[0]
            ok_class = True
            for (Hval, c6v, c7v) in reps[1:]:
                dH = Hval - base[0]
                # pick representative c6,c7 (these are determined by H within class iff law holds)
                dc6 = c6v[0] - base[1][0]
                dc7 = c7v[0] - base[2][0]
                if dH != 4 * dc6 + 2 * dc7:
                    ok_class = False
                    all_ok = False
            if len(set(t[0] for hs in hmap.values() for t in hs)) >= 2:
                c6_varies += 1
            if shown < 8:
                print(f"    sig{sig}: |skel|={len(skels)}  H={sorted(hmap)}  "
                      f"c6={sorted(set(t[0] for hs in hmap.values() for t in hs))} "
                      f"c7={sorted(set(t[1] for hs in hmap.values() for t in hs))}  "
                      f"law={'OK' if ok_class else 'FAIL'}", flush=True)
                shown += 1
        print(f"\n  within-class law Delta H=4Dc6+2Dc7 holds on all split classes: {all_ok}", flush=True)
        print(f"  split classes where c6 ALSO varies (THM-500 prediction): {c6_varies}/{len(split)}", flush=True)

    elif which == "8":
        n = 8
        print(f"=== OCF NON-SPECTRAL DEFECT  n={n}, {NS} random tournaments ===", flush=True)
        bad = [0, 0]
        for _ in range(NS):
            A = random_tournament(n, rng)
            r = analyze8(A)
            for i, b in enumerate(r["ok"]):
                if not b:
                    bad[i] += 1
        names = ["OCF H=1+2a1+4a2", "CLOSED H=skel+2c7+4c6+4c8+4Q44"]
        for nm, b in zip(names, bad):
            print(f"  [{ 'OK ' if b==0 else 'FAIL' }] {nm:34s}: {NS-b}/{NS} hold", flush=True)


if __name__ == "__main__":
    main()
