#!/usr/bin/env python3
"""
lrc14_pairsum_cert_adversarial_klein_S209.py

HYP-5731(c) adversarial phase: try to BREAK the pair-sum certificate system
(C1 gcd-exact union ledger [mac-mini HYP-5730] ∪ C4 Hunter ledger [klein-S209])
on k>=8 covering 13-sets with mid-band members.

Hill-climb objective (minimize): #certified rulers, tie-break minimize maxLM.
Moves: replace one cluster member with a random covering-preserving speed
(keeping >= s_min mid-band members and k >= 8). Multiple restarts, caps V.

Also tracks: does exact witness supply (max_q LM) ever approach 0? (If the
adversary can kill all certificates while supply stays fat, the GAP is in the
certificates; if supply itself thins, we've found genuinely hard instances.)
"""

import random
from math import gcd

random.seed(99)

QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def midband_count(S):
    V = max(S)
    return sum(1 for v in S if 14 * v > V and 14 * v < 9 * V)


def ruler_cert_and_LM(S, q):
    """(certified_by_C1_or_C4, LM) exact, integers only."""
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    reps = []
    seen = set()
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        reps.append(r)
    if any(r == 0 for r in reps):
        return False, 0  # dead ruler
    bad = bytearray(q)
    bad[0] = 1
    Bs = []
    for r in reps:
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        Bl = set()
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    p = p0 + t * qq
                    bad[p] = 1
                    if p:
                        Bl.add(p)
        Bs.append(Bl)
    LM = q - sum(bad)
    n = len(Bs)
    tot = sum(len(b) for b in Bs)
    c1 = tot < (q - 1)   # sum over classes of |B\{0}| < q-1
    if c1:
        return True, LM
    # Hunter: subtract max spanning tree of pairwise overlaps
    w = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            w[a][b] = w[b][a] = len(Bs[a] & Bs[b])
    in_tree = [False] * n
    in_tree[0] = True
    best = [w[0][i] for i in range(n)]
    tw = 0
    for _ in range(n - 1):
        j = max((i for i in range(n) if not in_tree[i]), key=lambda i: best[i])
        tw += best[j]
        in_tree[j] = True
        for i in range(n):
            if not in_tree[i] and w[j][i] > best[i]:
                best[i] = w[j][i]
    c4 = (tot - tw) < (q - 1)
    return c4, LM


def score(S):
    """(#certified rulers, maxLM) over all pair-sum rulers."""
    rulers = sorted({a + b for i, a in enumerate(S) for b in S[i:]})
    ncert = 0
    maxlm = 0
    for q in rulers:
        c, lm = ruler_cert_and_LM(S, q)
        ncert += 1 if c else 0
        maxlm = max(maxlm, lm)
    return ncert, maxlm


def gen_start(V):
    P = random.choice([(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13)])
    k = 13 - len(P)
    L = {V}
    missed = [q for q in QS if not any(p % q == 0 for p in P)]
    for q in missed:
        if any(u % q == 0 for u in L):
            continue
        lo, hi = -(-14 // q), V // q
        if lo > hi:
            return None
        L.add(q * random.randint(lo, hi))
    while len(L) < k:
        L.add(random.randint(14, V))
    S = sorted(set(P) | L)
    if len(S) == 13 and is_covering(S) and midband_count(S) >= 2:
        return S
    return None


def main():
    print("=" * 78)
    print("Adversarial stress of C1∪C4 pair-sum certificates on mid-band covering sets")
    print("=" * 78)
    global_min_cert = None
    global_min_lm = None
    for V in (120, 200, 280):
        print(f"\n-- cap V={V} --")
        for restart in range(6):
            S = None
            while S is None:
                S = gen_start(V)
            nc, ml = score(S)
            for step in range(120):
                # mutate: swap one non-essential member
                idx = random.randrange(13)
                cand = list(S)
                newv = random.randint(14, V)
                cand[idx] = newv
                cand = sorted(set(cand))
                if len(cand) != 13 or not is_covering(cand):
                    continue
                if midband_count(cand) < 2 or max(cand) != V:
                    continue
                nc2, ml2 = score(cand)
                if (nc2, ml2) < (nc, ml):
                    S, nc, ml = cand, nc2, ml2
            print(f"  restart {restart}: min #certified={nc}, maxLM={ml}, "
                  f"s_mid={midband_count(S)}  S={S}")
            if global_min_cert is None or (nc, ml) < (global_min_cert, global_min_lm):
                global_min_cert, global_min_lm = nc, ml
                worst = S
    print(f"\nGLOBAL adversarial minimum: #certified={global_min_cert}, maxLM={global_min_lm}")
    print(f"  worst instance: {worst}")
    print("  (certificates BROKEN iff #certified reaches 0 while the set is covering)")
    print("DONE.")


if __name__ == '__main__':
    main()
