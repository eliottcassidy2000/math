#!/usr/bin/env python3
"""Verifier (lane thinorbit, S93).  Independent computation of the base-2 port-chain block matrices.

Method (different from gen_blocks.py / gen_blocks_fast.py, which SOLVE the matrix from three glued base tournaments):
every entry is computed from closed path-sum formulas inside the block B only.
  Row h:  (H(B), 2 q_o(B,w), 2 q_e(B,w)),  q_par(M,w) = sum over paths Q of B[M] ending at w, |Q| = par (mod 2), of H(B[M] - Q)
          (one mixed odd cycle through w -> u, weight 2, OCF).
  Rows E', O' (paths from the new port u2 of even / odd order, weighted by H of the rest):
      a path that stays in B, vertex set S:   w in S: H(B - S) h;   w not in S: H(B-S) h + 2 q_o(B-S) E + 2 q_e(B-S) O
      a path Q P R (Q: u2 -> w in B, P from u in T, R a path of B - Q or empty): H(T - P) H(B - Q - R).
Then (a) compare the sets of distinct matrices with the lane's files, (b) brute-force check (plain DFS path counting,
no DP, no formula) of the transfer law on random glued tournaments.
Usage: v_gen.py KMAX
"""
import sys
import random
import numpy as np
from numba import njit

D = "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit/"


@njit(cache=True)
def all_blocks(k):
    npairs = k * (k - 1) // 2
    nt = 1 << npairs
    full = (1 << k) - 1
    nS = 1 << k
    out = np.zeros((nt * k * k, 9), np.int64)
    pc = np.zeros(nS, np.int64)
    for S in range(nS):
        c = 0
        for v in range(k):
            if (S >> v) & 1:
                c += 1
        pc[S] = c
    adj = np.zeros((k, k), np.int64)
    P = np.zeros((nS, k, k), np.int64)
    Hs = np.zeros(nS, np.int64)
    endAt = np.zeros((nS, k), np.int64)
    startAt = np.zeros((nS, k), np.int64)
    q = np.zeros((2, nS), np.int64)
    idx = 0
    for bits in range(nt):
        t = 0
        for i in range(k):
            for j in range(i + 1, k):
                if (bits >> t) & 1:
                    adj[i, j] = 1
                    adj[j, i] = 0
                else:
                    adj[j, i] = 1
                    adj[i, j] = 0
                t += 1
        P[:, :, :] = 0
        for a in range(k):
            P[1 << a, a, a] = 1
        for S in range(1, nS):
            for a in range(k):
                for b in range(k):
                    c = P[S, a, b]
                    if c == 0:
                        continue
                    for tt in range(k):
                        if adj[b, tt] == 1 and ((S >> tt) & 1) == 0:
                            P[S | (1 << tt), a, tt] += c
        Hs[0] = 1
        for S in range(1, nS):
            tot = 0
            for a in range(k):
                sa = 0
                for b in range(k):
                    sa += P[S, a, b]
                startAt[S, a] = sa
                tot += sa
            for b in range(k):
                sb = 0
                for a in range(k):
                    sb += P[S, a, b]
                endAt[S, b] = sb
            Hs[S] = tot
        for w in range(k):
            for M in range(nS):
                q[0, M] = 0
                q[1, M] = 0
                if ((M >> w) & 1) == 0:
                    continue
                S = M
                while S > 0:
                    if (S >> w) & 1:
                        q[pc[S] & 1, M] += endAt[S, w] * Hs[M ^ S]
                    S = (S - 1) & M
            for u2 in range(k):
                out[idx, 0] = Hs[full]
                out[idx, 1] = 2 * q[1, full]
                out[idx, 2] = 2 * q[0, full]
                for par in range(2):
                    ch = 0
                    cE = 0
                    cO = 0
                    for S in range(1, nS):
                        if ((S >> u2) & 1) == 0:
                            continue
                        rest = full ^ S
                        if (pc[S] & 1) == par:
                            sa = startAt[S, u2]
                            if sa != 0:
                                ch += sa * Hs[rest]
                                if ((S >> w) & 1) == 0:
                                    cE += sa * 2 * q[1, rest]
                                    cO += sa * 2 * q[0, rest]
                        if (S >> w) & 1:
                            c = P[S, u2, w]
                            if c != 0:
                                S2 = rest
                                while True:
                                    val = c * Hs[S2] * Hs[rest ^ S2]
                                    if ((pc[S] + pc[S2]) & 1) == par:
                                        cE += val
                                    else:
                                        cO += val
                                    if S2 == 0:
                                        break
                                    S2 = (S2 - 1) & rest
                    row = 3 if par == 0 else 6
                    out[idx, row] = ch
                    out[idx, row + 1] = cE
                    out[idx, row + 2] = cO
                idx += 1
    return out


def read_set(path):
    s = set()
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 9:
                s.add(tuple(int(x) for x in p[:9]))
    return s


# ---------- brute force (plain DFS, no DP) ----------
def ham_count(adj, verts):
    verts = list(verts)
    if not verts:
        return 1
    cnt = 0

    def rec(v, used, left):
        nonlocal cnt
        if left == 0:
            cnt += 1
            return
        for t in verts:
            if t not in used and adj[v][t]:
                used.add(t)
                rec(t, used, left - 1)
                used.discard(t)

    for s in verts:
        rec(s, {s}, len(verts) - 1)
    return cnt


def state_brute(adj, n, port):
    allv = set(range(n))
    h = ham_count(adj, allv)
    E = O = 0
    path = [port]

    def rec():
        nonlocal E, O
        val = ham_count(adj, allv - set(path))
        if len(path) % 2 == 0:
            E += val
        else:
            O += val
        v = path[-1]
        for t in range(n):
            if t not in path and adj[v][t]:
                path.append(t)
                rec()
                path.pop()

    rec()
    return (h, E, O)


def rand_tournament(n, rng):
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def tournament_from_bits(k, bits):
    adj = [[0] * k for _ in range(k)]
    t = 0
    for i in range(k):
        for j in range(i + 1, k):
            if (bits >> t) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            t += 1
    return adj


def main():
    kmax = int(sys.argv[1])
    mats = set()
    per_order = {}
    nblocks = 0
    arrs = {}
    for k in range(1, kmax + 1):
        arr = all_blocks(k)
        arrs[k] = arr
        nblocks += arr.shape[0]
        u = np.unique(arr, axis=0)
        s = set(tuple(int(x) for x in r) for r in u)
        per_order[k] = s
        mats |= s
        neg = int((arr < 0).sum())
        print(f"order {k}: {arr.shape[0]} marked blocks, {len(s)} distinct matrices of this order, "
              f"cumulative marked {nblocks}, cumulative distinct {len(mats)}, negative entries {neg}", flush=True)
    cum = set()
    for k in range(1, kmax + 1):
        cum |= per_order[k]
        fn = D + f"blocks2_k{k}.txt"
        try:
            lane = read_set(fn)
        except FileNotFoundError:
            continue
        print(f"compare order <= {k}: mine {len(cum)}, lane {len(lane)}, equal as sets: {cum == lane}")
    if kmax >= 6:
        lane6 = read_set(D + "fast_k6_only.txt")
        print(f"order exactly 6: mine {len(per_order[6])}, lane fast_k6_only {len(lane6)}, equal: {per_order[6] == lane6}")
    # brute-force check of the transfer law with my formula matrices
    rng = random.Random(93)
    fails = 0
    ntest = 0
    for trial in range(60):
        nT = rng.choice([1, 2, 3, 3, 4, 4])
        adjT = rand_tournament(nT, rng)
        u = rng.randrange(nT)
        k = rng.choice([1, 2, 3, 3, 4]) if nT == 4 else rng.choice([1, 2, 3, 4, 4])
        k = min(k, kmax)
        npairs = k * (k - 1) // 2
        bits = rng.randrange(1 << npairs)
        w = rng.randrange(k)
        u2 = rng.randrange(k)
        N = arrs[k][(bits * k + w) * k + u2]
        adjB = tournament_from_bits(k, bits)
        n = nT + k
        adj = [[0] * n for _ in range(n)]
        for i in range(nT):
            for j in range(nT):
                adj[i][j] = adjT[i][j]
        for i in range(k):
            for j in range(k):
                adj[nT + i][nT + j] = adjB[i][j]
        for i in range(nT):
            for j in range(k):
                adj[i][nT + j] = 1
        adj[u][nT + w] = 0
        adj[nT + w][u] = 1
        s = state_brute(adjT, nT, u)
        s2 = state_brute(adj, n, nT + u2)
        pred = tuple(int(N[3 * i] * s[0] + N[3 * i + 1] * s[1] + N[3 * i + 2] * s[2]) for i in range(3))
        ntest += 1
        if pred != s2:
            fails += 1
            print("BRUTE MISMATCH", nT, adjT, u, k, bits, w, u2, s, s2, pred)
    print(f"brute-force transfer checks (plain DFS path counting): {ntest} random glued tournaments, failures {fails}")
    # save my matrices
    for k in range(1, kmax + 1):
        cumk = set()
        for j in range(1, k + 1):
            cumk |= per_order[j]
        with open(D + f"verify/v_blocks2_k{k}.txt", "w") as f:
            for M in sorted(cumk):
                f.write(" ".join(map(str, M)) + "\n")


if __name__ == "__main__":
    main()
