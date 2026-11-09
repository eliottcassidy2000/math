#!/usr/bin/env python3
"""Block matrices of the port-chain transfer for I(Omega, x) (THM-5070).

A port chain is built from a single vertex u (the port) by repeatedly gluing a block (B, w, u'):
  T' = T (+) B (ordinal sum, every vertex of T beats every vertex of B), except that the single arc
  between the port u of T and w in B is reversed (w -> u);  the new port is u' in B.

State of (T, u): (h, E, O) with
  h = I(Omega(T), x),
  E = sum over paths P in T starting at u with |V(P)| even of I(Omega(T - V(P)), x),
  O = the same over paths with |V(P)| odd (the one-vertex path {u} included).
The new state is N(B, w, u') (h, E, O) with N a nonnegative integer 3x3 matrix (proof in PROOF.md).

Usage: chain_blocks.py X KMAX OUTFILE   (writes one matrix per line: 9 integers, rows h', E', O')
"""
import itertools
import sys


def tournaments(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << len(pairs)):
        adj = [[0] * k for _ in range(k)]
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj


def odd_cycle_sets(k, adj):
    """Vertex masks of odd directed cycles (each vertex set listed once per distinct cycle)."""
    out = []
    for size in range(3, k + 1, 2):
        for sub in itertools.combinations(range(k), size):
            s = sub[0]
            mask = 0
            for v in sub:
                mask |= 1 << v
            for perm in itertools.permutations(sub[1:]):
                seq = (s,) + perm
                if all(adj[seq[i]][seq[(i + 1) % size]] for i in range(size)):
                    out.append(mask)
    return out


def make_Hx(k, adj, x):
    cycles = odd_cycle_sets(k, adj)
    memo = {}

    def Hx(mask):
        if mask in memo:
            return memo[mask]
        cyc = [c for c in cycles if c & ~mask == 0]

        def rec(i, used):
            if i == len(cyc):
                return 1
            r = rec(i + 1, used)
            if not cyc[i] & used:
                r += x * rec(i + 1, used | cyc[i])
            return r

        memo[mask] = rec(0, 0)
        return memo[mask]

    return Hx


def paths_from(k, adj, start, allowed):
    out = []

    def dfs(v, mask, ln):
        out.append((mask, v, ln))
        for t in range(k):
            if allowed >> t & 1 and not mask >> t & 1 and adj[v][t]:
                dfs(t, mask | 1 << t, ln + 1)

    if allowed >> start & 1:
        dfs(start, 1 << start, 1)
    return out


def all_paths(k, adj, allowed):
    out = []
    for s in range(k):
        if allowed >> s & 1:
            out.extend(paths_from(k, adj, s, allowed))
    return out


def block_matrix(k, adj, w, up, x, Hx):
    full = (1 << k) - 1
    eta = Hx(full)

    def q_counts(R):
        """(q_even, q_odd): paths inside R ending at w, weighted by Hx(R - path)."""
        qe = qo = 0
        for (m, last, ln) in all_paths(k, adj, R):
            if last == w:
                wt = Hx(R & ~m)
                if ln % 2:
                    qo += wt
                else:
                    qe += wt
        return qe, qo

    qe, qo = q_counts(full)
    # rows: h', E', O' ; columns: h, E, O
    row_h = [eta, x * qo, x * qe]
    row = {0: [0, 0, 0], 1: [0, 0, 0]}  # parity 0 -> E', 1 -> O'
    for (m, last, ln) in paths_from(k, adj, up, full):
        par = ln % 2
        if not m >> w & 1:
            R = full & ~m
            re_, ro_ = q_counts(R)
            row[par][0] += Hx(R)
            row[par][1] += x * ro_   # T-path even + Q odd -> new odd cycle
            row[par][2] += x * re_
        else:
            row[par][0] += Hx(full & ~m)
            if last == w:
                rest = full & ~m
                conts = [(0, 0)] + [(m2, ln2) for (m2, l2, ln2) in all_paths(k, adj, rest)]
                for (m2, ln2) in conts:
                    wt = Hx(rest & ~m2)
                    pb = (ln + ln2) % 2
                    # T-part even keeps parity pb; T-part odd flips it
                    row[pb][1] += wt
                    row[1 - pb][2] += wt
    return tuple(row_h + row[0] + row[1])


def main():
    x = int(sys.argv[1])
    kmax = int(sys.argv[2])
    outfile = sys.argv[3]
    mats = set()
    reps = {}
    for k in range(1, kmax + 1):
        for adj in tournaments(k):
            Hx = make_Hx(k, adj, x)
            for w in range(k):
                for up in range(k):
                    M = block_matrix(k, adj, w, up, x, Hx)
                    mats.add(M)
                    if M not in reps:
                        reps[M] = (k, "".join(str(adj[i][j]) for i in range(k) for j in range(k)), w, up)
        print("order", k, "distinct matrices so far", len(mats), file=sys.stderr)
    with open(outfile, "w") as f:
        for M in sorted(mats):
            assert all(isinstance(v, int) and v >= 0 for v in M)
            f.write(" ".join(map(str, M)) + "\n")
    with open(outfile + ".reps", "w") as f:
        for M in sorted(mats):
            k, rows, w, up = reps[M]
            f.write(" ".join(map(str, M)) + f" | {k} {rows} {w} {up}\n")
    print(len(mats), "matrices written to", outfile)


if __name__ == "__main__":
    main()
