#!/usr/bin/env python3
"""AUDIT: random port chains, transfer prediction versus brute force (THM-5070 (1)), all bases x = 2..8.

Blocks and their matrices come from a1_blocks.py (blocks_all.txt, every block class of order <= 5, with its
black-box matrix in each base). Each chain is built as an explicit tournament (old vertices beat new block vertices
except w -> u), and a1_brute (C) recomputes (h, E, O) for every base plus the direct Hamiltonian path count.
Checks: predicted (h, E, O) == brute force for x = 2..8, and OCF: H == h_2.
Usage: a1_chains.py BLOCKS_ALL BRUTE_EXE TRIALS NMAX KMAX SEED
"""
import random
import subprocess
import sys


def main():
    path, exe, trials, nmax, kmax, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), \
        int(sys.argv[5]), int(sys.argv[6])
    blocks = []
    for line in open(path):
        parts = line.split("|")
        k, rows, w, up = parts[0].split()
        mats = [list(map(int, p.split())) for p in parts[1:]]
        if int(k) <= kmax:
            blocks.append((int(k), rows, int(w), int(up), mats))
    rng = random.Random(seed)
    inputs, preds = [], []
    for t in range(trials):
        adj, port = [[0]], 0
        states = [(1, 0, 1)] * 7
        target = rng.randint(max(2, nmax - 6), nmax)
        while len(adj) < target:
            n = len(adj)
            fit = [b for b in blocks if b[0] <= target - n]
            k, rows, w, up, mats = fit[rng.randrange(len(fit))]
            for r in adj:
                r.extend([1] * k)
            for i in range(k):
                adj.append([0] * n + [int(rows[i * k + j]) for j in range(k)])
            adj[port][n + w] = 0
            adj[n + w][port] = 1
            port = n + up
            new = []
            for xi, (h, E, O) in enumerate(states):
                M = mats[xi]
                new.append((M[0] * h + M[1] * E + M[2] * O, M[3] * h + M[4] * E + M[5] * O,
                            M[6] * h + M[7] * E + M[8] * O))
            states = new
        n = len(adj)
        assert all(adj[i][j] + adj[j][i] == 1 for i in range(n) for j in range(n) if i != j)
        inputs.append(f"{n} {port} " + "".join(str(b) for r in adj for b in r))
        preds.append((n, states))
    out = subprocess.run([exe], input="\n".join(inputs) + "\n", capture_output=True, text=True).stdout.split("\n")
    fails = 0
    orders = {}
    for (n, states), line in zip(preds, out):
        vals = list(map(int, line.split()))
        H, rest = vals[0], vals[1:]
        brute = [tuple(rest[3 * i:3 * i + 3]) for i in range(7)]
        ok = brute == list(states) and H == brute[0][0]
        orders[n] = orders.get(n, 0) + 1
        if not ok:
            fails += 1
            if fails <= 5:
                print("FAIL", n, states, brute, H)
    print(f"{trials} random chains, blocks of order <= {kmax}, orders <= {nmax} (distribution {sorted(orders.items())}):"
          f" {fails} failures over (h, E, O) in bases 2..8 and OCF")


if __name__ == "__main__":
    main()
