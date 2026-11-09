#!/usr/bin/env python3
"""Independent check of the port-chain transfer (THM-5070).

Builds random port chains from the representative blocks written by chain_blocks.py (file MATS.reps),
predicts (h, E, O) by the matrices, and compares h with an independent packing count of the actual
tournament (packing_count.c: Hamiltonian-cycle DP per vertex subset + set-partition DP), for x = 2 and 3.
Also checks OCF (direct Hamiltonian path count = I(Omega, 2)) on every chain.

Usage: chain_check.py REPS2 REPS3 PACKING_COUNT_BINARY TRIALS MAXN SEED
"""
import random
import subprocess
import sys


def load(path):
    blocks = []
    for line in open(path):
        mat, rep = line.split("|")
        M = list(map(int, mat.split()))
        k, rows, w, up = rep.split()
        blocks.append((M, int(k), rows, int(w), int(up)))
    return blocks


def build(word, blocks):
    """Return adjacency rows of the chain and the predicted state (h, E, O)."""
    adj = [[0]]
    port = 0
    state = (1, 0, 1)
    for idx in word:
        M, k, rows, w, up = blocks[idx]
        base = len(adj)
        for r in adj:
            r.extend([1] * k)              # old vertices beat the new block
        for i in range(k):
            adj.append([0] * base + [int(rows[i * k + j]) for j in range(k)])
        adj[port][base + w] = 0             # reverse the port-w arc
        adj[base + w][port] = 1
        port = base + up
        h, E, O = state
        state = (M[0] * h + M[1] * E + M[2] * O, M[3] * h + M[4] * E + M[5] * O, M[6] * h + M[7] * E + M[8] * O)
    return adj, state


def main():
    reps = {2: load(sys.argv[1]), 3: load(sys.argv[2])}
    prog, trials, maxn, seed = sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    rng = random.Random(seed)
    lines, preds = [], []
    for t in range(trials):
        x = 2 if t % 2 == 0 else 3
        blocks = reps[x]
        word, n = [], 1
        while True:
            idx = rng.randrange(len(blocks))
            if n + blocks[idx][1] > maxn:
                break
            word.append(idx)
            n += blocks[idx][1]
        adj, state = build(word, blocks)
        lines.append(f"{len(adj)} " + "".join(str(b) for r in adj for b in r))
        preds.append((x, state[0], len(adj)))
    out = subprocess.run([prog], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout.split("\n")
    fails = 0
    for (x, hpred, n), line in zip(preds, out):
        ham, i2, i3 = map(int, line.split())
        ok = (ham == i2) and ((i2 if x == 2 else i3) == hpred)
        if not ok:
            fails += 1
            if fails <= 5:
                print("FAIL", x, n, hpred, ham, i2, i3)
    print(f"{trials} random chains (orders <= {maxn}): {fails} failures")


if __name__ == "__main__":
    main()
