#!/usr/bin/env python3
"""spine_grammar_validate_20261104.py -- opus-2026-11-04-S79

Independent check of the grammar engine: random words of G(k, w) are decoded into explicit tournaments (chain order,
backward arcs from the slots, apex arcs) and their Hamiltonian-path counts are computed by a subset dynamic programme
here (pure Python, exact), then compared with the engine's word mode.
Usage: python3 spine_grammar_validate_20261104.py <grammar.exe> k w trials maxlen [seed]
"""
import sys, random, subprocess

def letters_for(c, w, k):
    L = []
    for T in range(1 << c):
        t = bin(T).count('1')
        for b in range(0, w - (c - t) + 1):
            for A in range(1 << k):
                L.append((T, b, A))
    return L

def random_word(k, w, maxlen, rng):
    """returns (letter indices, arcs as (tail, head) over chain vertices 0..m-1, apex bits per vertex)"""
    pending = []          # list of (head, group) in creation order
    groups = 0
    word = []; arcs = []; apexbits = []
    m = rng.randint(1, maxlen)
    for v in range(m):
        c = len(pending)
        L = letters_for(c, w, k)
        last = (v == m - 1)
        cands = []
        for li, (T, b, A) in enumerate(L):
            # at most one slot per group
            gs = [pending[i][1] for i in range(c) if T >> i & 1]
            if len(gs) != len(set(gs)): continue
            if last and (c - bin(T).count('1') + b) != 0: continue
            cands.append(li)
        if not cands:
            return None
        li = rng.choice(cands); T, b, A = L[li]
        word.append(li)
        newpending = []
        for i in range(c):
            if T >> i & 1: arcs.append((v, pending[i][0]))
            else: newpending.append(pending[i])
        for j in range(b): newpending.append((v, groups))
        if b: groups += 1
        pending = newpending
        apexbits.append(A)
    if pending: return None
    return word, arcs, apexbits

def hamilton(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w2 in out[v]:
                if not (S >> w2 & 1): dp[S | (1 << w2)][w2] += c
    return sum(dp[full][v] for v in range(n))

def build(k, arcs, apexbits):
    m = len(apexbits); n = m + k
    adj = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m): adj[i][j] = 1
    for (tail, head) in arcs:
        adj[head][tail] = 0; adj[tail][head] = 1
    for i in range(m):
        for a in range(k):
            x = m + a
            if apexbits[i] >> a & 1: adj[x][i] = 1
            else: adj[i][x] = 1
    if k == 2: adj[m][m + 1] = 1
    # sanity: tournament
    for i in range(n):
        for j in range(i + 1, n):
            assert adj[i][j] + adj[j][i] == 1, (i, j)
    out = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    return n, out

def main():
    exe, k, w, trials, maxlen = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 20261104
    rng = random.Random(seed); done = 0; bad = 0; tried = 0
    while done < trials and tried < 50 * trials:
        tried += 1
        r = random_word(k, w, maxlen, rng)
        if r is None: continue
        word, arcs, apexbits = r
        n, out = build(k, arcs, apexbits)
        if n > 15: continue
        h = hamilton(n, out)
        res = subprocess.run([exe, str(k), str(w), '1000', '64', '1', '0', ','.join(map(str, word))], capture_output=True, text=True)
        line = [l for l in res.stdout.splitlines() if l.startswith('H = ')]
        eng = int(line[0][4:]) if line else None
        done += 1
        if eng != h:
            bad += 1
            if bad <= 5: print('MISMATCH word', word, 'arcs', arcs, 'apex', apexbits, 'engine', eng, 'brute', h)
    print(f'G({k},{w}): {done} random words (chain length <= {maxlen}, order <= 15) compared, {bad} mismatches')

main()
