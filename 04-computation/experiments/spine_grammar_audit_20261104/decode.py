#!/usr/bin/env python3
"""decode.py -- auditor's own decoder of G(k,w) words from the letter semantics of the S79 note (section 2),
plus random-word comparison against the engine's word mode and the witness check of a census file.
Usage:
  python3 -I decode.py random k w N maxlen seed      (compare: my Python decode + hp batch DP, hp word mode, engine word mode)
  python3 -I decode.py witnesses k w censusfile       (decode every witness, recount H with hp batch, histogram of arcs)
"""
import sys, random, subprocess, re, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
HP = os.path.join(HERE, 'hp.exe')
ENGINE = 'C:/Users/Eliott/AppData/Local/Temp/s79/grammar4.exe'

def letters(c, w, k):
    L = []
    for T in range(1 << c):
        t = bin(T).count('1')
        for b in range(0, w - (c - t) + 1):
            for A in range(1 << k):
                L.append((T, b, A))
    return L

def decode(k, w, word, require_complete=True):
    """pending slots in creation order as (head, group); returns (arcs tail>head over chain vertices, apexbits) or None"""
    pending = []; groups = 0; arcs = []; apex = []
    for t, li in enumerate(word):
        c = len(pending); L = letters(c, w, k)
        if li < 0 or li >= len(L): return None
        T, b, A = L[li]
        gs = [pending[i][1] for i in range(c) if T >> i & 1]
        if len(gs) != len(set(gs)): return None
        newp = []
        for i in range(c):
            if T >> i & 1: arcs.append((t, pending[i][0]))
            else: newp.append(pending[i])
        for j in range(b): newp.append((t, groups))
        if b: groups += 1
        pending = newp; apex.append(A)
    if require_complete and pending: return None
    return arcs, apex

def adjacency(k, m, arcs, apex):
    """vertices: apexes 0..k-1 (0 -> 1), chain vertex t -> k + t; returns n and the n*n 0/1 string (row i col j = i -> j)"""
    n = m + k
    adj = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m): adj[k + i][k + j] = 1
    for (tail, head) in arcs:
        adj[k + head][k + tail] = 0; adj[k + tail][k + head] = 1
    for t in range(m):
        for a in range(k):
            if apex[t] >> a & 1: adj[a][k + t] = 1
            else: adj[k + t][a] = 1
    if k == 2: adj[0][1] = 1
    for i in range(n):
        for j in range(i + 1, n):
            assert adj[i][j] + adj[j][i] == 1
    return n, ''.join(str(adj[i][j]) for i in range(n) for j in range(n))

def batch_counts(items):
    """items: list of (n, adjstring); returns list of H via hp.exe batch"""
    if not items: return []
    inp = '\n'.join(f'{n} {s}' for n, s in items) + '\n'
    res = subprocess.run([HP, 'batch'], input=inp, capture_output=True, text=True)
    return [int(x) for x in res.stdout.split()]

def engine_word(k, w, word):
    res = subprocess.run([ENGINE, str(k), str(w), '1000', '64', '1', '0', ','.join(map(str, word))], capture_output=True, text=True)
    for l in res.stdout.splitlines():
        if l.startswith('H = '): return int(l[4:])
    return res.stdout.strip()[:40]

def hp_word(k, w, word):
    res = subprocess.run([HP, 'word', str(k), str(w), ','.join(map(str, word))], capture_output=True, text=True)
    for l in res.stdout.splitlines():
        if l.startswith('H = '): return int(l[4:])
    return res.stdout.strip()[:40]

def random_word(k, w, maxlen, rng):
    pending = []; groups = 0; word = []
    m = rng.randint(1, maxlen)
    for v in range(m):
        c = len(pending); L = letters(c, w, k); last = (v == m - 1)
        cands = []
        for li, (T, b, A) in enumerate(L):
            gs = [pending[i][1] for i in range(c) if T >> i & 1]
            if len(gs) != len(set(gs)): continue
            if last and c - bin(T).count('1') + b != 0: continue
            cands.append(li)
        if not cands: return None
        li = rng.choice(cands); T, b, A = L[li]; word.append(li)
        newp = [pending[i] for i in range(c) if not (T >> i & 1)]
        for j in range(b): newp.append((v, groups))
        if b: groups += 1
        pending = newp
    if pending: return None
    return word

def main():
    mode = sys.argv[1]
    if mode == 'random':
        k, w, N, maxlen, seed = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        rng = random.Random(seed); words = []
        while len(words) < N:
            wd = random_word(k, w, maxlen, rng)
            if wd is None: continue
            dec = decode(k, w, wd)
            if dec is None: continue
            arcs, apex = dec
            if len(wd) + k > 15: continue
            words.append((wd, arcs, apex))
        items = [adjacency(k, len(wd), arcs, apex) for (wd, arcs, apex) in words]
        mine = batch_counts(items)
        bad = 0
        for (wd, arcs, apex), h in zip(words, mine):
            e = engine_word(k, w, wd); h2 = hp_word(k, w, wd)
            ok = (e == h == h2)
            if not ok: bad += 1
            print(f"word {','.join(map(str, wd))} m={len(wd)} arcs={arcs} apex={''.join(map(str, apex))} myDP={h} hpWord={h2} engine={e} {'OK' if ok else 'MISMATCH'}")
        print(f"G({k},{w}): {len(words)} random words compared, {bad} mismatches")
    elif mode == 'witnesses':
        # H is multiplicative over the strong components (a Hamiltonian path visits the condensation in order), so each
        # witness is verified as the product of the subset-DP counts of its strong components (singletons contribute 1).
        k, w, path = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
        rows = []
        for line in open(path):
            m = re.match(r'\s+(\d+) = ([\d,]*)', line)
            if not m: continue
            N = int(m.group(1)); wd = [int(x) for x in m.group(2).split(',') if x != '']
            rows.append((N, wd))
        items = []; decs = []; comps_of = []
        for N, wd in rows:
            if not wd: decs.append(None); comps_of.append([]); continue
            dec = decode(k, w, wd)
            if dec is None: print('UNDECODABLE', N, wd); decs.append(None); comps_of.append([]); continue
            arcs, apex = dec; decs.append((arcs, apex))
            n, s = adjacency(k, len(wd), arcs, apex)
            adj = [[int(s[i * n + j]) for j in range(n)] for i in range(n)]
            comps = sccs(n, adj); idxs = []
            for comp in comps:
                if len(comp) == 1: continue
                if len(comp) > 20: idxs.append(None); continue
                sub = ''.join(str(adj[i][j]) for i in comp for j in comp)
                idxs.append(len(items)); items.append((len(comp), sub))
            comps_of.append(idxs)
        hs = batch_counts(items)
        bad = 0; unver = 0; byarcs = collections.Counter(); bym = collections.Counter(); maxm = 0; maxcomp = 0
        for (N, wd), dec, idxs in zip(rows, decs, comps_of):
            if dec is None:
                if wd: bad += 1
                continue
            arcs, apex = dec
            if any(i is None for i in idxs): unver += 1; print('UNVERIFIED (component > 20)', N); continue
            h = 1
            for i in idxs: h *= hs[i]
            if h != N: bad += 1; print('MISMATCH', N, wd, 'arcs', arcs, 'apex', apex, 'myH', h)
            byarcs[len(arcs)] += 1; bym[len(wd)] += 1; maxm = max(maxm, len(wd))
        print(f"{len(rows)} witnesses in {os.path.basename(path)}: {bad} mismatches, {unver} unverified; values by number of reversed arcs {sorted(byarcs.items())}; by chain length {sorted(bym.items())}; max chain length {maxm}")

def sccs(n, adj):
    """Tarjan; returns the strong components in condensation (topological) order"""
    index = [None] * n; low = [0] * n; onst = [False] * n; st = []; out = []; counter = [0]
    sys.setrecursionlimit(10000)
    def strong(v):
        index[v] = low[v] = counter[0]; counter[0] += 1; st.append(v); onst[v] = True
        for u in range(n):
            if not adj[v][u]: continue
            if index[u] is None: strong(u); low[v] = min(low[v], low[u])
            elif onst[u]: low[v] = min(low[v], index[u])
        if low[v] == index[v]:
            comp = []
            while True:
                u = st.pop(); onst[u] = False; comp.append(u)
                if u == v: break
            out.append(sorted(comp))
    for v in range(n):
        if index[v] is None: strong(v)
    return out

if __name__ == '__main__':
    main()
