"""Black native moves from regular tournaments to 'singleton + strong core'.
Native move tau_Q(T) = T^op xor E(Q): keep the arcs of the directed Hamiltonian path Q, reverse every other pair.
Q is BLUE in T iff its positional reversal q_i -> q_{n-1-i} is an anti-automorphism of T; BLACK otherwise.
Target: a vertex a of score 0 or n-1 whose deletion leaves a strong tournament.
Lower bound (proved): from a regular tournament, m=(n-1)/2 moves are necessary, since one move changes |s(v)-m| by at most 1.
Search: a must be an endpoint of every path (start if s(a)<=m, end if s(a)>=m), so the branching is over such black paths."""
import sys, random, time
sys.setrecursionlimit(10000)
def parse(line, n):
    out = [0] * n; p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[p] == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
            p += 1
    return out
def score(out, v): return bin(out[v]).count('1')
def is_strong(out, verts):
    verts = list(verts); k = len(verts); mask = sum(1 << v for v in verts)
    sc = sorted(bin(out[v] & mask).count('1') for v in verts)
    cum = 0
    for i in range(k - 1):
        cum += sc[i]
        if cum == (i + 1) * i // 2: return False
    return True
def move(out, n, q):
    full = (1 << n) - 1
    new = [full & ~(1 << i) & ~out[i] for i in range(n)]
    for i in range(n - 1):
        u, v = q[i], q[i + 1]; new[u] |= 1 << v; new[v] &= ~(1 << u)
    return new
def is_blue(out, n, q):
    for i in range(n):
        for j in range(i + 1, n):
            if ((out[q[i]] >> q[j]) & 1) != ((out[q[n - 1 - j]] >> q[n - 1 - i]) & 1): return False
    return True
def paths_from(out, n, a, rng=None):
    res = []; full = (1 << n) - 1
    def rec(v, used, path):
        if used == full: res.append(list(path)); return
        nx = out[v] & ~used; cand = [w for w in range(n) if (nx >> w) & 1]
        if rng: rng.shuffle(cand)
        for w in cand:
            path.append(w); rec(w, used | (1 << w), path); path.pop()
    rec(a, 1 << a, [a]); return res
def paths_to(out, n, a, rng=None):
    full = (1 << n) - 1; conv = [full & ~(1 << i) & ~out[i] for i in range(n)]
    return [p[::-1] for p in paths_from(conv, n, a, rng)]
def search(out, n, a, m, cap, rng, count_all=False):
    """returns (witness word or None, number of witnesses found if count_all)"""
    found = [None]; cnt = [0]
    def rec(T, depth, word):
        s = score(T, a); d = s - m
        if depth == m:
            if abs(d) == m and is_strong(T, [v for v in range(n) if v != a]):
                cnt[0] += 1
                if found[0] is None: found[0] = list(word)
                return not count_all
            return False
        if abs(d) != depth: return False
        cands = []
        if d <= 0: cands += paths_from(T, n, a, rng)
        if d >= 0 and (d > 0 or True): cands += (paths_to(T, n, a, rng) if d > 0 else paths_to(T, n, a, rng))
        if rng: rng.shuffle(cands)
        if cap: cands = cands[:cap]
        for q in cands:
            if is_blue(T, n, q): continue
            word.append(q)
            if rec(move(T, n, q), depth + 1, word): return True
            word.pop()
        return False
    rec(out, 0, []); return found[0], cnt[0]
def run(n, lines, cap, seed, count_all, label):
    m = (n - 1) // 2; rng = random.Random(seed)
    for idx, line in enumerate(lines):
        out = parse(line, n) if isinstance(line, str) else line
        if any(score(out, v) != m for v in range(n)): continue
        t0 = time.time(); best = None; tot = 0; vert = None
        for a in range(n):
            w, c = search(out, n, a, m, cap, rng, count_all)
            tot += c
            if w and best is None: best, vert = w, a
            if best and not count_all: break
        print(f"{label} n={n} class#{idx} regular: witness_len_{m}={'YES' if best else 'NO(cap)' if cap else 'NO'} a={vert} witnesses={tot if count_all else '-'} word={best}", flush=True)
        print(f"  ({label} class#{idx} {time.time()-t0:.1f}s)", file=sys.stderr)
if __name__ == '__main__':
    n = int(sys.argv[1]); mode = sys.argv[2]
    if mode == 'gentourng':
        cap = int(sys.argv[3]); count_all = sys.argv[4] == 'all'
        lines = [l.strip() for l in sys.stdin if len(l.strip()) == n * (n - 1) // 2]
        run(n, lines, cap, 1, count_all, 'gentourng')
    elif mode == 'paley':
        q = n; qr = {(x * x) % q for x in range(1, q)}
        out = [sum(1 << j for j in range(q) if j != i and ((j - i) % q) in qr) for i in range(q)]
        cap = int(sys.argv[3]); run(n, [out], cap, 2, False, 'paley')
    elif mode == 'randreg':
        q = n; qr = {(x * x) % q for x in range(1, q)}
        out = [sum(1 << j for j in range(q) if j != i and ((j - i) % q) in qr) for i in range(q)]
        rng = random.Random(5); samples = []
        for s in range(int(sys.argv[4])):
            for _ in range(200):  # random 3-cycle reversals keep regularity
                u, v, w = rng.sample(range(n), 3)
                if (out[u] >> v) & 1 and (out[v] >> w) & 1 and (out[w] >> u) & 1:
                    out[u] &= ~(1 << v); out[v] |= 1 << u; out[v] &= ~(1 << w); out[w] |= 1 << v; out[w] &= ~(1 << u); out[u] |= 1 << w
            samples.append(list(out))
        run(n, samples, int(sys.argv[3]), 3, False, 'randreg')
