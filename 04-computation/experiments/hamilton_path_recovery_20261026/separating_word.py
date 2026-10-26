# Constructive one-guard separation (THM-4807 (d)): for labelled tournaments T != T' on one vertex set, build an
# explicit vertex word that is a Hamilton path of exactly one of them.
#   (alpha) some arc traced by a Hamilton path P of T is reversed in T'      -> P (legal in T only);
#   (beta)  otherwise some untraced arc u -> v of T inside a strong component B differs
#           -> concatenate Hamilton paths of T'[C] for the components C before/after B, and for B the word
#              P'u, where uP' is a Hamilton path of T[B] from u to v (chord lemma); legal in T' only;
#   (gamma) otherwise only arcs between non-consecutive components differ: take the differing arc u -> v
#           (u in C_i, v in C_j) with i minimal and, for that i, j maximal; the word
#           Q_1..Q_(i-1) Q_(i+1)..Q_(j-1) [Q_j ending at v] [Q_i starting at u] Q_(j+1)..Q_k is legal in T' only.
# Exhaustive over all ordered pairs of labelled tournaments for n <= 5; random pairs (biased to small or untraced
# differences) for n = 6..9.
import sys, random, itertools
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import ham_paths
from ham_recovery import components

def sub(adj, keep):
    idx = {v: i for i, v in enumerate(keep)}
    return tuple(sum(1 << idx[w] for w in keep if (adj[v] >> w) & 1) for v in keep)

def hp_in(adj, verts, start=None, end=None):
    """a Hamilton path of adj[verts] (optionally with prescribed start / end); None if none exists"""
    verts = list(verts)
    if len(verts) == 1:
        v = verts[0]
        return [v] if start in (None, v) and end in (None, v) else None
    s = sub(adj, verts)
    for P in ham_paths(s):
        Q = [verts[i] for i in P]
        if (start is None or Q[0] == start) and (end is None or Q[-1] == end): return Q
    return None

def legal(adj, w):
    return len(set(w)) == len(adj) == len(w) and all((adj[a] >> b) & 1 for a, b in zip(w, w[1:]))

def separating_word(T, T2):
    n = len(T)
    HP = ham_paths(T)
    for P in HP:                                                   # (alpha)
        if any(not ((T2[a] >> b) & 1) for a, b in zip(P, P[1:])): return list(P), 'alpha'
    traced = {(a, b) for P in HP for a, b in zip(P, P[1:])}
    comps = components(T); ci = {v: i for i, c in enumerate(comps) for v in c}
    diff = [(u, v) for u in range(n) for v in range(n) if (T[u] >> v) & 1 and not ((T2[u] >> v) & 1)]
    intra = [(u, v) for u, v in diff if ci[u] == ci[v]]
    if intra:                                                      # (beta)
        u, v = intra[0]; B = comps[ci[u]]
        Bm = [x for x in B if x != u]
        cB = components(sub(T, Bm)); cB = [[Bm[i] for i in c] for c in cB]
        c0 = next(x for x in cB[0] if (T[u] >> x) & 1)
        Pp = []
        for t, C in enumerate(cB):
            Pp += hp_in(T, C, start=c0 if t == 0 else None, end=v if t == len(cB) - 1 else None)
        word = []
        for t, C in enumerate(comps):
            word += (Pp + [u]) if t == ci[u] else hp_in(T2, C)
        return word, 'beta'
    i = min(ci[u] for u, v in diff)                                # (gamma)
    j = max(ci[v] for u, v in diff if ci[u] == i)
    u, v = next((u, v) for u, v in diff if ci[u] == i and ci[v] == j)
    Q = {t: hp_in(T, C) for t, C in enumerate(comps)}
    word = []
    for t in range(i): word += Q[t]
    for t in range(i + 1, j): word += Q[t]
    word += hp_in(T, comps[j], end=v) + hp_in(T, comps[i], start=u)
    for t in range(j + 1, len(comps)): word += Q[t]
    return word, 'gamma'

def all_labelled(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
        yield tuple(rows)

from collections import Counter
nmax_exh = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for n in range(2, nmax_exh + 1):
    L = list(all_labelled(n)); bad = 0; cases = Counter()
    for T in L:
        for T2 in L:
            if T2 == T: continue
            w, case = separating_word(T, T2); cases[case] += 1
            bad += legal(T, w) == legal(T2, w)
    print(f"labelled n = {n}: {len(L)*(len(L)-1)} ordered pairs, cases {dict(cases)}, failures {bad}", flush=True)

rng = random.Random(7)
for n in range(6, 10):
    bad = 0; cases = Counter(); trials = 3000 if n <= 7 else 800
    for _ in range(trials):
        rows = [0] * n
        # biased generator: a random transitive-ish backbone with a few random reversals, so that many arcs are untraced
        order = list(range(n)); rng.shuffle(order)
        for a in range(n):
            for b in range(a + 1, n):
                x, y = order[a], order[b]
                if rng.random() < 0.15: rows[y] |= 1 << x
                else: rows[x] |= 1 << y
        T = tuple(rows)
        HP = ham_paths(T); traced = {(a, b) for P in HP for a, b in zip(P, P[1:])}
        arcs = [(u, v) for u in range(n) for v in range(n) if (T[u] >> v) & 1]
        untr = [e for e in arcs if e not in traced]
        pool = untr if (untr and rng.random() < 0.8) else arcs
        D = rng.sample(pool, rng.randint(1, min(4, len(pool))))
        T2 = list(T)
        for u, v in D: T2[u] &= ~(1 << v); T2[v] |= 1 << u
        T2 = tuple(T2)
        w, case = separating_word(T, T2); cases[case] += 1
        bad += legal(T, w) == legal(T2, w)
    print(f"random n = {n}: {trials} pairs, cases {dict(cases)}, failures {bad}", flush=True)
