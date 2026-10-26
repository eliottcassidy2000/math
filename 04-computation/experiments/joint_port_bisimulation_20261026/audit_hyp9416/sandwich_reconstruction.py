# Constructive check of the audit's reconstruction lemma (which PROVES HYP-9416 (1) for every n >= 3):
#   for a strong tournament C on n >= 3 vertices and s != t, the unlabelled sandwich G = e*h*e determines [C; s, t].
#   In G: d+(a) = 1, d-(b) = 1; d+(s) = d+_C(s), d-(t) = d-_C(t); every other C-vertex has d+, d- >= 2 (strongness
#   gives d+_C, d-_C >= 1).  Out-degree-1 vertices are a and possibly s (a -> s); a is the one whose out-neighbour has
#   out-degree 1, or the only one.  Symmetrically for b.  Then s = N+(a), t = N-(b), C = G - {a, b}.
# The script applies the recovery rule to a randomly relabelled G for every strong state n = 3..7 and checks that the
# recovered rooted class equals the original one.  Imports only the audit's tcanon.py.
import sys, os, random, subprocess
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, rooted_key, all_tournaments, is_strong, relabel

def glue(h, k):
    (oh, sh, th), (ok, sk, tk) = h, k
    n1, n2 = len(oh), len(ok); out = list(oh) + [o << n1 for o in ok]
    for u in range(n1):
        for w in range(n2):
            if u == th and w == sk: out[u] |= 1 << (n1 + w)
            else: out[n1 + w] |= 1 << u
    return (tuple(out), sh, n1 + tk)
def classes_of(n):
    if n <= 6:
        seen = {}
        for out in all_tournaments(n):
            c, aut, order = canon(out)
            if c not in seen: seen[c] = relabel(out, order)
        return list(seen.values())
    reps = []
    for s in subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split():
        out = [0] * n; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if s[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        reps.append(tuple(out))
    return reps

def recover(G):
    N = len(G)
    indeg = [sum((G[u] >> v) & 1 for u in range(N)) for v in range(N)]
    outdeg = [G[v].bit_count() for v in range(N)]
    A = [v for v in range(N) if outdeg[v] == 1]
    B = [v for v in range(N) if indeg[v] == 1]
    def nout(v): return (G[v] & -G[v]).bit_length() - 1
    def nin(v): return next(u for u in range(N) if (G[u] >> v) & 1)
    a = A[0] if len(A) == 1 else next(v for v in A if outdeg[nout(v)] == 1)
    b = B[0] if len(B) == 1 else next(v for v in B if indeg[nin(v)] == 1)
    assert len(A) <= 2 and len(B) <= 2
    s, t = nout(a), nin(b)
    keep = [v for v in range(N) if v not in (a, b)]
    pos = {v: i for i, v in enumerate(keep)}
    C = tuple(sum(1 << pos[w] for w in keep if (G[v] >> w) & 1) for v in keep)
    return C, pos[s], pos[t]

rnd = random.Random(9416); e = ((0,), 0, 0)
for n in range(3, 8):
    ok = 0; tot = 0
    for C in classes_of(n):
        if not is_strong(C): continue
        for s in range(n):
            for t in range(n):
                if s == t: continue
                G = glue(glue(e, (C, s, t)), e)[0]
                perm = list(range(n + 2)); rnd.shuffle(perm)
                C2, s2, t2 = recover(relabel(G, perm))
                tot += 1; ok += rooted_key(C2, s2, t2) == rooted_key(C, s, t)
    print(f"n = {n}: recovered the rooted class from a randomly relabelled e*h*e in {ok} of {tot} labelled strong (C, s, t)")
    assert ok == tot
