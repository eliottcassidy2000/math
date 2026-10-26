# Audit of HYP-9416 (2): anatomy of the one-sided collision pairs and of the 289/253 validation.
#   * R3 targets C, D (rooted-lines README section 4): one-sided/joint classes, strongness, H(e*h*e), and the number of
#     Hamilton paths from s to t (the size of the V1 successor multiset).
#   * For every one-sided collision pair at n = 3..7: is it separated by hst = H(C; s -> t) alone, at which V1
#     refinement round of the bare `free` observation it splits, and whether the pair is a converse-type pair.
# Imports only the audit's own modules.  Usage: collision_anatomy.py
import sys, os, itertools, subprocess
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, rooted_key, onesided_key, free_key, all_tournaments, is_strong, ham_paths, partner, relabel

def glue(h, k):
    (oh, sh, th), (ok, sk, tk) = h, k
    n1, n2 = len(oh), len(ok); out = list(oh) + [o << n1 for o in ok]
    for u in range(n1):
        for w in range(n2):
            if u == th and w == sk: out[u] |= 1 << (n1 + w)
            else: out[n1 + w] |= 1 << u
    return (tuple(out), sh, n1 + tk)
e = ((0,), 0, 0)
C = (26, 12, 41, 48, 38, 3); D = (38, 20, 56, 19, 33, 10)
for name, X in (('C', C), ('D', D)):
    for u in range(6):
        for v in range(u + 1, 6): assert ((X[u] >> v) & 1) ^ ((X[v] >> u) & 1)
    sand = len(ham_paths(glue(glue(e, (X, 0, 5)), e)[0]))
    hst = sum(1 for p in ham_paths(X) if p[0] == 0 and p[-1] == 5)
    print(f"{name}: strong {is_strong(X)}, H(e*h*e) = {sand}, Hamilton paths 0 -> 5: {hst}")
print("entrance-rooted equal:", onesided_key(C, 0) == onesided_key(D, 0), "; exit-rooted equal:",
      onesided_key(C, 5) == onesided_key(D, 5), "; joint equal:", rooted_key(C, 0, 5) == rooted_key(D, 0, 5),
      "; free equal:", free_key(C) == free_key(D))

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

for n in range(3, 8):
    rep = {}; succ = {}; hst = {}
    for Cc in classes_of(n):
        P = ham_paths(Cc)
        by = defaultdict(list)
        for p in P: by[(p[0], p[-1])].append(rooted_key(partner(Cc, p), p[0], p[-1]))
        for s in range(n):
            for t in range(n):
                if s != t:
                    x = rooted_key(Cc, s, t)
                    if x not in rep: rep[x] = (Cc, s, t); succ[x] = by.get((s, t), []); hst[x] = len(succ[x])
    states = sorted(rep)
    free = {x: free_key(rep[x][0]) for x in states}
    # V1 refinement of the bare class observation, recording the round at which each pair splits
    def norm(lab):
        keys = {}; return {x: keys.setdefault(lab[x], len(keys)) for x in states}
    history = [norm(free)]
    while True:
        cell = history[-1]
        new = norm({x: (cell[x], tuple(sorted(Counter(cell[y] for y in succ[x]).items()))) for x in states})
        if len(set(new.values())) == len(set(cell.values())): break
        history.append(new)
    S = [x for x in states if is_strong(rep[x][0])]
    fib = defaultdict(list)
    for x in S:
        Cc, s, t = rep[x]; fib[(onesided_key(Cc, s), onesided_key(Cc, t))].append(x)
    pairs = [(a, b) for v in fib.values() if len(v) > 1 for a, b in itertools.combinations(v, 2)]
    split_round = Counter()
    for a, b in pairs:
        r = next((i for i, c in enumerate(history) if c[a] != c[b]), None)
        split_round[r] += 1
    hsep = sum(1 for a, b in pairs if hst[a] != hst[b])
    # which pairs are 'dead' (no s -> t path) on one or both sides
    dead = Counter((hst[a] == 0) + (hst[b] == 0) for a, b in pairs)
    print(f"n = {n}: {len(pairs)} one-sided collision pairs; separated by hst = #(s->t Hamilton paths) alone: {hsep}; "
          f"V1 refinement round of the bare class observation at which the pair splits: {dict(sorted(split_round.items(), key=str))}; "
          f"pairs with 0/1/2 dead members: {dict(sorted(dead.items()))}")
