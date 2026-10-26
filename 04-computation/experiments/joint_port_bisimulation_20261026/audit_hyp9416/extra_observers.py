# Audit of the note's section 1 "Answer" bullets: V1 closure for weaker attachment observers.
#   Hx  = H(h*e) (scalar, exit side), Hn = H(e*h) (scalar, entrance side), Hxn = (H(h*e), H(e*h)),
#   Fx  = free class of h*e only, Fn = free class of e*h only, and sand = H(e*h*e) (two-sided scalar, for reference).
# Reports V1 and V2 quotient cells on strong states.  Imports only the audit's own modules.  Usage: extra_observers.py n
import sys, os, subprocess
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, rooted_key, free_key, all_tournaments, is_strong, ham_paths, partner, relabel

n = int(sys.argv[1])
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
def hcount(out):
    m = len(out); dp = [[0] * m for _ in range(1 << m)]
    for v in range(m): dp[1 << v][v] = 1
    for mask in range(1, 1 << m):
        row = dp[mask]
        for v in range(m):
            c = row[v]
            if not c: continue
            nxt = out[v] & ~mask
            while nxt:
                w = nxt & -nxt; nxt ^= w; dp[mask | w][w.bit_length() - 1] += c
    return sum(dp[(1 << m) - 1])

rep = {}; s1 = {}; s2 = {}
for C in classes_of(n):
    P = ham_paths(C); tg = [rooted_key(partner(C, p), p[0], p[-1]) for p in P]
    by = defaultdict(list)
    for p, x in zip(P, tg): by[(p[0], p[-1])].append(x)
    for s in range(n):
        for t in range(n):
            if s != t:
                x = rooted_key(C, s, t)
                if x not in rep: rep[x] = (C, s, t); s1[x] = by.get((s, t), []); s2[x] = tg
states = sorted(rep); S = [x for x in states if is_strong(rep[x][0])]
e = ((0,), 0, 0)
Hx = {x: hcount(glue(rep[x], e)[0]) for x in states}
Hn = {x: hcount(glue(e, rep[x])[0]) for x in states}
Fx = {x: free_key(glue(rep[x], e)[0]) for x in states}
Fn = {x: free_key(glue(e, rep[x])[0]) for x in states}
sand = {x: hcount(glue(glue(e, rep[x]), e)[0]) for x in states}
def norm(lab):
    keys = {}; return {x: keys.setdefault(lab[x], len(keys)) for x in states}
def quotient(succ, init):
    cell = norm(init)
    while True:
        new = norm({x: (cell[x], tuple(sorted(Counter(cell[y] for y in succ[x]).items()))) for x in states})
        if len(set(new.values())) == len(set(cell.values())): return cell
        cell = new
obs = {'Hx': Hx, 'Hn': Hn, 'Hxn': {x: (Hx[x], Hn[x]) for x in states}, 'Fx': Fx, 'Fn': Fn, 'sand': sand}
line = []
for name, o in obs.items():
    c1 = quotient(s1, o); c2 = quotient(s2, o)
    line.append(f"{name}: static {len({o[x] for x in S})}, V1 {len({c1[x] for x in S})}, V2 {len({c2[x] for x in S})}")
print(f"n = {n} ({len(S)} strong states), strong cells: " + "; ".join(line))
