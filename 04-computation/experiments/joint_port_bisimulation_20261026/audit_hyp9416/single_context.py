# Audit of HYP-9416 (1): is a SINGLE two-sided context already injective?  For each of the 9 two-sided contexts
# k*h*k' (k, k' rooted of order <= 2), count distinct free classes over strong states and over all states, n = 3..7.
# Expected by the reconstruction lemma (audit report): for strong C with n >= 3 the free class of e*h*e determines
# [C; s, t], because in e*h*e the new entrance vertex a is the out-degree-1 vertex whose out-neighbour is not of
# out-degree >= 2 ... (see report), and the new exit vertex b symmetrically.  Usage: single_context.py
import sys, os, subprocess
from collections import defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from tcanon import canon, rooted_key, free_key, all_tournaments, is_strong, relabel

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

e = ((0,), 0, 0); f = ((2, 0), 0, 1); g = ((0, 1), 0, 1)       # singleton; s -> t; t -> s
ctx = {'e': e, 'A2(s->t)': f, 'A2(t->s)': g}
for n in range(3, 8):
    rep = {}
    for C in classes_of(n):
        for s in range(n):
            for t in range(n):
                if s != t: rep.setdefault(rooted_key(C, s, t), (C, s, t))
    states = list(rep); S = [x for x in states if is_strong(rep[x][0])]
    line = []
    for a, k in ctx.items():
        for b, k2 in ctx.items():
            fc = {x: free_key(glue(glue(k, rep[x]), k2)[0]) for x in states}
            line.append(f"{a}*h*{b}: {len({fc[x] for x in S})}/{len(S)} strong, {len(set(fc.values()))}/{len(states)} all")
    print(f"n = {n}: " + "; ".join(line), flush=True)
