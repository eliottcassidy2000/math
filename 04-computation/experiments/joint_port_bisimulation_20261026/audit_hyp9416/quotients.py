# Independent audit of HYP-9416, checks B and C: fresh-re-marking quotients of rooted two-port states.
#   state h = [C; s, t];  V1: successors [C^op xor E(P); s, t] over Hamilton paths P of C from s to t (labelled paths);
#   V2: successors [C^op xor E(P); P_0, P_last] over all Hamilton paths of C.
#   quotient: coarsest partition refining an observation O in which equivalent states have equal multisets of
#   successor cells (also a normalised-distribution variant, "dist", which compares successor-cell frequencies and
#   gives states with no successor a separate dead mark).
# Observations: free; sand = H(e*h*e); glue1 = free classes of h*k and k*h over the 11 rooted k of order <= 3;
# glue2 = free classes of k*h*k' over the 9 pairs of rooted k, k' of order <= 2.  Adversarial extras:
# glue_e = (free class of h*e, free class of e*h) (one context only); hst = number of Hamilton paths from s to t;
# side = (entrance-rooted class, exit-rooted class).
# Imports only the audit's own tcanon.py.  Usage: quotients.py n
import sys, time, json, itertools, subprocess
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from tcanon import canon, rooted_key, onesided_key, free_key, all_tournaments, is_strong, ham_paths, partner, relabel

n = int(sys.argv[1]); t0 = time.time()

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
    cs = [canon(o) for o in reps]
    assert len({c[0] for c in cs}) == len(reps) and sum(factorial(n) // c[1] for c in cs) == 2 ** (n * (n - 1) // 2)
    return reps

def glue(h, k):
    """hinge product h*k: every cross arc from k to h except exit(h) -> entrance(k); ports (entrance h, exit k)"""
    (oh, sh, th), (ok, sk, tk) = h, k
    n1, n2 = len(oh), len(ok); out = list(oh) + [o << n1 for o in ok]
    for u in range(n1):
        for w in range(n2):
            if u == th and w == sk: out[u] |= 1 << (n1 + w)
            else: out[n1 + w] |= 1 << u
    return (tuple(out), sh, n1 + tk)

def rooted_contexts(maxn):
    ctx = [((0,), 0, 0)]
    for m in range(2, maxn + 1):
        seen = {}
        for out in classes_of(m):
            for s in range(m):
                for t in range(m):
                    if s != t: seen.setdefault(rooted_key(out, s, t), (out, s, t))
        ctx += [seen[k] for k in sorted(seen)]
    return ctx

def hp_batch(outs):
    """number of Hamilton paths, by a numpy subset DP vectorised over a batch of tournaments of equal order"""
    res = []
    for b0 in range(0, len(outs), 1500):
        chunk = outs[b0:b0 + 1500]; B = len(chunk); N = len(chunk[0])
        A = np.zeros((B, N, N), dtype=np.int64)
        for bi, out in enumerate(chunk):
            for u in range(N):
                for w in range(N):
                    if (out[u] >> w) & 1: A[bi, u, w] = 1
        dp = np.zeros((1 << N, B, N), dtype=np.int64)
        for v in range(N): dp[1 << v, :, v] = 1
        for mask in range(1, 1 << N):
            cur = dp[mask]
            if not cur.any(): continue
            for w in range(N):
                if (mask >> w) & 1: continue
                dp[mask | (1 << w), :, w] += np.einsum('bv,bv->b', cur, A[:, :, w])
        res += [int(x) for x in dp[(1 << N) - 1].sum(axis=1)]
    return res

classes = classes_of(n)
state_rep = {}; succ1 = {}; succ2 = {}; hst = {}
npaths = 0
for C in classes:
    P = ham_paths(C); npaths += len(P)
    tg = [rooted_key(partner(C, p), p[0], p[-1]) for p in P]
    allsucc = sorted(tg)
    by = defaultdict(list)
    for p, x in zip(P, tg): by[(p[0], p[-1])].append(x)
    for s in range(n):
        for t in range(n):
            if s == t: continue
            key = rooted_key(C, s, t); lst = sorted(by.get((s, t), []))
            if key in state_rep:
                assert succ1[key] == lst and succ2[key] == allsucc     # well defined on rooted classes
            else:
                state_rep[key] = (C, s, t); succ1[key] = lst; succ2[key] = allsucc; hst[key] = len(lst)
states = sorted(state_rep)
assert all(x in state_rep for x in states for y in (succ1[x] + succ2[x]) for x in [y])
strong = {x: is_strong(state_rep[x][0]) for x in states}
S = [x for x in states if strong[x]]
t1 = time.time()
print(f"n = {n}: {len(classes)} classes, {len(states)} rooted two-port states ({len(S)} strong); "
      f"{npaths} labelled Hamilton paths over class representatives; V1 transitions {sum(map(len, succ1.values()))}, "
      f"V2 transitions {sum(map(len, succ2.values()))}  [{t1 - t0:.0f}s]", flush=True)

free = {x: free_key(state_rep[x][0]) for x in states}
side = {x: (onesided_key(state_rep[x][0], state_rep[x][1]), onesided_key(state_rep[x][0], state_rep[x][2])) for x in states}
ctx3 = rooted_contexts(3); ctx2 = [c for c in ctx3 if len(c[0]) <= 2]
assert len(ctx3) == 11 and len(ctx2) == 3
glue1 = {}; glue_e = {}; glue2 = {}
for x in states:
    h = state_rep[x]
    right = tuple(free_key(glue(h, k)[0]) for k in ctx3)
    left = tuple(free_key(glue(k, h)[0]) for k in ctx3)
    glue1[x] = (right, left)
    glue_e[x] = (right[0], left[0])
    glue2[x] = tuple(free_key(glue(glue(k, h), k2)[0]) for k in ctx2 for k2 in ctx2)
e = ctx3[0]
sandH = hp_batch([glue(glue(e, state_rep[x]), e)[0] for x in states])
sand = dict(zip(states, sandH))
t2 = time.time()
print(f"  observations computed  [{t2 - t1:.0f}s]", flush=True)

def norm(lab):
    keys = {}; return {x: keys.setdefault(lab[x], len(keys)) for x in states}

def quotient(succ, init, mode='mult'):
    cell = norm(init); rounds = 0
    while True:
        sig = {}
        for x in states:
            c = Counter(cell[y] for y in succ[x])
            if mode == 'mult': sig[x] = (cell[x], tuple(sorted(c.items())))
            else:
                tot = len(succ[x])
                sig[x] = (cell[x], 'dead') if tot == 0 else (cell[x], tuple(sorted((k, Fraction(m, tot)) for k, m in c.items())))
        new = norm(sig); rounds += 1
        if len(set(new.values())) == len(set(cell.values())): return cell, rounds
        cell = new

def counts(cell):
    return len(set(cell.values())), len({cell[x] for x in S})

obs = {'free': free, 'sand': sand, 'glue1': glue1, 'glue_e': glue_e, 'side': side, 'hst': hst,
       'free+hst': {x: (free[x], hst[x]) for x in states}}
res = {'n': n, 'states': len(states), 'strong': len(S), 'static': {}, 'V1': {}, 'V2': {}, 'V1dist': {}}
for name in ('free', 'sand', 'glue1', 'glue_e', 'side', 'hst', 'free+hst'):
    res['static'][name] = counts(norm(obs[name]))
res['static']['glue2'] = (len(set(glue2.values())), len({glue2[x] for x in S}))
cells = {}
for dyn, succ, mode in (('V1', succ1, 'mult'), ('V2', succ2, 'mult'), ('V1dist', succ1, 'dist')):
    for name in ('free', 'sand', 'glue1', 'glue_e', 'side'):
        if dyn == 'V2' and name in ('glue_e', 'side'): continue
        cell, rounds = quotient(succ, obs[name], mode)
        cells[(dyn, name)] = cell
        res[dyn][name] = counts(cell) + (rounds,)
print("  static observation: cells (all states, strong states)")
for k, v in res['static'].items(): print(f"    {k:9s} {v}")
for dyn in ('V1', 'V2', 'V1dist'):
    print(f"  {dyn} quotient: (cells all, cells strong, refinement rounds)")
    for k, v in res[dyn].items(): print(f"    {k:9s} {v}")

# one-sided collisions among strong states
fib = defaultdict(list)
for x in S: fib[side[x]].append(x)
coll = [v for v in fib.values() if len(v) > 1]
pairs = [(a, b) for v in coll for a, b in itertools.combinations(v, 2)]
def sep(f): return sum(1 for a, b in pairs if f(a) != f(b))
report = {
    'fibres': len(coll), 'pairs': len(pairs), 'fibre_sizes': sorted(Counter(len(v) for v in coll).items()),
    'glue1_static': sep(lambda x: glue1[x]), 'glue_e_static': sep(lambda x: glue_e[x]),
    'free_static': sep(lambda x: free[x]),
    'glue2_static': sep(lambda x: glue2[x]), 'sand_static': sep(lambda x: sand[x]), 'hst_static': sep(lambda x: hst[x]),
    'sand_or_hst_static': sep(lambda x: (sand[x], hst[x])),
    'V1_glue1': sep(lambda x: cells[('V1', 'glue1')][x]), 'V1_free': sep(lambda x: cells[('V1', 'free')][x]),
    'V1_glue_e': sep(lambda x: cells[('V1', 'glue_e')][x]), 'V1_side': sep(lambda x: cells[('V1', 'side')][x]),
    'V1dist_glue1': sep(lambda x: cells[('V1dist', 'glue1')][x]), 'V1dist_side': sep(lambda x: cells[('V1dist', 'side')][x]),
    'V1dist_free': sep(lambda x: cells[('V1dist', 'free')][x]),
    'V2_glue1': sep(lambda x: cells[('V2', 'glue1')][x]),
}
res['collisions'] = report
print("  one-sided collisions among strong states:", report)
# glue1 is a function of the one-sided classes (check C, computational form): h*k depends on the exit-rooted class
# only, k*h on the entrance-rooted class only
f_right = {}; f_left = {}; ok = True
for x in states:
    ok &= f_right.setdefault(side[x][1], glue1[x][0]) == glue1[x][0]
    ok &= f_left.setdefault(side[x][0], glue1[x][1]) == glue1[x][1]
res['glue1_onesided'] = ok
print("  check C: h*k depends only on the exit-rooted class and k*h only on the entrance-rooted class, all states:", ok)
# successors of strong states: how many are nonstrong; states with no s->t path
res['strong_dead'] = sum(1 for x in S if hst[x] == 0)
res['strong_succ_nonstrong'] = sum(1 for x in S for y in succ1[x] if not strong[y])
print(f"  strong states with no s->t Hamilton path: {res['strong_dead']}; V1 arcs from strong to nonstrong states: "
      f"{res['strong_succ_nonstrong']}")
json.dump(res, open(f'quotients_n{n}.json', 'w'), indent=1, default=str)
print(f"  [{time.time() - t0:.0f}s]")
