"""a4_misc.py (AUDIT) - small checks of section 4:
 (a) Pi_2(C(B, O; TT_t)) = 2 [u 3^t + r2 2^t + s] with constants depending on the pair only (fit on t = 1, 2, 3,
     test on t = 4, 5, 6); u = 0 exactly when O is empty;
 (b) the 'fans': H(C(B, O; X => TT_t)) = r H(X) + M 2^t Pi_2(X) on random pairs and modules."""
import random, sys, os, itertools
from fractions import Fraction
from aud_common import *
out = []
def P(s):
    print(s); out.append(s); sys.stdout.flush()
rng = random.Random(99)
def closing(B, O, X):
    nb, nx = len(B), len(X)
    N = nb + nx + 1
    o = [set() for _ in range(N)]
    v = nb + nx
    for a in range(nb):
        o[a] |= set(B[a]); o[a] |= set(range(nb, nb + nx))
        if a in O: o[v].add(a)
        else: o[a].add(v)
    for a in range(nx):
        o[nb + a] |= {nb + u for u in X[a]}; o[nb + a].add(v)
    return o
def join(X, Y):
    n = len(X)
    return [set(s) | set(range(n, n + len(Y))) for s in X] + [set(n + u for u in s) for s in Y]
def pair_state(B, O):
    n = len(B); h = W = U = 0
    for p in itertools.permutations(range(n)):
        arcs = [p[i + 1] in B[p[i]] for i in range(n - 1)]
        if all(arcs): h += 1
        for c in range(0, n):
            if all(arcs[:max(c - 1, 0)]) and all(arcs[c:]):
                if p[c] in O:
                    U += 1
                    if c >= 1 and p[c - 1] in O: W += 1
    return h - W, U
bad = 0; cases = 0; uzero = 0; emptyO = 0
for trial in range(25):
    nb = rng.randint(1, 4)
    B = random_tournament(nb, rng)
    O = set(v for v in range(nb) if rng.random() < 0.5) if trial > 3 else set()
    vals = []
    for t in range(1, 7):
        T = closing(B, O, transitive(t))
        pk = P_tuple(T, 2)
        vals.append((pk[2] + 2 * pk[1]) // 2)
    # solve u 3^t + r 2^t + s = vals for t = 1, 2, 3
    import fractions
    M = [[Fraction(3 ** t), Fraction(2 ** t), Fraction(1), Fraction(vals[t - 1])] for t in (1, 2, 3)]
    for i in range(3):
        piv = next(r for r in range(i, 3) if M[r][i] != 0)
        M[i], M[piv] = M[piv], M[i]
        M[i] = [x / M[i][i] for x in M[i]]
        for r in range(3):
            if r != i:
                M[r] = [a - M[r][i] * b for a, b in zip(M[r], M[i])]
    u, r2, s = M[0][3], M[1][3], M[2][3]
    ok = all(u * 3 ** t + r2 * 2 ** t + s == vals[t - 1] for t in (4, 5, 6)) and u.denominator == 1 and r2.denominator == 1 and s.denominator == 1
    cases += 1; bad += not ok
    if not O: emptyO += 1
    if u == 0: uzero += 1
    if (u == 0) != (not O): P("   unexpected: u = %s with O = %s" % (u, O))
P("(a) Pi_2(C(B,O;TT_t))/2 = u 3^t + r2 2^t + s (integers, fitted on t <= 3, tested on t = 4, 5, 6): pairs %d, failures %d; u = 0 in %d cases = the %d cases with O empty" % (cases, bad, uzero, emptyO))
jobs = []
for trial in range(60):
    nb = rng.randint(1, 4)
    B = random_tournament(nb, rng)
    O = set(v for v in range(nb) if rng.random() < 0.6)
    r, M_ = pair_state(B, O)
    X = random_tournament(rng.randint(1, 4), rng)
    t = rng.randint(0, 3)
    pk = P_tuple(X, 2)
    Xt = join(X, transitive(t)) if t else X
    jobs.append((closing(B, O, Xt), r * pk[1] + M_ * 2 ** t * (pk[2] + 2 * pk[1])))
Hs = ham_c([j[0] for j in jobs])
P("(b) fans: H(C(B,O; X => TT_t)) = r H(X) + M 2^t Pi_2(X): cases %d, failures %d" % (len(jobs), sum(1 for (T, pred), h in zip(jobs, Hs) if h != pred)))
open(os.path.join(HERE, "a4_misc.out"), "w").write("\n".join(out) + "\n")
