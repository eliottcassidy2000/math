"""Auditor: independent re-check of HYP-9369 section 4 (decoration rule) with SESSION30's decorated engine
(session30-nonregular/decorated_engine.py: own core enumeration, own symbolic kernels, own contraction) instead of the
S54 lane's leaf_engine.py.  My own claw predicate; my own closed forms (Y(u), -u^T G u).
Tests on: almost-regular classes from the S54 files, random almost-regular-free tournaments (uniform random), the
transitive tournament.  Usage: python claw_check.py > claw_check.out"""
import os, sys, random, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
S30 = os.path.join(HERE, '..', '..', '..', '..', '05-knowledge', 'results', 'session30_nonregular_loops_20261022',
                   'session30-nonregular')
sys.path.insert(0, os.path.abspath(S30))
sys.path.insert(0, HERE)
from decorated_engine import cores, symbolic_kernel, evaluate      # Session30 engine
import tourn


def my_claw(ed):
    """centre c (no loop) with exactly three distinct non-loop neighbours; 4 vertices; every neighbour is a leaf
    (degree 1) or has exactly one loop and degree 3; no other edges."""
    V = 1 + max(max(a, b) for a, b in ed)
    if V != 4:
        return False
    deg = [0] * V; loops = [0] * V
    for a, b in ed:
        deg[a] += 1; deg[b] += 1
        if a == b:
            loops[a] += 1
    for c in range(V):
        nb = [b if a == c else a for a, b in ed if (a == c) != (b == c)]
        if loops[c] == 0 and len(nb) == 3 and len(set(nb)) == 3 and deg[c] == 3:
            if all(deg[x] == 1 or (deg[x] == 3 and loops[x] == 1) for x in nb):
                # every edge is a bridge c-x or a loop at an end
                if all((a == c) != (b == c) or (a == b and a != c) for a, b in ed):
                    return True
    return False


def two_arm(ed):
    """bridge skeleton: the path tree (2 leaves), the lollipop (1 leaf), the dumbbell (no leaf) -- here identified as:
    two 'ends' joined by one edge, each end a leaf or a vertex with exactly one loop."""
    V = 1 + max(max(a, b) for a, b in ed)
    if V != 2:
        return False
    deg = [0] * V; loops = [0] * V
    for a, b in ed:
        deg[a] += 1; deg[b] += 1
        if a == b:
            loops[a] += 1
    bridges = [(a, b) for a, b in ed if a != b]
    return len(bridges) == 1 and all(deg[x] == 1 or (deg[x] == 3 and loops[x] == 1) for x in range(V))


def mats(M):
    n = M.shape[0]; I = np.eye(n); M2 = M @ M
    G = -M2 @ np.linalg.inv(I + M2); H = M2 @ np.linalg.inv(I - M2)
    Qp = M @ np.linalg.inv(I + M2)
    d = np.diag(H - G)
    return G, H, Qp, d


CACHE = {}


def claw_cores():
    if 'claw' not in CACHE:
        out = {}
        for (e, l) in ((2, 0), (1, 1), (0, 2), (-1, 3)):
            sel = [(ed, a) for ed, a in cores(e, l) if my_claw(ed)]
            out[(e, l)] = [(ed, a, symbolic_kernel(ed, a)) for ed, a in sel]
        CACHE['claw'] = out
        out2 = {}
        for (e, l) in ((1, 0), (0, 1), (-1, 2)):
            sel = [(ed, a) for ed, a in cores(e, l) if two_arm(ed)]
            out2[(e, l)] = [(ed, a, symbolic_kernel(ed, a)) for ed, a in sel]
        CACHE['two'] = out2
    return CACHE['claw'], CACHE['two']


def check(name, A):
    n = A.shape[0]
    M = (A - A.T).astype(float) / n
    G, H, Qp, d = mats(M)
    one = np.ones(n); u = d + one
    claw, two = claw_cores()
    parts = {k: sum(evaluate(ed, sy, M) for ed, a, sy in v) for k, v in claw.items()}
    tot = sum(parts.values())
    Yu = -(2 / 3) * np.sum((G @ u) ** 3 + 3 * (G @ u) * (Qp @ u) ** 2)
    Yu_wrongQ = -(2 / 3) * np.sum((G @ u) ** 3 + 3 * (G @ u) * (M @ u) ** 2)      # Q read as Session30's Q = zM
    parts2 = {k: sum(evaluate(ed, sy, M) for ed, a, sy in v) for k, v in two.items()}
    sq = -u @ G @ u
    print(f'{name} (n={n}, ||S1||^2={int(np.sum(((A - A.T).sum(1)) ** 2))}):')
    print(f'   claw cores per class {", ".join(f"{k}: {len(claw[k])} core(s) {v:+.6e}" for k, v in parts.items())}')
    print(f'   claw sum {tot:+.12e}  Y(u) {Yu:+.12e}  diff {tot - Yu:+.1e}  | with Q := M instead: diff {tot - Yu_wrongQ:+.1e}')
    print(f'   two-arm cores {", ".join(f"{k}: {len(two[k])} {v:+.6e}" for k, v in parts2.items())}; '
          f'sum {sum(parts2.values()):+.12e}  -u^T G u {sq:+.12e}  diff {sum(parts2.values()) - sq:+.1e}')
    return tot - Yu, sum(parts2.values()) - sq


def main():
    claw, two = claw_cores()
    for k, v in claw.items():
        print(f'class {k}: claw cores selected by my predicate: {[(ed, a) for ed, a, _ in v]}')
    for k, v in two.items():
        print(f'class {k}: two-arm cores: {[(ed, a) for ed, a, _ in v]}')
    # also: does the S54 predicate is_claw select the same cores?  (re-implemented from decoration_rule.py's logic
    # on the S54 engine's core list; the S54 engine indexes leaves first)
    UP = os.path.join(HERE, '..')
    rng = random.Random(11)
    worst = [0.0, 0.0]
    tests = []
    for fn, n in (('ar8_classes.txt', 8), ('ar10_classes.txt', 10)):
        lines = open(os.path.join(UP, fn)).read().split('\n')
        for idx in (0, 7):
            tests.append((f'{fn}[{idx}]', tourn.from_code(lines[idx].split()[0])))
    for n in (7, 9, 12):
        code = ''.join(rng.choice('01') for _ in range(n * (n - 1) // 2))
        tests.append((f'uniform random {n}', tourn.from_code(code)))
    tests.append(('transitive 6', tourn.from_code('1' * 15)))
    tests.append(('C_3[K_2] x-ish: P_7 minus vertex', tourn.delete(tourn.circulant(7, [1, 2, 4]), 0)))
    for name, A in tests:
        a, b = check(name, A)
        worst = [max(worst[0], abs(a)), max(worst[1], abs(b))]
    print(f'max |claw - Y(u)| = {worst[0]:.2e}; max |two-arm - (-u^T G u)| = {worst[1]:.2e}')


if __name__ == '__main__':
    main()
