"""
H2 redone HONESTLY: is -log measS7 truly "irreducibly aggregate" or secretly low-order?
mac-mini-2026-06-20-S7 (followup; first run overfit: 56 pairwise params on 36 rows -> fake R^2=1).

We now use OVERDETERMINED data (rows >> params) and a TRAIN/TEST split, and climb the
interaction order: 1-body (additive), 2-body (pairwise), 3-body. The "free energy is -log of a
sum, not a sum of local terms" claim (the Gibbs interpretation of HYP-2704 / mac-mini-S6
"irreducibly aggregate") predicts that NO finite-order log-linear model fits exactly, and that
TEST R^2 stays < 1 and does not jump to ~1 at any low order.

Universe offsets 0..W-1 (0 forced = the e=0 clock). Subset size k.
"""
from itertools import combinations
from fractions import Fraction as F
import math
import random


def measS7_exact(E):
    E = list(E)
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        ae = abs(e)
        for j in range(ae):
            for kk in range(7):
                val = F(7 * j + kk, 7 * ae)
                if F(0) < val < F(1):
                    bps.add(val)
    pts = sorted(bps)
    total = F(0)
    for i in range(len(pts) - 1):
        lo, hi = pts[i], pts[i + 1]
        mid = (lo + hi) / 2
        cols = set()
        for e in E:
            cols.add(int(7 * ((e * mid) % 1)))
        if len(cols) == 7:
            total += (hi - lo)
    return total


def lstsq_solve(A, y, ridge=1e-10):
    cols = len(A[0])
    AtA = [[sum(A[r][i] * A[r][j] for r in range(len(A))) for j in range(cols)] for i in range(cols)]
    Aty = [sum(A[r][i] * y[r] for r in range(len(A))) for i in range(cols)]
    for i in range(cols):
        AtA[i][i] += ridge
    M = [AtA[i][:] + [Aty[i]] for i in range(cols)]
    for c in range(cols):
        p = max(range(c, cols), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        if abs(pv) < 1e-14:
            continue
        for j in range(c, cols + 1):
            M[c][j] /= pv
        for r in range(cols):
            if r != c and abs(M[r][c]) > 1e-300:
                f = M[r][c]
                for j in range(c, cols + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][cols] for i in range(cols)]


def r2(A, y, beta):
    pred = [sum(A[r][i] * beta[i] for i in range(len(beta))) for r in range(len(A))]
    ybar = sum(y) / len(y)
    ss_res = sum((y[r] - pred[r]) ** 2 for r in range(len(A)))
    ss_tot = sum((yy - ybar) ** 2 for yy in y)
    return 1 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0, ss_res


def features(s, offs, order):
    f = [1.0]
    f += [1.0 if e in s else 0.0 for e in offs]
    if order >= 2:
        f += [1.0 if (a in s and b in s) else 0.0 for a, b in combinations(offs, 2)]
    if order >= 3:
        f += [1.0 if (a in s and b in s and c in s) else 0.0 for a, b, c in combinations(offs, 3)]
    return f


def run(W, k, seed=1):
    offs = list(range(W))
    subs = [s for s in combinations(range(W), k) if 0 in s]
    data = []
    for s in subs:
        mv = measS7_exact(s)
        if mv > 0:
            data.append((set(s), math.log(float(mv))))
    random.seed(seed)
    random.shuffle(data)
    ntr = int(0.7 * len(data))
    train, test = data[:ntr], data[ntr:]
    print(f"  W={W} k={k}: {len(data)} positive subsets; train={len(train)} test={len(test)}")
    for order, name in [(1, "1-body additive"), (2, "2-body pairwise"), (3, "3-body triple ")]:
        ncols = len(features(set([0]), offs, order))
        if ncols > 0.8 * len(train):
            print(f"    {name}: SKIP ({ncols} params vs {len(train)} train rows -- underdetermined)")
            continue
        Atr = [features(s, offs, order) for s, _ in train]
        ytr = [v for _, v in train]
        beta = lstsq_solve(Atr, ytr)
        r2tr, _ = r2(Atr, ytr, beta)
        Ate = [features(s, offs, order) for s, _ in test]
        yte = [v for _, v in test]
        r2te, _ = r2(Ate, yte, beta)
        print(f"    {name}: params={ncols:3d}  train R^2={r2tr:.5f}  TEST R^2={r2te:.5f}")
    return data


print("=" * 78)
print("HONEST interaction-order climb for -log measS7 (overdetermined, train/test)")
print("=" * 78)
print("If LRC free energy were a SUM of local terms, some finite order -> TEST R^2 = 1.")
print("Aggregate (HYP-2704/S6) predicts TEST R^2 stays < 1 and doesn't snap to 1.\n")
run(12, 8)
print()
run(13, 7)
print()
print("INTERPRETATION:")
print("  -log measS7 = -log( a single sum-over-x measure ) = a free energy. A free energy is")
print("  generically NOT a finite-order interaction polynomial in the spins (offsets). If TEST")
print("  R^2 never reaches 1, that IS the precise meaning of 'irreducibly aggregate': the")
print("  obstruction is the log-of-a-sum (Gibbs free energy), not any term-by-term inequality.")
