"""Exact cubic norm controls and a separately labelled finite open-bound scout.

The theorem covers common-edge cycle families and recursive substitutions
with constant cross-block amplitudes. General weighted tournaments are
only a declared finite scout, not a proved universal bound.
"""
from itertools import combinations, product
from math import comb
import random

GATES = 0


def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def data(n, orientation, weights):
    a = [[0] * n for _ in range(n)]
    for (i, j), bit, weight in zip(combinations(range(n), 2), orientation, weights):
        a[i if bit else j][j if bit else i] = weight
    return a


def energy(a):
    return sum(x*x for row in a for x in row)


def cubic(a):
    cycles = transitive = 0
    for i, j, k in combinations(range(len(a)), 3):
        c = a[i][j]*a[j][k]*a[k][i] + a[i][k]*a[k][j]*a[j][i]
        total = (a[i][j]+a[j][i])*(a[i][k]+a[k][i])*(a[j][k]+a[k][j])
        cycles += c
        transitive += total-c
    return 3*cycles-transitive, cycles


def matrix_cubic(a):
    n = len(a)
    twice = sum((a[i][j]+a[j][i])*(a[j][k]-a[k][j])*(a[k][i]-a[i][k])
                for i in range(n) for j in range(n) for k in range(n))
    check(twice % 2 == 0, "integer cubic trace")
    return twice//2


def substitute(blocks, q):
    offsets = [0]
    for block in blocks:
        offsets.append(offsets[-1]+len(block))
    a = [[0]*offsets[-1] for _ in range(offsets[-1])]
    for i, b in enumerate(blocks):
        for u, row in enumerate(b):
            for v, val in enumerate(row):
                a[offsets[i]+u][offsets[i]+v] = val
    for i, j in combinations(range(len(blocks)), 2):
        for u in range(offsets[i], offsets[i+1]):
            for v in range(offsets[j], offsets[j+1]):
                a[u][v], a[v][u] = q[i][j], q[j][i]
    return a


def audit_substitution(blocks, q):
    a = substitute(blocks, q)
    sizes = list(map(len, blocks))
    internal = sum(cubic(b)[0] for b in blocks)
    cross = cross_cycles = 0
    for i, j, k in combinations(range(len(blocks)), 3):
        c = q[i][j]*q[j][k]*q[k][i] + q[i][k]*q[k][j]*q[j][i]
        t = (q[i][j]+q[j][i])*(q[i][k]+q[k][i])*(q[j][k]+q[k][j])
        cross += sizes[i]*sizes[j]*sizes[k]*(4*c-t)
        cross_cycles += sizes[i]*sizes[j]*sizes[k]*c
    tax = sum(sum(map(sum,b))*sum(sizes[j]*(q[i][j]+q[j][i])**2
                                   for j in range(len(blocks)) if i!=j)
              for i,b in enumerate(blocks))
    check(cubic(a)[0] == internal+cross-tax, "substitution cubic identity")
    check(cubic(a)[1] == sum(cubic(b)[1] for b in blocks)+cross_cycles,
          "substitution cyclic identity")
    check(tax >= 0, "nonnegative boundary tax")
    check(energy(a) == sum(energy(b) for b in blocks)
          + sum(sizes[i]*sizes[j]*(q[i][j]+q[j][i])**2
                for i,j in combinations(range(len(blocks)),2)), "energy sectors")
    f = cubic(a)[0]
    check(f <= 0 or 3*f*f <= energy(a)**3, "proved grammar norm bound")
    check(27*cubic(a)[1]**2 <= energy(a)**3, "proved grammar cyclic norm bound")
    return a


def main():
    count = 0
    for n in range(2, 5):
        m = comb(n, 2)
        for orientation in product((0,1), repeat=m):
            for weights in product((0,1,2), repeat=m):
                a = data(n, orientation, weights)
                f,c = cubic(a)
                check(f == matrix_cubic(a), "matrix/triple independent equality")
                check(27*c*c <= energy(a)**3, "common-edge cyclic bound")
                check(f <= 0 or 3*f*f <= energy(a)**3, "small-order production bound")
                count += 1

    # Every cycle through one common edge; equality permits several returns.
    common = 0
    for returns in range(1, 7):
        for aweight in (1,2,3):
            a = [[0]*(returns+2) for _ in range(returns+2)]
            a[0][1] = aweight
            for k in range(2, returns+2):
                a[1][k], a[k][0] = k-1, returns+3-k
            f,c = cubic(a)
            check(27*c*c <= energy(a)**3 and f == 3*c, "common-edge all-order control")
            common += 1
    # Rational equality with four returns, all five energy sectors explicit.
    equality = [[0]*6 for _ in range(6)]
    equality[0][1] = 2
    for k in range(2,6):
        equality[1][k] = equality[k][0] = 1
    f,c = cubic(equality)
    check(f == 24 and energy(equality) == 12, "six-vertex equality values")
    check(3*f*f == energy(equality)**3, "six-vertex exact equality")

    base = [[[0]], [[0,1],[0,0]], [[0,1,0],[0,0,1],[1,0,0]]]
    bank = []
    for indices in product(range(3), repeat=3):
        for weights in product((1,2), repeat=3):
            q = data(3, (1,0,1), weights)  # 0->1->2->0
            bank.append(audit_substitution([base[i] for i in indices],q))
    rng = random.Random(20260908)
    for _ in range(18):
        blocks = [rng.choice(bank) for _ in range(3)]
        q = data(3,(1,0,1),[rng.randrange(1,4) for _ in range(3)])
        audit_substitution(blocks,q)

    # The general bound is unproved here. Record finite outcomes; do not call this
    # sample a universal proof or make it a dependency of the proved cases.
    scout = 0
    failures_f, failures_c = [], []
    for n in range(5,11):
        m = comb(n,2)
        for _ in range(1000):
            a = data(n,[rng.randrange(2) for _ in range(m)],
                     [rng.randrange(5) for _ in range(m)])
            f,c = cubic(a)
            e = energy(a)
            if f>0 and 3*f*f>e**3:
                failures_f.append((n,a,f,e))
            if 27*c*c>e**3:
                failures_c.append((n,a,c,e))
            scout += 1
    print("PROVED-SCOPE exact controls; general weighted bound is unproved here")
    print(f"all labelled n=2..4 carriers with weights 0,1,2: {count}")
    print(f"common-edge controls: {common}; equality F=24,E=12 on six vertices")
    print(f"recursive constant-contact substitution controls: {len(bank)+18}")
    print(f"always-active proved-scope gates: {GATES}")
    print(f"FINITE-EXACT OPEN-BOUND SCOUT: seed20260908, n5..10, weights0..4, {scout} draws")
    print(f"production counterexamples: {len(failures_f)}; stronger cycle-count counterexamples: {len(failures_c)}")
    if failures_f:
        print("FIRST PRODUCTION HOSTILE:",failures_f[0])
    if failures_c:
        print("FIRST CYCLE-COUNT HOSTILE:",failures_c[0])
    print("No general weighted-tournament theorem or PDE evolution claim follows from the scout.")


if __name__ == "__main__":
    main()
