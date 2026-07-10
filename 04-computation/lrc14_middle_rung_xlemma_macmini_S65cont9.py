from itertools import combinations
from math import gcd
from functools import reduce
viol = tot = 0
tightX = []
Xdist = {}
for D in range(12, 22):
    h = D - 12
    for mid in combinations(range(1, D), 11):
        A = (0,) + mid + (D,)
        if reduce(gcd, A[1:]) != 1: continue
        tot += 1
        S = set(); F = set()
        for i in range(13):
            for j in range(i, 13):
                (S if i < j else F).add(A[i] + A[j])
        AA = S | F
        B = len(S)
        X = len(AA) - B - 2                    # monad identity check below
        # identity: extreme doubles 0,2D always outside S; middle doubles outside S counted by X
        Xdirect = sum(1 for i in range(1, 12) if 2*A[i] not in S)
        assert X == Xdirect, (A, X, Xdirect)
        missing = (2*D + 1) - len(AA)
        if missing > h - X:
            viol += 1
            if viol < 5: print("VIOLATION", A, "missing", missing, "h", h, "X", X)
        if X > 0 and missing == h - X:
            tightX.append((D, A, X, missing))
        Xdist[X] = Xdist.get(X, 0) + 1
print(f"exhausted {tot} gcd-1 13-sets, D in [12,21]: violations of [missing <= h - X]: {viol}")
print("X distribution:", dict(sorted(Xdist.items())))
print(f"tight cases with X > 0: {len(tightX)}; examples:")
for D, A, X, m in tightX[:6]:
    print(f"  D={D} X={X} missing={m} A={A}")
