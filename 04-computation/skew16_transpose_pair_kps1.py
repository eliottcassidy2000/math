# skew16_transpose_pair_kps1.py
# kind-pasteur-2026-06-09-S1, Branch D follow-up #2: confirm the transpose-pair structure.
# Known so far (hall_class_pin_kps1):
#   skew_16 ~ had.16.3 (explicit witness), skew_16 NOT ~ had.16.4 (exhaustive),
#   skew_16 NOT ~ skew_16^T (exhaustive).
# Now: skew_16^T ~ had.16.4 ?  had.16.3^T ~ had.16.4 ?  -> {16.3,16.4} = THE transpose pair.
# Also check self-transpose-equivalence of had.16.0/1/2 (completes transpose action on the 5 classes).

import os, time, itertools
import numpy as np

OUTPATH = os.path.join("05-knowledge", "results", "skew16_transpose_pair_kps1.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def w(self, *args):
        s = " ".join(str(a) for a in args)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()

T = Tee(OUTPATH)
W = T.w

def skew_tower(order):
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < order:
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    return S

def load_sloane(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if len(line) == 16 and set(line) <= {"+", "-"}:
                rows.append([1 if c == "+" else -1 for c in line])
            elif len(line) == 16 and set(line) <= {"0", "1"}:
                rows.append([1 if c == "0" else -1 for c in line])
    assert len(rows) == 16
    return np.array(rows, dtype=np.int64)

def colsig(rows_list, n):
    cols = []
    for c in range(n):
        v = tuple(r[c] for r in rows_list)
        nv = tuple(-x for x in v)
        cols.append(v if v <= nv else nv)
    cols.sort()
    return cols

def monomial_equiv(A, B, node_cap=20_000_000, time_cap=600.0):
    n = A.shape[0]
    Arows = [tuple(int(v) for v in A[i]) for i in range(n)]
    Brows = [tuple(int(v) for v in B[i]) for i in range(n)]
    sigB = [None] * (n + 1)
    for k in range(1, n + 1):
        sigB[k] = colsig(Brows[:k], n)
    nodes = [0]
    t0 = time.time()
    used = [False] * n
    chosen = []
    assign = []
    capped = [False]

    def bt(k):
        if k == n:
            return True
        for i in range(n):
            if used[i]:
                continue
            for s in (1, -1):
                nodes[0] += 1
                if nodes[0] > node_cap or time.time() - t0 > time_cap:
                    capped[0] = True
                    return False
                row = Arows[i] if s == 1 else tuple(-x for x in Arows[i])
                chosen.append(row)
                if colsig(chosen, n) == sigB[k + 1]:
                    used[i] = True
                    assign.append((i, s))
                    if bt(k + 1):
                        return True
                    used[i] = False
                    assign.pop()
                chosen.pop()
                if capped[0]:
                    return False
        return False

    ok = bt(0)
    if not ok:
        return ("CAP" if capped[0] else "NONE"), None, nodes[0], time.time() - t0
    P = np.array(chosen, dtype=np.int64)
    usedc = [False] * n
    tau = [None] * n
    delta = [None] * n
    for c in range(n):
        target = B[:, c]
        found = False
        for c2 in range(n):
            if usedc[c2]:
                continue
            if np.array_equal(P[:, c2], target):
                tau[c] = c2; delta[c] = 1; usedc[c2] = True; found = True; break
            if np.array_equal(P[:, c2], -target):
                tau[c] = c2; delta[c] = -1; usedc[c2] = True; found = True; break
        if not found:
            return "NONE", None, nodes[0], time.time() - t0
    sigma = [a[0] for a in assign]
    eps = [a[1] for a in assign]
    chk = np.empty((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            chk[i, j] = eps[i] * delta[j] * A[sigma[i], tau[j]]
    assert np.array_equal(chk, B), "witness verification FAILED"
    return "FOUND", (sigma, eps, tau, delta), nodes[0], time.time() - t0

S16 = skew_tower(16)
reps = {}
for name in ("had.16.0", "had.16.1", "had.16.2", "had.16.3", "had.16.4"):
    reps[name] = load_sloane(os.path.join("05-knowledge", "results", f"sloane_{name}_kps1.txt"))

W("Transpose-pair structure of the five order-16 classes + skew tower placement")
W("=" * 78)
jobs = [
    ("skew_16^T", S16.T.copy(), "had.16.4", reps["had.16.4"]),
    ("had.16.3^T", reps["had.16.3"].T.copy(), "had.16.4", reps["had.16.4"]),
    ("had.16.0^T", reps["had.16.0"].T.copy(), "had.16.0", reps["had.16.0"]),
    ("had.16.1^T", reps["had.16.1"].T.copy(), "had.16.1", reps["had.16.1"]),
    ("had.16.2^T", reps["had.16.2"].T.copy(), "had.16.2", reps["had.16.2"]),
]
for nameA, A, nameB, B in jobs:
    st, wit, nd, dt = monomial_equiv(A, B)
    W(f"{nameA} ~ {nameB}: {st}  nodes={nd}  time={dt:.1f}s")
    if st == "FOUND":
        sigma, eps, tau, delta = wit
        W(f"   row perm: {sigma}  row signs: {eps}")
        W(f"   col perm: {tau}  col signs: {delta}")

W("\nDONE.")
T.f.close()
