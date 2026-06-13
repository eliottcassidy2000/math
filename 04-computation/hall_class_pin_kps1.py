# hall_class_pin_kps1.py
# kind-pasteur-2026-06-09-S1, Branch D follow-up: pin the Hadamard equivalence class of the
# skew tower at order 16 against Sloane's library representatives had.16.0 .. had.16.4
# (the 5 classes of order 16) plus had.16.syl (Sylvester construction).
# Files downloaded to 05-knowledge/results/sloane_had.16.*_kps1.txt
#
# Method: invariants (quadruple profile, SNF over Z of H and (H+J)/2, GF(2) rank of (H+J)/2),
# then explicit monomial-equivalence witness search against the matching representative.

import sys, os, time, itertools, re
import numpy as np

OUTPATH = os.path.join("05-knowledge", "results", "hall_class_pin_kps1.out")

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

# ------------------------------------------------ towers
def skew_tower(order):
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < order:
        n = S.shape[0]
        I = np.eye(n, dtype=np.int64)
        S = np.block([[S, S], [S - 2 * I, 2 * I - S]])
    return S

def sylvester(order):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < order:
        H = np.block([[H, H], [H, -H]])
    return H

# ------------------------------------------------ load Sloane files
def load_sloane(path):
    rows = []
    aut = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if len(line) == 16 and set(line) <= {"+", "-"}:
                rows.append([1 if c == "+" else -1 for c in line])
            elif len(line) == 16 and set(line) <= {"0", "1"}:
                rows.append([1 if c == "0" else -1 for c in line])
            m = re.search(r"order\s*=\s*(\d+)", line)
            if m:
                aut = int(m.group(1))
    assert len(rows) == 16, f"{path}: got {len(rows)} rows"
    return np.array(rows, dtype=np.int64), aut

# ------------------------------------------------ invariants
def profile(H):
    n = H.shape[0]
    pairs = list(itertools.combinations(range(n), 2))
    PV = {p: H[p[0]] * H[p[1]] for p in pairs}
    counts = {}
    for (i, j, k, l) in itertools.combinations(range(n), 4):
        s = abs(int(np.dot(PV[(i, j)], PV[(k, l)])))
        counts[s] = counts.get(s, 0) + 1
    return dict(sorted(counts.items()))

def smith_normal_form_z(Ain):
    A = [[int(x) for x in row] for row in Ain]
    m, n = len(A), len(A[0])
    t = 0
    while t < min(m, n):
        piv = None; best = None
        for i in range(t, m):
            for j in range(t, n):
                a = A[i][j]
                if a != 0 and (best is None or abs(a) < best):
                    best = abs(a); piv = (i, j)
        if piv is None:
            break
        i0, j0 = piv
        A[t], A[i0] = A[i0], A[t]
        for r in range(m):
            A[r][t], A[r][j0] = A[r][j0], A[r][t]
        while True:
            again = False
            for i in range(t + 1, m):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    for j in range(t, n):
                        A[i][j] -= q * A[t][j]
                    if A[i][t] != 0:
                        A[t], A[i] = A[i], A[t]
                        again = True
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    for i in range(t, m):
                        A[i][j] -= q * A[i][t]
                    if A[t][j] != 0:
                        for i in range(t, m):
                            A[i][t], A[i][j] = A[i][j], A[i][t]
                        again = True
            if again:
                continue
            colc = all(A[i][t] == 0 for i in range(t + 1, m))
            rowc = all(A[t][j] == 0 for j in range(t + 1, n))
            if colc and rowc:
                d = A[t][t]
                bad = None
                for i in range(t + 1, m):
                    for j in range(t + 1, n):
                        if A[i][j] % d != 0:
                            bad = (i, j); break
                    if bad:
                        break
                if bad:
                    i, _ = bad
                    for jj in range(t, n):
                        A[t][jj] += A[i][jj]
                    continue
                break
        t += 1
    return [abs(A[i][i]) for i in range(min(m, n)) if A[i][i] != 0]

def snf_compact(d):
    out = []
    for v in d:
        if out and out[-1][0] == v:
            out[-1][1] += 1
        else:
            out.append([v, 1])
    return " ".join(f"{v}^{c}" if c > 1 else f"{v}" for v, c in out)

def gf2_rank(M):
    rows = []
    for r in M:
        x = 0
        for v in r:
            x = (x << 1) | (int(v) & 1)
        rows.append(x)
    rank = 0
    for b in range(M.shape[1] - 1, -1, -1):
        piv = None
        for idx, r in enumerate(rows):
            if (r >> b) & 1:
                piv = idx; break
        if piv is None:
            continue
        pr = rows.pop(piv)
        rows = [r ^ pr if (r >> b) & 1 else r for r in rows]
        rank += 1
    return rank

# ------------------------------------------------ witness search
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

# ================================================================ main
W("=" * 78)
W("HALL CLASS PINNING: skew_16 vs Sloane had.16.0..4 (5 classes) + had.16.syl")
W("=" * 78)

S16 = skew_tower(16)
H16 = sylvester(16)
J = np.ones((16, 16), dtype=np.int64)

mats = {}
for name in ("had.16.0", "had.16.1", "had.16.2", "had.16.3", "had.16.4", "had.16.syl"):
    path = os.path.join("05-knowledge", "results", f"sloane_{name}_kps1.txt")
    M, aut = load_sloane(path)
    n = 16
    assert np.array_equal(M @ M.T, n * np.eye(n, dtype=np.int64)), f"{name} not Hadamard!"
    mats[name] = (M, aut)
    W(f"loaded {name}: valid Hadamard, |Aut| = {aut}")

mats["SKEW_TOWER_16"] = (S16, None)
mats["SYLVESTER_TOWER_16"] = (H16, None)

W("\nInvariant table:")
W(f"{'matrix':22s} {'profile (|sum|:count)':38s} {'SNF((H+J)/2)':28s} {'gf2rk':5s}")
inv = {}
for name, (M, aut) in mats.items():
    p = profile(M)
    B = (M + J) // 2
    s = smith_normal_form_z(B)
    g = gf2_rank(B % 2)
    sh = smith_normal_form_z(M)
    inv[name] = (tuple(sorted(p.items())), tuple(s), g, tuple(sh))
    W(f"{name:22s} {str(p):38s} {snf_compact(s):28s} {g:5d}   SNF(H)={snf_compact(sh)}")

# which Sloane reps share ALL invariants with skew_16 / sylvester_16?
W("\nInvariant matching:")
for probe in ("SKEW_TOWER_16", "SYLVESTER_TOWER_16"):
    cands = [nm for nm in mats if nm.startswith("had") and inv[nm] == inv[probe]]
    W(f"  {probe}: invariant-matching Sloane reps: {cands}")

# explicit witness searches
W("\nWitness searches (monomial equivalence, exact):")
jobs = []
for probe in ("SYLVESTER_TOWER_16", "SKEW_TOWER_16"):
    for nm in ("had.16.0", "had.16.1", "had.16.2", "had.16.3", "had.16.4"):
        if inv[nm] == inv[probe]:
            jobs.append((probe, nm))
# also sanity: had.16.syl vs sylvester tower
jobs.append(("SYLVESTER_TOWER_16", "had.16.syl"))

for probe, nm in jobs:
    A = mats[probe][0]
    Bm = mats[nm][0]
    st, wit, nd, dt = monomial_equiv(A, Bm)
    W(f"  {probe} ~ {nm}: {st}  nodes={nd}  time={dt:.1f}s")
    if st == "FOUND":
        sigma, eps, tau, delta = wit
        W(f"     row perm:  {sigma}")
        W(f"     row signs: {eps}")
        W(f"     col perm:  {tau}")
        W(f"     col signs: {delta}")

# cross check: is skew_16 equivalent to its own transpose? (classes can split under transpose)
st, wit, nd, dt = monomial_equiv(S16, S16.T.copy())
W(f"\n  SKEW_TOWER_16 ~ its transpose: {st}  nodes={nd}  time={dt:.1f}s")

W("\nDONE.")
T.f.close()
