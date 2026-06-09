# verify_D_equivalence_kps1.py
# kind-pasteur-2026-06-09-S1 -- ADVERSARIAL VERIFICATION of Branch D Hadamard-equivalence claims
# (hadamard_tower_equivalence_kps1 / hall_class_pin_kps1 / skew16_transpose_pair_kps1).
#
# Everything here is INDEPENDENTLY recomputed with fresh code:
#   A. Bridge: skew_tower recursion == iterated D_skew from the 1-vertex tournament (orders 2..32),
#      plus known-good anchors H(C3)=3, H(D(C3))=45, H(regular T5)=15.
#   B. Sloane library authenticity: saved had.16.* files vs matrices re-fetched from
#      neilsloane.com/hadamard/ on 2026-06-09 (hardcoded below), Hadamard validity,
#      had.16.4 == had.16.3^T entrywise, had.16.0 == had.16.syl == Sylvester tower verbatim.
#   C. Independent invariants: quadruple profiles via DIRECT 4-row products (not the pair-dot
#      trick), GF(2) ranks of (H+J)/2 via numpy elimination, orders 16 and 32.
#   D. Entrywise re-verification of every printed witness (order-8 vs Sylvester;
#      skew_16~had.16.3; skew_16^T~had.16.4; had.16.1 and had.16.2 self-transpose).
#   E. Fresh exhaustive monomial-equivalence search (bitmask implementation, different code,
#      different traversal order): skew_16~had.16.3 (expect FOUND), skew_16 vs had.16.4
#      (expect NONE), skew_16 vs skew_16^T (expect NONE) + scramble self-tests of the checker.
#   F. Walsh claims: every entry of W^T S == 2 (mod 4) at m=8,16,32 (and the scope boundary at
#      m=2,4); (W^T S W)/m is a +-1 Hadamard matrix; Corr recursion; butterfly block identity.
#   G. Own 2-adic SNF (modular arithmetic, no coefficient blowup) for S/H at 16 and 32 and for
#      (H+J)/2, with slogdet check that no odd invariant-factor part is being missed.

import os, sys, time, itertools
import numpy as np

OUT = os.path.join("05-knowledge", "results", "verify_D_equivalence_kps1.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def w(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.f.write(s + "\n")
        self.f.flush()

T = Tee(OUT)
W = T.w
VERDICTS = []
def verdict(name, ok, detail=""):
    VERDICTS.append((name, ok, detail))
    W(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

# ---------------------------------------------------------------- towers (fresh code)
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

# ================================================================ PART A: bridge to Branch D
W("=" * 78)
W("PART A: skew_tower recursion vs ACTUAL Branch-D doubling D_skew (core library)")
W("=" * 78)
sys.path.insert(0, "04-computation")
from skew_doubling_core_kps1 import D_skew, H_count, M_of

A1 = np.zeros((1, 1), dtype=np.int64)   # the 1-vertex tournament
A = A1
ok_all = True
for level in range(5):
    A, Mp = D_skew(A)
    n = A.shape[0]
    S_from_D = Mp + np.eye(n, dtype=np.int64)
    same = np.array_equal(S_from_D, skew_tower(n))
    ok_all = ok_all and same
    had = np.array_equal(S_from_D @ S_from_D.T, n * np.eye(n, dtype=np.int64))
    skw = np.array_equal(S_from_D + S_from_D.T, 2 * np.eye(n, dtype=np.int64))
    W(f"  order {n:3d}: S(D^{level+1}(T1)) == skew_tower({n}): {same} | Hadamard: {had} | skew-type: {skw}")
    ok_all = ok_all and had and skw
verdict("A1 skew_tower == iterated D_skew, Hadamard+skew at 2,4,8,16,32", ok_all)

# anchors
C3 = np.zeros((3, 3), dtype=np.int64)
for i in range(3):
    C3[i, (i + 1) % 3] = 1
hC3 = H_count(C3)
AD, _ = D_skew(C3)
hDC3 = H_count(AD)
R5 = np.zeros((5, 5), dtype=np.int64)
for i in range(5):
    R5[i, (i + 1) % 5] = 1
    R5[i, (i + 2) % 5] = 1
hR5 = H_count(R5)
verdict("A2 anchors H(C3)=3, H(D(C3))=45, H(regularT5)=15",
        hC3 == 3 and hDC3 == 45 and hR5 == 15, f"got {hC3},{hDC3},{hR5}")

# ================================================================ PART B: Sloane authenticity
W("\n" + "=" * 78)
W("PART B: Sloane order-16 library -- authenticity and structure")
W("=" * 78)

# Matrices re-fetched from http://neilsloane.com/hadamard/had.16.{0..4}.txt on 2026-06-09
WEB = {
"had.16.0": """++++++++++++++++
+-+-+-+-+-+-+-+-
++--++--++--++--
+--++--++--++--+
++++----++++----
+-+--+-++-+--+-+
++----++++----++
+--+-++-+--+-++-
++++++++--------
+-+-+-+--+-+-+-+
++--++----++--++
+--++--+-++--++-
++++--------++++
+-+--+-+-+-++-+-
++----++--++++--
+--+-++--++-+--+""",
"had.16.1": """++++++++++++++++
+-+-+-+-+-+-+-+-
++--++--++--++--
+--++--++--++--+
++++----++++----
+-+--+-++-+--+-+
++----++++----++
+--+-++-+--+-++-
++++++++--------
+-+-+--+-+-+-++-
++--++----++--++
+--++-+--++--+-+
++++--------++++
+-+--++--+-++--+
++----++--++++--
+--+-+-+-++-+-+-""",
"had.16.2": """++++++++++++++++
+-+-+-+-+-+-+-+-
++--++--++--++--
+--++--++--++--+
++++----++++----
+-+--+-++-+--+-+
++----++++----++
+--+-++-+--+-++-
++++++++--------
++++--------++++
++--+-+---++-+-+
++---+-+--+++-+-
+-+-+--+-+-+-++-
+-+--++--+-++--+
+--+++---++---++
+--+--++-++-++--""",
"had.16.3": """++++++++++++++++
+-+-+-+-+-+-+-+-
++--++--++--++--
+--++--++--++--+
++++----++++----
+-+--+-++-+--+-+
++----++++----++
+--+-++-+--+-++-
++++++++--------
+++-+------+-+++
++-+---+--+-+++-
++---++---+++--+
+-++-+---+--+-++
+-+---++-+-+++--
+--++-+--++--+-+
+---++-+-+++--+-""",
"had.16.4": """++++++++++++++++
+-+-+-+-++++----
++--++--++--++--
+--++--++-+-+-+-
++++----++----++
+-+--+-++--++--+
++----+++--+-++-
+--+-++-+-+--+-+
++++++++--------
+-+-+-+-----++++
++--++----++--++
+--++--+-+-+-+-+
++++------++++--
+-+--+-+-++--++-
++----++-++-+--+
+--+-++--+-++-+-""",
}

def parse_pm(txt):
    rows = [[1 if c == "+" else -1 for c in line] for line in txt.strip().split("\n")]
    return np.array(rows, dtype=np.int64)

def load_saved(name):
    path = os.path.join("05-knowledge", "results", f"sloane_{name}_kps1.txt")
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

mats = {}
ok_auth = True
for nm in ("had.16.0", "had.16.1", "had.16.2", "had.16.3", "had.16.4"):
    saved = load_saved(nm)
    web = parse_pm(WEB[nm])
    same = np.array_equal(saved, web)
    ok_auth = ok_auth and same
    mats[nm] = saved
    W(f"  {nm}: saved file == web-refetched matrix: {same}")
verdict("B1 all five saved Sloane representatives authentic (== web refetch)", ok_auth)

syl_file = load_saved("had.16.syl")
ok_had = True
for nm, M in list(mats.items()) + [("had.16.syl", syl_file)]:
    ok_had = ok_had and np.array_equal(M @ M.T, 16 * np.eye(16, dtype=np.int64))
verdict("B2 all six saved matrices are valid Hadamard of order 16", ok_had)

H16 = sylvester(16)
S16 = skew_tower(16)
verdict("B3 had.16.4 == had.16.3^T entrywise (the literal transpose pair)",
        np.array_equal(mats["had.16.4"], mats["had.16.3"].T))
verdict("B4 sylvester tower == had.16.0 verbatim", np.array_equal(H16, mats["had.16.0"]))
verdict("B5 sylvester tower == had.16.syl verbatim", np.array_equal(H16, syl_file))
verdict("B6 had.16.0 is symmetric (self-transpose with identity witness)",
        np.array_equal(mats["had.16.0"], mats["had.16.0"].T))

# ================================================================ PART C: independent invariants
W("\n" + "=" * 78)
W("PART C: independent profile + GF(2)-rank recomputation (direct 4-row products)")
W("=" * 78)

def profile_direct(H):
    n = H.shape[0]
    rows = [H[i].astype(np.int64) for i in range(n)]
    counts = {}
    for (i, j, k, l) in itertools.combinations(range(n), 4):
        s = abs(int(np.sum(rows[i] * rows[j] * rows[k] * rows[l])))
        counts[s] = counts.get(s, 0) + 1
    return dict(sorted(counts.items()))

def gf2_rank_np(M):
    A = (np.asarray(M) % 2).astype(np.uint8).copy()
    r = 0
    nr, nc = A.shape
    for c in range(nc):
        piv = None
        for i in range(r, nr):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(nr):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r

CLAIMED_PROF = {
    "sylvester_16": {0: 1680, 16: 140},
    "skew_16":      {0: 1344, 8: 448, 16: 28},
    "had.16.0":     {0: 1680, 16: 140},
    "had.16.1":     {0: 1488, 8: 256, 16: 76},
    "had.16.2":     {0: 1392, 8: 384, 16: 44},
    "had.16.3":     {0: 1344, 8: 448, 16: 28},
    "had.16.4":     {0: 1344, 8: 448, 16: 28},
}
CLAIMED_GF2 = {"sylvester_16": 5, "skew_16": 8, "had.16.0": 5, "had.16.1": 6,
               "had.16.2": 7, "had.16.3": 8, "had.16.4": 8}

J16 = np.ones((16, 16), dtype=np.int64)
allmats = dict(mats)
allmats["sylvester_16"] = H16
allmats["skew_16"] = S16
ok_prof = ok_gf2 = True
for nm, M in allmats.items():
    p = profile_direct(M)
    g = gf2_rank_np((M + J16) // 2)
    okp = (p == CLAIMED_PROF[nm])
    okg = (g == CLAIMED_GF2[nm])
    ok_prof = ok_prof and okp
    ok_gf2 = ok_gf2 and okg
    W(f"  {nm:14s} profile={p}  (claimed match: {okp})   gf2rank((H+J)/2)={g} (claimed {CLAIMED_GF2[nm]}: {okg})")
verdict("C1 all order-16 quadruple profiles match claims (incl. Hall counts 140/76/44/28/28)", ok_prof)
verdict("C2 all order-16 GF(2) ranks of (H+J)/2 match claims (5/6/7/8/8; skew=8, syl=5)", ok_gf2)

# order 32
S32 = skew_tower(32)
H32 = sylvester(32)
t0 = time.time()
p32s = profile_direct(S32)
p32h = profile_direct(H32)
W(f"  skew_32      profile={p32s}")
W(f"  sylvester_32 profile={p32h}   ({time.time()-t0:.1f}s)")
verdict("C3 skew_32 profile == {0:26656, 8:7056, 16:1568, 24:624, 32:56}",
        p32s == {0: 26656, 8: 7056, 16: 1568, 24: 624, 32: 56})
verdict("C4 sylvester_32 profile == {0:34720, 32:1240} (and 1240 == C(32,3)-theory 32*31*30/24)",
        p32h == {0: 34720, 32: 1240} and 32 * 31 * 30 // 24 == 1240)
J32 = np.ones((32, 32), dtype=np.int64)
g32s = gf2_rank_np((S32 + J32) // 2)
g32h = gf2_rank_np((H32 + J32) // 2)
verdict("C5 GF(2) rank (H+J)/2 at 32: skew=16, sylvester=6", g32s == 16 and g32h == 6,
        f"got {g32s},{g32h}")

# ================================================================ PART D: witness re-verification
W("\n" + "=" * 78)
W("PART D: entrywise re-verification of every printed witness")
W("=" * 78)

def check_witness(Amat, Bmat, sigma, eps, tau, delta):
    n = Amat.shape[0]
    for i in range(n):
        for j in range(n):
            if Bmat[i, j] != eps[i] * delta[j] * Amat[sigma[i], tau[j]]:
                return False
    return True

S8 = skew_tower(8)
H8 = sylvester(8)
verdict("D1 order-8 witness: sylvester8[i,j]=eps[i]*delta[j]*skew8[sigma[i],tau[j]]",
        check_witness(S8, H8,
                      [0, 1, 2, 3, 4, 7, 5, 6], [1, 1, 1, -1, 1, -1, -1, 1],
                      [1, 2, 3, 4, 5, 6, 7, 0], [1] * 8))

verdict("D2 witness skew_16 ~ had.16.3 (from hall_class_pin out)",
        check_witness(S16, mats["had.16.3"],
                      [0, 1, 2, 3, 7, 4, 6, 5, 8, 13, 15, 14, 12, 9, 11, 10],
                      [1, 1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1],
                      [5, 2, 7, 4, 1, 6, 3, 8, 13, 10, 15, 12, 9, 14, 11, 0],
                      [1] * 16))

verdict("D3 witness skew_16^T ~ had.16.4 (from skew16_transpose_pair out)",
        check_witness(S16.T.copy(), mats["had.16.4"],
                      [0, 1, 2, 7, 3, 6, 13, 4, 8, 9, 10, 15, 11, 14, 5, 12],
                      [1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1],
                      [5, 1, 3, 4, 0, 6, 7, 2, 13, 9, 15, 11, 10, 14, 12, 8],
                      [-1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]))

verdict("D4 witness had.16.1^T ~ had.16.1 (self-transpose)",
        check_witness(mats["had.16.1"].T.copy(), mats["had.16.1"],
                      [0, 2, 1, 3, 8, 10, 9, 11, 4, 6, 5, 7, 12, 14, 13, 15], [1] * 16,
                      [0, 2, 1, 3, 8, 10, 9, 11, 4, 6, 5, 7, 12, 14, 13, 15], [1] * 16))

verdict("D5 witness had.16.2^T ~ had.16.2 (self-transpose)",
        check_witness(mats["had.16.2"].T.copy(), mats["had.16.2"],
                      [0, 1, 2, 3, 8, 9, 10, 11, 4, 12, 6, 14, 5, 13, 7, 15], [1] * 16,
                      [0, 1, 2, 3, 8, 12, 10, 14, 4, 5, 6, 7, 9, 13, 11, 15], [1] * 16))

# ================================================================ PART E: fresh exhaustive search
W("\n" + "=" * 78)
W("PART E: independent exhaustive monomial-equivalence search (bitmask, fresh code)")
W("=" * 78)

def equiv_search(Amat, Bmat, node_cap=60_000_000, time_cap=1800.0):
    """Exhaustive search for signed row perm + signed col perm with
    B[i,j] = eps[i]*delta[j]*A[sigma[i],tau[j]].  Status: FOUND / NONE / CAP.
    Traversal order deliberately DIFFERENT from the sibling script (source rows
    descending, negative sign first); for a NONE proof the full pruned tree is
    traversed either way, so the node count is order-invariant and comparable."""
    n = Amat.shape[0]
    full = (1 << n) - 1

    def mask_of(row):
        x = 0
        for v in row:
            x = (x << 1) | (1 if v > 0 else 0)
        return x

    def bits_of(m):
        return tuple((m >> (n - 1 - j)) & 1 for j in range(n))

    Am = [mask_of(Amat[i]) for i in range(n)]
    Abits = [(bits_of(m), bits_of(m ^ full)) for m in Am]
    Bbits = [bits_of(mask_of(Bmat[i])) for i in range(n)]

    sigB = []
    colv = [0] * n
    for k in range(n):
        colv = [(v << 1) | b for v, b in zip(colv, Bbits[k])]
        mk = (1 << (k + 1)) - 1
        sigB.append(sorted(v if v <= (v ^ mk) else (v ^ mk) for v in colv))

    nodes = 0
    t0 = time.time()
    used = [False] * n
    assign = []
    status = ["NONE"]

    def dfs(k, colvals):
        nonlocal nodes
        mk = (1 << (k + 1)) - 1
        target = sigB[k]
        for i in range(n - 1, -1, -1):          # descending source order
            if used[i]:
                continue
            for s in (1, 0):                     # negative sign first
                nodes += 1
                if nodes > node_cap:
                    status[0] = "CAP"
                    return False
                if (nodes & 0xFFFFF) == 0 and time.time() - t0 > time_cap:
                    status[0] = "CAP"
                    return False
                bits = Abits[i][s]
                nv = [(v << 1) | b for v, b in zip(colvals, bits)]
                sig = sorted(x if x <= (x ^ mk) else (x ^ mk) for x in nv)
                if sig == target:
                    used[i] = True
                    assign.append((i, s))
                    if k + 1 == n or dfs(k + 1, nv):
                        return True
                    used[i] = False
                    assign.pop()
                if status[0] == "CAP":
                    return False
        return False

    found = dfs(0, [0] * n)
    dt = time.time() - t0
    if not found:
        return status[0], None, nodes, dt
    # reconstruct witness
    sigma = [a[0] for a in assign]
    eps = [1 if a[1] == 0 else -1 for a in assign]
    P = np.array([Amat[sigma[i]] * eps[i] for i in range(n)], dtype=np.int64)
    usedc = [False] * n
    tau = [None] * n
    delta = [None] * n
    for c in range(n):
        for c2 in range(n):
            if usedc[c2]:
                continue
            if np.array_equal(P[:, c2], Bmat[:, c]):
                tau[c], delta[c], usedc[c2] = c2, 1, True
                break
            if np.array_equal(P[:, c2], -Bmat[:, c]):
                tau[c], delta[c], usedc[c2] = c2, -1, True
                break
        assert tau[c] is not None, "column matching failed though signature matched"
    assert check_witness(Amat, Bmat, sigma, eps, tau, delta), "witness check failed"
    return "FOUND", (sigma, eps, tau, delta), nodes, dt

rng = np.random.default_rng(424242)
def scramble(Hm):
    n = Hm.shape[0]
    rp = rng.permutation(n)
    cp = rng.permutation(n)
    rs = rng.choice(np.array([-1, 1], dtype=np.int64), n)
    cs = rng.choice(np.array([-1, 1], dtype=np.int64), n)
    return (rs[:, None] * Hm[np.ix_(rp, cp)]) * cs[None, :]

# checker self-tests
st1, _, nd1, dt1 = equiv_search(scramble(S16), S16)
W(f"  selftest scramble(skew_16) ~ skew_16:        {st1}  nodes={nd1} time={dt1:.1f}s (expect FOUND)")
st2, _, nd2, dt2 = equiv_search(scramble(mats["had.16.3"]), S16)
W(f"  selftest scramble(had.16.3) ~ skew_16:       {st2}  nodes={nd2} time={dt2:.1f}s (expect FOUND)")
st3, _, nd3, dt3 = equiv_search(mats["had.16.0"], mats["had.16.1"])
W(f"  selftest had.16.0 vs had.16.1 (known noneq): {st3}  nodes={nd3} time={dt3:.1f}s (expect NONE)")
verdict("E0 checker self-tests (2x FOUND on scrambles, 1x NONE on known-inequivalent)",
        st1 == "FOUND" and st2 == "FOUND" and st3 == "NONE")

st, wit, nd, dt = equiv_search(S16, mats["had.16.3"])
W(f"  skew_16 ~ had.16.3:   {st}  nodes={nd}  time={dt:.1f}s")
verdict("E1 skew_16 ~ had.16.3 (independent search FOUND + witness verified)", st == "FOUND")

st, wit, nd, dt = equiv_search(S16, mats["had.16.4"])
W(f"  skew_16 ~ had.16.4:   {st}  nodes={nd}  time={dt:.1f}s   (sibling: NONE @ 10,188,512 nodes)")
verdict("E2 skew_16 NOT ~ had.16.4 (independent exhaustive NONE)", st == "NONE",
        f"nodes={nd}")

st, wit, nd, dt = equiv_search(S16, S16.T.copy())
W(f"  skew_16 ~ skew_16^T:  {st}  nodes={nd}  time={dt:.1f}s   (sibling: NONE @ 10,188,512 nodes)")
verdict("E3 skew_16 NOT ~ its own transpose (independent exhaustive NONE)", st == "NONE",
        f"nodes={nd}")

# ================================================================ PART F: Walsh-domain claims
W("\n" + "=" * 78)
W("PART F: Walsh-domain claims")
W("=" * 78)
ok_mod4 = True
for m in (2, 4, 8, 16, 32, 64):
    Wm = sylvester(m)
    Sm = skew_tower(m)
    N = Wm.T @ Sm
    all2mod4 = bool(np.all(N % 4 == 2))
    nz = int(np.count_nonzero(N))
    W(f"  m={m:3d}: all entries of W^T S == 2 (mod 4): {all2mod4}   nonzeros {nz}/{m*m}")
    if m in (8, 16, 32):
        ok_mod4 = ok_mod4 and all2mod4 and nz == m * m
verdict("F1 every entry of W^T S == 2 (mod 4) at m=8,16,32 (=> 100% dense, no zeros possible)",
        ok_mod4)

ok_conj = True
hists = {}
for m in (8, 16, 32):
    Wm = sylvester(m)
    Sm = skew_tower(m)
    G = Wm.T @ Sm @ Wm
    ok = bool(np.all(G % m == 0))
    G = G // m
    pm1 = bool(np.all(np.abs(G) == 1))
    had = np.array_equal(G @ G.T, m * np.eye(m, dtype=np.int64))
    vals, cts = np.unique(G, return_counts=True)
    hists[m] = {int(v): int(c) for v, c in zip(vals, cts)}
    ok_conj = ok_conj and ok and pm1 and had
    W(f"  m={m:3d}: (W^T S W)/m integral={ok}, +-1 entries={pm1}, Hadamard={had}, hist={hists[m]}")
verdict("F2 (W^T S W)/m is +-1 Hadamard at m=8,16,32 with claimed histograms",
        ok_conj and hists[8] == {-1: 28, 1: 36} and hists[16] == {-1: 120, 1: 136}
        and hists[32] == {-1: 496, 1: 528})

# Corr tower: define Corr := S - W and CHECK it satisfies the claimed recursion
ok_corr = True
for m in (2, 4, 8, 16):
    Cm = skew_tower(m) - sylvester(m)
    C2m = skew_tower(2 * m) - sylvester(2 * m)
    I = np.eye(m, dtype=np.int64)
    pred = np.block([[Cm, Cm], [Cm - 2 * I, 2 * I - Cm]])
    ok_corr = ok_corr and np.array_equal(C2m, pred)
verdict("F3 Corr := S - W satisfies Corr_2m = [[C,C],[C-2I,2I-C]] (orders 4..32)", ok_corr)

ok_bfly = True
rng2 = np.random.default_rng(7)
for m in (8, 64, 512):
    S2m = skew_tower(2 * m)
    Sm = skew_tower(m)
    x = rng2.integers(-99, 100, size=2 * m).astype(np.int64)
    a, b = x[:m], x[m:]
    lhs = S2m @ x
    rhs = np.concatenate([Sm @ (a + b), Sm @ (a - b) - 2 * (a - b)])
    ok_bfly = ok_bfly and np.array_equal(lhs, rhs)
verdict("F4 butterfly block identity S_2m x = [S_m(x1+x2); S_m(x1-x2)-2(x1-x2)] at m=8,64,512",
        ok_bfly)

# ================================================================ PART G: 2-adic SNF (own impl)
W("\n" + "=" * 78)
W("PART G: own 2-adic SNF (modular, K=40 bits) + odd-part exclusion via slogdet")
W("=" * 78)

def snf_2adic(Min, K=40):
    mod = 1 << K
    Awork = [[int(x) % mod for x in row] for row in Min]
    n = len(Awork)
    m = len(Awork[0])
    def v2(x):
        return K if x == 0 else (x & -x).bit_length() - 1
    facs = []
    t = 0
    while t < min(n, m):
        best, bi, bj = K, -1, -1
        for i in range(t, n):
            for j in range(t, m):
                v = v2(Awork[i][j])
                if v < best:
                    best, bi, bj = v, i, j
        if bi < 0 or best >= K:
            break
        Awork[t], Awork[bi] = Awork[bi], Awork[t]
        for r in range(n):
            Awork[r][t], Awork[r][bj] = Awork[r][bj], Awork[r][t]
        v = best
        u = (Awork[t][t] >> v) % mod
        uinv = pow(u, -1, mod)
        Awork[t] = [(x * uinv) % mod for x in Awork[t]]
        for i in range(t + 1, n):
            a = Awork[i][t]
            if a:
                f = (a >> v) % mod
                rowt = Awork[t]
                Awork[i] = [(x - f * y) % mod for x, y in zip(Awork[i], rowt)]
        for j in range(t + 1, m):
            Awork[t][j] = 0   # column ops touch only row t (column below pivot already 0)
        facs.append(1 << v)
        t += 1
    return facs

def compact(facs):
    out = []
    for v in facs:
        if out and out[-1][0] == v:
            out[-1][1] += 1
        else:
            out.append([v, 1])
    return " ".join(f"{v}^{c}" if c > 1 else f"{v}" for v, c in out)

CLAIMED_SNF = {
    "skew_16 H":  [1] + [2] * 7 + [8] * 7 + [16],
    "syl_16 H":   [1] + [2] * 4 + [4] * 6 + [8] * 4 + [16],
    "skew_16 B":  [1] * 8 + [4] * 7 + [8],
    "syl_16 B":   [1] * 5 + [2] * 6 + [4] * 4 + [8],
    "skew_32 H":  [1] + [2] * 15 + [16] * 15 + [32],
    "syl_32 H":   [1] + [2] * 5 + [4] * 10 + [8] * 10 + [16] * 5 + [32],
    "skew_32 B":  [1] * 16 + [8] * 15 + [16],
    "syl_32 B":   [1] * 6 + [2] * 10 + [4] * 10 + [8] * 5 + [16],
}
import math
ok_snf = True
jobs = [("skew_16 H", S16), ("syl_16 H", H16),
        ("skew_16 B", (S16 + J16) // 2), ("syl_16 B", (H16 + J16) // 2),
        ("skew_32 H", S32), ("syl_32 H", H32),
        ("skew_32 B", (S32 + J32) // 2), ("syl_32 B", (H32 + J32) // 2)]
for nm, M in jobs:
    facs = snf_2adic(M)
    sgn, logdet = np.linalg.slogdet(M.astype(np.float64))
    log2det = logdet / math.log(2.0)
    log2facs = sum(math.log2(f) for f in facs)
    odd_free = abs(log2det - log2facs) < 1e-6 and len(facs) == M.shape[0]
    okm = (facs == CLAIMED_SNF[nm]) and odd_free
    ok_snf = ok_snf and okm
    W(f"  {nm:10s}: 2-adic SNF = {compact(facs):28s} log2|det|={log2det:.6f} "
      f"(sum log2 facs={log2facs:.0f})  match claim+no odd part: {okm}")
verdict("G1 all eight SNFs match claims; |det| fully accounted (no odd invariant parts)", ok_snf)

# ================================================================ summary
W("\n" + "=" * 78)
W("SUMMARY")
W("=" * 78)
nfail = sum(1 for _, ok, _ in VERDICTS if not ok)
for nm, ok, det in VERDICTS:
    W(f"  [{'PASS' if ok else 'FAIL'}] {nm}")
W(f"\n{len(VERDICTS) - nfail}/{len(VERDICTS)} checks passed, {nfail} failed.")
W("DONE.")
T.f.close()
