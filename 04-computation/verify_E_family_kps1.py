#!/usr/bin/env python3
"""
verify_E_family_kps1.py — kind-pasteur-2026-06-09-S1 (ADVERSARIAL VERIFIER)

Independent recomputation of the branch-E claims made in:
  block_doubling_classification_kps1.py / .out
  regular7_doubling_Htest_kps1.py / .out
  consecutive_circulant_iso_kps1.py / .out

EVERYTHING here is implemented fresh (nothing imported from
skew_doubling_core_kps1):
  * H counted by a PREDECESSOR-direction bitmask DP (sibling used a
    successor-direction push DP), validated against raw permutation brute
    force at orders 6, 8 and 10 before being trusted at order 14;
  * fresh canonical form (full-matrix bytes key, sibling used triu+tril);
  * fresh doubling constructor via np.block (sibling used np.kron);
  * iso claims re-proved with a certificate-producing backtracker, every
    positive answer re-verified by an explicit permutation check, and the
    small non-iso claims settled by exhaustive permutation sweep.

Checks (numbered in the output):
  V0  anchors: H(C3)=3, H(regT5)=15, H(Paley T7)=189 (brute AND dp)
  V1  16-member validity (n=3,4 exhaustive), symmetry identities g1,g2,g3
      (n=3 exhaustive), orbit structure 4+8+4 with named reps separated
  V2  Kron square law and block-diag(M'^2) <=> a+d=0  (n=4 exhaustive)
  V3  skew-Hadamard 16-member test on bordered-C3 seed (expect exactly the
      8 trP=0 members) + towers to order 32 for two orbit-B members
  V4  H tables: n=3,4 (dp vs BRUTE at orders 6,8 vs published),
      n=5 (dp vs published incl. means/argmax; brute spot checks at order 10
      on all three regular-T5 doubles), within-orbit split structure
  V5  n=6 exhaustive max H (claim: 45) + the three C3 doubles
  V6  n=7 regular classes re-found with a fresh seed; all 9 H values at
      order 14 recomputed and compared to the published table
  V7  iso certificates: K2(C3)~SCblow(C3); D_skew(C3) vs both (exhaustive
      non-iso); K2(regT5)~SCblow(regT5); K2(U7)~SCblow(U7) at order 14

Output: 05-knowledge/results/verify_E_family_kps1.out
"""
import itertools
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUTPATH = os.path.join(ROOT, "05-knowledge", "results",
                       "verify_E_family_kps1.out")
out = open(OUTPATH, "w", encoding="utf-8")
FAILURES = []


def w(line=""):
    out.write(line + "\n")
    out.flush()
    print(line)


def check(label, ok):
    w("  [%s] %s" % ("PASS" if ok else "FAIL", label))
    if not ok:
        FAILURES.append(label)
    return ok


# ---------------- fresh H counters ----------------

def H_dp(A):
    """Hamiltonian path count, predecessor-direction subset DP."""
    n = A.shape[0]
    Al = A.tolist()
    inmask = []
    for v in range(n):
        m = 0
        for u in range(n):
            if Al[u][v]:
                m |= 1 << u
        inmask.append(m)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(3, full + 1):
        if mask & (mask - 1) == 0:
            continue
        row = dp[mask]
        mm = mask
        while mm:
            lb = mm & -mm
            mm ^= lb
            v = lb.bit_length() - 1
            prev = mask ^ lb
            preds = inmask[v] & prev
            t = 0
            rowp = dp[prev]
            while preds:
                pb = preds & -preds
                preds ^= pb
                t += rowp[pb.bit_length() - 1]
            row[v] = t
    return sum(dp[full])


def H_brute(A):
    """Raw permutation sweep with early abort (ground truth, orders <= 10)."""
    n = A.shape[0]
    Al = [tuple(bool(x) for x in row) for row in A.tolist()]
    c = 0
    for p in itertools.permutations(range(n)):
        prev = p[0]
        for k in range(1, n):
            q = p[k]
            if not Al[prev][q]:
                break
            prev = q
        else:
            c += 1
    return c


# ---------------- fresh constructions ----------------

def circulant(n, S):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if (j - i) % n in S:
                A[i, j] = 1
    return A


def member(A, s):
    """M' = [[aM, bM+beta I],[bM-beta I, dM]]; returns (adjacency, M')."""
    a, b, d, be = s
    n = A.shape[0]
    M = A - A.T
    I = np.eye(n, dtype=np.int64)
    Mp = np.block([[a * M, b * M + be * I], [b * M - be * I, d * M]])
    return (Mp > 0).astype(np.int64), Mp


REP = {"D_skew": (1, 1, -1, 1), "TK2": (1, 1, 1, 1), "SCblow": (1, -1, 1, 1)}
ALL16 = [(a, b, d, be) for a in (1, -1) for b in (1, -1)
         for d in (1, -1) for be in (1, -1)]


def is_tournament(Ad):
    m = Ad.shape[0]
    return np.array_equal(Ad + Ad.T,
                          np.ones((m, m), dtype=np.int64)
                          - np.eye(m, dtype=np.int64))


def canon_key(A):
    n = A.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        k = A[np.ix_(p, p)].tobytes()
        if best is None or k < best:
            best = k
    return best


def iso_classes_fresh(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen = {}
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        key = canon_key(A)
        if key not in seen:
            seen[key] = A.copy()
    return list(seen.values())


def is_DRT_fresh(A):
    n = A.shape[0]
    s = A.sum(axis=1)
    if len(set(int(x) for x in s)) != 1:
        return False
    vals = set()
    for u in range(n):
        for v in range(u + 1, n):
            vals.add(int(np.dot(A[u], A[v])))
    return len(vals) == 1


def is_SH(S):
    m = S.shape[0]
    I = np.eye(m, dtype=np.int64)
    return ((np.abs(S) == 1).all()
            and np.array_equal(S + S.T, 2 * I)
            and np.array_equal(S @ S.T, m * I))


# ---------------- fresh iso machinery ----------------

def find_iso(A, B):
    """Backtracking; returns mapping list u->v or None."""
    n = A.shape[0]
    sA = [int(x) for x in A.sum(axis=1)]
    sB = [int(x) for x in B.sum(axis=1)]
    if sorted(sA) != sorted(sB):
        return None
    cand = [[v for v in range(n) if sB[v] == sA[u]] for u in range(n)]
    order = sorted(range(n), key=lambda u: len(cand[u]))
    mapping = [-1] * n
    used = [False] * n

    def dfs(k):
        if k == n:
            return True
        u = order[k]
        for v in cand[u]:
            if used[v]:
                continue
            ok = True
            for kk in range(k):
                uu = order[kk]
                vv = mapping[uu]
                if A[u, uu] != B[v, vv] or A[uu, u] != B[vv, v]:
                    ok = False
                    break
            if ok:
                mapping[u] = v
                used[v] = True
                if dfs(k + 1):
                    return True
                mapping[u] = -1
                used[v] = False
        return False

    return mapping if dfs(0) else None


def mapping_valid(A, B, mp):
    p = np.array(mp)
    return np.array_equal(B[np.ix_(p, p)], A)


def iso_exhaustive(A, B):
    n = A.shape[0]
    for p in itertools.permutations(range(n)):
        if np.array_equal(A[np.ix_(p, p)], B):
            return True
    return False


# ---------------- orbit machinery (fresh) ----------------

def g1(s):
    return (s[2], s[1], s[0], -s[3])


def g2(s):
    return (-s[0], -s[1], -s[2], -s[3])


def g3(s):
    return (-s[0], -s[1], -s[2], s[3])


def orbits_fresh():
    seen = set()
    orbs = []
    for s in ALL16:
        if s in seen:
            continue
        orb = {s}
        frontier = [s]
        while frontier:
            t = frontier.pop()
            for g in (g1, g2, g3):
                u = g(t)
                if u not in orb:
                    orb.add(u)
                    frontier.append(u)
        seen |= orb
        orbs.append(orb)
    return orbs


# ================= published numbers (from sibling .out files) =============

PUB3 = sorted([(1, (0, 1, 2), 13, 1, 41), (3, (1, 1, 1), 45, 45, 45)])
PUB4 = sorted([(1, (0, 1, 2, 3), 95, 1, 629), (3, (1, 1, 1, 3), 189, 45, 633),
               (5, (1, 1, 2, 2), 523, 393, 653), (3, (0, 2, 2, 2), 333, 45, 633)])
PUB5 = sorted([
    (1, (0, 1, 2, 3, 4), 1033, 1, 14937),
    (3, (1, 1, 1, 3, 4), 1809, 45, 14961),
    (5, (1, 1, 2, 2, 4), 2817, 393, 15109),
    (3, (0, 2, 2, 2, 4), 1809, 45, 14973),
    (9, (1, 1, 2, 3, 3), 8137, 3225, 15313),
    (5, (0, 2, 2, 3, 3), 5289, 393, 15109),
    (9, (1, 1, 2, 3, 3), 8145, 2785, 15201),
    (3, (0, 1, 3, 3, 3), 3561, 45, 14961),
    (11, (1, 2, 2, 2, 3), 11017, 6069, 15461),
    (15, (1, 2, 2, 2, 3), 12129, 10773, 15261),
    (13, (1, 2, 2, 2, 3), 11625, 8097, 15305),
    (15, (2, 2, 2, 2, 2), 15505, 15565, 15565)])
PUB7 = {171: (24577581, 24392325, 24485541),
        175: (24540117, 24855901, 24855901),
        189: (24793755, 24589929, 24453597)}
PUB5_MEANS = (82876, 47436, 182156)  # column sums D, K2, SC over 12 classes


def main():
    t_start = time.time()
    w("=== verify_E_family_kps1 — adversarial independent verification ===")
    w("")

    # ---------------- V0 anchors ----------------
    w("--- V0: anchors (brute force AND fresh dp) ---")
    C3 = circulant(3, {1})
    regT5 = circulant(5, {1, 2})
    P7 = circulant(7, {1, 2, 4})
    U7 = circulant(7, {1, 2, 3})
    check("H(C3) = 3 brute=%d dp=%d" % (H_brute(C3), H_dp(C3)),
          H_brute(C3) == 3 == H_dp(C3))
    check("H(regT5) = 15 brute=%d dp=%d" % (H_brute(regT5), H_dp(regT5)),
          H_brute(regT5) == 15 == H_dp(regT5))
    check("H(Paley T7) = 189 brute=%d dp=%d" % (H_brute(P7), H_dp(P7)),
          H_brute(P7) == 189 == H_dp(P7))
    check("U7 = C7({1,2,3}): H=%d (expect 175), DRT=%s (expect False)"
          % (H_dp(U7), is_DRT_fresh(U7)),
          H_dp(U7) == 175 and not is_DRT_fresh(U7))
    check("P7 DRT = True", is_DRT_fresh(P7))
    w("")

    # ---------------- V1 validity / symmetries / orbits ----------------
    w("--- V1: 16-member validity, symmetry identities, orbits ---")
    for n in (3, 4):
        bad = 0
        tot = 0
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for bits in range(1 << len(pairs)):
            A = np.zeros((n, n), dtype=np.int64)
            for k, (i, j) in enumerate(pairs):
                if (bits >> k) & 1:
                    A[i, j] = 1
                else:
                    A[j, i] = 1
            tot += 1
            for s in ALL16:
                if not is_tournament(member(A, s)[0]):
                    bad += 1
        check("validity n=%d: %d tournaments x 16, invalid = %d" % (n, tot, bad),
              bad == 0)
    # symmetry identities, exhaustive n=3
    SW = np.zeros((6, 6), dtype=np.int64)
    SW[:3, 3:] = np.eye(3, dtype=np.int64)
    SW[3:, :3] = np.eye(3, dtype=np.int64)
    ok1 = ok2 = ok3 = True
    pairs3 = [(i, j) for i in range(3) for j in range(i + 1, 3)]
    for bits in range(1 << 3):
        A = np.zeros((3, 3), dtype=np.int64)
        for k, (i, j) in enumerate(pairs3):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        for s in ALL16:
            Mp = member(A, s)[1]
            ok1 &= np.array_equal(SW @ Mp @ SW, member(A, g1(s))[1])
            ok2 &= np.array_equal(-Mp, member(A, g2(s))[1])
            ok3 &= np.array_equal(member(A.T.copy(), s)[1], member(A, g3(s))[1])
    check("g1 swap identity (exh n=3)", ok1)
    check("g2 negate identity (exh n=3)", ok2)
    check("g3 pre-op identity (exh n=3)", ok3)
    orbs = orbits_fresh()
    sizes = sorted(len(o) for o in orbs)
    check("3 orbits, sizes %s (expect [4,4,8])" % sizes,
          len(orbs) == 3 and sizes == [4, 4, 8])
    orb_of = {}
    for k, o in enumerate(orbs):
        for s in o:
            orb_of[s] = k
    named = {nm: orb_of[s] for nm, s in REP.items()}
    check("named reps in pairwise distinct orbits %s" % named,
          len(set(named.values())) == 3)
    check("D_skew rep in the size-8 orbit",
          len(orbs[named["D_skew"]]) == 8)
    dso = orbs[named["D_skew"]]
    check("D_skew orbit = exactly the a+d=0 members",
          all((s[0] + s[2] == 0) == (s in dso) for s in ALL16))
    w("")

    # ---------------- V2 Kron square law / block-diag ----------------
    w("--- V2: square law + block-diagonality (exhaustive n=4) ---")
    okk = True
    okbd = True
    okform = True
    pairs4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    I4 = np.eye(4, dtype=np.int64)
    for bits in range(1 << 6):
        A = np.zeros((4, 4), dtype=np.int64)
        for k, (i, j) in enumerate(pairs4):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        M = A - A.T
        M2 = M @ M
        for s in ALL16:
            a, b, d, be = s
            Mp = member(A, s)[1]
            sq = Mp @ Mp
            P = np.array([[a, b], [b, d]], dtype=np.int64)
            Q = be * np.array([[0, 1], [-1, 0]], dtype=np.int64)
            law = (np.kron(P @ P, M2) + np.kron(P @ Q + Q @ P, M)
                   - np.kron(np.eye(2, dtype=np.int64), I4))
            okk &= np.array_equal(sq, law)
            off = sq[:4, 4:]
            okbd &= (not off.any()) == (a + d == 0)
            okform &= np.array_equal(off, (a + d) * (b * M2 + be * M))
    check("Kron square law M'^2 = P^2(x)M^2 + {P,Q}(x)M - I (64x16)", okk)
    check("block-diag(M'^2) <=> a+d = 0 (64x16)", okbd)
    check("off-diag block = (a+d)(b M^2 + beta M) (64x16)", okform)
    w("")

    # ---------------- V3 skew-Hadamard ----------------
    w("--- V3: skew-Hadamard preservation (bordered C3 seed) ---")
    A4 = np.zeros((4, 4), dtype=np.int64)
    A4[0, 1:] = 1
    A4[1:, 1:] = C3
    S4 = (A4 - A4.T) + np.eye(4, dtype=np.int64)
    check("seed order-4 skew-Hadamard", is_SH(S4))
    preservers = []
    for s in ALL16:
        Mp = member(A4, s)[1]
        Sp = Mp + np.eye(8, dtype=np.int64)
        if is_SH(Sp):
            preservers.append(s)
    check("preservers = 8/16 (got %d)" % len(preservers), len(preservers) == 8)
    check("preservers are EXACTLY the a+d=0 set",
          all(s[0] + s[2] == 0 for s in preservers)
          and len(preservers) == 8)
    for s in ((1, 1, -1, 1), (-1, 1, 1, 1)):
        T = A4
        oks = []
        for _ in range(3):
            Ad, Mp = member(T, s)
            Sp = Mp + np.eye(Mp.shape[0], dtype=np.int64)
            oks.append(is_SH(Sp))
            T = Ad
        check("tower %s: orders 8,16,32 SH = %s" % (str(s), oks), all(oks))
    w("")

    # ---------------- V4 H tables ----------------
    w("--- V4: H tables n=3,4 (dp vs brute vs published) ---")
    for n, pub in ((3, PUB3), (4, PUB4)):
        rows = []
        ok_bv = True
        for A in iso_classes_fresh(n):
            hT = H_dp(A)
            trip = []
            for nm in ("D_skew", "TK2", "SCblow"):
                Ad = member(A, REP[nm])[0]
                hd = H_dp(Ad)
                hb = H_brute(Ad)
                ok_bv &= (hd == hb)
                trip.append(hd)
            rows.append((hT, tuple(sorted(int(x) for x in A.sum(axis=1))),
                         trip[0], trip[1], trip[2]))
        check("n=%d: dp == brute on all %d doubles (order %d)"
              % (n, 3 * len(rows), 2 * n), ok_bv)
        check("n=%d: table matches published" % n, sorted(rows) == pub)
        for r in sorted(rows):
            w("      %s" % str(r))

    w("")
    w("--- V4b: full n=5 table (order-10 dp; brute spot check on regular class) ---")
    cls5 = iso_classes_fresh(5)
    check("n=5 has 12 iso classes (got %d)" % len(cls5), len(cls5) == 12)
    rows5 = []
    sums = [0, 0, 0]
    argmax_ok = True
    t0 = time.time()
    Hall5 = {}  # (classidx, s) -> H for all 16 members (for orbit splits)
    for i, A in enumerate(cls5):
        hT = H_dp(A)
        trip = []
        for j, nm in enumerate(("D_skew", "TK2", "SCblow")):
            hd = H_dp(member(A, REP[nm])[0])
            trip.append(hd)
            sums[j] += hd
        for s in ALL16:
            Hall5[(i, s)] = H_dp(member(A, s)[0])
        rows5.append((hT, tuple(sorted(int(x) for x in A.sum(axis=1))),
                      trip[0], trip[1], trip[2]))
        argmax_ok &= (trip[2] == max(trip))
    check("n=5 table matches published (12 rows, all 3 columns)",
          sorted(rows5) == PUB5)
    check("column sums (D,K2,SC) = %s expect %s" % (tuple(sums), PUB5_MEANS),
          tuple(sums) == PUB5_MEANS)
    w("      => means %.1f %.1f %.1f (published 6906.3 / 3953.0 / 15179.7)"
      % (sums[0] / 12, sums[1] / 12, sums[2] / 12))
    check("SCblow is argmax (with ties) in all 12 classes", argmax_ok)
    # brute spot checks at order 10 on the regular class
    Areg = [A for A in cls5
            if tuple(sorted(int(x) for x in A.sum(axis=1))) == (2,) * 5][0]
    for nm, expect in (("D_skew", 15505), ("TK2", 15565), ("SCblow", 15565)):
        Ad = member(Areg, REP[nm])[0]
        tb = time.time()
        hb = H_brute(Ad)
        check("order-10 BRUTE %s(regT5) = %d expect %d (%.1fs)"
              % (nm, hb, expect, time.time() - tb), hb == expect)
    # within-orbit split structure
    sc5 = [canon_key(c) == canon_key(np.ascontiguousarray(c.T))
           for c in cls5]
    split_report = []
    orbit_ok = True
    for k, orb in enumerate(orbs):
        nm = [n_ for n_, s_ in REP.items() if s_ in orb][0]
        for i in range(12):
            vals = sorted(set(Hall5[(i, s)] for s in orb))
            if len(vals) > 1:
                split_report.append((nm, i, vals))
                if nm != "D_skew" or sc5[i]:
                    orbit_ok = False
    check("only D_skew orbit splits, and only on non-SC classes", orbit_ok)
    w("      splits: %s" % split_report)
    dsplit_vals = sorted(tuple(v) for _, _, v in split_report)
    check("split values are {1809,3561} x2 and {2817,5289} x2",
          dsplit_vals == sorted([(1809, 3561), (1809, 3561),
                                 (2817, 5289), (2817, 5289)]))
    w("      (n=5 block done in %.1fs)" % (time.time() - t0))
    w("")

    # ---------------- V5 n=6 exhaustive max ----------------
    w("--- V5: n=6 exhaustive max H over all 32768 labeled tournaments ---")
    t0 = time.time()
    pairs6 = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    mx = 0
    nmax = 0
    for bits in range(1 << 15):
        A = np.zeros((6, 6), dtype=np.int64)
        for k, (i, j) in enumerate(pairs6):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        h = H_dp(A)
        if h > mx:
            mx, nmax = h, 1
        elif h == mx:
            nmax += 1
    check("max H at n=6 = %d (expect 45); attained by %d labeled tournaments"
          % (mx, nmax), mx == 45)
    w("      (%.1fs)" % (time.time() - t0))
    trips = {}
    for nm in ("D_skew", "TK2", "SCblow"):
        Ad = member(C3, REP[nm])[0]
        trips[nm] = Ad
        sc = tuple(sorted(int(x) for x in Ad.sum(axis=1)))
        check("%s(C3): H=%d (expect 45), scores=%s" % (nm, H_dp(Ad), sc),
              H_dp(Ad) == 45 and sc == (2, 2, 2, 3, 3, 3))
    w("")

    # ---------------- V6 n=7 regular classes at order 14 ----------------
    w("--- V6: the three regular n=7 classes, doubles at order 14 ---")
    rng = np.random.default_rng(424242)  # different seed from sibling
    found = {}
    tries = 0
    while len(found) < 3 and tries < 500000:
        tries += 1
        A = np.zeros((7, 7), dtype=np.int64)
        for i in range(7):
            for j in range(i + 1, 7):
                if rng.integers(2):
                    A[i, j] = 1
                else:
                    A[j, i] = 1
        if tuple(int(x) for x in A.sum(axis=1)) != (3,) * 7:
            continue
        key = canon_key(A)
        if key not in found:
            found[key] = A.copy()
    check("found exactly 3 distinct regular n=7 classes in %d tries" % tries,
          len(found) == 3)
    kU, kP = canon_key(U7), canon_key(P7)
    ids = {}
    for key, A in found.items():
        h = H_dp(A)
        tag = "U7" if key == kU else ("Paley" if key == kP else "third")
        ids[h] = (tag, A)
    check("classes have H = 171, 175, 189 (got %s)" % sorted(ids),
          sorted(ids) == [171, 175, 189])
    check("H=175 class is U7; H=189 class is Paley",
          ids.get(175, ("", 0))[0] == "U7"
          and ids.get(189, ("", 0))[0] == "Paley")
    all7_ok = True
    for h in sorted(ids):
        tag, A = ids[h]
        got = []
        for nm in ("D_skew", "TK2", "SCblow"):
            tb = time.time()
            got.append(H_dp(member(A, REP[nm])[0]))
            dt = time.time() - tb
        exp = PUB7[h]
        ok = tuple(got) == exp
        all7_ok &= ok
        w("    H(T)=%d (%s): D=%d K2=%d SC=%d  expect %s  [%s] (last dp %.1fs)"
          % (h, tag, got[0], got[1], got[2], exp, "PASS" if ok else "FAIL", dt))
        if not ok:
            FAILURES.append("n=7 H mismatch at H(T)=%d" % h)
    check("all 9 order-14 H values match published", all7_ok)
    check("U7 coincidence: K2 == SC == 24855901",
          PUB7[175][1] == PUB7[175][2] == 24855901)
    w("")

    # ---------------- V7 iso certificates ----------------
    w("--- V7: isomorphism claims with explicit certificates ---")
    AK3, AS3, AD3 = (member(C3, REP[nm])[0] for nm in ("TK2", "SCblow", "D_skew"))
    mp = find_iso(AK3, AS3)
    check("K2(C3) ~ SCblow(C3): certificate found and verified",
          mp is not None and mapping_valid(AK3, AS3, mp))
    check("D_skew(C3) !~ K2(C3) (exhaustive 720 perms)",
          not iso_exhaustive(AD3, AK3))
    check("D_skew(C3) !~ SCblow(C3) (exhaustive 720 perms)",
          not iso_exhaustive(AD3, AS3))
    AK5 = member(regT5, REP["TK2"])[0]
    AS5 = member(regT5, REP["SCblow"])[0]
    mp = find_iso(AK5, AS5)
    check("K2(regT5) ~ SCblow(regT5) at order 10: certificate verified",
          mp is not None and mapping_valid(AK5, AS5, mp))
    AK7 = member(U7, REP["TK2"])[0]
    AS7 = member(U7, REP["SCblow"])[0]
    t0 = time.time()
    mp = find_iso(AK7, AS7)
    ok = mp is not None and mapping_valid(AK7, AS7, mp)
    check("K2(U7) ~ SCblow(U7) at order 14: certificate verified (%.2fs)"
          % (time.time() - t0), ok)
    if mp is not None:
        w("      certificate permutation: %s" % mp)
    # Paley split is already non-iso by H (24589929 != 24453597, recomputed in V6)
    w("")

    # ---------------- summary ----------------
    w("=== SUMMARY ===")
    if FAILURES:
        w("FAILURES (%d):" % len(FAILURES))
        for f in FAILURES:
            w("  - " + f)
    else:
        w("ALL CHECKS PASSED — no discrepancy found with the branch-E claims.")
    w("total time %.1fs" % (time.time() - t_start))
    out.close()


if __name__ == "__main__":
    main()
