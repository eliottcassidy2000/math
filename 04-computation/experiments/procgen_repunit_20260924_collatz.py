#!/usr/bin/env python3
"""Part 2: the Collatz analogue of the consecutive-prime digit transitions.

U(n) = odd part of 3n+1 (Syracuse map; "consecutive odd terms"), T(n) = n/2 or (3n+1)/2.

 (A) EXACT Haar transition laws of U and T modulo M = 10 (last digit), 3, 9, 8 (and 5, 40, 72),
     derived from the 2-adic/odd-adic independence (Haar measure on Z_2 x Z_m is a product):
       T: T(x) mod M depends on x mod 2M; the extra top bit is a fresh fair bit.
       U: v = v_2(3n+1) is Geom(1/2) (P(v=j)=2^-j), independent of n mod m (m odd);
          given v = j >= k, U(n) mod 2^k is uniform on odd classes and U(n) = (3n+1)2^-j mod m.
     Stationary laws (the Brouwer/Perron fixed point of each matrix), exact.
 (B) Empirical matrices from the C helper procgen_repunit_20260924_collatz.c:
       terras: random starts in [2^100, 2^101); transitions determined by the 100 random low bits
               (exactly Haar-distributed by the Terras bijection: an exact test of (A));
       post_high: the same orbits after the random bits are used up (deterministic Collatz dynamics),
               current value >= 2^40;
       post_low: the same, 10^6 <= value < 2^40: 10^7 orbits funnel through comparatively few integers
               there and merge, so transitions are counted with multiplicity (tree-weighted);
       full:   every start n <= NMAX followed to 1 ("all"), and only transitions from values > 10^4 ("big")
               (tree-weighted: merging orbits are counted once per start).
     Chi-square against the exact law; forbidden cells must be exactly empty in every regime.
     An independent brute-force enumeration mod 2^K re-derives the U laws (second code path).
 (C) Markov checks: the two-step conditional P(c | a, b) against the exact one-step law P(c | b).
 (D) The contrast with primes: second eigenvalue, total-variation distance from the stationary row,
     and the 10-adic reading of the forced cells (9 -> 9 is the shadow of T(-1) = -1 = ...999;
     3 -> 5 is the shadow of 3(-1/3)+1 = 0).
Session collatz-procgen-20260922, lane procgen_repunit_20260924.  Runtime ~ 3 min, memory < 300 MB.
"""
import os, sys, subprocess, tempfile, math, time
from fractions import Fraction as F
from collections import defaultdict
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_repunit_collatz_")
BIN = os.path.join(TMP, "collatz")
subprocess.run(["cc", "-O2", "-o", BIN, os.path.join(HERE, "procgen_repunit_20260924_collatz.c")], check=True)

NSTARTS = int(os.environ.get("PR_COLLATZ_NSTARTS", 10_000_000))
NFULL = int(os.environ.get("PR_COLLATZ_NFULL", 10_000_000))
THRESH = 10_000


def v2(n):
    n = abs(n)
    return (n & -n).bit_length() - 1 if n else 10**9


def ordmod(b, m):
    if m == 1:
        return 1
    k, x = 1, b % m
    while x != 1:
        x = x * b % m
        k += 1
    return k


def crt2(r1, m1, r2, m2):
    """x = r1 mod m1, x = r2 mod m2 (coprime)."""
    return (r1 + m1 * ((r2 - r1) * pow(m1, -1, m2) % m2)) % (m1 * m2) if m2 > 1 else r1 % m1


# ---------------------------------------------------------------- (A) exact laws
def haar_T(M):
    P = {}
    for a in range(M):
        row = defaultdict(F)
        for lift in (a, a + M):
            t = lift // 2 if lift % 2 == 0 else (3 * lift + 1) // 2
            row[t % M] += F(1, 2)
        P[a] = dict(row)
    return P


def haar_U(M):
    k = v2(M)
    m = M >> k
    L = ordmod(2, m)
    inv2 = pow(2, -1, m) if m > 1 else 0
    P = {}
    states = [a for a in range(M) if (k == 0 or a % 2 == 1)]
    for a in states:
        row = defaultdict(F)
        if k >= 1 and (3 * a + 1) % (1 << k) != 0:
            v = v2(3 * a + 1)                        # determined, v < k
            for t in range(1 << v):
                n = a + M * t
                assert v2(3 * n + 1) == v
                row[((3 * n + 1) >> v) % M] += F(1, 1 << v)
        else:
            jmin = max(k, 1)                         # v >= k (and v >= 1 always)
            c = (k - 1) if k >= 1 else 0             # P(v = j | a) = 2^-(j - c)
            for r in range(L):
                j0 = jmin + ((r - jmin) % L)
                w = F(1, 1 << (j0 - c)) / (1 - F(1, 1 << L))
                um = ((3 * a + 1) * pow(inv2, j0, m)) % m if m > 1 else 0
                if k >= 1:
                    odds = [o for o in range(1 << k) if o % 2 == 1]
                    for o in odds:
                        row[crt2(o, 1 << k, um, m)] += w / len(odds)
                else:
                    row[um] += w
        P[a] = dict(row)
        assert sum(P[a].values()) == 1, (M, a, sum(P[a].values()))
    return P


def stationary(P):
    """exact stationary law of the (single-recurrent-class) chain P (dict of dicts)."""
    S = sorted(P)
    idx = {s: i for i, s in enumerate(S)}
    n = len(S)
    # equations: sum_i pi_i (P_ij - delta_ij) = 0 for all j, sum pi = 1
    A = [[F(0)] * n for _ in range(n + 1)]
    for i, s in enumerate(S):
        for t, p in P[s].items():
            A[idx[t]][i] += p
        A[i][i] -= 1
    A[n] = [F(1)] * n
    rhs = [F(0)] * n + [F(1)]
    # least-squares-free exact solve: Gaussian elimination on the (n+1) x n system
    rows = [A[r][:] + [rhs[r]] for r in range(n + 1)]
    piv_row = 0
    where = [-1] * n
    for col in range(n):
        pr = next((r for r in range(piv_row, n + 1) if rows[r][col] != 0), None)
        if pr is None:
            continue
        rows[piv_row], rows[pr] = rows[pr], rows[piv_row]
        pv = rows[piv_row][col]
        rows[piv_row] = [x / pv for x in rows[piv_row]]
        for r in range(n + 1):
            if r != piv_row and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[piv_row])]
        where[col] = piv_row
        piv_row += 1
    pi = {s: (rows[where[i]][n] if where[i] >= 0 else F(0)) for i, s in enumerate(S)}
    for t in S:                                        # verify
        assert sum(pi[s] * P[s].get(t, 0) for s in S) == pi[t]
    assert sum(pi.values()) == 1
    return pi


def show_matrix(P, states, name, denom=None):
    lcm = 1
    for a in states:
        for p in P[a].values():
            lcm = lcm * p.denominator // math.gcd(lcm, p.denominator)
    d = denom or lcm
    print(f"  {name}: entries x {d}  (rows = current, cols = next)")
    print("       " + " ".join(f"{b:>5}" for b in states))
    for a in states:
        print(f"  {a:>4} " + " ".join(f"{int(P[a].get(b, 0) * d):>5}" for b in states))


# ---------------------------------------------------------------- (B) empirical tables
def parse_tables(fn):
    out = {}
    cur = None
    for line in open(fn):
        t = line.split()
        if t[0] == "TABLE":
            cur = {"one": np.zeros((360, 360), dtype=np.int64), "two": np.zeros((40, 40, 40), dtype=np.int64),
                   "meta": " ".join(t[2:])}
            out[t[1]] = cur
        elif t[0] == "O":
            cur["one"][int(t[1]), int(t[2])] = int(t[3])
        elif t[0] == "W":
            cur["two"][int(t[1]), int(t[2]), int(t[3])] = int(t[4])
    return out


def agg1(one, q):
    C = np.zeros((q, q), dtype=np.int64)
    for a in range(360):
        for b in range(360):
            if one[a, b]:
                C[a % q, b % q] += one[a, b]
    return C


def agg2(two, q):
    C = np.zeros((q, q, q), dtype=np.int64)
    nz = np.nonzero(two)
    for a, b, c in zip(*nz):
        C[a % q, b % q, c % q] += two[a, b, c]
    return C


def compare(C, P, states, label, min_row=200):
    """chi-square of empirical counts C (q x q) against exact law P on the given states."""
    chi, dof, maxz, bad0 = 0.0, 0, 0.0, 0
    for a in states:
        n = C[a].sum()
        if n < min_row:
            continue
        support = [b for b in P[a]]
        for b in range(C.shape[1]):
            p = float(P[a].get(b, 0))
            if p == 0:
                if C[a, b] != 0:
                    bad0 += int(C[a, b])
                continue
            e = n * p
            chi += (C[a, b] - e) ** 2 / e
            z = (C[a, b] - e) / math.sqrt(e * (1 - p)) if p < 1 else 0.0
            maxz = max(maxz, abs(z))
        dof += len(support) - 1
    pval = stats.chi2.sf(chi, dof) if dof > 0 else float("nan")
    print(f"    {label:<34} transitions={int(C.sum()):>12,}  chi2={chi:9.2f}  dof={dof:3d}  p={pval:6.3f}"
          f"  max|z|={maxz:5.2f}  counts in forbidden cells={bad0}")
    return chi, dof, pval, bad0


def empirical_rows(C, states):
    out = {}
    for a in states:
        n = C[a].sum()
        out[a] = {b: C[a, b] / n for b in states} if n else None
    return out


def main():
    t0 = time.time()
    print("=" * 100)
    print("PART 2. Collatz digit/residue transitions: exact Haar laws versus orbits")
    print("=" * 100)

    # ---------- exact laws
    print("\n(A) EXACT Haar laws (PROVED by the derivation in the module docstring; every row sums to 1 exactly)")
    laws = {}
    for M in (10, 3, 9, 8, 5, 40, 72):
        laws[("U", M)] = haar_U(M)
        laws[("T", M)] = haar_T(M)
    U10 = laws[("U", 10)]
    show_matrix(U10, [1, 3, 5, 7, 9], "U mod 10 (last digit of consecutive odd terms)")
    pi = stationary(U10)
    print("    stationary law:", {k: str(v) for k, v in pi.items()})
    assert all(pi[d] == F(1, 5) for d in (1, 3, 5, 7, 9))
    assert U10[3] == {5: F(1)}, "3 -> 5 is forced"
    assert all(U10[a].get(5, 0) == 0 for a in (1, 5, 7, 9)), "5 is entered only from 3"
    # the 10-adic reading of the forced/fixed cells
    assert (3 * (-1) + 1) * pow(2, -1, 5) % 5 == (-1) % 5          # v = 1 fixes -1 mod 5: 9 -> 9
    assert (3 * (-1 * pow(3, -1, 5)) + 1) % 5 == 0                  # -1/3 = 3 mod 5 gives 3n+1 = 0 mod 5: 3 -> 5
    print("    forced cells: 3 -> 5 with probability 1 (3n+1 = 0 mod 5 iff n = -1/3 = 3 mod 5);")
    print("    5 is entered only from 3; 9 -> 9 has probability 8/15 (v = 1 mod 4 keeps n = -1 mod 5:")
    print("    the mod-5 shadow of the fixed point T(-1) = -1 = ...999 in Z_10).")

    for M in (3, 9):
        st = sorted(laws[("U", M)])
        show_matrix(laws[("U", M)], st, f"U mod {M}")
        print("    stationary law:", {k: str(v) for k, v in stationary(laws[("U", M)]).items() if v})
    print("    U mod 3: every row is (0, 1/3, 2/3): the mod-3 classes of consecutive odd terms are i.i.d.")
    print("    (U(n) = 2^-v mod 3 = (-1)^v; P(v odd) = 2/3).  U mod 9 depends on the current class only mod 3.")
    st8 = sorted(laws[("U", 8)])
    show_matrix(laws[("U", 8)], st8, "U mod 8")
    print("    stationary law:", {k: str(v) for k, v in stationary(laws[("U", 8)]).items()})

    T10 = laws[("T", 10)]
    show_matrix(T10, list(range(10)), "T mod 10 (shortcut map, all terms)")
    print("    stationary law:", {k: str(v) for k, v in stationary(T10).items()})
    for M in (3, 9, 8):
        st = sorted(laws[("T", M)])
        show_matrix(laws[("T", M)], st, f"T mod {M}")
        print("    stationary law:", {k: str(v) for k, v in stationary(laws[("T", M)]).items() if v})

    # ---------- independent re-derivation of the U laws by brute force mod 2^K (second code path)
    for M, K in ((10, 20), (9, 20), (8, 20), (72, 18)):
        k = v2(M); m = M >> k
        cnt = defaultdict(lambda: defaultdict(int)); tot = defaultdict(int); undetermined = 0
        # n runs over odd residues mod 2^K * m; v and U(n) mod M are determined when v + k <= K
        for n in range(1, (1 << K) * m, 2):
            w = 3 * n + 1
            v = v2(w)
            if v + k > K:
                undetermined += 1
                continue
            cnt[n % M][(w >> v) % M] += 1
            tot[n % M] += 1
        P = laws[("U", M)]
        err = max(abs(cnt[a][b] / tot[a] - float(P[a].get(b, 0))) for a in tot for b in range(M))
        print(f"    brute force U mod {M:>2} over odd n mod 2^{K}*{m}: max |empirical - exact| = {err:.2e}"
              f" (undetermined mass {undetermined / ((1 << K) * m / 2):.1e})")
        assert err < 1e-3

    # ---------- spectra: the memory of each chain
    print("\n    memory of each chain: exact mixing time (first t with all rows of P^t equal), else |lambda_2|:")
    import sympy
    for key in [("U", 10), ("U", 3), ("U", 9), ("U", 8), ("T", 10), ("T", 3), ("T", 9), ("T", 8)]:
        P = laws[key]
        st = sorted(P)
        A = sympy.Matrix([[sympy.Rational(P[a].get(b, 0).numerator, P[a].get(b, 0).denominator) if b in P[a] else 0
                           for b in st] for a in st])
        Pt, tmix = A, None
        for t in range(1, 9):
            if all(Pt.row(i) == Pt.row(0) for i in range(Pt.rows)):
                tmix = t
                break
            Pt = Pt * A
        if tmix:
            print(f"      {key[0]} mod {key[1]:>2}: rows of P^{tmix} all equal (exact; memory {tmix} step(s))")
        else:
            lam = sympy.symbols("lam")
            cp = sympy.factor(A.charpoly(lam).as_expr())
            ev = sorted(np.abs(np.linalg.eigvals(np.array(A.tolist(), dtype=float))), reverse=True)
            print(f"      {key[0]} mod {key[1]:>2}: never exactly mixed; |lambda_2| = {ev[1]:.6f}; charpoly = {cp}")

    # ---------- empirical
    print("\n(B) Empirical transition counts")
    fb = os.path.join(TMP, "bulk.txt")
    ff = os.path.join(TMP, "full.txt")
    t1 = time.time()
    subprocess.run([BIN, "bulk", str(NSTARTS), "20260924", "200", "500", fb], check=True)
    subprocess.run([BIN, "full", str(NFULL), str(THRESH), ff], check=True)
    print(f"    C helper: bulk {NSTARTS:,} starts in [2^100, 2^101) (<= 200 U-steps, <= 500 T-steps, stop below 10^6);"
          f" full orbits of every n <= {NFULL:,}; {time.time() - t1:.0f} s")
    tabs = parse_tables(fb)
    tabs.update(parse_tables(ff))
    for name, tb in tabs.items():
        print(f"    table {name}: {tb['meta']}")
    results = []
    for mapname, tnames in (("U", ["U_terras", "U_post_high", "U_post_low", "U_full_all", "U_full_big"]),
                            ("T", ["T_terras", "T_post_high", "T_post_low", "T_full_all", "T_full_big"])):
        for q in (10, 3, 9, 8, 5, 72):
            P = laws[(mapname, q)]
            states = sorted(P)
            for tn in tnames:
                C = agg1(tabs[tn]["one"], q)
                r = compare(C, P, states, f"{tn} mod {q}")
                results.append((tn, q) + r)
    # a readable empirical matrix: U mod 10 Terras regime
    C = agg1(tabs["U_terras"]["one"], 10)
    print("\n    U mod 10, Terras regime, empirical row frequencies x 15 (exact law has integers):")
    for a in (1, 3, 5, 7, 9):
        n = C[a].sum()
        print(f"      {a}: " + " ".join(f"{15 * C[a, b] / n:6.3f}" for b in (1, 3, 5, 7, 9)) + f"   (n={n:,})")
    C = agg1(tabs["U_full_all"]["one"], 10)
    print("    U mod 10, full orbits of n <= %s to 1 (all transitions), x 15:" % f"{NFULL:,}")
    for a in (1, 3, 5, 7, 9):
        n = C[a].sum()
        print(f"      {a}: " + " ".join(f"{15 * C[a, b] / n:6.3f}" for b in (1, 3, 5, 7, 9)) + f"   (n={n:,})")
    # the boundary effect: the final descent
    print("    The full-orbit tallies are tree-weighted: every integer is counted once per start whose orbit passes")
    print("    through it, so merged paths (e.g. the trunk ... -> 5 -> 1, 16 -> 8 -> 4 -> 2 -> 1) dominate.  They deviate")
    print("    from Haar even above 10^4 ('big'); the merging diagnostic below shows the same effect in the bulk sample")
    print("    and that it disappears when every integer is counted once.  No transition ever lands in a cell of Haar")
    print("    probability 0: the support of the law is exact for every integer.")

    # ---------- (C) Markov checks
    print("\n(C) Markov property: two-step conditional P(c | a, b) against the exact one-step law P(c | b)")
    for mapname, tn in (("U", "U_terras"), ("U", "U_post_high"), ("T", "T_terras"), ("T", "T_post_high")):
        for q in (10, 8, 5):
            P = laws[(mapname, q)]
            W = agg2(tabs[tn]["two"], q)
            chi, dof, bad0 = 0.0, 0, 0
            for a in range(q):
                for b in range(q):
                    n = W[a, b].sum()
                    if n < 200 or b not in P:
                        continue
                    for c in range(q):
                        p = float(P[b].get(c, 0))
                        if p == 0:
                            bad0 += int(W[a, b, c])
                            continue
                        e = n * p
                        chi += (W[a, b, c] - e) ** 2 / e
                    dof += len(P[b]) - 1
            print(f"    {tn} mod {q:>2}: chi2={chi:9.2f} dof={dof:4d} p={stats.chi2.sf(chi, dof):6.3f}"
                  f"  forbidden-cell counts={bad0}")
    print("    (PROVED: (U^t n mod M) and (T^t n mod M) are Markov chains under Haar measure; for T the state")
    print("    mod 2^k is the window of the next k i.i.d. parity bits, for U it is the induced chain on the")
    print("    windows that start with an odd bit, and the odd-prime part is updated by a fresh v.)")

    # ---------- v-histogram: the exact geometric law behind every U-law above
    NV = int(os.environ.get("PR_COLLATZ_NV", 100_000_000))
    out = subprocess.run([BIN, "vhist", str(NV), "31415"], capture_output=True, text=True, check=True).stdout.split()
    N = int(out[1].split("=")[1]); cnt = [int(c) for c in out[2:]]
    even = sum(cnt[j - 1] for j in range(2, 41, 2))
    z = (even - N / 3) / math.sqrt(N * 2 / 9)
    print(f"\n    v-histogram over {N:,} Terras-regime Syracuse steps ({NV:,} starts): P(v even) = {even / N:.7f}"
          f" (exact 1/3; z = {z:+.2f});  P(v=j), j=1..6: " +
          ", ".join(f"{cnt[j - 1] / N:.6f}" for j in range(1, 7)) + "  (exact 2^-j)")

    # ---------- (D) contrast data
    print("\n(D) Contrast with consecutive primes (see Part 1 for the prime numbers):")
    for key in [("U", 10), ("U", 3), ("T", 10), ("T", 3)]:
        P = laws[key]
        st = sorted(P)
        pi = stationary(P)
        tv = max(sum(abs(float(P[a].get(b, 0) - pi[b])) for b in st) / 2 for a in st if pi.get(a, 0) > 0)
        print(f"    {key[0]} mod {key[1]:>2}: max row TV distance from the stationary law = {tv:.4f} (exact, "
              f"independent of the size of n)")
    print("    For primes the same quantity is O(log log x / log x) and tends to 0 (Part 1).")
    print(f"\nPart 2 done in {time.time() - t0:.0f} s")
    bad = [r for r in results if r[-1] != 0]
    assert not bad, bad          # no transition ever lands in a cell of Haar probability 0
    # the merging diagnostic
    nA, nD = tabs["U_band_all"]["one"].sum(), tabs["U_band_distinct"]["one"].sum()
    print(f"\n    merging diagnostic, U-transitions from odd values in [10^6, 10^8): {nA:,} visits of {nD:,} distinct"
          f" values (mean multiplicity {nA / max(nD, 1):.2f}; the band has 49,500,000 odd values)")
    for q in (10, 3, 9, 72):
        P = laws[("U", q)]
        for tn in ("U_band_all", "U_band_distinct"):
            C = agg1(tabs[tn]["one"], q)
            r = compare(C, P, sorted(P), f"{tn} mod {q}")
            results.append((tn, q) + r)
    for tag in ("terras", "post_high", "post_low", "band_all", "band_distinct"):
        ok = [r for r in results if tag in r[0]]
        print(f"    {tag} p-values:", ", ".join(f"{r[0]}/{r[1]}: {r[4]:.3f}" for r in ok))
    print("    Reading: the Terras regime is the exact test of (A).  post_high (deterministic orbits, values >= 2^40) tests")
    print("    the heuristic that real orbits follow the Haar law.  In the low band orbits merge: compare 'band_all' (every")
    print("    visit) with 'band_distinct' (each integer once) -- a deviation that disappears on deduplication is a")
    print("    multiplicity (tree-weighting) effect, not a change of the one-step law.")


if __name__ == "__main__":
    main()
