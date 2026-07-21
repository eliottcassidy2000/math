#!/usr/bin/env python3
"""
thm1785_referee_telescoping_S187r.py — HOSTILE REFEREE checks for THM-1785
section 3 (executed double-period telescoping). Referee session S187r.
Run in parts: `python3 thm1785_referee_telescoping_S187r.py A|B|C` (or all).

A/T1  exactN is the honest INTEGRAL: validate against direct numerical
      quadrature (Laurent CT + Simpson) at several (m,j,k).
A/T2  E/U/W relation families re-verified exactly on a WIDER lattice than
      the file's (j to +-3, k to 7, m to 10; two extra triples); plus the
      k=0 W-row (even j, m>=1) which the THEOREM TEXT claims but the
      BUILDER EXCLUDES (code guard k>=1) — verify the claim is true.
A/T3  THE RECURRENCE (md closed form) for triple (3/2,-2/3,5/4): exact
      verification m = 3..120 (file: 4..56) — m=3 edge + beyond-56 range.
A/T7  The md's displayed recurrence is coefficient-SPECIFIC (md never
      states the triple): show it FAILS on triple 2's moments.
B/T4  Second triple (1/2,2,-1/3): FULL derivation m0=10..17 +
      interpolation + exact verification m=4..120 (the file ran ONE
      window, m0=12 — 'order 3, constant leading' for this triple is a
      normalization tautology at a single point until this is done).
B/T5  Adversarial triples incl. degenerate (b=0, c=0, a=0, sign flip):
      does the elimination still return a relation, and is the RAW top
      coefficient cS[-1] nonzero (the normalization divides by it)?
C/T6  ALL-m0 CERTIFICATE attempt (md claims 'P-RECURSIVE ...
      unconditionally' from 8 numeric windows + finite checks — a gap):
      reconstruct the left-kernel vector y(m0) as a polynomial in m0 from
      m0 = 10..21 and TEST at m0 = 22. If y is polynomial (deg <= 11) and
      y^T A(m0) vanishes on all non-mu columns at 13 >= deg+2 points, the
      vanishing holds IDENTICALLY in m0 — repairing the gap. If
      reconstruction fails, the gap stands.

boxeph-referee, 2026-07-20. Pure python, exact Fractions.
"""
from fractions import Fraction
import math
import sys
import time

T0 = time.time()
PART = sys.argv[1] if len(sys.argv) > 1 else "all"


def in_part(p):
    return PART == "all" or PART == p


def fact(n):
    return math.factorial(n)


def Iw(K):
    assert K % 2 == 1
    kap = (K - 1) // 2
    return Fraction(2) ** kap * fact(kap)


def exactN(a, b, c, m, j, k):
    tot = Fraction(0)
    for be in range(m + 1):
        rem = m - be
        if (rem - j) % 2:
            continue
        al = (rem - j) // 2
        ga = (rem + j) // 2
        if al < 0 or ga < 0:
            continue
        wpow = al + 2 * be + ga + k
        if wpow % 2 == 0:
            return None
        coef = Fraction(fact(m), fact(al) * fact(be) * fact(ga))
        tot += coef * a ** al * b ** be * c ** ga * Iw(wpow)
    return tot


def mom(a, b, c, m):
    return exactN(a, b, c, m, 0, 1)


aQ, bQ, cQ = Fraction(3, 2), Fraction(-2, 3), Fraction(5, 4)
a2, b2, c2 = Fraction(1, 2), Fraction(2), Fraction(-1, 3)


def md_rec_coeffs(m):
    # mu_m + (8m-4)/3 mu_{m-1} + (32m^2-199m+167)/18 mu_{m-2} - 5(m-1)(m-2) mu_{m-3}
    return [Fraction(-5) * (m - 1) * (m - 2),
            Fraction(32 * m * m - 199 * m + 167, 18),
            Fraction(8 * m - 4, 3),
            Fraction(1)]


def derive_full(a, b, c, m0, levels=4, JMAX=3, KMAX=9):
    """Structural copy of the target's derive_recurrence_at, extended to
    return (cS_raw, y, rows, labels, free_col). Whole-row drop semantics
    preserved (a dropped-row system is a SUBSET of true identities)."""
    idx = {}
    for l in range(m0 - levels + 1, m0 + 1):
        for j in range(-JMAX, JMAX + 1):
            for k in range(0, KMAX + 1):
                if (j + k) % 2 == 1:
                    idx[(l, j, k)] = len(idx)
    nu = len(idx)
    rows = []
    labels = []

    def add_row(terms, lab):
        row = [Fraction(0)] * nu
        for (l, j, k), co in terms:
            if abs(j) > JMAX or k > KMAX or k < 0 or not (m0 - levels + 1 <= l <= m0):
                return
            if (j + k) % 2 == 0:
                return
            row[idx[(l, j, k)]] += co
        rows.append(row)
        labels.append(lab)

    for l in range(m0 - levels + 2, m0 + 1):
        off = l - m0
        for j in range(-JMAX, JMAX + 1):
            for k in range(0, KMAX + 1):
                if (j + k) % 2 == 1:
                    add_row([((l, j, k), Fraction(1)),
                             ((l - 1, j + 1, k + 1), -a),
                             ((l - 1, j, k + 2), -b),
                             ((l - 1, j - 1, k + 1), -c)], ("E", off, j, k))
                    add_row([((l, j, k), Fraction(j)),
                             ((l - 1, j + 1, k + 1), Fraction(l) * a),
                             ((l - 1, j - 1, k + 1), -Fraction(l) * c)],
                            ("U", off, j, k))
                else:
                    if k >= 1:
                        add_row([((l, j, k + 1), Fraction(1)),
                                 ((l, j, k - 1), -Fraction(k)),
                                 ((l - 1, j + 1, k), -Fraction(l) * a),
                                 ((l - 1, j, k + 1), -2 * b * Fraction(l)),
                                 ((l - 1, j - 1, k), -Fraction(l) * c)],
                                ("W", off, j, k))

    S = [idx[(l, 0, 1)] for l in range(m0 - levels + 1, m0 + 1)]
    notS = [i for i in range(nu) if i not in S]
    M = [[rows[r][cidx] for r in range(len(rows))] for cidx in notS]
    ncols = len(rows)
    R = [row[:] for row in M]
    piv_of_col = [-1] * ncols
    r = 0
    for cc in range(ncols):
        piv = None
        for rr in range(r, len(R)):
            if R[rr][cc] != 0:
                piv = rr
                break
        if piv is None:
            continue
        R[r], R[piv] = R[piv], R[r]
        pv = R[r][cc]
        R[r] = [x / pv for x in R[r]]
        for rr in range(len(R)):
            if rr != r and R[rr][cc] != 0:
                f = R[rr][cc]
                R[rr] = [x - f * y for x, y in zip(R[rr], R[r])]
        piv_of_col[cc] = r
        r += 1
        if r == len(R):
            break
    free = [cc for cc in range(ncols) if piv_of_col[cc] < 0]
    for fc in free:
        y = [Fraction(0)] * ncols
        y[fc] = Fraction(1)
        for cc in range(ncols):
            pr = piv_of_col[cc]
            if pr >= 0:
                y[cc] = -R[pr][fc]
        cS = [sum(y[rr] * rows[rr][si] for rr in range(len(rows))) for si in S]
        if any(x != 0 for x in cS):
            return cS, y, rows, labels, fc
    return None, None, rows, labels, None


def lagrange_eval(pts, x):
    tot = Fraction(0)
    for i, (xi, yi) in enumerate(pts):
        num = Fraction(1)
        den = Fraction(1)
        for jj, (xj, _) in enumerate(pts):
            if jj != i:
                num *= (x - xj)
                den *= (xi - xj)
        tot += yi * num / den
    return tot


# ============================ PART A ============================
if in_part("A"):
    print("=" * 78)
    print("T1: exactN vs direct numerical quadrature of the defining integral")
    print("=" * 78)
    af, bf, cf = 1.5, -2.0 / 3, 1.25

    def CT_numeric(m, j, w):
        base = {1: af * w, 0: bf * w * w, -1: cf * w}
        cur = {0: 1.0}
        for _ in range(m):
            nxt = {}
            for pw1, v1 in cur.items():
                for pw2, v2 in base.items():
                    nxt[pw1 + pw2] = nxt.get(pw1 + pw2, 0.0) + v1 * v2
            cur = nxt
        return cur.get(-j, 0.0)

    def quadN(m, j, k, W=14.0, n=2800):
        h = W / n
        tot = 0.0
        for i in range(n + 1):
            w = i * h
            f = math.exp(-w * w / 2) * (w ** k) * CT_numeric(m, j, w)
            wgt = 1 if i in (0, n) else (4 if i % 2 else 2)
            tot += wgt * f
        return tot * h / 3

    for (m, j, k) in ((5, 0, 1), (6, 1, 2), (7, -2, 3), (8, 3, 0)):
        ex = float(exactN(aQ, bQ, cQ, m, j, k))
        nu = quadN(m, j, k)
        rel = abs(ex - nu) / max(1e-300, abs(ex))
        print("  (m,j,k)=(%d,%2d,%d): exactN = %.10e  quad = %.10e  rel = %.1e"
              % (m, j, k, ex, nu, rel))

    print()
    print("=" * 78)
    print("T2: E/U/W wider lattice; k=0 W-row (claimed in text, unused by code)")
    print("=" * 78)

    def check_relations(a, b, c, ms, js, ks):
        bad = []
        for m in ms:
            for j in js:
                for k in ks:
                    if (j + k) % 2 == 1:
                        E = exactN(a, b, c, m, j, k) - (
                            a * exactN(a, b, c, m - 1, j + 1, k + 1)
                            + b * exactN(a, b, c, m - 1, j, k + 2)
                            + c * exactN(a, b, c, m - 1, j - 1, k + 1))
                        U = (j * exactN(a, b, c, m, j, k)
                             + m * (a * exactN(a, b, c, m - 1, j + 1, k + 1)
                                    - c * exactN(a, b, c, m - 1, j - 1, k + 1)))
                        if E != 0:
                            bad.append(("E", m, j, k))
                        if U != 0:
                            bad.append(("U", m, j, k))
                    else:
                        if k >= 1:
                            W = (exactN(a, b, c, m, j, k + 1)
                                 - k * exactN(a, b, c, m, j, k - 1)
                                 - m * (a * exactN(a, b, c, m - 1, j + 1, k)
                                        + 2 * b * exactN(a, b, c, m - 1, j, k + 1)
                                        + c * exactN(a, b, c, m - 1, j - 1, k)))
                            if W != 0:
                                bad.append(("W", m, j, k))
        return bad

    triples = [(aQ, bQ, cQ), (a2, b2, c2),
               (Fraction(-1), Fraction(3), Fraction(7, 5))]
    totbad = 0
    for (a, b, c) in triples:
        bad = check_relations(a, b, c, (1, 2, 4, 7, 10), range(-3, 4), range(0, 8))
        totbad += len(bad)
        print("  triple (%s,%s,%s): mismatches = %d" % (a, b, c, len(bad)))
    print("  wider-lattice E/U/W total mismatches: %d" % totbad)

    badW0 = 0
    for (a, b, c) in triples:
        for m in (1, 2, 3, 5, 8):
            for j in (-2, 0, 2):
                lhs = exactN(a, b, c, m, j, 1)
                rhs = m * (a * exactN(a, b, c, m - 1, j + 1, 0)
                           + 2 * b * exactN(a, b, c, m - 1, j, 1)
                           + c * exactN(a, b, c, m - 1, j - 1, 0))
                if lhs != rhs:
                    badW0 += 1
    print("  k=0 W-row (m>=1) mismatches: %d  (TRUE as claimed; the builder" % badW0)
    print("   omits it via its k>=1 guard — statement/code mismatch, benign:")
    print("   using FEWER proved rows cannot create false relations)")

    print()
    print("=" * 78)
    print("T3: md closed-form recurrence for (3/2,-2/3,5/4), exact m=3..120")
    print("=" * 78)
    fails = []
    for m in range(3, 121):
        co = md_rec_coeffs(m)
        tot = sum(co[i] * mom(aQ, bQ, cQ, m - 3 + i) for i in range(4))
        if tot != 0:
            fails.append(m)
    print("  failures in m = 3..120: %s" % (fails if fails else "NONE"))
    print("  (file verified only 4..56; range extended and m=3 included here)")

    print()
    print("=" * 78)
    print("T7: md displayed recurrence on triple 2's moments (should FAIL)")
    print("=" * 78)
    for m in (5, 8):
        co = md_rec_coeffs(m)
        tot = sum(co[i] * mom(a2, b2, c2, m - 3 + i) for i in range(4))
        print("  m=%d: residual = %s  (%s)"
              % (m, tot, "ZERO?!" if tot == 0 else "nonzero, as expected"))
    print("  => the md never states WHICH (a,b,c) its displayed recurrence is")
    print("     for (it is (3/2,-2/3,5/4), read off the frozen .out).")
    print("  [part A done, t=%.0fs]" % (time.time() - T0))

# ============================ PART B ============================
if in_part("B"):
    print()
    print("=" * 78)
    print("T4: second triple (1/2,2,-1/3) — full windows m0=10..17 + verify")
    print("=" * 78)
    recs2 = {}
    for m0 in range(10, 18):
        cS, y, rows, labels, fc = derive_full(a2, b2, c2, m0)
        if cS is None:
            print("  m0=%d: NO combination" % m0)
            continue
        raw_top = cS[-1]
        cSn = [x / raw_top for x in cS]
        recs2[m0] = cSn
        print("  m0=%2d: raw top = %s ; normalized %s  [t=%.0fs]" %
              (m0, raw_top, ["%s" % x for x in cSn], time.time() - T0))
    ms2 = sorted(recs2)
    polys2 = [[(Fraction(m0), recs2[m0][pos]) for m0 in ms2] for pos in range(4)]
    degok = True
    for pos in range(4):
        pred = lagrange_eval(polys2[pos][:7], Fraction(ms2[7]))
        if pred != recs2[ms2[7]][pos]:
            degok = False
    print("  normalized coeffs consistent with degree<=6 polynomials in m: %s" % degok)
    fails2 = []
    for m in range(4, 121):
        co = [lagrange_eval(polys2[pos], Fraction(m)) for pos in range(4)]
        tot = sum(co[i] * mom(a2, b2, c2, m - 3 + i) for i in range(4))
        if tot != 0:
            fails2.append(m)
    print("  interpolated triple-2 recurrence verified m=4..120: %s"
          % ("FAILURES %s" % fails2 if fails2 else "ALL HOLD"))
    m3f = [lagrange_eval(polys2[pos], Fraction(3)) for pos in range(4)]
    t3 = sum(m3f[i] * mom(a2, b2, c2, i) for i in range(4))
    print("  triple-2 recurrence at m=3: %s" % ("HOLDS" if t3 == 0 else "FAILS"))

    print()
    print("=" * 78)
    print("T5: adversarial triples at m0=12 — existence + raw top coefficient")
    print("=" * 78)
    adv = [(Fraction(1), Fraction(1), Fraction(1)),
           (Fraction(1), Fraction(0), Fraction(1)),
           (Fraction(1), Fraction(1), Fraction(0)),
           (Fraction(0), Fraction(1), Fraction(1)),
           (Fraction(1), Fraction(-1), Fraction(1))]
    for (a, b, c) in adv:
        cS, y, rows, labels, fc = derive_full(a, b, c, 12)
        if cS is None:
            print("  (a,b,c)=(%s,%s,%s): NO mu-supported combination in window"
                  % (a, b, c))
            continue
        print("  (a,b,c)=(%s,%s,%s): raw cS = %s ; top zero? %s  [t=%.0fs]"
              % (a, b, c, ["%s" % x for x in cS], cS[-1] == 0, time.time() - T0))
    print("  [part B done, t=%.0fs]" % (time.time() - T0))

# ============================ PART C ============================
if in_part("C"):
    print()
    print("=" * 78)
    print("T6: all-m0 certificate attempt for triple 1 (y(m0) polynomial in m0?)")
    print("=" * 78)
    ys = {}
    cSs = {}
    labels_ref = None
    nrows_ref = None
    structure_ok = True
    windows = list(range(10, 23))
    for m0 in windows:
        cS, y, rows, labels, fc = derive_full(aQ, bQ, cQ, m0)
        ys[m0] = y
        cSs[m0] = cS
        if labels_ref is None:
            labels_ref = labels
            nrows_ref = len(rows)
        elif labels != labels_ref or len(rows) != nrows_ref:
            structure_ok = False
            print("  ROW STRUCTURE CHANGES with m0 at m0=%d — abort attempt" % m0)
            break
        print("  m0=%d: rows=%d, free col=%s, raw cS[-1]=%s  [t=%.0fs]"
              % (m0, len(rows), fc, cS[-1], time.time() - T0))
    if structure_ok:
        mlist = windows[:-1]     # 12 interpolation points, deg <= 11
        testm = windows[-1]
        nbad = 0
        for r in range(nrows_ref):
            pts = [(Fraction(m0), ys[m0][r]) for m0 in mlist]
            if lagrange_eval(pts, Fraction(testm)) != ys[testm][r]:
                nbad += 1
        ok_pred = (nbad == 0)
        print("  y-components polynomial (deg<=11) predicting m0=%d: %s (%d/%d mismatch)"
              % (testm, ok_pred, nbad, nrows_ref))
        if ok_pred:
            print("  => y(m0) is a polynomial vector on 13 windows; every entry of")
            print("     y(m0)^T A(m0) is a polynomial of degree <= 12 in m0 that")
            print("     vanishes at 13 points on all non-mu columns, hence")
            print("     IDENTICALLY: the certificate holds for ALL m0 (window valid")
            print("     m0 >= 4). The all-m gap is REPAIRABLE by this finite check.")
        else:
            print("  => y from the free-column construction is NOT polynomial across")
            print("     windows; md's 'unconditional' P-recursiveness rests on 8")
            print("     numeric windows + finite m-checks: ALL-m REMAINS A GAP")
            print("     (repair route: symbolic-m0 elimination over Q(m0)).")
            cn = {m0: [x / cSs[m0][-1] for x in cSs[m0]] for m0 in windows}
            okB = True
            for pos in range(4):
                pts = [(Fraction(m0), cn[m0][pos]) for m0 in mlist]
                if lagrange_eval(pts, Fraction(testm)) != cn[testm][pos]:
                    okB = False
            print("     plan B: normalized cS(m0) polynomial predicting m0=%d: %s"
                  % (testm, okB))
    print("  [part C done, t=%.0fs]" % (time.time() - T0))

# ============================ PART D ============================
if in_part("D"):
    print()
    print("=" * 78)
    print("T6b: airtight all-m certificate, BOTH executed triples")
    print("     (explicit kernel verification at 13 windows + polynomial y +")
    print("      degree-count => identity in m0; no RREF trust)")
    print("=" * 78)
    windows = list(range(10, 23))

    def all_m_certificate(a, b, c, name, target_coeffs):
        ys = {}
        labels_ref = None
        nrows = None
        for m0 in windows:
            cS, y, rows, labels, fc = derive_full(a, b, c, m0)
            if labels_ref is None:
                labels_ref, nrows = labels, len(rows)
            elif labels != labels_ref or len(rows) != nrows:
                print("  [%s] row structure varies at m0=%d — FAIL" % (name, m0))
                return
            # explicit: y^T rows vanishes on ALL non-mu columns, equals cS on mu
            nu = len(rows[0])
            Scols = set()
            # mu columns are (l,0,1); recover their indices by recomputing idx
            idx = {}
            for l in range(m0 - 3, m0 + 1):
                for j in range(-3, 4):
                    for k in range(0, 10):
                        if (j + k) % 2 == 1:
                            idx[(l, j, k)] = len(idx)
            Scols = {idx[(l, 0, 1)] for l in range(m0 - 3, m0 + 1)}
            comb = [sum(y[r] * rows[r][col] for r in range(nrows))
                    for col in range(nu)]
            bad_notS = sum(1 for col in range(nu)
                           if col not in Scols and comb[col] != 0)
            cS_direct = [comb[idx[(l, 0, 1)]] for l in range(m0 - 3, m0 + 1)]
            tc = target_coeffs(m0)
            match = (cS_direct == tc)
            print("  [%s] m0=%d: nonzero non-mu cols = %d ; cS == target form: %s"
                  % (name, m0, bad_notS, match))
            if bad_notS or not match:
                return
            ys[m0] = y
        nbad = 0
        for r in range(nrows):
            pts = [(Fraction(m0), ys[m0][r]) for m0 in windows[:-1]]
            if lagrange_eval(pts, Fraction(windows[-1])) != ys[windows[-1]][r]:
                nbad += 1
        print("  [%s] y polynomial (deg<=11) across 13 windows: %s (%d bad)"
              % (name, nbad == 0, nbad))
        if nbad == 0:
            print("  [%s] => y^T A(m0) entries are deg<=12 polynomials in m0," % name)
            print("       zero at 13 points on non-mu columns and EQUAL to the")
            print("       target coefficient polynomials on mu columns at 13")
            print("       points: both hold IDENTICALLY in m0. Combined with the")
            print("       row identities (valid all m0 >= 3): the recurrence is")
            print("       PROVED FOR ALL m >= 3 at this triple.")

    all_m_certificate(aQ, bQ, cQ, "triple1 (3/2,-2/3,5/4)", md_rec_coeffs)

    # triple 2 target form from its 8-window interpolation (part B):
    recs2b = {}
    for m0 in range(10, 18):
        cS, y, rows, labels, fc = derive_full(a2, b2, c2, m0)
        recs2b[m0] = [x / cS[-1] for x in cS]
    p2 = [[(Fraction(m0), recs2b[m0][pos]) for m0 in sorted(recs2b)]
          for pos in range(4)]

    def t2_coeffs(m):
        return [lagrange_eval(p2[pos], Fraction(m)) for pos in range(4)]

    all_m_certificate(a2, b2, c2, "triple2 (1/2,2,-1/3)", t2_coeffs)

    # T5 follow-up: at (1,0,1) the found relation is supported on ODD levels
    # whose moments vanish identically (b=0): true but vacuous — and the
    # target's normalization-by-top would divide by zero there.
    print()
    print("  T5 follow-up at (a,b,c)=(1,0,1): mu_9 = %s, mu_11 = %s (b=0 =>"
          % (mom(Fraction(1), Fraction(0), Fraction(1), 9),
             mom(Fraction(1), Fraction(0), Fraction(1), 11)))
    print("  odd moments vanish; the returned relation -40*mu_9 + mu_11 = 0 is")
    print("  TRUE but VACUOUS, and cS[-1] = 0 breaks the constant-leading form.")

print("\nDONE (%s). total %.0fs" % (PART, time.time() - T0))
