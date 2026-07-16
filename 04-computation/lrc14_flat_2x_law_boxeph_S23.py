#!/usr/bin/env python3
"""
THE FLAT = 2x GOOD LAW  (boxeph-2026-07-16-S23)

Resolves the "2x identity" backlogged by mac-mini-S113 / handed around by kps cont.28:

    lambda(F) = 6617/97020 = 2 x 6617/194040 = 2 x m({1..12}; 1/14)

where F = flat set of klein's FT pair-energy budget
    S2(t) = sum_{m=1}^{12} (13-m) * (1/7 - ||m t||)_+   (floor 6/7, LEM-020/MISTAKE-152)
and m({1..12}; 1/14) = measure of the deep-well good set G = {t : ||mt|| >= 1/14, m<=12}
(THM-805 / THM-853(II): = H_12/91).

MECHANISM PROVED HERE (exact verification of every step):
 (1) F and G both live exactly on the 12 two-gap Farey-12 cells (a/i, b/j), i+j = 13.
     Three-gap cells (i+j >= 14) carry ZERO flat measure: by the successor rule the
     walk visits m -> m+i mod 13-ish only in the 2-gap case; in EVERY cell some
     alpha-gap is adjacent to a beta-gap, so "all gaps <= 1/7 & adjacent 2-sums >= 1/7"
     pins alpha+beta = 1/7 (a single t) unless the cell is two-gap.
 (2) On the two-gap cell with left denominator i < right denominator j (i+j=13):
     the 13 orbit gaps take two values alpha = i t - a (j copies), beta = b - j t
     (i copies), j*alpha + i*beta = 1; the circular gap word is the rotation coding
     w_r = [r*i mod 13 in {0..j-1}]: minority letter never self-adjacent, majority
     self-adjacent iff strict majority.  Hence
        F cap cell = { alpha in [1/14, (7-i)/(7(13-2i))] }   (closed)
        G cap cell = { alpha in [1/14, ...], beta >= 1/14 }
     with t-lengths  |F| = 1/(14 min(i,j) |j-i|),  |G| = 1/(14 i j),  and cl(G) c F.
 (3) Aggregate: sum over the 12 sites. The identity
        sum_{i unit mod 13, i<13/2} 1/(i(13-2i)) = 2 sum 1/(i(13-i)) = 2 H_12 / 13
     is the even/odd harmonic split of (Z/13)^* under u = 2v (the x(-2) permutation).
     => lambda(F) = 2 H_12/91 = 2 lambda(G).  The overhang F \ G has measure EXACTLY
     lambda(G): the flat set = good set + equal-measure overhang.
 (4) GENERAL N (points {t..Nt}, kernel threshold kappa = 2/(N+1), good threshold
     1/(N+1)): the law lambda(F_N) = 2 m({1..N-1}; 1/(N+1)) holds for EVERY ODD N
     and FAILS for every even N (even units of Z/N are empty).  LRC(n): N = n-1,
     so the 2x law is an EVEN-n phenomenon; n = 14 qualifies.

Pure Python / fractions. No numpy.
"""

from fractions import Fraction as Fr
from math import gcd

# ---------------------------------------------------------------- basics

def norm_dist(x):
    """||x|| = distance to nearest integer, exact."""
    f = x - int(x)
    if f < 0:
        f += 1
    return min(f, 1 - f)

def S2(t, N):
    """pair-energy budget: sum_{m=1}^{N-1} (N-m) (kappa - ||mt||)_+, kappa = 2/(N+1)."""
    kappa = Fr(2, N + 1)
    s = Fr(0)
    for m in range(1, N):
        d = norm_dist(m * t)
        if d < kappa:
            s += (N - m) * (kappa - d)
    return s

def flat_set(N):
    """Exact flat set of S2 on [0,1]: maximal runs where S2 == floor = (N-1)/(N+1).
    S2 is piecewise linear with kinks in {c/(2m)} u {(c +- kappa)/m}. Returns list of
    closed intervals [x,y] with y > x, plus the exact measure."""
    kappa = Fr(2, N + 1)
    kinks = set([Fr(0), Fr(1)])
    for m in range(1, N):
        for c in range(0, 2 * m + 1):
            kinks.add(Fr(c, 2 * m))
        for c in range(0, m + 1):
            for s in (+1, -1):
                x = Fr(c + s * kappa, m)
                if 0 <= x <= 1:
                    kinks.add(x)
    grid = sorted(kinks)
    floor = Fr(N - 1, N + 1)
    vals = [S2(x, N) for x in grid]
    assert min(vals) == floor, "floor violated!"
    ivs = []
    cur = None
    for idx in range(len(grid) - 1):
        # segment [grid[idx], grid[idx+1]] is flat iff both endpoints at floor
        if vals[idx] == floor and vals[idx + 1] == floor:
            if cur is None:
                cur = [grid[idx], grid[idx + 1]]
            else:
                cur[1] = grid[idx + 1]
        else:
            if cur is not None:
                ivs.append(tuple(cur)); cur = None
    if cur is not None:
        ivs.append(tuple(cur))
    meas = sum(y - x for x, y in ivs)
    return ivs, meas

def good_set(N):
    """G = {t : ||mt|| >= 1/(N+1), m=1..N-1}: exact intervals + measure."""
    lam = Fr(1, N + 1)
    bad = []
    for m in range(1, N):
        for c in range(0, m + 1):
            lo, hi = Fr(c, m) - lam / m, Fr(c, m) + lam / m
            bad.append((max(lo, Fr(0)), min(hi, Fr(1))))
    bad.sort()
    merged = []
    for lo, hi in bad:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    good = []
    prev = Fr(0)
    for lo, hi in merged:
        if lo > prev:
            good.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1:
        good.append((prev, Fr(1)))
    meas = sum(y - x for x, y in good)
    return good, meas

def harmonic(n):
    return sum(Fr(1, k) for k in range(1, n + 1))

# ---------------------------------------------------------------- N = 13 deep dive

def farey_adjacencies(n):
    """adjacent pairs (a/i, b/j) in Farey_n on [0,1]."""
    seq = sorted(set(Fr(a, q) for q in range(1, n + 1) for a in range(0, q + 1)))
    return [(seq[k], seq[k + 1]) for k in range(len(seq) - 1)]

def orbit_word(t, N):
    """sorted orbit {mt}, m=1..N: gap values in circular order + gap 'letters'.
    Returns (gaps list in circular order, set of distinct values)."""
    pts = sorted(set((m * t) - int(m * t) for m in range(1, N + 1)))
    gaps = [pts[k + 1] - pts[k] for k in range(len(pts) - 1)]
    gaps.append(pts[0] + 1 - pts[-1])
    return gaps

def run_N13():
    N = 13
    kappa = Fr(2, N + 1)          # 1/7
    lam = Fr(1, N + 1)            # 1/14
    print("=" * 72)
    print("PART 1 -- N = 13: exact censuses")
    F, mF = flat_set(N)
    G, mG = good_set(N)
    H12 = harmonic(12)
    print(f"flat set: {len(F)} intervals, measure {mF}  (claim 6617/97020: {mF == Fr(6617,97020)})")
    print(f"good set: {len(G)} intervals, measure {mG}  (claim 6617/194040: {mG == Fr(6617,194040)})")
    print(f"H_12/91 == mG: {H12/91 == mG};   mF == 2*mG: {mF == 2*mG}")

    print()
    print("PART 2 -- per-cell verification on Farey-12 sites with i+j = 13")
    sites = []
    for L, R in farey_adjacencies(12):
        i, j = L.denominator, R.denominator
        if i + j == 13:
            sites.append((L, R, i, j))
    print(f"sites found: {len(sites)} (expect 12)")
    okF = okG = okC = 0
    for (L, R, i, j) in sites:
        a = L.numerator
        # predicted F interval
        if i < j:
            lo = L + Fr(1, 14 * i)
            hi = L + Fr(7 - i, 7 * i * (13 - 2 * i))
        else:
            # mirror of the (j', i') = (j, i) site: F ends at R, starts symmetric
            hi = R - Fr(1, 14 * j)
            lo = R - Fr(7 - j, 7 * j * (13 - 2 * j))
        predF = (lo, hi)
        predLenF = Fr(1, 14 * min(i, j) * abs(j - i))
        # find the actual F interval inside this cell
        actF = [(x, y) for x, y in F if x >= L and y <= R]
        # allow endpoint-touching (flat interval may end exactly at a Farey point)
        match = len(actF) == 1 and actF[0][0] == predF[0] and actF[0][1] == predF[1]
        lenok = len(actF) == 1 and (actF[0][1] - actF[0][0]) == predLenF
        okF += match and lenok
        # predicted G interval (closed)
        gl, gh = L + lam / i, R - lam / j
        actG = [(x, y) for x, y in G if x >= L and y <= R]
        gok = len(actG) == 1 and actG[0] == (gl, gh) and (gh - gl) == Fr(1, 14 * i * j)
        okG += gok
        # containment cl(G) subset F
        cont = len(actF) == 1 and len(actG) == 1 and actF[0][0] <= actG[0][0] and actG[0][1] <= actF[0][1]
        okC += cont
        if not (match and lenok and gok and cont):
            print(f"  MISMATCH at site ({L},{R}) i={i} j={j}: F={actF} pred={predF}, G={actG}")
    print(f"F intervals match prediction: {okF}/12; G match: {okG}/12; cl(G) c F: {okC}/12")
    # total F mass ON sites = all of F?
    onsite = sum(min(y, R) - max(x, L) for x, y in F for (L, R, i, j) in sites if x < R and y > L)
    print(f"flat mass on i+j=13 sites: {onsite} == total {mF}: {onsite == mF}  (=> zero 3-gap flat mass)")

    print()
    print("PART 3 -- gap-word structure per cell (orbit of 13 points)")
    twog = threeg = 0
    ok_two = ok_ab = 0
    for (L, R) in farey_adjacencies(12):
        i, j = L.denominator, R.denominator
        t = (L + R) / 2
        if t == Fr(L.numerator + R.numerator, i + j):
            t = L + (R - L) * Fr(4, 9)   # avoid the mediant (alpha == beta there)
        gaps = orbit_word(t, N)
        vals = sorted(set(gaps))
        alpha = i * t - L.numerator      # ||i t||, the left-approach gap value
        beta = R.numerator - j * t       # ||j t||, the right-approach gap value
        assert alpha > 0 and beta > 0 and alpha != beta
        adj = set()
        for k in range(len(gaps)):
            adj.add(frozenset({gaps[k], gaps[(k + 1) % len(gaps)]}))
        if len(vals) == 2:
            twog += 1
            assert i + j == 13, f"two-gap cell with i+j={i+j}?"
            assert set([alpha, beta]) == set(vals)
            # counts: alpha has 13-i = j copies, beta has 13-j = i copies
            cok = (gaps.count(alpha) == j and gaps.count(beta) == i)
            maj, mino = (alpha, beta) if j > i else (beta, alpha)
            maj_self = frozenset({maj}) in adj
            min_self = frozenset({mino}) in adj
            ab = frozenset({alpha, beta}) in adj
            if maj_self and (not min_self) and ab and cok:
                ok_two += 1
            else:
                print(f"  word anomaly at 2-gap cell ({L},{R}): majself={maj_self} minself={min_self} ab={ab} counts_ok={cok}")
        else:
            threeg += 1
            assert vals == sorted([alpha, beta, alpha + beta]), f"cell ({L},{R}): vals {vals}"
            # alpha-beta adjacency claim (the load-bearing Lemma B)
            if frozenset({alpha, beta}) in adj:
                ok_ab += 1
            else:
                print(f"  NO alpha-beta adjacency in 3-gap cell ({L},{R}) i={i} j={j}")
            # multiplicities (13-i, 13-j, i+j-13) attached to (alpha, beta, alpha+beta)
            c = [gaps.count(alpha), gaps.count(beta), gaps.count(alpha + beta)]
            assert c == [13 - i, 13 - j, i + j - 13], f"multiplicity fail ({L},{R}): {c}"
    print(f"two-gap cells: {twog} (expect 12), word structure ok: {ok_two}/12")
    print(f"three-gap cells: {threeg}, alpha-beta adjacency present: {ok_ab}/{threeg}")

    print()
    print("PART 4 -- overhang anatomy: F = cl(G) + equal-measure overhang")
    over = mF - mG
    print(f"overhang measure = {over} == lambda(G) = {mG}: {over == mG}")
    # per-site ratio law |F|/|G| = max(i,j)/|i-j|
    print("per-site ratio |F cap cell| / |G cap cell| = max(i,j)/|j-i|:")
    for (L, R, i, j) in sites:
        fint = [(x, y) for x, y in F if x >= L and y <= R][0]
        gint = [(x, y) for x, y in G if x >= L and y <= R][0]
        ratio = (fint[1] - fint[0]) / (gint[1] - gint[0])
        assert ratio == Fr(max(i, j), abs(j - i))
    print("  all 12 sites: ratio == max(i,j)/|j-i|  [7, 8/3, 9/5, 10/7, 11/9, 12/11 x2]  OK")

def run_generalN():
    print()
    print("=" * 72)
    print("PART 5 -- the general-N law: flat = 2 x good iff N odd")
    print(f"{'N':>3} {'n=N+1':>6} {'flat':>22} {'good':>22} {'ratio':>10} {'=2?':>4}")
    for N in range(3, 18):
        F, mF = flat_set(N)
        G, mG = good_set(N)
        ratio = mF / mG if mG else None
        print(f"{N:>3} {N+1:>6} {str(mF):>22} {str(mG):>22} {str(ratio):>10} {'YES' if ratio == 2 else 'no':>4}")

def run_harmonic_identities():
    print()
    print("=" * 72)
    print("PART 6 -- the unit-group harmonic identities (the mechanism)")
    print("odd N: sum_{i unit < N/2} 1/(i(N-2i)) == 2 sum_{i unit < N/2} 1/(i(N-i)) == 2 H-mass")
    for N in range(3, 26):
        units = [i for i in range(1, (N + 1) // 2) if gcd(i, N) == 1 and N - 2 * i != 0]
        lhs = sum(Fr(1, i * (N - 2 * i)) for i in units)
        rhs = 2 * sum(Fr(1, i * (N - i)) for i in units)
        eq = lhs == rhs
        par = "odd " if N % 2 else "EVEN"
        print(f"  N={N:>2} ({par}): LHS={str(lhs):>18} RHS={str(rhs):>18} equal={eq}")
        if N % 2 == 1:
            assert eq, f"odd-N identity FAILS at N={N}"
        else:
            assert not eq, f"even-N identity unexpectedly holds at N={N}"

def run_harmonic_law_correction():
    print()
    print("=" * 72)
    print("PART 7 -- independent per-site re-derivation of THM-819's prime criterion")
    print("(THM-819, kps cont.9, already proved this: m({1..k};1/(k+2)) =")
    print(" (2/((k+1)(k+2))) sum_{units} 1/u; = H_k/C(k+2,2) iff k+1 prime.")
    print(" This part re-verifies it through n = 28 via the Farey-cell route,")
    print(" including the odd-composite-coprime-to-6 station n = 26.)")
    print("per-site closed form (all N): m({1..N-1}; 1/(N+1)) =")
    print("    (2/(N+1)) * sum_{i in (Z/N)^*, i<N/2} 1/(i(N-i))")
    for N in range(3, 28):
        n = N + 1
        G, mG = good_set(N)
        pred = Fr(2, N + 1) * sum(Fr(1, i * (N - i)) for i in range(1, (N + 1) // 2)
                                  if gcd(i, N) == 1)
        law = harmonic(n - 2) / (Fr(n * (n - 1), 2))
        isprime = N > 1 and all(N % p for p in range(2, int(N ** 0.5) + 1))
        status = []
        assert mG == pred, f"universal closed form FAILS at N={N}: {mG} vs {pred}"
        if n % 2 == 0:
            holds = (mG == law)
            assert holds == isprime, f"prime criterion fails at n={n}"
            status.append(f"harmonic law {'HOLDS' if holds else 'fails'} (N={N} {'prime' if isprime else 'composite'})")
            if n == 26:
                status.append("<< the predicted HYP-6885 counterexample (26 = 2 mod 6)")
        print(f"  n={n:>2}: m={str(mG):>20} closed-form OK; {' '.join(status)}")

def run_doubling_law():
    """PART 8 (S24 addendum, joint with death-star HYP-7013): THE DOUBLING LAW.
    For every ODD N: F_N = 2*G_N mod 1 EXACTLY, componentwise; the induced site
    permutation is u -> u * 2^{-1} on ordered (Z/N)^* labels (label u = left
    denominator i; u even: image label u/2; u odd: image label (N+u)/2).
    Verified endpoint-exactly per site. Even N: the label map needs 2
    invertible mod N -- fails, and F != 2G (checked)."""
    print()
    print("=" * 72)
    print("PART 8 -- THE DOUBLING LAW: F = 2G mod 1 (odd N), site map u -> u/2")
    for N in range(3, 18):
        F, mF = flat_set(N)
        G, mG = good_set(N)
        # double each G component, fold mod 1, merge
        d = []
        for lo, hi in G:
            a, b = 2 * lo, 2 * hi
            if b <= 1:
                d.append((a, b))
            elif a >= 1:
                d.append((a - 1, b - 1))
            else:
                d.append((a, Fr(1))); d.append((Fr(0), b - 1))
        d.sort()
        m = []
        for lo, hi in d:
            if m and lo <= m[-1][1]:
                m[-1][1] = max(m[-1][1], hi)
            else:
                m.append([lo, hi])
        m = [(a, b) for a, b in m]
        eq = (m == F)
        if N % 2 == 1:
            assert eq, f"doubling law FAILS at odd N={N}"
            # per-site label map: G component at site with left denom i maps to
            # F component at left denom i/2 (i even) or (N+i)/2 (i odd)
            ok_sites = 0
            adj = farey_adjacencies(N - 1)
            sites = [(L, R) for (L, R) in adj if L.denominator + R.denominator == N]
            kap = Fr(2, N + 1)

            def Fpiece(L, R):
                """the proved per-cell flat interval at cell (L, R), i+j=N"""
                i, j = L.denominator, R.denominator
                if i < j:
                    return (L + kap / (2 * i),
                            L + (1 - i * kap) / (i * (j - i)))
                else:
                    return (R - (1 - j * kap) / (j * (i - j)),
                            R - kap / (2 * j))

            piece_at_label = {}
            for (L, R) in sites:
                piece_at_label[L.denominator] = Fpiece(L, R)
            for (L, R) in sites:
                i, j = L.denominator, R.denominator
                lam = Fr(1, N + 1)
                gcomp = (L + lam / i, R - lam / j)
                x, y = 2 * gcomp[0], 2 * gcomp[1]
                k = int(x)
                x, y = x - k, y - k
                target_label = i // 2 if i % 2 == 0 else (N + i) // 2
                okc = piece_at_label.get(target_label) == (x, y)
                ok_sites += okc
                if not okc:
                    print(f"    site i={i} -> label {target_label}: doubled "
                          f"({x},{y}) vs piece {piece_at_label.get(target_label)}")
            # the label map must be a bijection on labels
            labels = sorted(piece_at_label)
            image = sorted((i // 2 if i % 2 == 0 else (N + i) // 2) for i in labels)
            print(f"  N={N:>2} (odd):  F == 2G: {eq};  per-site doubled-G == F-piece "
                  f"at label u*2^-1: {ok_sites}/{len(sites)}; label map bijective: "
                  f"{image == labels}")
            assert ok_sites == len(sites) and image == labels
        else:
            print(f"  N={N:>2} (EVEN): F == 2G: {eq} (must be False)")
            assert not eq

if __name__ == "__main__":
    run_N13()
    run_generalN()
    run_harmonic_identities()
    run_harmonic_law_correction()
    run_doubling_law()
    print()
    print("ALL CHECKS PASSED")
