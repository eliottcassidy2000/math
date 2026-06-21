"""
lrc14_paley_simultaneous_moment_floor_angle6.py

ANGLE 6 follow-up. The 1st test REFUTED "min sum|lambda|^4 <=> max H" (non-monotone at p=11),
but CONFIRMED flat-spectrum = QR difference set = global argmin of sum|lambda|^4.

NEW SHARPER HYPOTHESIS (falsifiable):

  Paley does NOT win by minimizing a single moment; it wins because the FLAT spectrum
  SIMULTANEOUSLY MINIMIZES *every* even moment  M_{2k} = sum_j |lambda_j|^{2k}  for all k>=2,
  among all tournament connection sets on Z_p.  (Schur-convexity: the constant vector
  majorization-minimizes every symmetric convex function under a fixed sum constraint.
  Parseval fixes M_2 = sum|lambda|^2 = constant across ALL tournament sets, so the |lambda_j|^2
  vector lies on a fixed-sum simplex; flat = the majorization-MINIMAL point => minimizes EVERY
  Schur-convex M_{2k} simultaneously.)

  Since H = 1 + 2N + (alpha-corrections) and the cycle counts N, c_k are spectral moment
  expressions where the *dominant* contributions are Schur-CONCAVE in |lambda|^2 (i.e. the
  e_k elementary-symmetric basis, which is MAXIMIZED at the flat point), Paley maximizes H.

  TEST 1: Is M_{2k}(Paley) = min over all sets, for k = 2,3,4,5 simultaneously? (p=7,11,19)
  TEST 2: Is the |lambda|^2 vector at Paley MAJORIZED BY every other set's |lambda|^2 vector?
          (the master Schur statement -- if true, EVERY Schur-convex functional is min at Paley
           and EVERY Schur-concave functional is MAX at Paley, no case analysis needed.)
  TEST 3: Express H in the elementary-symmetric basis e_k(|lambda|^2) and check that, although
          the single-moment (power-sum) story fails, the e_k story (H increasing in each e_k,
          all e_k max at Paley) explains H-max. Concretely: is H a function with all-nonneg
          coefficients in (e_2,...,e_m) MINUS small terms, and does Paley dominate every set in
          e_k for ALL k (so any positive-combination must be maxed at Paley)?
  TEST 4 (the crux of the wall): even though H is NOT monotone in M_4, is H monotone in
          MAJORIZATION ORDER?  i.e. if set A's |lambda|^2 vector majorizes set B's, is H(A) <= H(B)?
          If YES, then "H is Schur-concave in |lambda|^2" and Paley (the majorization-minimum)
          is the unique max -- a GLOBAL proof. If NO, find the counterexample pair (the true
          obstruction).

Run: python3 04-computation/lrc14_paley_simultaneous_moment_floor_angle6.py
"""
import sys
from fractions import Fraction
from itertools import combinations

sys.stdout.reconfigure(line_buffering=True)


def qr_set(p):
    return frozenset((a * a) % p for a in range(1, p))


def tournament_connection_sets(p):
    m = (p - 1) // 2
    pairs = [(t, (p - t) % p) for t in range(1, m + 1)]
    sets = []
    for bits in range(1 << m):
        S = set()
        for i, (a, b) in enumerate(pairs):
            S.add(a if (bits >> i & 1) else b)
        sets.append(frozenset(S))
    return sets


def lambda_sq_vector_exact(S, p):
    """|lambda_k|^2 for k=1..p-1, computed EXACTLY.
       lambda_k = sum_{s in S} omega^{ks}, omega = e^{2pi i /p}.
       |lambda_k|^2 = sum_{a,b in S} cos(2pi k (a-b)/p).
       We use the autocorrelation: |lambda_k|^2 = sum_t c(t) * cos(2pi k t / p).
       For EXACT majorization we instead use the integer power-sum trick:
       M_{2k} = sum_j |lambda_j|^2k is an integer (symmetric in roots of unity).
       For the |lambda|^2 VECTOR (needed for majorization), |lambda_k|^2 are algebraic;
       we compute them as high-precision floats but ALSO compute integer moments exactly.
    """
    import cmath, math
    omega = cmath.exp(2j * math.pi / p)
    vec = []
    for k in range(1, p):
        lam = sum(omega ** (k * s) for s in S)
        vec.append(abs(lam) ** 2)
    return vec


def moment_2k_exact(S, p, k):
    """M_{2k} = sum_{j=1}^{p-1} |lambda_j|^{2k}  EXACTLY (integer).
       |lambda_j|^{2k} for the FULL set j=0..p-1 minus j=0 term.
       Use: sum_{j=0}^{p-1} |lambda_j|^{2k} = p * (number of solutions to a closed walk).
       Concretely |lambda_j|^2 = sum_t c(t) omega^{jt}, so
       sum_j |lambda_j|^{2k} = sum_j (sum_t c(t) omega^{jt})^k.
       Expanding: = p * sum_{t_1+...+t_k = 0 mod p} c(t_1)...c(t_k).
       Then subtract the j=0 term |lambda_0|^{2k} = |S|^{2k}.
    """
    Sset = set(S)
    c = [0] * p
    for a in S:
        for b in S:
            c[(a - b) % p] += 1
    # number-theoretic transform style: count weighted (k)-tuples summing to 0 mod p
    # DP over convolution of c with itself k times (mod p index), weighted.
    # poly[t] = c(t); convolve k times under cyclic add.
    poly = c[:]  # k=1
    for _ in range(k - 1):
        new = [0] * p
        for t1 in range(p):
            if poly[t1] == 0:
                continue
            pt = poly[t1]
            for t2 in range(p):
                if c[t2] == 0:
                    continue
                new[(t1 + t2) % p] += pt * c[t2]
        poly = new
    total_with_zero = p * poly[0]  # sum_{j=0}^{p-1} |lambda_j|^{2k}
    lam0_2k = len(S) ** (2 * k)
    return total_with_zero - lam0_2k


def majorizes(u, v, tol=1e-7):
    """Does vector u majorize v? (both same length, equal sum).
       u majorizes v iff sorted-desc partial sums of u >= those of v, equal total.
       Returns 'u>v' if u majorizes v, 'v>u' if v majorizes u, 'eq', or 'incomp'."""
    us = sorted(u, reverse=True)
    vs = sorted(v, reverse=True)
    su = 0.0
    sv = 0.0
    u_ge = True
    v_ge = True
    for i in range(len(us)):
        su += us[i]
        sv += vs[i]
        if su < sv - tol:
            u_ge = False
        if sv < su - tol:
            v_ge = False
    if u_ge and v_ge:
        return 'eq'
    if u_ge:
        return 'u>v'
    if v_ge:
        return 'v>u'
    return 'incomp'


def ham_paths_circulant(S, p):
    n = p
    adj = [0] * n
    for i in range(n):
        mask = 0
        for j in range(n):
            if i != j and (j - i) % n in S:
                mask |= (1 << j)
        adj[i] = mask
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            cnt = dp[mask][v]
            if cnt == 0:
                continue
            avail = adj[v] & ~mask
            mm = avail
            while mm:
                u = (mm & -mm).bit_length() - 1
                dp[mask | (1 << u)][u] += cnt
                mm &= mm - 1
    full = (1 << n) - 1
    return sum(dp[full])


def main():
    for p in [7, 11, 19]:
        print("=" * 72)
        print(f"p = {p}")
        QR = qr_set(p)
        sets = tournament_connection_sets(p)

        # TEST 1: simultaneous moment floor
        print("  TEST 1: is M_{2k}(Paley) = min over all sets, for k=2..5?")
        for k in range(2, 6):
            mvals = [(moment_2k_exact(S, p, k), S) for S in sets]
            mmin = min(m for m, _ in mvals)
            argmin = [S for m, S in mvals if m == mmin]
            paley_min = (moment_2k_exact(QR, p, k) == mmin)
            print(f"    k={k}: M_min={mmin}, #argmin={len(argmin)}, "
                  f"Paley achieves min: {paley_min}, "
                  f"argmin==2sets: {len(argmin)==2}")

        # TEST 2: majorization -- Paley |lambda|^2 vector majorized by EVERY other set
        print("  TEST 2: is Paley's |lambda|^2 vector MAJORIZED BY every other set?")
        paley_vec = lambda_sq_vector_exact(QR, p)
        all_majorize_paley = True
        incomp_count = 0
        examples = []
        for S in sets:
            if S == QR:
                continue
            v = lambda_sq_vector_exact(S, p)
            rel = majorizes(v, paley_vec)
            if rel == 'v>u' or rel == 'eq':
                pass  # v majorizes paley -> good (paley is minimal)
            else:
                all_majorize_paley = False
                if rel == 'incomp':
                    incomp_count += 1
                examples.append((sorted(S), rel))
        print(f"    Every set majorizes Paley (Paley = majorization MIN)? {all_majorize_paley}")
        print(f"    # incomparable to Paley: {incomp_count}")
        if examples[:3]:
            print(f"    first non-majorizing examples: {examples[:3]}")

        # TEST 4: H monotone in majorization order? (the crux)
        if p <= 17:
            print("  TEST 4: is H Schur-CONCAVE in |lambda|^2 (A majorizes B => H(A) <= H(B))?")
            data = []
            for S in sets:
                v = lambda_sq_vector_exact(S, p)
                H = ham_paths_circulant(S, p)
                data.append((S, v, H))
            violations = 0
            checked = 0
            viol_examples = []
            for i in range(len(data)):
                for j in range(len(data)):
                    if i == j:
                        continue
                    Si, vi, Hi = data[i]
                    Sj, vj, Hj = data[j]
                    rel = majorizes(vi, vj)
                    if rel == 'u>v':  # vi majorizes vj => expect Hi <= Hj
                        checked += 1
                        if Hi > Hj:
                            violations += 1
                            if len(viol_examples) < 3:
                                viol_examples.append(
                                    (sorted(Si), Hi, sorted(Sj), Hj))
            print(f"    comparable ordered pairs checked: {checked}, "
                  f"Schur-concavity VIOLATIONS: {violations}")
            if viol_examples:
                print(f"    violation examples (A maj B but H(A)>H(B)): {viol_examples}")
            print(f"    => H Schur-concave in |lambda|^2 on the circulant family: "
                  f"{violations == 0}")

            # Also: is Paley strictly the unique global max H?
            maxH = max(H for _, _, H in data)
            argmax = [sorted(S) for S, _, H in data if H == maxH]
            print(f"    max H = {maxH}, argmax sets = {argmax}, "
                  f"Paley is unique-up-to-complement max: "
                  f"{set(map(tuple,argmax)) == {tuple(sorted(QR)), tuple(sorted(frozenset(t for t in range(1,p) if t not in QR)))}}")
        print()


if __name__ == "__main__":
    main()
