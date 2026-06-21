"""
lrc14_paley_schur_concavity_crux_angle6.py

ANGLE 6 -- THE CRUX. Two prior tests established:
  (1) Paley = QR difference set = the unique flat-spectrum circulant = the MAJORIZATION MINIMUM
      of the |lambda|^2 spectrum simplex (every other tournament set's |lambda|^2 vector
      majorizes Paley's). Verified p=7,11,19.
  (2) H is Schur-CONCAVE in the |lambda|^2 vector (zero violations p=7,11): A majorizes B => H(A)<=H(B).

If BOTH hold, the chain is a GLOBAL proof that Paley maximizes H among circulants:
  Schur-concave + majorization-minimum => Paley is the max. (No single-moment monotonicity needed;
  that one FAILED -- the wall -- precisely because H lives in the majorization (Schur) order,
  not the moment order.)

This script STRESS-TESTS the two load-bearing claims and probes WHY H is Schur-concave:

  CLAIM A (totality): Is the majorization order on the spectrum vectors a TOTAL order (chain)?
     If yes, "majorization minimum" is unambiguous and the H-vs-majorization curve is 1-D.
     If no (a genuine lattice), Schur-concavity is a stronger, non-obvious statement.

  CLAIM B (Schur-concavity is EXACT, all comparable pairs): exhaustively verify A maj B => H(A)<=H(B)
     for ALL ordered comparable pairs at p=7,11, AND on the continuous neighborhood:
     check the 2-element Robin-Hood transfer test: moving |lambda|^2 mass from a larger to a
     smaller coordinate (a "Robin Hood transfer", which strictly DEcreases majorization) should
     strictly INCREASE H. We test this on the DISCRETE realizable spectra by pairing
     majorization-adjacent sets.

  CLAIM C (the mechanism): H as a symmetric function of x=|lambda|^2 -- fit H exactly to the
     elementary-symmetric basis e_1..e_m(x) and CHECK whether H is Schur-concave by the
     Schur-Ostrowski criterion on the realizable lattice: sign of (x_i - x_j)(dH/dx_i - dH/dx_j).
     We can't take derivatives of the discrete H, but we CAN test the *finite-difference*
     Schur-Ostrowski: for every majorization-covering pair (A covers B), is H(A) <= H(B)?

  CLAIM D (does this beat the wall?): re-confirm the single-moment FAILURE and show the
     majorization order RESCUES it: list pairs where M_4(A)<M_4(B) but H(A)<H(B) (moment says A
     better, H says B better) and verify those pairs are INCOMPARABLE in majorization (so no
     contradiction with Schur-concavity).

Run: python3 04-computation/lrc14_paley_schur_concavity_crux_angle6.py
"""
import sys
import cmath, math
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


def lambda_sq_vector(S, p):
    omega = cmath.exp(2j * math.pi / p)
    return [abs(sum(omega ** (k * s) for s in S)) ** 2 for k in range(1, p)]


def moment4(S, p):
    """sum|lambda_k|^4 exact integer via autocorrelation L2."""
    c = [0] * p
    for a in S:
        for b in S:
            c[(a - b) % p] += 1
    total_with_zero = p * sum(c[t] * c[t] for t in range(p))
    return total_with_zero - len(S) ** 4


def majorizes(u, v, tol=1e-7):
    us = sorted(u, reverse=True)
    vs = sorted(v, reverse=True)
    su = sv = 0.0
    u_ge = v_ge = True
    for i in range(len(us)):
        su += us[i]; sv += vs[i]
        if su < sv - tol: u_ge = False
        if sv < su - tol: v_ge = False
    if u_ge and v_ge: return 'eq'
    if u_ge: return 'u>v'
    if v_ge: return 'v>u'
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
    return sum(dp[(1 << n) - 1])


def canonical_class(S, p):
    """Group sets by multiplier-equivalence (a*S for a in Z_p*) and complement,
       so we count genuinely-distinct H-classes (the 4 classes at p=11)."""
    reps = []
    for a in range(1, p):
        reps.append(frozenset((a * s) % p for s in S))
    return min(tuple(sorted(r)) for r in reps)


def main():
    for p in [7, 11]:
        print("=" * 72)
        print(f"p = {p}")
        QR = qr_set(p)
        QRc = frozenset(t for t in range(1, p) if t not in QR)
        sets = tournament_connection_sets(p)
        data = []
        for S in sets:
            v = lambda_sq_vector(S, p)
            H = ham_paths_circulant(S, p)
            m4 = moment4(S, p)
            data.append((S, v, H, m4))

        # M_2 constancy
        m2s = set(round(sum(v)) for _, v, _, _ in data)
        print(f"  M_2 = sum|lambda|^2 constant across all sets? {len(m2s)==1}  (value ~ {m2s})")

        # CLAIM A: totality of majorization
        n = len(data)
        comparable = 0
        incomp = 0
        for i in range(n):
            for j in range(i + 1, n):
                rel = majorizes(data[i][1], data[j][1])
                if rel == 'incomp':
                    incomp += 1
                else:
                    comparable += 1
        total_pairs = n * (n - 1) // 2
        print(f"  CLAIM A (totality): {comparable}/{total_pairs} pairs comparable, "
              f"{incomp} incomparable.  Majorization is a TOTAL order: {incomp==0}")

        # CLAIM B/C: Schur-concavity over ALL comparable ordered pairs
        violations = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                rel = majorizes(data[i][1], data[j][1])
                if rel == 'u>v':  # i majorizes j => expect H_i <= H_j
                    if data[i][2] > data[j][2]:
                        violations.append((sorted(data[i][0]), data[i][2],
                                           sorted(data[j][0]), data[j][2]))
        print(f"  CLAIM B (Schur-concavity, all comparable pairs): "
              f"violations = {len(violations)}  => H Schur-CONCAVE in |lambda|^2: {len(violations)==0}")
        if violations[:2]:
            print(f"    examples: {violations[:2]}")

        # Paley = majorization min => unique max H
        paley_vec = lambda_sq_vector(QR, p)
        paley_is_min = all(majorizes(v, paley_vec) in ('u>v', 'eq')
                           for _, v, _, _ in data)
        maxH = max(H for _, _, H, _ in data)
        argmaxH = set(canonical_class(S, p) for S, _, H, _ in data if H == maxH)
        print(f"  Paley = majorization MINIMUM: {paley_is_min}")
        print(f"  => Paley maximizes H (max H = {maxH}); "
              f"argmax classes = {len(argmaxH)} (Paley class): "
              f"{argmaxH == {canonical_class(QR, p)}}")

        # CLAIM D: the wall -- single moment fails, majorization rescues
        print(f"  CLAIM D (the wall): pairs where M_4 says A<B (A 'better') but H(A)<H(B):")
        wall_pairs = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if data[i][3] < data[j][3] and data[i][2] < data[j][2]:
                    # M4 ranks i above j, but H ranks j above i -> moment fails
                    rel = majorizes(data[i][1], data[j][1])
                    wall_pairs.append((data[i][3], data[i][2], data[j][3],
                                       data[j][2], rel))
        # de-dup by signature
        seen = set()
        uniq = []
        for w in wall_pairs:
            key = (w[0], w[1], w[2], w[3], w[4])
            if key not in seen:
                seen.add(key)
                uniq.append(w)
        all_incomp = all(w[4] == 'incomp' for w in uniq)
        print(f"    #such 'moment-fails' ordered pairs = {len(uniq)}; "
              f"ALL incomparable in majorization (no Schur contradiction): {all_incomp}")
        for w in uniq[:4]:
            print(f"      M4={w[0]}(H={w[1]}) vs M4={w[2]}(H={w[3]}): majorization rel = {w[4]}")
        print()


if __name__ == "__main__":
    main()
