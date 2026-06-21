"""
lrc14_paley_difference_set_flatspectrum_angle6.py

ANGLE 6 (tournament Paley H-max via graph symmetry / spectral).

OUT-OF-BOX HYPOTHESIS (HYP candidate, falsifiable):

  For a CIRCULANT tournament on Z_p (p prime, p = 3 mod 4), with connection set S
  (|S| = (p-1)/2, S disjoint from -S, S u (-S) = Z_p \ {0}):

  (A) FLAT SPECTRUM  <=>  S is a (p, (p-1)/2, (p-3)/4) DIFFERENCE SET
      i.e. the autocorrelation  c(t) = #{(a,b) in SxS : a - b = t}  is CONSTANT = (p-3)/4
      for all t != 0.  This is the Paley / QR difference set, unique up to multiplier.

  (B) Sum |lambda_k|^4  has a GLOBAL Boolean floor over all tournament connection sets,
      and the floor is achieved EXACTLY at the difference set (flat spectrum).
      Equivalently:  sum_t c(t)^2  is minimized exactly at the constant-autocorrelation set.

  (C) min sum|lambda|^4  <=>  max H  (so the difference-set / Paley set is the H-max),
      and this is a STRICT global statement, not merely the local Hessian (THM-134).

KEY IDENTITY exploited (Parseval / Wiener-Khinchin on Z_p):
  The eigenvalues are lambda_k = sum_{s in S} omega^{ks}.  Then
      sum_k |lambda_k|^4  =  p * sum_t |c(t)|^2     (c = autocorrelation of indicator 1_S)
  where c(t) = sum_a 1_S(a) 1_S(a - t) = #{ pairs in S with difference t }.
  So minimizing sum|lambda|^4 == minimizing the L2 norm of the autocorrelation
  == making the autocorrelation as FLAT as possible == difference set.

This is a HARD Boolean optimization with a clean combinatorial floor. If the floor is a
difference set and difference set <=> H-max, we have a *global* (not local) symmetry/spectral
characterization of Paley H-max -- the tournament side of the wall, where Doyle-Holt failed
(HYP-2748: Doyle-Holt is the OPPOSITE, half-arc / non-self-converse side).

Run: python3 04-computation/lrc14_paley_difference_set_flatspectrum_angle6.py
"""
import sys
from fractions import Fraction
from itertools import combinations

sys.stdout.reconfigure(line_buffering=True)


def qr_set(p):
    return frozenset((a * a) % p for a in range(1, p))


def tournament_connection_sets(p):
    """All valid tournament connection sets S on Z_p:
       for each nonzero pair {t, -t} pick exactly one. (p-1)/2 binary choices.
       Returns list of frozensets S with |S| = (p-1)/2, S cap (-S) = empty."""
    m = (p - 1) // 2
    pairs = [(t, (p - t) % p) for t in range(1, m + 1)]  # {t, p-t}
    sets = []
    for bits in range(1 << m):
        S = set()
        for i, (a, b) in enumerate(pairs):
            if bits >> i & 1:
                S.add(a)
            else:
                S.add(b)
        sets.append(frozenset(S))
    return sets


def autocorr_l2_sq(S, p):
    """sum_{t != 0} c(t)^2 where c(t) = #{(a,b) in SxS : (a-b) mod p = t}.
       Exact integer."""
    Sset = set(S)
    c = [0] * p
    for a in S:
        for b in S:
            c[(a - b) % p] += 1
    # c[0] = |S| always. Off-zero autocorrelation L2:
    return sum(c[t] * c[t] for t in range(1, p)), c


def autocorr_is_flat(c, p):
    """Off-zero autocorrelation constant?"""
    vals = set(c[t] for t in range(1, p))
    return len(vals) == 1, vals


def sum_lambda4_via_autocorr(S, p):
    """sum_k |lambda_k|^4 = p * sum_t |c(t)|^2 over ALL t (including t=0).
       Verify the Wiener-Khinchin identity numerically-exactly via integer c.
       Returns exact Fraction times p? Actually sum_k|lambda_k|^4 = p * sum_t c(t)^2 (all t).
    """
    Sset = set(S)
    c = [0] * p
    for a in S:
        for b in S:
            c[(a - b) % p] += 1
    full = sum(c[t] * c[t] for t in range(p))  # includes t=0 term |S|^2
    return p * full, c


def ham_paths_circulant(S, p):
    """Count directed Hamiltonian paths in the circulant tournament (Held-Karp DP).
       H = number of HP. (We use HP count as the H proxy = A038375 quantity.)"""
    n = p
    adj = [0] * n  # adj[i] bitmask of out-neighbors
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
    for p in [7, 11, 19, 23]:
        print("=" * 70)
        print(f"p = {p}  (p mod 4 = {p % 4})")
        QR = qr_set(p)
        sets = tournament_connection_sets(p)
        print(f"  # tournament connection sets = {len(sets)} = 2^{(p-1)//2}")

        # find QR among them (it's a valid tournament set iff p = 3 mod 4)
        results = []
        flat_sets = []
        for S in sets:
            l2off, c = autocorr_l2_sq(S, p)
            flat, vals = autocorr_is_flat(c, p)
            results.append((l2off, S, flat, vals, c))
            if flat:
                flat_sets.append((S, vals))

        results.sort(key=lambda r: r[0])
        min_l2 = results[0][0]
        # which sets achieve the min L2 autocorrelation?
        argmin = [r for r in results if r[0] == min_l2]

        print(f"  min off-zero autocorr L2 (sum c(t)^2) = {min_l2}")
        print(f"  # sets achieving min = {len(argmin)}")
        print(f"  # flat-spectrum (constant off-zero autocorr) sets = {len(flat_sets)}")

        # Difference set parameter check
        lam_target = (p - 3) // 4
        is_para_diffset = (p % 4 == 3)
        print(f"  expected difference-set lambda = (p-3)/4 = {Fraction(p-3,4)}  "
              f"(integer iff p=3 mod4: {is_para_diffset})")

        # Are the flat sets exactly QR and its complement (= -QR)?
        QRc = frozenset(t for t in range(1, p) if t not in QR)
        flat_set_frozen = set(S for (S, v) in flat_sets)
        print(f"  QR_p = {sorted(QR)}")
        if p % 4 == 3:
            print(f"  QR is a valid tournament set: "
                  f"{QR in set(sets)}   (QR cap -QR empty: {len(QR & frozenset((-x)%p for x in QR))==0})")
        print(f"  flat sets == {{QR, -QR}}? "
              f"{flat_set_frozen == {QR, QRc}}   (flat sets: {[sorted(S) for S in flat_set_frozen]})")

        # (B) does flat == argmin(sum|lambda|^4)?
        print(f"  argmin(L2 autocorr) sets: {[sorted(r[1]) for r in argmin]}")
        argmin_frozen = set(r[1] for r in argmin)
        print(f"  argmin == flat sets? {argmin_frozen == flat_set_frozen}")

        # (C) does flat == H-max?  (only feasible to brute-H for p <= 17 or so)
        if p <= 17:
            print("  computing H (HP count) for all sets ...")
            Hs = []
            for (l2off, S, flat, vals, c) in results:
                H = ham_paths_circulant(S, p)
                Hs.append((H, l2off, S, flat))
            Hs.sort(key=lambda x: -x[0])
            maxH = Hs[0][0]
            argmaxH = [x for x in Hs if x[0] == maxH]
            print(f"  max H = {maxH}, # argmax = {len(argmaxH)}, "
                  f"argmax sets: {[sorted(x[2]) for x in argmaxH]}")
            # is argmax H == flat sets?
            argmaxH_frozen = set(x[2] for x in argmaxH)
            print(f"  argmax(H) == flat sets? {argmaxH_frozen == flat_set_frozen}")
            # monotone: does lower L2-autocorr always mean higher H?
            # sort by L2 ascending and check H descending
            by_l2 = sorted([(x[1], x[0]) for x in Hs])  # (l2, H)
            mono = all(by_l2[i][1] >= by_l2[i+1][1] for i in range(len(by_l2)-1))
            print(f"  STRICT monotone (L2 autocorr up => H down, over ALL sets)? {mono}")
            if not mono:
                # show first violation
                for i in range(len(by_l2)-1):
                    if by_l2[i][1] < by_l2[i+1][1]:
                        print(f"    VIOLATION: L2={by_l2[i][0]} H={by_l2[i][1]} then "
                              f"L2={by_l2[i+1][0]} H={by_l2[i+1][1]}")
                        break
        else:
            print("  (skipping brute H for p>17)")
        print()


if __name__ == "__main__":
    main()
