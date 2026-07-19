"""
opus-2026-07-19-S401: REFUTATION of Sungkawichai-Trakulthongchai Conjecture 7.1
(arXiv:2604.23906, verbatim: "Let k+1 be a positive integer. There exists a
constant D such that, for any integer d >= D, every non-tight speed tuple
v in Z_{>0}^k with gcd(v)=1 has a witness time in (1/d)Z.")

THE COUNTEREXAMPLE (three lines, self-contained -- no dependence on their
paper beyond the conjecture's statement):

  For any d with (k+1) NOT dividing d, the cluster V_d = {d+1, d+2, ..., d+k}:
  (i)  is coprime (consecutive integers) and NON-TIGHT: at t* = 1/(2d+k+1)
       the positions (d+i)/(2d+k+1) lie in a band symmetric about 1/2 with
       min distance (d+1)/(2d+k+1) = 1/(k+1) + d(k-1)/((k+1)(2d+k+1)),
       STRICTLY above 1/(k+1) (margin -> (k-1)/(2(k+1)));
  (ii) has NO witness time in (1/d)Z: at any j/d, since d | d we have
       ||(d+i) j/d|| = ||i j/d||  -- the cluster's distances at the d-grid
       EQUAL the AP {1..k}'s -- and no j/d with (k+1) not| d is an AP witness:
  (iii) LEMMA (pigeonhole + quantization): if min_{1<=i<=k} ||i t|| >= 1/(k+1)
       at t = j/d, then among the k+1 points {0, t, 2t, ..., kt} two lie
       within 1/(k+1) (k+1 circular gaps sum to 1), so ||m t|| <= 1/(k+1)
       for some 1 <= m <= k; with the hypothesis, ||m j/d|| = 1/(k+1)
       exactly; but ||m j/d|| is a multiple of 1/d' (d' = d/gcd(j,d)),
       so (k+1) | d' | d -- contradiction.
  Hence for every D there is a BAD d >= D (any d with (k+1) not| d) and a
  non-tight coprime family with no witness in (1/d)Z.  Conjecture 7.1 is
  FALSE for every k >= 2.  QED.

  Mechanism, in repo language: the d-grid restriction is TRANSLATION-BLIND
  when d | translation -- the cluster is the AP translated by d, M jumps from
  1/(k+1) to ~1/2, but the (1/d)Z-restricted problem cannot see the
  translation (the repo's translation-invariance lore, fifth appearance,
  now striking a literature conjecture); and the AP's witness set is the
  measure-zero {s/(k+1)} (their own Remark 3.2), which the d-grid misses
  entirely unless (k+1) | d.

This script verifies everything EXACTLY for k = 12 and k = 13 over a ladder
of d values including large ones, plus the band-witness margins.
"""
from fractions import Fraction
from math import gcd

def ap_equal_check(k, d, samples):
    """||(d+i) j/d|| == ||i j/d|| for all i, sampled j."""
    for j in samples:
        for i in range(1, k+1):
            a = ((d+i)*j) % d
            b = (i*j) % d
            if min(a, d-a) != min(b, d-b):
                return False
    return True

def no_witness_in_grid(k, d):
    """check exhaustively: no j in [0,d) has min_i ||i j/d|| >= 1/(k+1);
    equivalently (k+1)*min_residue < d for every j."""
    for j in range(0, d):
        ok = True
        for i in range(1, k+1):
            r = (i*j) % d
            r = min(r, d-r)
            if (k+1)*r < d:
                ok = False
                break
        if ok:
            return j  # found a witness -- would falsify the refutation
    return None

def band_margin(k, d):
    """strict non-tightness witness value at t* = 1/(2d+k+1)."""
    q = 2*d + k + 1
    vals = [Fraction(min(d+i, q-(d+i)), q) for i in range(1, k+1)]
    return min(vals) - Fraction(1, k+1)

if __name__ == "__main__":
    for k in (12, 13):
        print(f"=== k = {k} (threshold 1/{k+1}) ===")
        # d values: small, medium, large; all with (k+1) not| d
        ds = [x for x in (25, 100, 997, 10007, 100003) if x % (k+1) != 0]
        for d in ds:
            V = [d+i for i in range(1, k+1)]
            g = 0
            for v in V: g = gcd(g, v)
            m = band_margin(k, d)
            same = ap_equal_check(k, d, [1, 2, d//3, d//2, d-1])
            w = no_witness_in_grid(k, d)
            print(f"  d = {d}: gcd = {g}; band margin over 1/{k+1} = {m} "
                  f"(> 0: {m > 0}); cluster==AP at d-grid: {same}; "
                  f"witness j/{d} found: {w}  (None = BAD d, as claimed)")
        # control: d a multiple of (k+1) SHOULD have a witness (j/d = s/(k+1))
        dctrl = (k+1) * 100
        w = no_witness_in_grid(k, dctrl)
        print(f"  control d = {dctrl} = (k+1)*100: witness found at j = {w} "
              f"(expected j with j/{dctrl} = s/{k+1})")
    print("\nCONCLUSION: for every D pick any d >= D with (k+1) not| d;")
    print("V_d = {d+1,...,d+k} is coprime, non-tight (strict band witness),")
    print("and has NO witness in (1/d)Z.  Conjecture 7.1 (uniform D over v)")
    print("is FALSE for every k >= 2.")
