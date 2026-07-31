#!/usr/bin/env python3
"""eq(27) decode: it is a LOGIT gate, its weight is NOT pinned, and the weight's
primes all come from the bias data.  Independent of death-star's decode script.

Fragment under audit (external agent, relayed 2026-07-31):

    t_A = 389/2181,        (1+t_A)/(1-t_A) = 1285/896,
    t_B = 5872957/11821757,(1+t_B)/(1-t_B) = 8847357/2974400,
    (2457/6592)*log(8847357/2974400) - log(1285/896) > 1/25,

certified by the artanh sandwich
    2(t+t^3/3+t^5/5) <= log((1+t)/(1-t)) <= 2(t+t^3/3+t^5/(5(1-t^2))),
with claimed exact slack
    391926968594914200867482400554891567498742649630277
  / 82738859282193417287303438726081463937219800938169600.

WHAT THIS SCRIPT ESTABLISHES
  1. DECODE.  With p+q=1 and t=p-q=2p-1, (1+t)/(1-t) = p/q exactly.  So each log
     is a LOGIT (log-odds) and eq(27) is
         (2457/6592) * logit(p_B) - logit(p_A) > 1/25,
         p_A = 1285/2181,  p_B = 8847357/11821757.
     This answers "what functional of (p_A,p_B) are the two logs".
  2. EXACT RE-CERTIFICATION.  The claimed slack is reproduced as an exact
     Fraction, byte-identical.  Sandwich orientation is checked (the minus term
     uses the UPPER bound, the plus term the LOWER bound).
  3. THE WEIGHT IS NOT PINNED.  The admissible set of weights is the half-line
     alpha > alpha_min = (r_A + 1/25)/r_B, computed to 80 digits.  3/8, 2/5,
     1/2, 37/100, 41/110 and 7/19 all certify; 43/117, 18/49, 11/30, 4/11 do
     not.  So 2457/6592 is an OUTPUT of the construction, not a consequence of
     the inequality, and the simplest certified replacement is
         3*logit(p_B) > 8*logit(p_A) + 8/25.
  4. STRADDLE PERTURBATION TEST (the open question in the lane-D decode).
     The "capacity straddle" sigma_B < alpha < sigma_A holds on a WIDE window;
     intersected with the certificate window it is still wide, so the straddle
     is NOT decisive evidence for alpha = 2457/6592.  Reported honestly.
  5. MULTIPLICATIVE INDEPENDENCE.  257 divides only the A-side and 2949119 only
     the B-side, so r_A/r_B is irrational and alpha*r_B - r_A never vanishes for
     rational alpha.  The linear form always sits above a floor; only its SIZE
     is an irrationality-measure question.
  6. PRIME FINGERPRINT (new).  Every prime of the weight divides one of the four
     bias integers, and in particular 103 | 6592 AND 103 | 5872957 (the
     numerator of t_B).  Reported with an explicit coincidence heuristic.

Reproduce: python3 04-computation/amm12592_eq27_logit_decode_and_weight_window_klein.py
"""

from fractions import Fraction as Fr
from sympy import factorint
from mpmath import mp, mpf, log, floor

mp.dps = 90

tA = Fr(389, 2181)
tB = Fr(5872957, 11821757)
W = Fr(2457, 6592)
GAP = Fr(1, 25)
CLAIM = Fr(391926968594914200867482400554891567498742649630277,
           82738859282193417287303438726081463937219800938169600)


def art_lo(t, m=3):
    return 2 * sum(t ** (2 * i + 1) / (2 * i + 1) for i in range(m))


def art_hi(t, m=3):
    return 2 * (sum(t ** (2 * i + 1) / (2 * i + 1) for i in range(m - 1))
                + t ** (2 * m - 1) / ((2 * m - 1) * (1 - t ** 2)))


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


def main():
    rule("1. DECODE: the two logs are LOGITS of two biases")
    pA, pB = (1 + tA) / 2, (1 + tB) / 2
    qA, qB = 1 - pA, 1 - pB
    print(f"  t = p - q = 2p - 1  =>  (1+t)/(1-t) = p/q  (the odds).")
    print(f"  p_A = {pA}   q_A = {qA}   odds p_A/q_A = {pA / qA}"
          f"   matches 1285/896: {pA / qA == Fr(1285, 896)}")
    print(f"  p_B = {pB}   q_B = {qB}   odds p_B/q_B = {pB / qB}"
          f"   matches 8847357/2974400: {pB / qB == Fr(8847357, 2974400)}")
    print("  => eq(27) IS   (2457/6592)*logit(p_B) - logit(p_A) > 1/25,"
          "   logit(p) = log(p/q).")
    ok1 = pA / qA == Fr(1285, 896) and pB / qB == Fr(8847357, 2974400)

    print()
    rule("2. EXACT RE-CERTIFICATION of the claimed slack")
    slack = W * art_lo(tB) - art_hi(tA) - GAP
    print(f"  computed = {slack}")
    print(f"  claimed  = {CLAIM}")
    ok2 = slack == CLAIM
    print(f"  BYTE-EXACT MATCH: {ok2}   (> 0: {slack > 0}, = {float(slack):.17f})")
    print("  sandwich orientation: minus term uses UPPER(t_A), plus term uses"
          " LOWER(t_B) -- sound direction.")

    print()
    rule("3. THE WEIGHT IS NOT PINNED: the admissible set is a half-line")
    rA = log(mpf(1285) / 896)
    rB = log(mpf(8847357) / 2974400)
    amin = (rA + mpf(1) / 25) / rB
    print(f"  r_A = logit(p_A) = {mp.nstr(rA, 30)}")
    print(f"  r_B = logit(p_B) = {mp.nstr(rB, 30)}")
    print(f"  alpha_min = (r_A + 1/25)/r_B = {mp.nstr(amin, 50)}")
    print("  admissible weights = { alpha : alpha > alpha_min }.  Exact tests:")
    cands = [Fr(2457, 6592), Fr(41, 110), Fr(3, 8), Fr(2, 5), Fr(1, 2),
             Fr(37, 100), Fr(7, 19), Fr(43, 117), Fr(18, 49), Fr(11, 30), Fr(4, 11)]
    ok3 = True
    for a in cands:
        v = a * art_lo(tB) - art_hi(tA) - GAP
        pred = mpf(a.numerator) / a.denominator > amin
        print(f"    alpha = {str(a):11s} = {float(a):.8f}  certified: {str(v > 0):5s}"
              f"  (half-line predicts {str(pred):5s})  margin {float(v):+.3e}")
        # the sandwich is slightly lossy, so certified => predicted, not iff
        if v > 0 and not pred:
            ok3 = False
    print("  SIMPLEST CERTIFIED REPLACEMENT:  3*logit(p_B) > 8*logit(p_A) + 8/25")
    print(f"    exact slack at alpha = 3/8: {Fr(3, 8) * art_lo(tB) - art_hi(tA) - GAP}")
    print("  => 2457/6592 is an OUTPUT of the construction, not forced by (27).")

    print()
    rule("4. STRADDLE PERTURBATION TEST (was an open question)")
    sA = Fr(896, 1285)
    sB = Fr(2974400, 8847357)
    print(f"  dual ray slopes sigma = q/p:  sigma_A = {sA} = {float(sA):.8f},"
          f"  sigma_B = {sB} = {float(sB):.8f}")
    print(f"  'capacity straddle' asks sigma_B < alpha < sigma_A, i.e."
          f" alpha in ({float(sB):.6f}, {float(sA):.6f})")
    lo = max(float(sB), float(amin))
    print(f"  intersect with the certificate window alpha > {float(amin):.6f}:"
          f"  alpha in ({lo:.6f}, {float(sA):.6f})")
    print(f"  window width = {float(sA) - lo:.6f}, which contains 3/8, 2/5, 37/100,"
          f" 7/19, 41/110 and 2457/6592 alike.")
    print("  VERDICT: the straddle is NOT decisive.  Do not cite it as evidence"
          " for the specific weight.")

    print()
    rule("5. MULTIPLICATIVE INDEPENDENCE => the linear form never vanishes")
    fa = {**factorint(1285), **factorint(896)}
    fb = {**factorint(8847357), **factorint(2974400)}
    onlyA = sorted(set(fa) - set(fb))
    onlyB = sorted(set(fb) - set(fa))
    print(f"  1285*896   primes: {sorted(set(fa))}")
    print(f"  8847357*2974400 primes: {sorted(set(fb))}")
    print(f"  primes on the A side only: {onlyA};  on the B side only: {onlyB}")
    ok5 = bool(onlyA) and bool(onlyB)
    print("  If n*r_A = m*r_B then (1285/896)^n = (8847357/2974400)^m; comparing the")
    print(f"  exponent of {onlyA[-1]} forces n = 0, hence m = 0.  So r_A/r_B is IRRATIONAL")
    print("  and alpha*r_B - r_A != 0 for every rational alpha: there is always a floor.")
    print(f"  PROVED: {ok5}")

    print()
    rule("6. PRIME FINGERPRINT of the weight (new observation)")
    print(f"  2457 = {factorint(2457)}      6592 = {factorint(6592)}")
    print(f"  1285 = {factorint(1285)}   896 = {factorint(896)}")
    print(f"  8847357 = {factorint(8847357)}   2974400 = {factorint(2974400)}")
    print(f"  5872957 (num t_B) = {factorint(5872957)}   11821757 = {factorint(11821757)}")
    bias_primes = set()
    for n in (1285, 896, 8847357, 2974400, 389, 2181, 5872957, 11821757):
        bias_primes |= set(factorint(n))
    wprimes = set(factorint(2457)) | set(factorint(6592))
    print(f"  primes of the weight: {sorted(wprimes)}")
    print(f"  all of them divide some bias integer: {wprimes <= bias_primes}")
    print("  STRIKING: 103 | 6592 AND 103 | 5872957 (the numerator of t_B).")
    print("  Heuristic: a random 7-digit integer is divisible by 103 with"
          " probability ~1/103 ~ 1%.")
    print("  So this is suggestive of a shared construction, NOT proof.  The other")
    print("  weight primes (2,3,7,13) are small enough to divide something by chance.")
    print("  ACTIONABLE: search constructions whose output ratio is 2457:6592 =")
    print("  (3^3*7*13):(2^6*103), with 103 inherited from t_B's numerator.")

    print()
    rule("SUMMARY")
    print(f"  decode={ok1}  exact-recertification={ok2}  half-line-consistent={ok3}"
          f"  mult-independence={ok5}")
    print("  C = 1 + gamma reading, if gamma = 2457/6592:  C = 9049/6592 ="
          f" {float(Fr(9049, 6592)):.15f}")


if __name__ == "__main__":
    main()
