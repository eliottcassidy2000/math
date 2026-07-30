#!/usr/bin/env python3
"""Exact referee: decode of the two-bias artanh certificate for the AMM 12592
minimal-C frontier (biased-coin critical-run extraction, THM-2160/THM-2225).

Fragment supplied to session death-star 2026-07-30 (source: proposer-side
notes on Monthly Problem 12592, equation numbering external to this repo):

    upper bound  U(t) = 2(t + t^3/3 + t^5/(5(1-t^2))) >= log((1+t)/(1-t))
    lower bound  L(t) = 2(t + t^3/3 + t^5/5)          <= log((1+t)/(1-t))
    t_A = 389/2181      (1+t_A)/(1-t_A) = 1285/896
    t_B = 5872957/11821757   (1+t_B)/(1-t_B) = 8847357/2974400
    claim: "right side of (27) minus 1/25 is at least F > 0"

This referee proves, in exact rational arithmetic:

  (1) the two artanh polynomial bounds hold at every 0<t<1 (tail-comparison),
  (2) (2457/6592)*L(t_B) - U(t_A) = 1/25 + F exactly, where F is the
      51/53-digit fraction from the fragment; hence the external equation (27)
      is the strict log inequality

        (2457/6592)*log(8847357/2974400) - log(1285/896) > 1/25,      (27)

      equivalently  (1285/896)^6592 * e^(6592/25) < (8847357/2974400)^2457.
  (3) the bias translation: t = p-q gives p_A = 1285/2181, p_B =
      8847357/11821757, so both logs are log-likelihood ratios log(p/q)
      (= 2*artanh(t), the rapidity of the bias).

Discovery route: PSLQ integer-relation on (L(t_B), U(t_A), 1, 1/25+F) found
(-2457, 6592, 0, 6592); verified exactly below without PSLQ.
"""

from fractions import Fraction as Fr


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def upper(t: Fr) -> Fr:
    return 2 * (t + t**3 / 3 + t**5 / (5 * (1 - t**2)))


def lower(t: Fr) -> Fr:
    return 2 * (t + t**3 / 3 + t**5 / 5)


T_A = Fr(389, 2181)
T_B = Fr(5872957, 11821757)
F = Fr(
    391926968594914200867482400554891567498742649630277,
    82738859282193417287303438726081463937219800938169600,
)


def main() -> None:
    # ratio simplifications from the fragment
    require((1 + T_A) / (1 - T_A) == Fr(1285, 896), "t_A ratio drift")
    require((1 + T_B) / (1 - T_B) == Fr(8847357, 2974400), "t_B ratio drift")

    # bias translation t = p - q, p + q = 1
    p_a, p_b = (1 + T_A) / 2, (1 + T_B) / 2
    require(p_a == Fr(1285, 2181), "p_A drift")
    require(p_b == Fr(8847357, 11821757), "p_B drift")
    require(p_a / (1 - p_a) == Fr(1285, 896), "p_A likelihood ratio drift")
    require(p_b / (1 - p_b) == Fr(8847357, 2974400), "p_B likelihood ratio drift")

    # exact identity: (2457/6592) L(t_B) - U(t_A) == 1/25 + F
    combo = Fr(2457, 6592) * lower(T_B) - upper(T_A)
    require(combo == Fr(1, 25) + F, "certificate identity drift")
    require(F > 0, "certificate margin not positive")

    # soundness of the two polynomial bounds for every 0<t<1:
    # log((1+t)/(1-t)) = 2 sum_{m>=0} t^(2m+1)/(2m+1).
    # lower: truncation of a positive series.
    # upper: tail sum 2*sum_{m>=3} t^(2m+1)/(2m+1) <= (2/7)*t^7/(1-t^2)
    #        and the replacement 5(1-t^2) term adds exactly
    #        2*t^5/5*(1/(1-t^2)-1) = (2/5) t^7/(1-t^2) >= (2/7) t^7/(1-t^2).
    # Both hold with strict inequality for 0<t<1; spot-check exactness margins
    # at the two certificate points with a 60-term partial sum plus tail bound.
    for t, name in ((T_A, "t_A"), (T_B, "t_B")):
        partial = 2 * sum(t ** (2 * m + 1) / (2 * m + 1) for m in range(60))
        tail_hi = 2 * t**121 / (121 * (1 - t**2))
        require(partial > lower(t), f"lower bound not strict at {name}")
        require(partial + tail_hi < upper(t), f"upper bound not strict at {name}")

    print("t_A =", T_A, " p_A =", p_a, " p_A/q_A = 1285/896")
    print("t_B =", T_B, " p_B =", p_b, " p_B/q_B = 8847357/2974400")
    print("identity: (2457/6592)*L(t_B) - U(t_A) - 1/25 == F  [exact]")
    print("hence (27): (2457/6592)*log(8847357/2974400) - log(1285/896) > 1/25")
    print("margin F =", float(F), "(exact fraction preserved in source)")
    print("status=CERTIFICATE_DECODED_VERIFIED_EXACT")


if __name__ == "__main__":
    main()
