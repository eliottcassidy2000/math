#!/usr/bin/env python3
r"""Compare the 2026 finite-sieve proof scale with repo certificate clocks.

Session S580.  Convention: paper k = moving speeds, repo total n = k + 1.
The prime lists and proof diagrams are from Sungkawichai-Trakulthongchai,
arXiv:2604.23906, Table 1 and Section 5.  The certificate side is a
methodology-size model, not a completed proof.  THM-398/HYP-2103 adds a
dodge-first interval filter: if G(S\{v}) has a component longer than 2/(n v),
then S is already loose.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb, exp, gcd, log, log10
from typing import Iterable


PAPER_PRIMES = {
    10: [
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
        269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
        349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
        431, 433, 439, 443, 449, 457, 461, 463, 467,
    ],
    11: [
        23, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263,
        269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
        349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
        431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
        503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
    ],
    12: [
        167, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
        241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
        317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
        401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
        479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
        571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
        647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
    ],
}


@dataclass(frozen=True)
class Case:
    total_n: int
    paper_k: int
    diagram: str
    final_l_if_exact: int
    paper_runtime: str
    note: str


CASES = [
    Case(11, 10, "x2,x2,x2,project; polynomial removes AP", 11, "~45 min", "odd prime n=11; c=11 lift avoided"),
    Case(12, 11, "x2,x2,x2,x2,x3,x3", 144, "~40 h", "composite n=12; small-factor lifts"),
    Case(13, 12, "x2,x2,x2,project; polynomial removes AP", 13, "~40 d", "odd prime n=13; c=13 lift avoided"),
]

VITALI_SAMPLE = [
    (6, 72.4, 331),
    (8, 78.7, 255),
    (10, 88.9, 133),
    (12, 91.5, 102),
    (14, 96.8, 39),
]


def logsumexp(values: Iterable[float]) -> float:
    values = list(values)
    m = max(values)
    return m + log(sum(exp(v - m) for v in values))


def euler_phi(n: int) -> int:
    out = 0
    for a in range(1, n + 1):
        if gcd(a, n) == 1:
            out += 1
    return out


def ln_B(k: int) -> float:
    return k * (((k - 1) * log(comb(k + 1, 2))) - log(k))


def fmt_log10(x: float) -> str:
    return f"{x:8.2f}"


def main() -> None:
    print("LRC paper finite-sieve scale vs certificate-clock scale (S580)")
    print("Convention: total n = paper k + 1")
    print()
    print("PAPER PROOF SCHEDULE")
    print("n  k  |P| p-range   ln(prod P)  ln(B_k)  runtime   diagram")
    for case in CASES:
        primes = PAPER_PRIMES[case.paper_k]
        print(
            f"{case.total_n:2d} {case.paper_k:2d} "
            f"{len(primes):4d} {min(primes):3d}-{max(primes):3d} "
            f"{sum(log(p) for p in primes):11.3f} {ln_B(case.paper_k):8.3f} "
            f"{case.paper_runtime:>8s}  {case.diagram}"
        )

    print()
    print("LOG10 WORK-SCALE COMPARISON")
    print(
        "n  raw sum p^k   exact final ansatz   paper heuristic   "
        "cert scan ops   raw->cert   final->cert"
    )
    for case in CASES:
        k = case.paper_k
        primes = PAPER_PRIMES[k]
        raw = logsumexp(k * log(p) for p in primes) / log(10)
        final = logsumexp(k * log(p * case.final_l_if_exact) for p in primes) / log(10)
        heuristic = logsumexp(((k + 1) / 2) * log(p) - log(k) - k * log(2) for p in primes) / log(10)

        C = 2 * case.total_n - 1
        unit_shells = euler_phi(C) // 2
        all_shells = case.total_n - 1
        nonunit_shells = all_shells - unit_shells
        d_obligations = case.total_n - 2
        n_obligations = case.total_n - 1
        obligations = d_obligations + unit_shells + n_obligations
        pair_clocks = comb(k, 2)
        transversals = 2**k
        cert_ops = transversals * (pair_clocks + obligations + 1)
        cert = log10(cert_ops)

        print(
            f"{case.total_n:2d} {fmt_log10(raw)} {fmt_log10(final)} "
            f"{fmt_log10(heuristic)} {fmt_log10(cert)} "
            f"{raw - cert:10.2f} {final - cert:12.2f}"
        )
        print(
            f"   C=2n-1={C}; unit/nonunit shells={unit_shells}/{nonunit_shells}; "
            f"D/U/N={d_obligations}/{unit_shells}/{n_obligations}; "
            f"pairs={pair_clocks}; shell transversals={transversals}"
        )

    print()
    print("KEY LIFT PENALTIES")
    for case in CASES:
        k = case.paper_k
        if case.total_n in (11, 13):
            penalty = (case.total_n / 2) ** k
            print(
                f"n={case.total_n}: direct c={case.total_n} lift costs "
                f"{penalty:.3e} times a c=2 lift for fixed input set "
                f"(paper avoids it by the odd-prime polynomial method)."
            )
        else:
            print(
                f"n={case.total_n}: paper uses composite lift chain {case.diagram}; "
                "certificate clock C=23 is prime, so the repo unit-shell residual is 0."
            )

    print()
    print("THM-398 DODGE-FIRST FILTER")
    print("n  dominance cut       all-short cap at v=n   n-clock safe gap   gap/cap")
    for case in CASES:
        n = case.total_n
        cap = 2 / (n * n)
        gap = (n - 2) / (n * n)
        print(
            f"{n:2d} v > {n - 1:2d}*max(other) "
            f"{cap:21.6f} {gap:18.6f} {gap / cap:9.2f}"
        )
    print(
        "   General interval criterion: for any chosen runner v, a component of "
        "G(S\\{v}) longer than 2/(n*v) certifies looseness."
    )
    print(
        "   For a multiple v=n*w, the only residual is all-short: every component "
        "of G(S\\{v}) has length <= 2/(n^2*w)."
    )

    print()
    print("S573 VITALI SAMPLE FOR THE INTERVAL FILTER")
    print("n  mult-of-n configs proved by B'   all-short loose residual")
    for n, pct, residual in VITALI_SAMPLE:
        print(f"{n:2d} {pct:27.1f}% {residual:26d}")
    print("   Sampled all-short tight cases: 0 for every listed n.")

    print()
    print("INTERPRETATION")
    print("- Paper exactness is prime-fiber exactness: certify J(k,p)=empty for many p.")
    print("- Repo clocks ask for a proof quotient first: n-clock, C=2n-1 shells, pair pinches, D/U/N owners.")
    print("- THM-398 inserts an interval oracle before the quotient scan: long safe components dodge any one runner.")
    print("- S573 identifies n|v as the Vitali handoff; the fallback is arc alignment, not ordinary tuple search.")
    print("- If the certificate quotient is complete, the expensive object is not p^k tuples but at most a few thousand shell/fiber rows plus O(k^2) pair tests.")
    print("- n=12 is the clearest rotated case: paper n=12 is composite, but C=23 is prime, so every antipodal shell is a unit shell.")


if __name__ == "__main__":
    main()
