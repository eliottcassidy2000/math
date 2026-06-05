#!/usr/bin/env python3
"""S644: Vieta carriers, perfect numbers, and the 14/21 aliquot shadow.

The previous S642 packet put Goldbach and Lemoine on the same prime-pair
carrier:

    E = p + q,   O = p + 2q.

This supplement follows the user's new prompt: the unique doubled/tripled-prime
crossing is 6, and perfect numbers should enter the long carrier story without
replacing the LRC proof target.

The useful object here is the triangular pair-count carrier

    A = C(n,2),       8A + 1 = (2n - 1)^2.

Thus the LRC shell clock C=2n-1 is literally the Vieta/discriminant square root
of the pair count.  Perfect numbers are fixed points of the aliquot map
s(A)=sigma(A)-A inside this triangular carrier.  The n=14 row is not perfect,
but it has a striking exact side shadow:

    A = C(14,2) = 91 = 7*13,   s(A) = 1+7+13 = 21.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt


LIMIT = 20_000


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    step = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step
    return True


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def fmt_factor(n: int) -> str:
    parts = []
    for p, e in factor(n).items():
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts) if parts else "1"


def sigma(n: int) -> int:
    total = 1
    for p, e in factor(n).items():
        total *= (p ** (e + 1) - 1) // (p - 1)
    return total


def aliquot(n: int) -> int:
    return sigma(n) - n


def choose2(n: int) -> int:
    return n * (n - 1) // 2


def square_root_or_none(n: int) -> int | None:
    r = isqrt(n)
    return r if r * r == n else None


def triangular_root_from_pair_count(a: int) -> tuple[int | None, int | None]:
    """Return (n, shell) if a=C(n,2), else (None, None)."""
    shell = square_root_or_none(8 * a + 1)
    if shell is None:
        return None, None
    if (shell + 1) % 2:
        return None, shell
    n = (shell + 1) // 2
    if choose2(n) != a:
        return None, shell
    return n, shell


@dataclass(frozen=True)
class Route:
    name: str
    preserves_lrc_clock: int
    exposes_fixed_point: int
    arithmetic_certainty: int
    proof_transfer: int
    scalar_risk: int


ROUTES = [
    Route("triangular_vieta_square_root_C_2n_minus_1", 5, 4, 5, 5, 1),
    Route("aliquot_shadow_C14_pair_count_to_21", 5, 5, 5, 5, 1),
    Route("even_perfect_numbers_as_triangular_fixed_controls", 4, 5, 5, 4, 1),
    Route("semiprime_n_2p_family_s_Cn2_equals_3p", 4, 4, 5, 4, 1),
    Route("doubled_tripled_prime_six_seam", 3, 5, 5, 4, 1),
    Route("goldbach_lemoine_diagonal_packet", 3, 4, 4, 4, 2),
    Route("raw_perfect_number_numerology", 1, 2, 3, 1, 5),
]


def route_score(route: Route) -> int:
    return (
        3 * route.preserves_lrc_clock
        + 2 * route.exposes_fixed_point
        + 2 * route.arithmetic_certainty
        + 2 * route.proof_transfer
        - 3 * route.scalar_risk
    )


def beats(a: Route, b: Route) -> bool:
    sa = route_score(a)
    sb = route_score(b)
    if sa != sb:
        return sa > sb
    return a.name < b.name


def directed_triangles(routes: list[Route]) -> int:
    total = 0
    for a in range(len(routes)):
        for b in range(a + 1, len(routes)):
            for c in range(b + 1, len(routes)):
                wins = []
                for x, y in ((a, b), (a, c), (b, c)):
                    wins.append((x, y) if beats(routes[x], routes[y]) else (y, x))
                out = Counter(x for x, _ in wins)
                if sorted(out.values()) == [1, 1, 1]:
                    total += 1
    return total


def hamiltonian_paths(routes: list[Route]) -> int:
    n = len(routes)
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(routes[last], routes[nxt]):
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_header() -> None:
    print("S644 Vieta-perfect aliquot carriers")
    print("===================================")
    print()
    print("Carrier identities")
    print("------------------")
    print("Triangular pair count: A = C(n,2)")
    print("Vieta/discriminant carrier: 8A + 1 = (2n - 1)^2")
    print("LRC shell clock: C = 2n - 1 = sqrt(8A + 1)")
    print("Aliquot map: s(A) = sigma(A) - A, the sum of proper divisors")
    print()


def print_six_seam() -> None:
    print("The unique doubled/tripled-prime crossing")
    print("-----------------------------------------")
    hits = []
    for n in range(1, LIMIT + 1):
        doubled_prime = n % 2 == 0 and is_prime(n // 2)
        tripled_prime = n % 3 == 0 and is_prime(n // 3)
        if doubled_prime and tripled_prime:
            hits.append(n)
    print(f"Search n<= {LIMIT}: {hits}")
    print("Proof: if 2p = 3q with p,q prime, then 3 divides p, so p=3.")
    print("Then 6 = 3q, hence q=2 and the common value is 6.")
    print("The same 6 is the unique distinct product-sum resonance:")
    print("  6 = 2*3 = 1+2+3 = 1*2*3")
    print(f"  sigma(6)={sigma(6)}, aliquot(6)={aliquot(6)}, 8*6+1={8*6+1}=7^2")
    print()


def print_perfect_triangular_controls() -> None:
    print("Even perfect numbers as triangular fixed controls")
    print("-------------------------------------------------")
    rows = []
    for p in range(2, 32):
        m = (1 << p) - 1
        if is_prime(p) and is_prime(m):
            n = 1 << p
            a = (1 << (p - 1)) * m
            shell = 2 * n - 1
            rows.append((p, m, n, a, shell, aliquot(a)))
    print("Euclid-Euler rows with Mersenne prime M=2^p-1:")
    print("  p   M        n=2^p     A=C(n,2)        shell      s(A)")
    for p, m, n, a, shell, s_a in rows[:7]:
        print(f"  {p:<3} {m:<8} {n:<9} {a:<15} {shell:<10} {s_a}")
    print("Here s(A)=A, so perfect numbers are fixed points of the aliquot")
    print("map inside the same triangular/Vieta carrier used by LRC clocks.")
    print()


def print_lrc_triangular_rows() -> None:
    print("LRC triangular rows and aliquot shadows")
    print("--------------------------------------")
    ns = [4, 6, 8, 10, 12, 14, 16, 18, 20, 21, 22, 32, 38, 62, 128]
    print("  n    A=C(n,2)   factor(A)        shell=2n-1   s(A)   class")
    for n in ns:
        a = choose2(n)
        s_a = aliquot(a)
        if s_a == a:
            cls = "perfect"
        elif s_a < a:
            cls = "deficient"
        else:
            cls = "abundant"
        marker = "  <== n=14 gives s(91)=21" if n == 14 else ""
        print(
            f"  {n:<4} {a:<10} {fmt_factor(a):<16} {2*n-1:<12} "
            f"{s_a:<6} {cls}{marker}"
        )
    print()
    print("Key row:")
    a14 = choose2(14)
    print(f"  C(14,2) = {a14} = {fmt_factor(a14)}")
    print(f"  s(91) = {aliquot(a14)}")
    print("  sqrt(8*91+1) = 27 = 2*14-1")
    print("So 21 is an exact divisor-sum shadow of the n=14 pair-count carrier.")
    print()


def print_semiprime_family() -> None:
    print("Semiprime family n=2p with 2p-1 prime")
    print("-------------------------------------")
    rows = []
    for p in range(2, 300):
        q = 2 * p - 1
        if is_prime(p) and is_prime(q):
            n = 2 * p
            a = choose2(n)
            rows.append((p, q, n, a, aliquot(a)))
    print("If p and q=2p-1 are prime, then C(2p,2)=p*q and")
    print("  s(C(2p,2)) = 1 + p + q = 3p.")
    print("This is a whole family where the triangular pair count casts a")
    print("tripled-prime aliquot shadow.")
    print()
    print("  p    q=2p-1   n=2p   A=C(n,2)   s(A)=3p   diagonal packet")
    for p, q, n, a, s_a in rows[:12]:
        marker = "  <== p=7 gives n=14, shadow 21" if p == 7 else ""
        print(f"  {p:<4} {q:<8} {n:<6} {a:<11} {s_a:<9} (2p,3p)=({2*p},{3*p}){marker}")
    print()


def print_vieta_examples() -> None:
    print("Vieta/discriminant examples")
    print("---------------------------")
    for a in [6, 15, 28, 91, 231, 496, 703, 8128]:
        n, shell = triangular_root_from_pair_count(a)
        print(
            f"  A={a:<5} factor={fmt_factor(a):<12} "
            f"8A+1={8*a+1:<7} sqrt={shell!s:<4} n={n!s:<4} s(A)={aliquot(a)}"
        )
    print()
    print("This is the Vieta carrier x^2+x-2A=0.  For A=C(n,2),")
    print("the positive root is x=n-1 and the discriminant root is 2n-1.")
    print()


def print_tournament_analysis() -> None:
    print("Tournament Analysis over proof lenses")
    print("-------------------------------------")
    routes = list(ROUTES)
    scores = {r.name: sum(1 for other in routes if r is not other and beats(r, other)) for r in routes}
    hist = Counter(scores.values())
    ranking = sorted(routes, key=lambda r: (-scores[r.name], r.name))
    print(f"vertices={len(routes)}")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={directed_triangles(routes)}")
    print(f"hamiltonian_paths={hamiltonian_paths(routes)}")
    print("tie Hamiltonian path:")
    for route in ranking:
        print(f"  score={scores[route.name]} {route.name}")
    print()
    print("Pairwise observable: which lens better preserves the hidden carrier")
    print("needed for LRC-style proof transfer: pair-count square root, aliquot")
    print("fixed/shadow data, prime 2/3/7 seams, and side-channel risk.")
    print("Vertices are proof lenses, not runners, primes, or perfect numbers.")


def main() -> None:
    print_header()
    print_six_seam()
    print_perfect_triangular_controls()
    print_lrc_triangular_rows()
    print_semiprime_family()
    print_vieta_examples()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
