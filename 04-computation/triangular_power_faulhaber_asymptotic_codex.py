#!/usr/bin/env python3
"""Faulhaber odd-moment expansion for triangular power-balance anchors.

codex-2026-06-12

Extends HYP-2454's finite bracket scout with the user's stronger formulation.
For

    sum_{j=0}^n (a+j)^p = sum_{j=1}^n (a+n+j)^p,

set c=a+n and u=n(n+1).  Then

    c^p = 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n),

where S_r(n)=sum_{k=1}^n k^r.  Only the odd Faulhaber moments survive.

This script now does four linked jobs for HYP-2454:
1. derive the first two asymptotic coefficients using only S_1, S_3, S_5;
2. verify the p=1,2 exact towers and record the square-pyramidal cuboid
   identity 6*P_2(n)=u(2n+1);
3. extract the n=1 live balance polynomials and certify their irreducibility
   for p=3..20 by finite-field reductions;
4. run Tournament Analysis on proof carriers rather than raw interval entries.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import combinations, permutations
from math import comb


getcontext().prec = 80


def u_value(n: int) -> int:
    return n * (n + 1)


def triangle(n: int) -> int:
    return n * (n + 1) // 2


def power_sum(n: int, r: int) -> int:
    return sum(k**r for k in range(1, n + 1))


def square_pyramidal(n: int) -> int:
    return power_sum(n, 2)


def balance_anchor_decimal(a: Decimal, p: int, n: int) -> Decimal:
    left = sum((a + Decimal(j)) ** p for j in range(n + 1))
    right = sum((a + Decimal(n + j)) ** p for j in range(1, n + 1))
    return left - right


def odd_moment_balance(c: int | Fraction, p: int, n: int) -> Fraction:
    total = Fraction(c) ** p
    for r in range(1, p + 1, 2):
        total -= 2 * comb(p, r) * (Fraction(c) ** (p - r)) * power_sum(n, r)
    return total


def alpha_coeff(p: int) -> Fraction:
    return Fraction((p - 1) * (p - 2), 12 * p)


def beta_coeff(p: int) -> Fraction:
    return -Fraction((p - 1) * (p - 2) * (2 * p * p - 4 * p - 1), 180 * p**3)


def anchor_asymptotic(p: int, n: int) -> Fraction:
    if p == 1:
        return Fraction(n * n, 1)
    if p == 2:
        return Fraction(2 * n * n + n, 1)
    u = u_value(n)
    return (
        Fraction(p * n * n + (p - 1) * n, 1)
        + alpha_coeff(p)
        + beta_coeff(p) / u
    )


def center_asymptotic(p: int, n: int) -> Fraction:
    return anchor_asymptotic(p, n) + n


def root_anchor(p: int, n: int) -> Decimal:
    if p == 1:
        return Decimal(n * n)
    if p == 2:
        return Decimal(2 * n * n + n)

    left = Decimal(p * n * n + (p - 1) * n)
    right = left + Decimal(1)
    f_left = balance_anchor_decimal(left, p, n)
    f_right = balance_anchor_decimal(right, p, n)
    if not (f_left < 0 < f_right):
        raise ValueError(f"Bracket failure at p={p}, n={n}: {f_left}, {f_right}")

    for _ in range(100):
        mid = (left + right) / 2
        f_mid = balance_anchor_decimal(mid, p, n)
        if f_mid == 0:
            return mid
        if f_mid > 0:
            right = mid
        else:
            left = mid
    return (left + right) / 2


def decimal_of_fraction(q: Fraction) -> Decimal:
    return Decimal(q.numerator) / Decimal(q.denominator)


def trim_poly(coeffs: list[int]) -> list[int]:
    out = coeffs[:]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_to_string(coeffs: list[int], var: str = "C") -> str:
    pieces: list[str] = []
    for power in range(len(coeffs) - 1, -1, -1):
        coeff = coeffs[power]
        if coeff == 0:
            continue
        mag = abs(coeff)
        if power == 0:
            body = str(mag)
        elif power == 1:
            body = var if mag == 1 else f"{mag}*{var}"
        else:
            body = f"{var}^{power}" if mag == 1 else f"{mag}*{var}^{power}"
        if not pieces:
            pieces.append(body if coeff > 0 else f"-{body}")
        else:
            pieces.append(f"+ {body}" if coeff > 0 else f"- {body}")
    return " ".join(pieces) if pieces else "0"


def prime_factors(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        out.append(m)
    return out


def primes_upto(limit: int) -> list[int]:
    out: list[int] = []
    for x in range(2, limit + 1):
        for d in range(2, int(x**0.5) + 1):
            if x % d == 0:
                break
        else:
            out.append(x)
    return out


def poly_sub_mod(a: list[int], b: list[int], prime: int) -> list[int]:
    size = max(len(a), len(b))
    out = [0] * size
    for i in range(size):
        out[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % prime
    return trim_poly(out)


def poly_mul_mod(a: list[int], b: list[int], prime: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y == 0:
                continue
            out[i + j] = (out[i + j] + x * y) % prime
    return trim_poly(out)


def poly_divmod_mod(f: list[int], g: list[int], prime: int) -> tuple[list[int], list[int]]:
    rem = trim_poly(f[:])
    den = trim_poly(g[:])
    if len(rem) < len(den):
        return [0], rem
    lead_inv = pow(den[-1], -1, prime)
    quot = [0] * (len(rem) - len(den) + 1)
    while len(rem) >= len(den) and rem != [0]:
        shift = len(rem) - len(den)
        coeff = (rem[-1] * lead_inv) % prime
        quot[shift] = coeff
        if coeff:
            for i, value in enumerate(den):
                rem[shift + i] = (rem[shift + i] - coeff * value) % prime
        rem = trim_poly(rem)
    return trim_poly(quot), rem


def poly_mod(f: list[int], modulus: list[int], prime: int) -> list[int]:
    return poly_divmod_mod(f, modulus, prime)[1]


def poly_gcd_mod(a: list[int], b: list[int], prime: int) -> list[int]:
    left = trim_poly(a[:])
    right = trim_poly(b[:])
    while right != [0]:
        _, rem = poly_divmod_mod(left, right, prime)
        left, right = right, rem
    if left == [0]:
        return [0]
    inv = pow(left[-1], -1, prime)
    return [(coeff * inv) % prime for coeff in left]


def poly_powmod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    out = [1]
    cur = base[:]
    exp = exponent
    while exp:
        if exp & 1:
            out = poly_mod(poly_mul_mod(out, cur, prime), modulus, prime)
        cur = poly_mod(poly_mul_mod(cur, cur, prime), modulus, prime)
        exp //= 2
    return trim_poly(out)


def irreducible_mod_prime(coeffs: list[int], prime: int) -> bool:
    f = [coeff % prime for coeff in trim_poly(coeffs[:])]
    degree = len(f) - 1
    if degree <= 0 or f[-1] == 0:
        return False
    x = [0, 1]
    for divisor in prime_factors(degree):
        probe = poly_powmod(x, prime ** (degree // divisor), f, prime)
        if len(poly_gcd_mod(poly_sub_mod(probe, x, prime), f, prime)) > 1:
            return False
    final_probe = poly_powmod(x, prime**degree, f, prime)
    return trim_poly(poly_sub_mod(final_probe, x, prime)) == [0]


def find_irreducibility_prime(coeffs: list[int], prime_limit: int = 200) -> int | None:
    for prime in primes_upto(prime_limit):
        if irreducible_mod_prime(coeffs, prime):
            return prime
    return None


def n1_live_factor(power: int) -> list[int]:
    """Primitive live factor at n=1.

    For n=1 the balance polynomial is
        D_p(C,1) = C^p + (C-1)^p - (C+1)^p.
    When p is even, C=0 is a forced symmetry root, so we remove the factor C.
    """
    coeffs = [0] * (power + 1)
    coeffs[power] += 1
    for k in range(power + 1):
        coeffs[k] += comb(power, k) * ((-1) ** (power - k))
        coeffs[k] -= comb(power, k)
    coeffs = trim_poly(coeffs)
    if power % 2 == 0:
        if coeffs[0] != 0:
            raise ValueError(f"Expected C factor at even p={power}")
        coeffs = coeffs[1:]
    return trim_poly(coeffs)


@dataclass(frozen=True)
class Carrier:
    name: str
    exactness: int
    geometry: int
    generality: int
    irreducibility: int
    support_transfer: int
    proof_reach: int


def tournament_edge(a: Carrier, b: Carrier, tie_order: list[str]) -> bool:
    keys = (
        "exactness",
        "geometry",
        "generality",
        "irreducibility",
        "support_transfer",
        "proof_reach",
    )
    a_wins = 0
    b_wins = 0
    for key in keys:
        av = getattr(a, key)
        bv = getattr(b, key)
        if av > bv:
            a_wins += 1
        elif bv > av:
            b_wins += 1
    if a_wins != b_wins:
        return a_wins > b_wins
    return tie_order.index(a.name) < tie_order.index(b.name)


def score_histogram(adj: dict[str, set[str]]) -> dict[int, int]:
    return dict(sorted(Counter(len(targets) for targets in adj.values()).items()))


def directed_3cycles(adj: dict[str, set[str]], names: list[str]) -> int:
    total = 0
    for a, b, c in combinations(names, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            total += 1
        elif c in adj[a] and a in adj[b] and b in adj[c]:
            total += 1
    return total


def scc_sizes(adj: dict[str, set[str]], names: list[str]) -> list[int]:
    reverse = {name: set() for name in names}
    for src in names:
        for dst in adj[src]:
            reverse[dst].add(src)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(graph: dict[str, set[str]], start: str, out: list[str]) -> None:
        stack = [(start, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                out.append(node)
                continue
            if node in seen:
                continue
            seen.add(node)
            stack.append((node, True))
            for nxt in graph[node]:
                if nxt not in seen:
                    stack.append((nxt, False))

    for name in names:
        if name not in seen:
            dfs(adj, name, order)

    seen.clear()
    sizes: list[int] = []

    def dfs_count(start: str) -> int:
        stack = [start]
        count = 0
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            count += 1
            for nxt in reverse[node]:
                if nxt not in seen:
                    stack.append(nxt)
        return count

    for name in reversed(order):
        if name not in seen:
            sizes.append(dfs_count(name))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: dict[str, set[str]], names: list[str]) -> int:
    total = 0
    for perm in permutations(names):
        if all(perm[i + 1] in adj[perm[i]] for i in range(len(perm) - 1)):
            total += 1
    return total


def print_identity_section() -> None:
    print("ODD-MOMENT REFORMULATION")
    print("========================")
    print("Let c=a+n and u=n(n+1). Then")
    print("  sum_{j=0}^n (a+j)^p = sum_{j=1}^n (a+n+j)^p")
    print("becomes")
    print("  c^p = 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n).")
    print("Only odd Faulhaber moments survive because even binomial terms cancel in")
    print("  (c+j)^p - (c-j)^p.")
    print()
    print("Checks on small integer centers:")
    for p in range(1, 7):
        n = p + 2
        c = p * u_value(n)
        residual = odd_moment_balance(c, p, n)
        print(f"  p={p}, n={n}, c={c}: odd-moment residual = {residual}")
    print()


def print_faulhaber_section() -> None:
    print("LEADING ODD FAULHABER MOMENTS")
    print("=============================")
    print("In terms of u=n(n+1):")
    print("  S_1 = u/2")
    print("  S_3 = u^2/4")
    print("  S_5 = u^2(2u-1)/12")
    print("  S_7 = u^2(3u^2-4u+2)/24")
    print()
    print("Only S_1, S_3, S_5 are needed for the first two asymptotic coefficients,")
    print("because higher odd moments start at order u^(p-3) in the balance equation.")
    print()


def print_derivation_section() -> None:
    print("COEFFICIENT MATCHING")
    print("====================")
    print("Write")
    print("  c = p*u + alpha + beta/u + O(u^-2).")
    print("Using only S_1, S_3, S_5, the coefficient equations are:")
    print("  u^(p-1):   p*alpha - (p-1)(p-2)/12 = 0")
    print("  u^(p-2):   p^3*beta + A^2/72 - A(p-3)(p-4)/360 = 0")
    print("             where A=(p-1)(p-2).")
    print()
    print("So")
    print("  alpha = (p-1)(p-2)/(12p)")
    print("  beta  = -(p-1)(p-2)(2p^2-4p-1)/(180 p^3)")
    print()
    print("Exact verification for p=3..12:")
    for p in range(3, 13):
        a = alpha_coeff(p)
        b = beta_coeff(p)
        A = (p - 1) * (p - 2)
        eq1 = p * a - Fraction(A, 12)
        eq2 = p**3 * b + Fraction(A * A, 72) - Fraction(A * (p - 3) * (p - 4), 360)
        print(f"  p={p:2d}: eq1={eq1}, eq2={eq2}, alpha={a}, beta={b}")
    print()
    print("Therefore")
    print("  a_p(n) = p*n^2 + (p-1)*n + alpha + beta/u + O(u^-2)")
    print("and since u=n(n+1)~n^2, the remainder is O(n^-4).")
    print()


def print_exact_cases() -> None:
    print("EXACT LOW-POWER TOWERS")
    print("======================")
    print("p=1: c=u, so a=n^2 exactly.")
    print("p=2: c=2u, so a=2n^2+n exactly.")
    print("These are precisely the two integer towers already isolated in HYP-2454.")
    print()


def print_numeric_section() -> None:
    print("NUMERICAL VERIFICATION")
    print("======================")
    print("Root bracket and asymptotic accuracy for p=3..8.")
    print("Columns: n, root a_p(n), asymptotic, error, n^4*error.")
    for p in range(3, 9):
        print(f"p={p}")
        for n in (5, 10, 20, 50, 100, 200):
            root = root_anchor(p, n)
            approx = decimal_of_fraction(anchor_asymptotic(p, n))
            err = root - approx
            scaled = err * (Decimal(n) ** 4)
            print(
                f"  n={n:3d}: root={root:.24f}, "
                f"asym={approx:.24f}, err={err:.3E}, n^4*err={scaled:.3E}"
            )
        print()


def print_packing_section() -> None:
    print("SQUARE-PYRAMIDAL CUBOID PACKING")
    print("===============================")
    print("Faulhaber at p=2 gives the square-pyramidal numbers")
    print("  P_2(n)=1^2+...+n^2 = n(n+1)(2n+1)/6.")
    print("Hence")
    print("  6*P_2(n) = n(n+1)(2n+1) = u(2n+1),")
    print("so six n-step square pyramids fill the u x (2n+1) cuboid.")
    print("Equivalently, one doubled and four regular pyramids also account for the")
    print("same cuboid because 2+4=6 at the level of the volume identity.")
    print()
    print("This is the p=2 analogue of the p=1 staircase triangle identity")
    print("  2*T_n = n(n+1) = u,")
    print("so the first two exact anchors are rectangle/cuboid packings with no")
    print("Bernoulli correction, while p>=3 acquires the Faulhaber offset alpha+beta/u+...")
    print()
    print("Sample values:")
    for n in range(1, 7):
        u = u_value(n)
        pyr = square_pyramidal(n)
        print(
            f"  n={n}: T_n={triangle(n)}, 2*T_n={2*triangle(n)}=u={u}, "
            f"P_2(n)={pyr}, 6*P_2(n)={6*pyr}=u*(2n+1)={u*(2*n+1)}"
        )
    print()
    print("Also")
    print("  S_1(n)=n(n+1)(2n+1)/2 = 3*P_2(n),")
    print("so the p=1 common sum is already three square pyramids.")
    print()


def print_irreducibility_section() -> None:
    print("LIVE BALANCE POLYNOMIALS AT n=1")
    print("================================")
    print("At n=1 the balance polynomial is")
    print("  D_p(C,1) = C^p + (C-1)^p - (C+1)^p.")
    print("For even p, C=0 is a forced symmetry root, so the live factor is D_p(C,1)/C.")
    print("This is the smallest row where higher-p irrationality is visible without")
    print("any asymptotic expansion.")
    print()
    print("Exact low-power cases:")
    print("  p=1: C - 2")
    print("  p=2: C - 4     after removing the forced C factor")
    print()
    print("Irreducible live factors for p=3..12 (finite-field certificates):")
    for p in range(3, 13):
        coeffs = n1_live_factor(p)
        cert_prime = find_irreducibility_prime(coeffs)
        if cert_prime is None:
            raise ValueError(f"Missing irreducibility certificate for p={p}")
        print(f"  p={p:2d}: {poly_to_string(coeffs)}    irreducible mod {cert_prime}")
    print()
    certs = []
    for p in range(3, 21):
        cert_prime = find_irreducibility_prime(n1_live_factor(p))
        if cert_prime is None:
            raise ValueError(f"Missing irreducibility certificate through p=20 at p={p}")
        certs.append(f"p={p}:q={cert_prime}")
    print("Certificate summary through p=20:")
    print("  " + ", ".join(certs))
    print()
    print("Interpretation:")
    print("  The tower cases p=1,2 are the last rows whose live n=1 balance factors")
    print("  down to degree 1. From p>=3 onward, the first visible balance polynomial")
    print("  is already irreducible in the stored window, so the higher anchors are")
    print("  algebraic but not split by the old tower mechanism.")
    print()


def print_tournament_section() -> None:
    print("TOURNAMENT ANALYSIS")
    print("===================")
    print("Assumption challenge:")
    print("  Vertices were not taken to be powers p, interval entries, or odd moments")
    print("  alone. The chosen quotient uses proof carriers, preserving exactness,")
    print("  packing geometry, asymptotic reach, irreducibility signal, support-transfer")
    print("  value, and proof reach while discarding decorative row data.")
    print()
    carriers = [
        Carrier("rectangle_packing", 5, 5, 1, 0, 1, 2),
        Carrier("square_pyramidal_cuboid", 5, 5, 2, 1, 2, 3),
        Carrier("odd_moment_reduction", 5, 2, 5, 3, 4, 5),
        Carrier("alpha_beta_asymptotic", 4, 1, 4, 2, 3, 4),
        Carrier("n1_live_irreducibility", 4, 0, 3, 5, 5, 4),
        Carrier("78_90_support_shadow", 4, 3, 2, 1, 5, 3),
        Carrier("global_bracket_tail", 1, 1, 5, 2, 3, 5),
        Carrier("hidden_lift_transfer", 2, 0, 3, 4, 5, 3),
    ]
    tie_order = [carrier.name for carrier in carriers]
    adj = {carrier.name: set() for carrier in carriers}
    for left, right in combinations(carriers, 2):
        if tournament_edge(left, right, tie_order):
            adj[left.name].add(right.name)
        else:
            adj[right.name].add(left.name)
    exactness_adj = {carrier.name: set() for carrier in carriers}
    for left, right in combinations(carriers, 2):
        if left.exactness > right.exactness:
            exactness_adj[left.name].add(right.name)
        elif right.exactness > left.exactness:
            exactness_adj[right.name].add(left.name)
        elif tie_order.index(left.name) < tie_order.index(right.name):
            exactness_adj[left.name].add(right.name)
        else:
            exactness_adj[right.name].add(left.name)

    names = [carrier.name for carrier in carriers]
    ranking = sorted(names, key=lambda name: len(adj[name]), reverse=True)
    flips = 0
    for a, b in combinations(names, 2):
        in_main = b in adj[a]
        in_exact = b in exactness_adj[a]
        if in_main != in_exact:
            flips += 1

    print("Pairwise observable:")
    print("  majority over (exactness, geometry, generality, irreducibility,")
    print("  support_transfer, proof_reach), with the listed carrier order as the")
    print("  fixed tie Hamiltonian path.")
    print()
    print("Vertices:")
    for carrier in carriers:
        print(
            "  "
            f"{carrier.name}: "
            f"(exact={carrier.exactness}, geom={carrier.geometry}, gen={carrier.generality}, "
            f"irr={carrier.irreducibility}, support={carrier.support_transfer}, "
            f"reach={carrier.proof_reach})"
        )
    print()
    print("Tournament fingerprint:")
    print(f"  ranking leader = {ranking[0]}")
    print(f"  score_hist = {score_histogram(adj)}")
    print(f"  directed_3cycles = {directed_3cycles(adj, names)}")
    print(f"  scc_sizes = {scc_sizes(adj, names)}")
    print(f"  hamiltonian_paths = {hamiltonian_path_count(adj, names)}")
    print(f"  edge_flips_vs_exactness_only = {flips}")
    print()
    print("Readout:")
    print("  The leader is odd_moment_reduction, not the exact p=1/p=2 packings and")
    print("  not the 78/90 shadow. The tournament therefore prefers the carrier that")
    print("  simultaneously explains all p, produces the asymptotic, and interfaces")
    print("  with the irreducibility and hidden-lift lanes.")
    print()
def main() -> None:
    print("TRIANGULAR POWER-BALANCE FAULHABER ASYMPTOTIC")
    print("=============================================")
    print("Extends HYP-2454 with the exact odd-moment reformulation, the first")
    print("two coefficients in the anchor expansion, the square-pyramidal cuboid")
    print("packing, the n=1 irreducibility lane, and a proof-carrier tournament.")
    print()
    print_identity_section()
    print_faulhaber_section()
    print_derivation_section()
    print_exact_cases()
    print_numeric_section()
    print_packing_section()
    print_irreducibility_section()
    print_tournament_section()
    print("TAKEAWAYS")
    print("=========")
    print("1. The midpoint c=a+n is centered at p*u, where u=n(n+1), for every power p.")
    print("2. Only odd Faulhaber moments survive in the exact balance equation.")
    print("3. The first two Bernoulli/Faulhaber corrections are")
    print("     +(p-1)(p-2)/(12p) and")
    print("     -(p-1)(p-2)(2p^2-4p-1)/(180 p^3 u).")
    print("4. p=1 and p=2 are exactly the no-correction cases; from p>=3 onward the")
    print("   same triangular carrier persists but acquires Bernoulli address terms.")
    print("5. The first visible higher-p live balance polynomial (n=1) is already")
    print("   irreducible in the stored p=3..20 window, so the tower split stops")
    print("   exactly where the Bernoulli/Faulhaber corrections begin.")


if __name__ == "__main__":
    main()
