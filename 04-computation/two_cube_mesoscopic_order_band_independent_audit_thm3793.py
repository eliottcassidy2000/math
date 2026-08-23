#!/usr/bin/env python3
"""Independent hostile audit of the two-cube mesoscopic product-order band.

This checker deliberately does not import the provisional primary.  It uses
target-by-target cube lookup for the singleton fibres, explicit tuple
enumeration for the collision union bound, and SciPy's Poisson CDF for the
asymptotic numerical controls.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import sys
from fractions import Fraction

from scipy.special import pdtr


sys.stdout.reconfigure(encoding="utf-8", newline="\n")
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def primes_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [n for n in range(2, limit + 1) if sieve[n]]


def valuation(n: int, p: int) -> int:
    answer = 0
    while n % p == 0:
        answer += 1
        n //= p
    return answer


def representations(target: int, coordinate_cap: int, cubes: dict[int, int]) -> list[tuple[int, int]]:
    """Find every positive u<v with u^3+v^3=target.

    For a target built from a row of sum d, target<d^3, so every competing
    coordinate is <d.  Therefore coordinate_cap=d-1 is a proved complete
    universe, rather than a construction-family lookup.
    """

    result: list[tuple[int, int]] = []
    for u in range(1, coordinate_cap + 1):
        v = cubes.get(target - u**3)
        if v is not None and u < v <= coordinate_cap:
            result.append((u, v))
    return result


def elementary_by_subsets(weights: list[Fraction], degree: int) -> Fraction:
    total = Fraction(0)
    for subset in itertools.combinations(weights, degree):
        product = Fraction(1)
        for weight in subset:
            product *= weight
        total += product
    return total


# ---------------------------------------------------------------------------
# 1. Global singleton fibres and cross-order disjointness on a complete bank.
# ---------------------------------------------------------------------------

bank = [p for p in primes_to(17) if p >= 5 and p % 3 == 2]
gate(bank == [5, 11, 17], "independent inert-prime sieve")
rows: list[tuple[int, int, tuple[int, ...]]] = []
for order in range(1, len(bank) + 1):
    for prime_set in itertools.combinations(bank, order):
        rows.append((order, math.prod(prime_set), prime_set))

max_d = max(d for _, d, _ in rows)
cube_lookup = {v**3: v for v in range(1, max_d + 1)}
seen: dict[int, tuple[int, int, int]] = {}
layer_counts = {order: 0 for order in range(1, len(bank) + 1)}
layer_surrogates = {order: Fraction(0) for order in range(1, len(bank) + 1)}
actual_layer_mass = {order: 0.0 for order in range(1, len(bank) + 1)}

for order, d, prime_set in rows:
    gate(d % 2 == 1, f"odd row d={d}")
    gate(math.prod(prime_set) == d and len(set(prime_set)) == order, f"squarefree row d={d}")
    row_count = 0
    for x in range(1, (d - 1) // 2 + 1):
        y = d - x
        m = x**3 + y**3
        g = math.gcd(x, y)
        s = d // g
        gate(x < y and x + y == d, f"ordered split d={d},x={x}")
        gate(m < d**3, f"strict height d={d},x={x}")
        gate(all(p % 3 == 2 and valuation(s, p) <= 2 for p in prime_set), f"singleton hypotheses d={d},x={x}")

        # Recheck the local valuation mechanism in primitive coordinates.
        X, Y = x // g, y // g
        primitive_cofactor = X * X - X * Y + Y * Y
        for p in prime_set:
            alpha = valuation(g, p)
            exponent = valuation(s, p)
            gate(primitive_cofactor % p != 0, f"inert cofactor d={d},x={x},p={p}")
            gate(valuation(m, p) == 3 * alpha + exponent, f"valuation invoice d={d},x={x},p={p}")

        fibre = representations(m, d - 1, cube_lookup)
        gate(fibre == [(x, y)], f"complete global singleton d={d},x={x}")
        gate(m not in seen, f"cross-order/cross-row collision m={m}")
        seen[m] = (order, x, y)
        layer_counts[order] += 1
        layer_surrogates[order] += Fraction(1, d * d)
        actual_layer_mass[order] += m ** (-2.0 / 3.0)
        row_count += 1
    gate(row_count == (d - 1) // 2, f"complete odd row d={d}")
    gate(Fraction(row_count, d * d) >= Fraction(2, 5 * d), f"row floor d={d}")

gate(len(seen) == 644, "constructed support size")
gate(set(order for order, _, _ in seen.values()) == {1, 2, 3}, "all orders represented")

weights = [Fraction(1, p) for p in bank]
elementary_floor_all = Fraction(0)
for J in range(1, 4):
    elementary_floor = sum(
        (Fraction(2, 5) * elementary_by_subsets(weights, j) for j in range(1, J + 1)),
        Fraction(0),
    )
    if J == 3:
        elementary_floor_all = elementary_floor
    surrogate = sum((layer_surrogates[j] for j in range(1, J + 1)), Fraction(0))
    actual = sum(actual_layer_mass[j] for j in range(1, J + 1))
    gate(surrogate >= elementary_floor, f"finite all-order surrogate J={J}")
    gate(actual > float(surrogate), f"strict weighted mass J={J}")
    for m, (order, _, _) in seen.items():
        if order <= J:
            gate(m < 17 ** (3 * J), f"all-order height J={J},m={m}")
gate(elementary_floor_all == Fraction(722, 4675), "exact finite elementary floor")


# ---------------------------------------------------------------------------
# 2. Ordered-tuple collision inequality by an implementation independent of
#    the coefficient recursion in the provisional primary.
# ---------------------------------------------------------------------------

tuple_bank = [p for p in primes_to(41) if p >= 5 and p % 3 == 2]
tuple_weights = [Fraction(1, p) for p in tuple_bank]
A_small = sum(tuple_weights, Fraction(0))
B_small = sum((w * w for w in tuple_weights), Fraction(0))
for j in range(2, 6):
    total_tuple_mass = Fraction(0)
    distinct_tuple_mass = Fraction(0)
    collision_mass = Fraction(0)
    for indices in itertools.product(range(len(tuple_weights)), repeat=j):
        weight = Fraction(1)
        for index in indices:
            weight *= tuple_weights[index]
        total_tuple_mass += weight
        if len(set(indices)) == j:
            distinct_tuple_mass += weight
        else:
            collision_mass += weight
    e_j = elementary_by_subsets(tuple_weights, j)
    union_budget = math.comb(j, 2) * B_small * A_small ** (j - 2)
    gate(total_tuple_mass == A_small**j, f"tuple total j={j}")
    gate(distinct_tuple_mass == math.factorial(j) * e_j, f"distinct tuple normalization j={j}")
    gate(total_tuple_mass == distinct_tuple_mass + collision_mass, f"tuple partition j={j}")
    gate(collision_mass <= union_budget, f"collision union bound j={j}")
    gate(e_j >= (A_small**j - union_budget) / math.factorial(j), f"elementary lower bound j={j}")

large_bank = [p for p in primes_to(101) if p >= 5 and p % 3 == 2]
large_weights = [Fraction(1, p) for p in large_bank]
A_large = sum(large_weights, Fraction(0))
B_large = sum((w * w for w in large_weights), Fraction(0))
gate(B_large < Fraction(1, 4), "finite square-sum control")
for j in range(2, len(large_weights) + 1):
    e_j = elementary_by_subsets(large_weights, j)
    lower = (A_large**j - math.comb(j, 2) * B_large * A_large ** (j - 2)) / math.factorial(j)
    gate(e_j >= lower, f"large-bank collision inequality j={j}")

# The one-sided proxy step is purely monotone once A>=a>=J.  Check it with
# exact rationals across a hostile grid in which A-a and B vary independently.
for a_int in (8, 13, 25, 50):
    for J in range(2, a_int + 1):
        for delta in (Fraction(0), Fraction(1, 11), Fraction(7, 5)):
            proxy_A = Fraction(a_int) + delta
            for proxy_B in (Fraction(0), Fraction(1, 7), Fraction(999, 4000)):
                for j in range(2, J + 1):
                    collision = math.comb(j, 2) * proxy_B * proxy_A ** (j - 2)
                    exact_proxy = (proxy_A**j - collision) / math.factorial(j)
                    crude_A = proxy_A**j * (1 - Fraction(j * j, 8) / proxy_A**2) / math.factorial(j)
                    floor_a = Fraction(7, 8) * Fraction(a_int**j, math.factorial(j))
                    gate(exact_proxy >= crude_A, f"B-only reduction a={a_int},J={J},j={j}")
                    gate(crude_A >= floor_a, f"A-only monotonicity a={a_int},J={J},j={j}")

# The sharpened cutoff lies slightly above a.  The same argument remains
# legal with the positive uniform factor 1-J^2/(8a^2), which tends to 7/8.
for a_int in (50, 100, 200):
    J = a_int + math.floor(a_int ** (2 / 3))
    uniform_factor = 1 - Fraction(J * J, 8 * a_int * a_int)
    gate(uniform_factor > 0, f"above-mean uniform factor a={a_int}")
    for delta in (Fraction(0), Fraction(1, 11), Fraction(7, 5)):
        proxy_A = Fraction(a_int) + delta
        for proxy_B in (Fraction(0), Fraction(1, 7), Fraction(999, 4000)):
            for j in range(2, J + 1):
                collision = math.comb(j, 2) * proxy_B * proxy_A ** (j - 2)
                exact_proxy = (proxy_A**j - collision) / math.factorial(j)
                strong_floor = uniform_factor * Fraction(a_int**j, math.factorial(j))
                gate(exact_proxy >= strong_floor, f"above-mean proxy a={a_int},J={J},j={j}")


# ---------------------------------------------------------------------------
# 3. Floor/window, Poisson normalization, height conversion, and constants.
# ---------------------------------------------------------------------------

C_P = 14 / math.e + 0.5 + (math.log(3) + 2 / 9) / 2
normal_band = 0.5 - 0.5 * (1 + math.erf(-1 / math.sqrt(2)))
c_star = (7 / 20) * math.exp(-C_P) * math.sqrt(2 / 3) * normal_band
c_sharp = (7 / 20) * math.exp(-C_P) * math.sqrt(2 / 3)
gate(abs(c_star - 1.772150772662629e-4) < 2e-18, "displayed c_star")
gate(c_sharp > 2.9 * c_star, "full-lower-tail strengthening")

candidate_samples: list[tuple[float, int, float, int, float, float, float]] = []
strong_samples: list[tuple[float, int, float, float, float, float]] = []
for L in (1_000_000.0, 10_000_000.0, 100_000_000.0, 1_000_000_000.0, 10_000_000_000.0):
    # Candidate floor and one-standard-deviation window.
    J = math.floor(L / 2 - math.log(L))
    a = 0.5 * (L - math.log(3 * J)) - C_P
    j0 = J - math.floor(math.sqrt(J))
    upper_z = (J - a) / math.sqrt(a)
    lower_z_inclusive = (j0 - 1 - a) / math.sqrt(a)
    band_probability = float(pdtr(J, a) - pdtr(j0 - 1, a))
    scaling = math.exp(-C_P) * math.sqrt(L / (3 * J))
    gate(2 <= j0 <= J <= a, f"candidate floor ordering L={L}")
    gate(1 - J * J / (8 * a * a) >= Fraction(7, 8), f"candidate collision factor L={L}")
    gate(0 < band_probability < 1, f"candidate Poisson probability L={L}")
    gate(scaling > 0, f"candidate height scaling L={L}")
    candidate_samples.append((L, J, a, j0, upper_z, lower_z_inclusive, band_probability))

    # A stronger legal order cutoff: it is above the proxy mean by much more
    # than sqrt(a), but differs from it by o(a).  Hence the entire artificial
    # Poisson lower tail is captured while the collision factor still tends
    # to 7/8.
    J_plus = math.floor(L / 2 - 0.5 * math.log(L) + L ** (2 / 3))
    a_plus = 0.5 * (L - math.log(3 * J_plus)) - C_P
    factor_plus = 1 - J_plus * J_plus / (8 * a_plus * a_plus)
    full_probability = float(pdtr(J_plus, a_plus) - pdtr(1, a_plus))
    scaled_plus = (2 / 5) * factor_plus * full_probability * math.exp(-C_P) * math.sqrt(L / (3 * J_plus))
    gate(J_plus >= 2 and a_plus > 0, f"strong cutoff positivity L={L}")
    gate(J_plus / a_plus < math.sqrt(8), f"strong collision positivity L={L}")
    gate(0 < factor_plus < Fraction(7, 8), f"strong factor approaches below L={L}")
    gate(0 < full_probability <= 1, f"strong Poisson lower tail L={L}")
    strong_samples.append((L, J_plus, a_plus, factor_plus, full_probability, scaled_plus))

gate(abs(candidate_samples[-1][4]) < 0.002, "candidate upper endpoint tends zero")
gate(abs(candidate_samples[-1][5] + 1) < 0.002, "candidate lower endpoint tends minus one")
gate(abs(candidate_samples[-1][6] - normal_band) < 3e-5, "candidate Poisson CLT")
gate(abs(strong_samples[-1][3] - 7 / 8) < 3e-4, "strong collision limit")
gate(strong_samples[-1][4] > 0.999999, "strong lower-tail probability")
gate(strong_samples[-1][5] > 0.995 * c_sharp, "strong scaled constant")

semantic = (
    "global singleton fibres are complete for every split in each tested squarefree inert row",
    "different prime products and different product orders cannot double count",
    "the ordered-tuple collision inequality has j! elementary-symmetric normalization",
    "the candidate cutoff gives standardized Poisson endpoints minus one and zero",
    "the Poisson law is an artificial exponential-series identity, not Euler-product Bernoulli normalization",
    "the stated c_star is valid, and a legal upper cutoff captures the full lower tail for a stronger constant",
)
semantic_hash = hashlib.sha256(repr(semantic).encode()).hexdigest()
print("MESOSCOPIC_INDEPENDENT_SINGLETON", f"rows={len(rows)};values={len(seen)};orders={sorted(layer_counts)}")
print("MESOSCOPIC_LAYER_COUNTS", ";".join(f"j={j}:{layer_counts[j]}" for j in sorted(layer_counts)))
print("MESOSCOPIC_SURROGATE_TOTAL", str(sum(layer_surrogates.values(), Fraction(0))))
print("MESOSCOPIC_ELEMENTARY_FLOOR", str(elementary_floor_all))
for L, J, a, j0, upper_z, lower_z, probability in candidate_samples:
    print(
        "MESOSCOPIC_CANDIDATE_WINDOW",
        f"L={int(L)};J={J};a={a:.12f};j0={j0};lower_z={lower_z:.9f};upper_z={upper_z:.9f};prob={probability:.12f}",
    )
for L, J, a, factor, probability, scaled in strong_samples:
    print(
        "MESOSCOPIC_STRONG_CUTOFF",
        f"L={int(L)};J={J};a={a:.12f};factor={factor:.12f};prob={probability:.12f};scaled={scaled:.12e}",
    )
print("MESOSCOPIC_C_STAR", f"{c_star:.15e}")
print("MESOSCOPIC_C_SHARP", f"{c_sharp:.15e}")
print("SEMANTIC_SHA256", semantic_hash)
print("GATES", GATES)
print("RESULT", "PASS")
