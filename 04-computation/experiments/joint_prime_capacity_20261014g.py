"""Exact reciprocal-capacity bounds for first-rise Giuga/Korselt sources.

No hypothetical source is assumed to exist. A finite prime sieve supplies
directed rational upper bounds, then a proved tail bound handles every larger
prime. Four marked support bits retain all ten normal-support branches.
"""
from fractions import Fraction
from itertools import combinations
from math import isqrt, prod
import json

CHECKS = 0
SCALE = 1 << 64
MARKED = (3, 7, 11, 19)


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(label)


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer domain")
    return n


def sieve(limit):
    integer(limit, 2)
    bits = bytearray(b"\x01") * (limit + 1)
    bits[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if bits[p]:
            bits[p*p::p] = b"\x00" * ((limit - p*p)//p + 1)
    return bits


def prime_trial(n):
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def normal(ps):
    return all((p-1) % q for p in ps for q in ps if p != q)


def patterns(marked=MARKED):
    return [tuple(p for i, p in enumerate(marked) if mask >> i & 1)
            for mask in range(1 << len(marked))
            if normal(tuple(p for i, p in enumerate(marked) if mask >> i & 1))]


def eligible(p, selected, marked=MARKED):
    # Below the largest marked prime there are no other primes 3 mod 4.
    return p in selected or (p not in marked and all((p-1) % q for q in selected))


def capacity_row(ps, cutoff, selected, marked=MARKED, index_floor=2,
                 index_ceiling=4, rank_class=3):
    if selected not in patterns(marked):
        raise ValueError("not a normal marked pattern")
    bank = [p for p in ps if p <= cutoff and eligible(p, selected, marked)]
    q = next(p for p in ps if p > cutoff and eligible(p, selected, marked))
    lo = sum(SCALE//p for p in bank)
    hi = sum((SCALE+p-1)//p for p in bank)
    check(0 < hi < index_floor*SCALE, "capacity below required index")
    count = len(bank)
    excluded = count + ((index_floor*SCALE-hi)*q)//SCALE
    rho = excluded + 1 + (rank_class-(excluded+1)) % 8
    def upper_numerator(r):
        return hi*q + max(0, r-count)*SCALE
    check(upper_numerator(excluded) <= index_floor*SCALE*q,
          "all ranks through excluded have sum at most required index")
    check(upper_numerator(excluded+1) > index_floor*SCALE*q,
          "arithmetic floor is sharp for this relaxed bound")
    for r in range(excluded+1, rho):
        check(upper_numerator(r) < index_ceiling*SCALE*q and r % 8 != rank_class,
              "index pinned at its minimum, then rank congruence excludes")
    check(rho % 8 == rank_class, "first permitted lattice rank")
    return {"selected_marked_primes": list(selected), "eligible_count": count,
            "index_floor": index_floor, "next_allowed_index": index_ceiling,
            "rank_class_at_minimum_index_mod8": rank_class,
            "largest_eligible_below_cutoff": bank[-1], "first_eligible_above_cutoff": q,
            "reciprocal_lower_numerator": lo, "reciprocal_upper_numerator": hi,
            "reciprocal_denominator": SCALE, "excluded_factor_count": excluded,
            "necessary_factor_count_with_lattice": rho}


def independent_controls():
    bits = sieve(20000)
    for n in range(2, 20001):
        check(bool(bits[n]) == prime_trial(n), "sieve versus independent trial primality")
    edges = [(p, q) for p, q in combinations(MARKED, 2)
             if (q-1) % p == 0 or (p-1) % q == 0]
    check(edges == [(3, 7), (3, 19)], "normal-support graph is P3 plus isolate")
    check(len(patterns()) == 10, "complete independent-set count")
    # Literal independent pair exclusions, not the arithmetic normal() helper.
    for mask in range(16):
        chosen = tuple(p for i, p in enumerate(MARKED) if mask >> i & 1)
        independent = not ((3 in chosen and 7 in chosen) or
                           (3 in chosen and 19 in chosen))
        check((chosen in patterns()) == independent, "all support patterns")
    small = [p for p in range(3, 2000, 4) if prime_trial(p)]
    for chosen in patterns():
        bank = [p for p in small if eligible(p, chosen)]
        exact = sum((Fraction(1, p) for p in bank), Fraction())
        lo = sum(SCALE//p for p in bank)
        hi = sum((SCALE+p-1)//p for p in bank)
        check(Fraction(lo, SCALE) <= exact <= Fraction(hi, SCALE), "directed rational sums")
    # Every small actual normal support enters its exact marked branch.
    tiny = [p for p in small if p < 100]
    for size in range(1, 5):
        for support in combinations(tiny, size):
            if normal(support):
                chosen = tuple(p for p in MARKED if p in support)
                check(chosen in patterns() and all(eligible(p, chosen) for p in support),
                      "complete branch assignment for normal support")
    # Necessary hypotheses: ordinary first-rise Carmichael 8911 is not Giuga.
    support = (7, 19, 67)
    n = prod(support)
    check(n == 8911 and n % 4 == 3, "first-rise hostile source")
    check(all((n-1) % (p-1) == 0 for p in support), "hostile satisfies full Korselt")
    index = sum((Fraction(1, p) for p in support), Fraction())-Fraction(1,n)
    check(index.denominator != 1, "hostile fails integral Giuga index")
    for value in (True, 2.0, -1):
        try:
            sieve(value)
        except ValueError:
            check(True, "malformed input rejected")
        else:
            check(False, "malformed input accepted")
    return {"trial_primality_interval": [2, 20000], "marked_graph_edges": edges,
            "marked_independent_sets": 10, "hostile_source": n,
            "hostile_nonintegral_index": str(index)}


def main():
    controls = independent_controls()
    cutoff = 10_000_000
    # The explicit margin is verified to contain an eligible successor in every
    # branch; next() would fail, rather than silently assume a prime-gap bound.
    bits = sieve(cutoff+1000)
    ps = [p for p in range(3, cutoff+1001, 4) if bits[p]]
    rows = [capacity_row(ps, cutoff, chosen) for chosen in patterns()]
    for row in rows:
        q = row['first_eligible_above_cutoff']
        check(prime_trial(q), "successor independently prime")
        chosen = tuple(row['selected_marked_primes'])
        for n in range(cutoff+1, q):
            if n % 4 == 3 and eligible(n, chosen):
                check(not prime_trial(n), "no skipped eligible prime after cutoff")
    raw = min(row['excluded_factor_count'] for row in rows)
    threshold = min(row['necessary_factor_count_with_lattice'] for row in rows)
    check(raw == 11833260 and threshold == 11833267, "reproducible global capacity")
    check(3**12 > 2**19, "exact lower bound for log2(3)")
    exponent_min = 19*threshold//12+1
    marked_mersenne = (7, 23, 31, 47)
    edges = [(p, q) for p, q in combinations(marked_mersenne, 2)
             if (q-1) % p == 0 or (p-1) % q == 0]
    check(edges == [(23, 47)] and len(patterns(marked_mersenne)) == 12,
          "odd-exponent Mersenne support graph")
    ps7 = [p for p in ps if p % 8 == 7]
    mersenne_rows = [capacity_row(ps7, cutoff, chosen, marked_mersenne, 4, 8, 5)
                     for chosen in patterns(marked_mersenne)]
    for row in mersenne_rows:
        q = row['first_eligible_above_cutoff']
        check(prime_trial(q), "Mersenne-branch successor independently prime")
        chosen = tuple(row['selected_marked_primes'])
        for n in range(cutoff+1, q):
            if n % 8 == 7 and eligible(n, chosen, marked_mersenne):
                check(not prime_trial(n), "no skipped Mersenne-branch prime")
    mersenne_threshold = min(row['necessary_factor_count_with_lattice'] for row in mersenne_rows)
    # All distinct factors are at least seven; 7^5 > 2^14.
    check(7**5 > 2**14, "exact lower bound for log2(7)")
    odd_exponent_min = 14*mersenne_threshold//5+1
    if odd_exponent_min % 2 == 0:
        odd_exponent_min += 1
    return {"status": "PROVED conditional capacity; FINITE-EXACT prime census; no ROOT claim",
            "checks": CHECKS, "cutoff": cutoff, "marked_primes": list(MARKED),
            "controls": controls, "branches": rows,
            "global_excluded_factor_count": raw,
            "global_necessary_factor_count_with_lattice": threshold,
            "necessary_Mersenne_exponent_if_joint": exponent_min,
            "odd_exponent_Mersenne_capacity": {
                "marked_primes": list(marked_mersenne), "marked_graph_edges": edges,
                "branches": mersenne_rows,
                "necessary_factor_count_with_lattice": mersenne_threshold,
                "necessary_odd_exponent": odd_exponent_min,
                "hypotheses": "Full joint target M_E, E odd>=3; primes7mod8, index0mod4"},
            "scope": "Full Giuga plus Korselt, N=3 mod4; every larger prime handled by tail bound"}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
