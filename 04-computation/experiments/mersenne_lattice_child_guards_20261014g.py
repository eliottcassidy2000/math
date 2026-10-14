"""Source-owned Mersenne arithmetic filters and a local/global descent hostile.

No function materializes the current Mersenne source or discovers a ROOT word.
The modular prefix reader is formal until the separate growth guard is checked.
"""
from dataclasses import dataclass, replace
from functools import lru_cache
from fractions import Fraction
from math import gcd, isqrt, lcm
import hashlib
import json

import mersenne_joint_guards_20261014c as native

EXPONENT = 99_708_991_741
RAY_WORD = (1, 1, 1, 3, 1, 1, 1, 5)
CHECKPOINT_DEPTH = 18_843
CHECKPOINT_COST = 37_416
FULL_COST = CHECKPOINT_COST + sum(RAY_WORD)
PHASE_BITS = FULL_COST - 2


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer outside domain")
    return n


def prime(n):
    integer(n)
    return n >= 2 and (n == 2 or n % 2 != 0 and
                       all(n % d for d in range(3, isqrt(n)+1, 2)))


def factor_distinct(n):
    integer(n, 1)
    answer = []
    for d in range(2, isqrt(n)+1):
        if n % d == 0:
            answer.append(d)
            while n % d == 0:
                n //= d
    if n > 1:
        answer.append(n)
    return tuple(answer)


def order_two(p):
    integer(p, 3)
    if not prime(p):
        raise ValueError("odd prime required")
    result = p-1
    for q in factor_distinct(p-1):
        while result % q == 0 and pow(2, result//q, p) == 1:
            result //= q
    return result


def exponent_fermat_remainder(exponent):
    """Nonzero excludes Carmichael; zero is only a necessary condition."""
    integer(exponent, 2)
    return (pow(2, exponent, exponent)-2) % exponent


def order_group_capacity(exponent, factors):
    """Strict reciprocal-support bound from a full unchanged exponent factorization.

    Its proof groups prime divisors of M_E by their exact order of two.
    It does not factor M_E or supply a Collatz receipt.
    """
    integer(exponent, 3)
    if exponent % 2 == 0 or type(factors) is not tuple or not factors:
        raise ValueError("odd exponent and exact full factorization required")
    value, divisors, last, multiplicity = 1, 1, 1, 0
    for row in factors:
        if type(row) is not tuple or len(row) != 2:
            raise ValueError("exact prime-power row required")
        p, a = row
        integer(p, 3)
        integer(a, 1)
        if p <= last or not prime(p):
            raise ValueError("canonical odd prime factors required")
        value *= p**a
        divisors *= a+1
        multiplicity += a
        last = p
    if value != exponent:
        raise ValueError("factorization changes or omits the source exponent")
    upper = Fraction(divisors-1,5)
    return {"divisor_count":divisors,"strict_support_upper":upper,
            "giuga_excluded":upper <= 1,
            "joint_excluded":upper <= 4,
            "both_signed_integral_indices_excluded":multiplicity >= 2 and upper+Fraction(1,7) <= 1}


@dataclass(frozen=True)
class ExclusionPhase:
    prime_divisor: int
    order: int
    residue: int
    modulus: int


def audit_exclusion(phase):
    if type(phase) is not ExclusionPhase:
        raise ValueError("exact exclusion phase required")
    for n in (phase.prime_divisor, phase.order, phase.modulus):
        integer(n, 1)
    integer(phase.residue)
    q, order = phase.prime_divisor, phase.order
    if (not prime(q) or q == 2 or order != order_two(q)
            or phase.modulus != lcm(q, order)
            or phase.residue >= phase.modulus or phase.residue % q
            or phase.residue % order == 1 % order):
        raise ValueError("invalid source-owned exclusion phase")
    return phase


def exclusion_phase(exponent, divisor):
    integer(exponent, 2)
    integer(divisor, 3)
    if not prime(divisor) or exponent % divisor:
        raise ValueError("supplied prime must divide the unchanged exponent")
    order = order_two(divisor)
    modulus = lcm(divisor, order)
    return audit_exclusion(ExclusionPhase(divisor, order, exponent % modulus, modulus))


def excludes_exponent(phase, exponent):
    audit_exclusion(phase)
    integer(exponent, 2)
    return exponent % phase.modulus == phase.residue


def joint_mod8_residue_allowed(index, factor_count):
    """Only the two proved residue constraints; not a Giuga/Korselt verifier."""
    integer(index, 1)
    integer(factor_count, 1)
    return (factor_count % 2 == 1
            and (factor_count-index-1) % 8 == 0
            and (factor_count+index-1) % 8 == 0)


def mersenne_prefix(exponent, steps, precision):
    """Exact valuations after the initial ones run, with explicit finite bits.

    This reader alone does not assert that an earlier actual state was not ROOT.
    """
    integer(exponent, 2)
    integer(steps)
    integer(precision, 2)
    bits = precision
    x = (2*pow(3, exponent-1, 1 << bits)-1) % (1 << bits)
    word = []
    for _ in range(steps):
        z = (3*x+1) % (1 << bits)
        if not z:
            raise ValueError("valuation unresolved at supplied precision")
        a = (z & -z).bit_length()-1
        if a >= bits:
            raise ValueError("odd endpoint unresolved at supplied precision")
        bits -= a
        x = z >> a
        word.append(a)
    return tuple(word), x, bits


def word_digest(word):
    native.word(word)
    return hashlib.sha256(",".join(map(str, word)).encode("ascii")).hexdigest()


@lru_cache(maxsize=1)
def _checkpoint_data():
    # A frozen finite reader, not a search for the next occurrence.
    full, _, bits = mersenne_prefix(EXPONENT, CHECKPOINT_DEPTH+len(RAY_WORD), FULL_COST+1)
    prefix = full[:-len(RAY_WORD)]
    if (full[-len(RAY_WORD):] != RAY_WORD or sum(prefix) != CHECKPOINT_COST
            or sum(full) != FULL_COST or bits != 1):
        raise ValueError("frozen checkpoint metadata mismatch")
    return prefix, full, native.carrier(prefix), native.carrier(full)


def checkpoint_guard(exponent):
    integer(exponent, 2)
    return exponent % (1 << PHASE_BITS) == EXPONENT


def checkpoint_residue(exponent, modulus, after_descent=False):
    """Exact actual checkpoint residue after a source-authenticated growth cut."""
    integer(exponent, 2)
    integer(modulus, 1)
    if type(after_descent) is not bool:
        raise ValueError("exact Boolean stage required")
    if not checkpoint_guard(exponent) or exponent-1 < 2*FULL_COST:
        raise ValueError("native source or first-hit growth guard failed")
    data = _checkpoint_data()
    p, q, b = data[3 if after_descent else 2]
    x = (2*pow(3, exponent-1, q*modulus)-1) % (q*modulus)
    numerator = p*x+b
    if numerator % q:
        raise ValueError("full-denominator division failed")
    return (numerator//q) % modulus


def j2_partner_candidates():
    """Complete gaps for this one fixed head and relation X=81Y+10."""
    p, q, b = native.carrier(RAY_WORD)
    rows = []
    for gap in (-1, 0, 1):
        cost = sum(RAY_WORD)-2*gap
        target = (b+10*p+(1 << cost)*(4**gap-1)//3 if gap >= 0
                  else b+10*p-q*(4**(-gap)-1)//3)
        rows.append((gap, cost, target, native.decode(len(RAY_WORD)+4, cost, target)))
    return tuple(rows)


def literal(source, word):
    integer(source, 1)
    if source % 2 == 0:
        raise ValueError("positive odd source required")
    native.word(word)
    states = [source]
    for a in word:
        if source == 1:
            raise ValueError("ROOT padding")
        z = 3*source+1
        actual = (z & -z).bit_length()-1
        if actual != a:
            raise ValueError("wrong actual valuation")
        source = z >> a
        states.append(source)
    return tuple(states)


def main():
    checks = 0
    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError("check failed: " + str(checks))
    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True)
        else:
            check(False)

    # Exact multiplicative order of 2 modulo M_E, and the necessary test.
    for e in range(2, 21):
        n = (1 << e)-1
        check(pow(2, e, n) == 1)
        check(all(pow(2, a, n) != 1 for a in range(1, e)))
        check((pow(2, n-1, n) == 1) == (exponent_fermat_remainder(e) == 0))
    check(exponent_fermat_remainder(11) == 0)
    check(2047 == 23*89 and 2046 % 88 != 0)  # Passing is not Korselt.
    check(prime(47779) and prime(2086879) and 47779*2086879 == EXPONENT)
    check(pow(2, EXPONENT, EXPONENT) == 13454057925)
    check(exponent_fermat_remainder(EXPONENT) == 13454057923)
    capacity = order_group_capacity(EXPONENT,((47779,1),(2086879,1)))
    check(capacity["divisor_count"] == 4 and capacity["strict_support_upper"] == Fraction(3,5))
    check(capacity["giuga_excluded"] and capacity["both_signed_integral_indices_excluded"])
    check(capacity["strict_support_upper"]+Fraction(1,7) == Fraction(26,35))
    for e in range(3,32,2):
        factors = []
        for p in factor_distinct(e):
            a, z = 0, e
            while z % p == 0:
                z //= p
                a += 1
            factors.append((p,a))
        bound = order_group_capacity(e,tuple(factors))
        n = (1 << e)-1
        support = factor_distinct(n)
        groups = {}
        for p in support:
            d = order_two(p)
            check(d > 1 and e % d == 0 and d % 2 == 1 and p >= 2*d+1)
            groups.setdefault(d,[]).append(p)
        for d, ps in groups.items():
            product = 1
            for p in ps:
                product *= p
            check(((1 << d)-1) % product == 0)
            check(sum((Fraction(1,p) for p in ps),Fraction()) < Fraction(1,5))
        total = sum((Fraction(1,p) for p in support),Fraction())
        check(total < bound["strict_support_upper"])
        if bound["both_signed_integral_indices_excluded"]:
            check(0 < total-Fraction(1,n) < 1 and 0 < total+Fraction(1,n) < 1)
    phase = exclusion_phase(EXPONENT, 47779)
    check(phase == ExclusionPhase(47779,15926,27377367,760928354))
    for t in range(32):
        e = phase.residue+t*phase.modulus
        check(excludes_exponent(phase,e))
        check(exponent_fermat_remainder(e) != 0)
    # Independent prime-order check of the conditional factor restriction.
    for p in range(3, 2000, 2):
        if prime(p):
            order = order_two(p)
            check(pow(2, order, p) == 1)
            check(all(pow(2, d, p) != 1 for d in range(1, order) if order % d == 0))
            if p % 4 == 3 and order % 2:
                check(p % 8 == 7)
    residue_rows = []
    for k in range(1,33):
        for r in range(1,33):
            actual = joint_mod8_residue_allowed(k,r)
            check(actual == (k % 4 == 0 and r % 8 == (k+1) % 8))
            if actual and k <= 8 and r <= 8:
                residue_rows.append([k,r])

    check(native.carrier(RAY_WORD) == (6561,16384,9425))
    check(native.native_cell(RAY_WORD) == (8911,32768))
    for j in range(1,len(RAY_WORD)):
        p,q,b = native.carrier(RAY_WORD[:j])
        check(p > q and b > 0)
    for t in range(64):
        n = 8911+32768*t
        states = literal(n,RAY_WORD)
        check(states[-1] == 3569+13122*t < n)
        check(all(x > n for x in states[1:-1]))
    for e in range(1,65):
        check((pow(2,e,32768)-1) % 32768 != 8911)
    partners = j2_partner_candidates()
    check(all(row[-1] is None for row in partners))
    check(all(row[2] < 3**11 for row in partners))
    check(14-2*2 < 12 and 9425+10*6561-16384*(4**2-1)//3 < 1)

    # Independent literal-reader controls on manageable sources, no ROOT padding.
    for e in range(3,25):
        word, _, _ = mersenne_prefix(e,4,40)
        x = 2*3**(e-1)-1
        y = x
        for a in word:
            z = 3*y+1
            check((z & -z).bit_length()-1 == a)
            y = z >> a
        p,q,b = native.carrier(word)
        check((p*x+b)//q == y and (p*x+b) % q == 0)

    prefix, full, carrier0, carrier1 = _checkpoint_data()
    p,q,b = carrier1
    start_mod = (2*pow(3,EXPONENT-1,2*q)-1) % (2*q)
    check((p*start_mod+b) % (2*q) == q)  # Independent full native cylinder.
    check(EXPONENT-1 >= 2*FULL_COST and EXPONENT < 1 << PHASE_BITS)
    check(checkpoint_guard(EXPONENT) and checkpoint_guard(EXPONENT+(1 << PHASE_BITS)))
    check(not checkpoint_guard(EXPONENT+2))
    combined_step = lcm(phase.modulus,1 << PHASE_BITS)
    check(combined_step == 380464177*(1 << PHASE_BITS))
    check(checkpoint_guard(EXPONENT+combined_step)
          and excludes_exponent(phase,EXPONENT+combined_step))
    check(checkpoint_residue(EXPONENT,32768) == 8911)
    for modulus in (1,2,3,5,7,16,81,32768):
        # Reduce at the full extra denominator, then use the eight-step ray map.
        z = checkpoint_residue(EXPONENT,16384*modulus)
        expected = (6561*z+9425)//16384 % modulus
        check((6561*z+9425) % 16384 == 0)
        check(checkpoint_residue(EXPONENT,modulus,True) == expected)
    for bits in range(3,31):
        modulus = 1 << bits
        order = 1 << (bits-2)
        check(pow(3,order,modulus) == 1 and pow(3,order//2,modulus) != 1)
    for bad in (True, 3.0, 0, -1):
        reject(lambda bad=bad: exponent_fermat_remainder(bad))
        reject(lambda bad=bad: checkpoint_guard(bad))
    for bad in (replace(phase,prime_divisor=47779.0), replace(phase,order=True),
                replace(phase,residue=phase.residue+1), replace(phase,modulus=phase.modulus+2)):
        reject(lambda bad=bad: audit_exclusion(bad))
    reject(lambda: exclusion_phase(11,3))
    reject(lambda: exclusion_phase(33,3))  # E=1 mod ord_3(2), no exclusion.
    reject(lambda: checkpoint_residue(EXPONENT+2,7))
    reject(lambda: checkpoint_residue(EXPONENT,7,1))
    reject(lambda: mersenne_prefix(EXPONENT,3,2))
    reject(lambda: literal(1,(2,)))
    for e, factors in ((EXPONENT,((47779,1),)), (EXPONENT,((47779.0,1),(2086879,1))),
                       (EXPONENT,((47779,True),(2086879,1))), (9,[(3,2)]),
                       (9,((9,1),)), (15,((5,1),(3,1))), (8,((2,3),))):
        reject(lambda e=e,factors=factors: order_group_capacity(e,factors))

    result = {
        "status": "PROVED scoped filters; FINITE-EXACT readers; giant ROOT remains OPEN",
        "checks": checks,
        "actual_minimum_exponent": EXPONENT,
        "fermat_exclusion": {"pow_2_E_mod_E":13454057925,"necessary_remainder":13454057923,
                              "exponent_prime_factors":[47779,2086879],"phase":phase.__dict__},
        "conditional_joint_mersenne": {"prime_residue_mod_8":7,"index_divisible_by":4,
                                        "admissible_small_index_factor_count_residues":residue_rows,
                                        "positive_joint_sources_asserted":0},
        "order_group_capacity": {"full_exponent_factorization":[[47779,1],[2086879,1]],
            "divisor_count":4,"strict_reciprocal_support_upper":"3/5",
            "strict_signed_index_upper":"26/35","actual_source_not_giuga":True,
            "actual_both_signed_indices_nonintegral":True,
            "universal_odd_exponent_giuga_necessary_divisor_count":7,
            "universal_odd_exponent_joint_necessary_divisor_count":22,
            "literal_support_control_exponents":[3,31,2]},
        "ray": {"word":RAY_WORD,"source_base":8911,"source_step":32768,
                "child_base":3569,"child_step":13122,"direct_mersenne_intersection":"empty",
                "j2_complete_signed_partner_candidates":partners},
        "current_source_checkpoint": {"post_ascent_depth":CHECKPOINT_DEPTH,
            "prefix_cost":CHECKPOINT_COST,"prefix_sha256_ascii_comma":word_digest(prefix),
            "next_word":RAY_WORD,"full_cost":FULL_COST,"exponent_phase_bits":PHASE_BITS,
            "initial_precision_bits":FULL_COST+1,"full_remaining_precision_bits":1,
            "local_descent":True,"paid_against_original_mersenne":False,
            "all_height_scope":"E=99708991741 modulo 2^37428, positive exponents",
            "source_materialized":False,"root_certificate_asserted":False},
        "universes": {"literal_order_exponents":[2,20],"prime_order_bound_exclusive":2000,
                      "residue_index_and_count":[1,32],"ray_parameters":[0,63],
                      "small_reader_exponents":[3,24]},
    }
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
