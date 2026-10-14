"""Exact order-group capacity certificates, without constructing M_E.

Universal profile checks are finite certificates for all exponent sizes under
the declared divisor-count cap. They do not assert a joint source or ROOT.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import prod
import hashlib
import json

import mersenne_lattice_child_guards_20261014g as prior

EXPONENT = prior.EXPONENT
ODD_PRIMES = (3,5,7,11,13,17,19,23,29,31,37)
PROFILE_CAP = 1919
SERIES_CAP = 24


def integer(n, minimum=0):
    return prior.integer(n, minimum)


def _unit_log_bounds(x, terms):
    # 1 <= x <= 2, so 0 <= z <= 1/3. Exact geometric remainder.
    z = (x-1)/(x+1)
    total, power = F(), z
    for j in range(terms):
        total += 2*power/(2*j+1)
        power *= z*z
    return total, total+2*power/((2*terms+1)*(1-z*z))


def log_bounds(value, terms=16):
    """Rational lower/upper bounds for natural log of an exact number >=1."""
    integer(terms,1)
    if type(value) is int:
        integer(value,1)
        value = F(value)
    elif type(value) is not F or value < 1:
        raise ValueError("exact rational at least one required")
    k = value.numerator.bit_length()-value.denominator.bit_length()
    if value < 2**k:
        k -= 1
    x = value/(1 << k)
    lower, upper = _unit_log_bounds(x,terms)
    l2,u2 = _unit_log_bounds(F(2),terms)
    return k*l2+lower, k*u2+upper


def log2_bounds(n, terms=16):
    integer(n,2)
    integer(terms,1)
    k = n.bit_length()-1
    low, high = _unit_log_bounds(F(n,1 << k),terms)
    l2,u2 = _unit_log_bounds(F(2),terms)
    return k+low/u2, k+high/l2


def order_floor(d, minimum_prime):
    """Exact largest r with minimum_prime**r < 2**d, using no such powers."""
    integer(d,3)
    integer(minimum_prime,3)
    if minimum_prime % 2 == 0:
        raise ValueError("odd minimum required for strict logarithmic boundary")
    terms = 4
    while True:
        low, high = log2_bounds(minimum_prime,terms)
        a, b = F(d)/high, F(d)/low
        rlow = a.numerator//a.denominator
        rhigh = b.numerator//b.denominator
        if rlow == rhigh:
            return rlow, terms
        terms *= 2


def order_group(d, joint=False):
    """Three bounds for a marked odd order class. Joint is an extra premise."""
    integer(d,3)
    if d % 2 == 0 or type(joint) is not bool:
        raise ValueError("odd order and exact Boolean premise required")
    multiplier = 6 if joint and d % 4 == 1 else 2
    low_prime = multiplier*d+1
    step = (8 if joint else 2)*d
    lowlog,_ = log2_bounds(low_prime,24)
    analytic_upper = F(d,low_prime)/lowlog
    count, terms = order_floor(d,low_prime)
    if count == 0:
        spaced = F()
    else:
        ratio = F(low_prime+step*(count-1),low_prime)
        _, upper_log = log_bounds(ratio,24)
        spaced = F(1,low_prime)+upper_log/step
    return {"order":d,"minimum_prime":low_prime,"step":step,
            "exact_count_cap":count,"count_log_terms":terms,
            "strict_analytic_upper":analytic_upper,
            "integer_count_upper":F(count,low_prime),"spaced_upper":spaced}


def source_budget(exponent, factors, joint=False, divisor_cap=4096):
    """Full source factorization required; never replace the actual exponent."""
    integer(divisor_cap,1)
    if type(joint) is not bool:
        raise ValueError("exact Boolean premise required")
    old = prior.order_group_capacity(exponent,factors)
    if old["divisor_count"] > divisor_cap:
        raise ValueError("explicit divisor enumeration cap exceeded")
    divisors = [1]
    for p,a in factors:
        divisors = [d*p**j for d in divisors for j in range(a+1)]
    rows = tuple(order_group(d,joint) for d in sorted(divisors) if d > 1)
    return rows, sum((row["spaced_upper"] for row in rows),F())


@lru_cache(maxsize=None)
def _geometric(p,a,capped):
    if capped and a > SERIES_CAP:
        return F(p,p-1), F(p,(p-1)**2)
    power = p**a
    return (F(p*power-1,power*(p-1)),
            F(p*power-(a+1)*p+a,power*(p-1)**2))


@lru_cache(maxsize=None)
def _prime_logs(p):
    return log_bounds(p,16)


def profile_bound(profile, upper=True, capped=True):
    """Canonical exponent shape on the first odd primes; H envelope bound."""
    if type(profile) is not tuple or type(upper) is not bool or type(capped) is not bool:
        raise ValueError("exact tuple and Boolean flags required")
    if len(profile) > len(ODD_PRIMES):
        raise ValueError("retained prime list too short")
    for a in profile:
        integer(a,1)
    if tuple(sorted(profile,reverse=True)) != profile:
        raise ValueError("canonical nonincreasing multiplicities required")
    if capped and not upper:
        raise ValueError("infinite geometric caps are upper bounds only")
    gs, js, logs = [], [], []
    for p,a in zip(ODD_PRIMES,profile):
        g,j = _geometric(p,a,capped)
        gs.append(g)
        js.append(j)
        logs.append(_prime_logs(p)[int(upper)])
    R = prod(gs)
    return (R-1+sum((R*js[i]*logs[i]/gs[i] for i in range(len(gs))),F()))/2


def profiles(cap, maximum=None):
    integer(cap,1)
    if maximum is not None:
        integer(maximum,1)
    yield ()
    for a in range(1,min(cap-1,maximum if maximum is not None else cap-1)+1):
        for tail in profiles(cap//(a+1),a):
            yield (a,)+tail


def profile_certificate(cap=PROFILE_CAP):
    integer(cap,1)
    if cap >= 1 << (len(ODD_PRIMES)+1):
        raise ValueError("divisor cap exceeds the retained prime list")
    rows, best, winner = 0, F(), ()
    digest = hashlib.sha256()
    for aa in profiles(cap):
        bound = profile_bound(aa)
        rows += 1
        digest.update((str(aa)+":"+str(bound)+"\n").encode("ascii"))
        if bound > best:
            best,winner = bound,aa
    return {"cap":cap,"profile_count":rows,"upper":best,"maximizing_envelope_profile":winner,
            "profile_stream_sha256":digest.hexdigest(),"certifies_joint_exclusion":best <= 4}


def infinite_multiplicity_bound(distinct_primes, upper=True):
    """All exponent multiplicities unbounded, on the smallest odd bases."""
    integer(distinct_primes,1)
    if distinct_primes > len(ODD_PRIMES) or type(upper) is not bool:
        raise ValueError("retained prime count and exact Boolean required")
    ps = ODD_PRIMES[:distinct_primes]
    R = prod(F(p,p-1) for p in ps)
    return (R-1+R*sum((_prime_logs(p)[int(upper)]/(p-1) for p in ps),F()))/2


def paley_edge(p,x,y):
    integer(p,3)
    integer(x); integer(y)
    if p % 4 != 3 or not prior.prime(p) or x >= p or y >= p or x == y:
        raise ValueError("distinct vertices of a prime Paley tournament required")
    return pow((y-x)%p,(p-1)//2,p) == 1


def fraction_json(value):
    if type(value) is F:
        return {"numerator":value.numerator,"denominator":value.denominator}
    raise TypeError("unsupported JSON value")


def main():
    checks = 0
    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError("check failed: "+str(checks))
    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True)
        else:
            check(False)

    check(all(prior.prime(p) for p in ODD_PRIMES))
    check(tuple(p for p in range(3,38,2) if prior.prime(p)) == ODD_PRIMES)
    # Independent defining-power checks of the integer packing count.
    for d in range(3,130,2):
        for joint in (False,True):
            row = order_group(d,joint)
            R,L = row["exact_count_cap"],row["minimum_prime"]
            check(L**R < 1 << d and L**(R+1) > 1 << d)
            check(row["integer_count_upper"] < row["strict_analytic_upper"])
            candidates = [L+row["step"]*j for j in range(R)]
            check(sum((F(1,p) for p in candidates),F()) <= row["spaced_upper"])
    # Prime support is actually factored only in this small independent universe.
    for e in range(3,32,2):
        support = prior.factor_distinct((1 << e)-1)
        groups = {}
        for p in support:
            groups.setdefault(prior.order_two(p),[]).append(p)
        for d, ps in groups.items():
            row = order_group(d)
            S = sum((F(1,p) for p in ps),F())
            check(S <= row["integer_count_upper"] and S <= row["spaced_upper"])
            if all(p % 8 == 7 for p in ps):
                jointrow = order_group(d,True)
                check(S <= jointrow["spaced_upper"])
    # Independent divisor enumeration for the closed geometric formula.
    for length in range(4):
        for aa in product(range(1,4),repeat=length):
            profile = tuple(sorted(aa,reverse=True))
            ds = [1]
            for p,a in zip(ODD_PRIMES,profile):
                ds = [d*p**j for d in ds for j in range(a+1)]
            lower = sum(((1+log_bounds(d,20)[0])/(2*d) for d in ds if d > 1),F())
            upper = sum(((1+log_bounds(d,20)[1])/(2*d) for d in ds if d > 1),F())
            check(profile_bound(profile,False,False) <= upper)
            check(profile_bound(profile,True,False) >= lower)
    # Canonicalization checks keep actual arithmetic and the relaxed cap distinct.
    for aa in profiles(24):
        check(prod(a+1 for a in aa) <= 24 and tuple(sorted(aa,reverse=True)) == aa)
    tiny = set(profiles(12))
    independent = {tuple(sorted(aa,reverse=True)) for length in range(4)
                   for aa in product(range(1,12),repeat=length) if prod(a+1 for a in aa) <= 12}
    check(tiny == independent)
    cert = profile_certificate()
    check(cert["profile_count"] == 16894)
    check(cert["certifies_joint_exclusion"] and cert["upper"] < F(3999,1000))
    check(cert["maximizing_envelope_profile"] == (3,2,1,1,1,1,1,1,1))
    first_unpaid = (4,2,1,1,1,1,1,1,1)
    check(prod(a+1 for a in first_unpaid) == 1920)
    check(profile_bound(first_unpaid,False,False) > 4)
    pure_upper = infinite_multiplicity_bound(1)
    seven_upper = infinite_multiplicity_bound(7)
    check(pure_upper < F(2,3))
    check(seven_upper < F(399,100))
    check(infinite_multiplicity_bound(8,False) > 4)
    # The source-owned, unchanged giant label uses just three nontrivial divisors.
    factors = ((47779,1),(2086879,1))
    rows, total = source_budget(EXPONENT,factors)
    jointrows, jointtotal = source_budget(EXPONENT,factors,True)
    check(tuple(row["order"] for row in rows) == (47779,2086879,EXPONENT))
    check(total < F(1,10000))
    check(total+F(1,EXPONENT) < F(1,10000))
    check(jointtotal < total)
    for row in rows+jointrows:
        check(row["spaced_upper"] <= row["integer_count_upper"])
    paley_rows = []
    for p in (7,23,31,47,71,3,11,19):
        preserving = p % 8 == 7
        for x in range(p):
            for y in range(p):
                if x != y:
                    check((paley_edge(p,2*x%p,2*y%p) == paley_edge(p,x,y)) == preserving)
        order = prior.order_two(p)
        check((order % 2 == 1) if preserving else (order % 2 == 0))
        seen = set()
        lengths = []
        for x in range(1,p):
            if x in seen:
                continue
            y, length = x,0
            while y not in seen:
                seen.add(y); length += 1; y = 2*y % p
            check(y == x and length == order)
            lengths.append(length)
        check(len(lengths) == (p-1)//order)
        paley_rows.append({"prime":p,"doubling_preserves_orientation":preserving,
                           "order_of_two":order,"nonzero_orbits":len(lengths)})
    for bad in (True,3.0,0,-1):
        reject(lambda bad=bad: order_group(bad))
        reject(lambda bad=bad: profile_certificate(bad))
    reject(lambda: source_budget(EXPONENT,((47779,1),)))
    reject(lambda: source_budget(EXPONENT,factors,1))
    reject(lambda: source_budget(EXPONENT,factors,divisor_cap=3))
    reject(lambda: order_group(4))
    reject(lambda: order_group(3,1))
    reject(lambda: order_floor(3,8))
    reject(lambda: profile_bound((1,2)))
    reject(lambda: profile_bound((True,)))
    reject(lambda: profile_bound([2,1]))
    reject(lambda: profile_bound((2,),False,True))
    reject(lambda: log_bounds(1.0))
    reject(lambda: paley_edge(7,True,2))
    reject(lambda: paley_edge(11,1,1))
    reject(lambda: profile_certificate(1 << (len(ODD_PRIMES)+1)))
    reject(lambda: infinite_multiplicity_bound(True))
    reject(lambda: infinite_multiplicity_bound(7,1))
    result = {"status":"PROVED necessary capacity bounds; no joint source or Collatz ROOT asserted",
              "checks":checks,"universal_joint_divisor_count_minimum":1920,
              "profile_certificate":cert,"geometric_series_cap":SERIES_CAP,
              "infinite_multiplicity_certificates":{
                  "prime_power_exponent_support_upper":pure_upper,
                  "seven_distinct_exponent_primes_support_upper":seven_upper,
                  "joint_distinct_exponent_prime_minimum":8,
                  "odd_prime_power_exponents_giuga_excluded":True},
              "relaxed_unpaid_profile_at_1920":first_unpaid,
              "unchanged_exponent":EXPONENT,"exponent_factorization":factors,
              "unconditional_source_groups":rows,"unconditional_source_support_upper":total,
              "conditional_joint_source_groups":jointrows,"conditional_joint_support_upper":jointtotal,
              "paley_doubling_controls":paley_rows,
              "universes":{"integer_count_orders":[3,129,2],"literal_support_exponents":[3,31,2],
                           "unbounded_exponent_profiles_tau_cap":PROFILE_CAP},
              "scope":"Profile cap excludes every exponent size under that cap. The 1920 row only fails this relaxed upper-bound test; it is not an actual joint source."}
    print(json.dumps(result,indent=2,sort_keys=True,default=fraction_json))


if __name__ == "__main__":
    main()
