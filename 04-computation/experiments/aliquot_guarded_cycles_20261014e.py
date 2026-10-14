"""Exact aliquot cycles in fixed coprime prime-extension templates.

No production orbit discovery and no primality heuristic. The finite census
bounds cores and length, not the prime variables determined by each template.
"""
from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, prod
import json


def nat(n, minimum=1):
    if type(n) is not int:
        raise TypeError("exact integer required")
    if n < minimum:
        raise ValueError("integer below required minimum")
    return n


def prime(n):
    nat(n, 0)
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


@lru_cache(None)
def _factor(n):
    answer = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                e += 1
                n //= p
            answer.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def factor(n):
    # Exact type validation precedes the private cache (bool/float aliases).
    nat(n)
    return _factor(n)


def audit_factors(n, factors):
    nat(n)
    if type(factors) is not tuple:
        raise TypeError("tuple factorization required")
    last, value, sigma = 1, 1, 1
    for row in factors:
        if type(row) is not tuple or len(row) != 2:
            raise TypeError("prime/exponent tuple required")
        p, e = row
        nat(p, 2)
        nat(e)
        if p <= last or not prime(p):
            raise ValueError("prime factors must be strictly increasing")
        value *= p ** e
        sigma *= (p ** (e + 1) - 1) // (p - 1)
        last = p
    if value != n:
        raise ValueError("factorization belongs to a different source")
    return sigma


def sigma(n):
    return prod((p ** (e + 1) - 1) // (p - 1) for p, e in factor(n))


def aliquot(n):
    return sigma(n) - n


def direct_sigma(n):
    """Independent divisor-pair engine used only in controls."""
    nat(n)
    return sum(d + n // d if d*d != n else d
               for d in range(1, isqrt(n) + 1) if n % d == 0)


@dataclass(frozen=True)
class FactorRow:
    source: int
    factors: tuple


@dataclass(frozen=True)
class Cycle:
    rows: tuple


def audit_cycle(receipt):
    if type(receipt) is not Cycle or type(receipt.rows) is not tuple or not receipt.rows:
        raise TypeError("nonempty exact cycle receipt required")
    ns, sigmas = [], []
    for row in receipt.rows:
        if type(row) is not FactorRow:
            raise TypeError("exact factor row required")
        sigmas.append(audit_factors(row.source, row.factors))
        ns.append(row.source)
    if len(set(ns)) != len(ns):
        raise ValueError("cycle is not primitive")
    for n, s, target in zip(ns, sigmas, ns[1:] + ns[:1]):
        if s - n != target:
            raise ValueError("false aliquot edge")
    return tuple(ns)


def make_cycle(nodes):
    if type(nodes) is not tuple or not nodes:
        raise TypeError("nonempty tuple nodes required")
    result = Cycle(tuple(FactorRow(n, factor(n)) for n in nodes))
    audit_cycle(result)
    return result


def core_word(cores):
    if type(cores) is not tuple or not cores:
        raise TypeError("nonempty tuple of cores required")
    for c in cores:
        nat(c, 2)
    return cores


def core_carrier(cores):
    core_word(cores)
    p, q, carry = 1, 1, 0
    for c, d in zip(cores, cores[1:] + cores[:1]):
        s = sigma(c)
        a = s - c
        p, q, carry = a*p, d*q, a*carry + s*q
    return p, q, carry


@dataclass(frozen=True)
class CoreCycle:
    cores: tuple
    primes: tuple
    cycle: Cycle


def audit_core(receipt):
    if type(receipt) is not CoreCycle:
        raise TypeError("exact core-cycle receipt required")
    cs = core_word(receipt.cores)
    ps = receipt.primes
    if type(ps) is not tuple or len(ps) != len(cs):
        raise TypeError("one exact prime per core required")
    for c, p in zip(cs, ps):
        nat(p, 2)
        if not prime(p) or gcd(c, p) != 1:
            raise ValueError("prime or coprimality guard failed")
    ns = audit_cycle(receipt.cycle)
    if ns != tuple(c*p for c,p in zip(cs,ps)):
        raise ValueError("receipt changes the marked sources")
    p, q, carry = core_carrier(cs)
    if q <= p or (q-p)*ps[0] != carry:
        raise ValueError("affine closure failed")
    return receipt


def analyze_core(cores):
    """Total finite decision for one template; a rejected template is not
    a statement about any source outside that template.
    """
    p, q, carry = core_carrier(cores)
    if q <= p:
        return "slope", None
    if carry % (q-p):
        return "nonintegral_start", None
    current = carry // (q-p)
    ps = []
    for c, d in zip(cores, cores[1:] + cores[:1]):
        if not prime(current) or gcd(c, current) != 1:
            return "native_guard", None
        ps.append(current)
        numerator = (sigma(c)-c)*current + sigma(c)
        if numerator % d:
            return "native_guard", None
        current = numerator // d
    nodes = tuple(c*p for c,p in zip(cores,ps))
    if len(set(nodes)) != len(nodes):
        return "repeated_nodes", None
    packet = CoreCycle(cores, tuple(ps), make_cycle(nodes))
    return "cycle", audit_core(packet)


def formal_primes(cores):
    """Rational candidate before native guards; never a cycle certificate."""
    p, q, carry = core_carrier(cores)
    if q <= p:
        return None
    value = Fraction(carry, q-p)
    result = []
    for c,d in zip(cores, cores[1:]+cores[:1]):
        result.append(value)
        value = ((sigma(c)-c)*value+sigma(c))/d
    return tuple(result)


def canonical(word):
    return min(word[i:]+word[:i] for i in range(len(word)))


KNOWN = (
    (6,), (28,), (220,284),
    (1264460,1547860,1727636,1305184),
    (12496,14288,15472,14536,14264),
    (14316,19116,31704,47616,83328,177792,295488,629072,589786,
     294896,358336,418904,366556,274924,275444,243760,376736,381028,
     285778,152990,122410,97946,48976,45946,22976,22744,19916,17716),
)


def main():
    checks = 0
    def check(value, label):
        nonlocal checks
        checks += 1
        if not value:
            raise RuntimeError(label)
    def rejects(f, label):
        try:
            f()
        except (TypeError, ValueError):
            check(True, label)
        else:
            check(False, label)

    for n in range(1, 4097):
        s = sigma(n)
        check(s == direct_sigma(n), "independent divisor engine")
        exceptional = isqrt(n)**2 == n or (n%2 == 0 and isqrt(n//2)**2 == n//2)
        check(bool(s%2) == exceptional, "square/twice-square parity")
        check(audit_factors(n, factor(n)) == s, "source factor receipt")

    # Fixed-cofactor, fixed-prime-exponent monotonicity controls.
    primes = tuple(p for p in range(2,44) if prime(p))
    for c in range(1,17):
        for e in range(1,5):
            values = []
            for p in primes:
                if gcd(c,p) == 1:
                    value = (sigma(c)-c)*p**e + sigma(c)*sum(p**j for j in range(e))
                    check(value == aliquot(c*p**e), "prime-power polynomial")
                    values.append(value)
            if c == e == 1:
                check(set(values) == {1}, "prime constant exception")
            else:
                check(all(x<y for x,y in zip(values,values[1:])), "strict template order")

    known_rows = []
    for nodes in KNOWN:
        receipt = make_cycle(nodes)
        check(audit_cycle(receipt) == nodes, "known primitive cycle")
        check(all(direct_sigma(n)-n == target for n,target in zip(nodes,nodes[1:]+nodes[:1])), "independent known edges")
        ratios = tuple(Fraction(aliquot(n),n) for n in nodes)
        abundancies = tuple(x+1 for x in ratios)
        check(prod(ratios) == 1, "cycle ratio product")
        check(sum(sigma(n) for n in nodes) == 2*sum(nodes), "weighted abundance balance")
        check(prod(abundancies) >= 2**len(nodes), "abundance geometric bound")
        if len(nodes)>1:
            check(aliquot(min(nodes))>min(nodes) and aliquot(max(nodes))<max(nodes), "cycle extremes")
            check(sum(abundancies)>2*len(nodes) and prod(abundancies)>2**len(nodes), "strict nontrivial abundance")
            g = gcd(*nodes)
            check(sigma(g)<2*g, "common divisor must be deficient")
            if len(nodes)%2:
                d = gcd(*(sigma(n) for n in nodes))
                check(g%(d//2 if d%4==2 else d)==0, "cited odd-cycle gcd guard")
        switches = [n for n in nodes if n%2 != aliquot(n)%2]
        check(len(switches)%2 == 0, "cycle parity switches")
        known_rows.append(dict(nodes=nodes,factors=[row.factors for row in receipt.rows],
                               abundance_signs=[(aliquot(n)>n)-(aliquot(n)<n) for n in nodes],
                               parity_switches=switches))

    counts, found, cycle_set = Counter(), [], set()
    for length in range(1,4):
        for cores in product(range(2,65), repeat=length):
            if canonical(cores) != cores:
                continue
            counts[(length,"templates")]+=1
            p,q,b = core_carrier(cores)
            sigs = tuple(direct_sigma(c) for c in cores)
            aa = tuple(s-c for s,c in zip(sigs,cores))
            dd = cores[1:]+cores[:1]
            explicit_b = sum(sigs[i]*prod(aa[i+1:])*prod(dd[:i]) for i in range(length))
            check((p,q,b)==(prod(aa),prod(cores),explicit_b), "independent carrier expansion")
            status, packet = analyze_core(cores)
            counts[(length,status)]+=1
            if packet is not None:
                ns = audit_cycle(packet.cycle)
                check(all(direct_sigma(n)-n == m for n,m in zip(ns,ns[1:]+ns[:1])), "native compiled edges")
                check(formal_primes(cores)==tuple(Fraction(v) for v in packet.primes), "rational candidate agrees")
                cycle_set.add(canonical(ns))
                found.append(dict(cores=cores,primes=packet.primes,nodes=ns,carrier=(p,q,b)))
    check(sum(v for (l,s),v in counts.items() if s=="templates")==85470, "exact core universe")
    check(counts[(3,"cycle")]==0, "restricted no3 result")
    check(cycle_set=={(6,),(28,),(496,),(8128,),(220,284),(5020,5564)}, "cycle support")

    # A longer actual cycle also fits the compiler when its cores are retained.
    five_cores = (176,304,16,184,8)
    status, five = analyze_core(five_cores)
    check(status=="cycle" and audit_cycle(five.cycle)==KNOWN[4], "five-cycle core reconstruction")
    check(analyze_core((20,4))[1].primes==(11,71), "amicable marked orientation")
    # Formal affine closure is not a native aliquot cycle: p=1 is not prime.
    ghost = (2,2,16)
    check(formal_primes(ghost)==(Fraction(23),Fraction(13),Fraction(1)), "formal closure hostile")
    check(analyze_core(ghost)==("native_guard",None), "ghost rejected")
    check(aliquot(46)==26 and aliquot(26)==16 and aliquot(16)==15, "first failed ghost edge")
    rejects(lambda: make_cycle((46,26,16)), "factor authenticated ghost rejection")
    # Same initial sign does not select a terminal or a cycle length.
    terminating=(12,16,15,9,4,3,1,0)
    check(all(aliquot(n)==m for n,m in zip(terminating,terminating[1:])), "abundant terminating source")
    check(aliquot(12)>12 and aliquot(220)>220, "abundant sign hostile")
    check(4*5*251==5020 and 4*13*107==5564, "fixed-shape amicable pair")
    check(5<13 and 251>107, "prime coordinates must cross")

    good = analyze_core((20,4))[1]
    for bad in (replace(good, cores=(20.0,4)), replace(good, primes=(11,True)),
                replace(good, primes=(11.0,71)), replace(good, cores=[20,4])):
        rejects(lambda bad=bad: audit_core(bad), "typed core receipt")
    for bad in (True,1.0,0,-1):
        rejects(lambda bad=bad: factor(bad), "factor cache type guard")
    for cs in ((),[2,4],(1,4),(True,4)):
        rejects(lambda cs=cs: analyze_core(cs), "invalid core template")
    rejects(lambda: audit_factors(6,((2,1),(3.0,1))), "factor prime type")
    rejects(lambda: audit_factors(6,((2,1),(5,1))), "factor source binding")
    rejects(lambda: make_cycle((6,6)), "primitive cycle requirement")

    result=dict(status="PROVED fixed-template iff; FINITE-EXACT census; general period3 existence OPEN",
                checks=checks,core_universe=dict(core_min=2,core_max=64,lengths=[1,2,3],
                    rotation_quotient=True,templates=85470,prime_variable_bound=None),
                rejection_counts={str(l):{s:v for (ll,s),v in sorted(counts.items()) if ll==l}
                                  for l in range(1,4)},compiled=found,distinct_cycles=sorted(cycle_set),
                known_cycle_controls=known_rows,five_cycle_cores=five_cores,
                formal_ghost=dict(cores=ghost,formal_primes=[23,13,1],formal_nodes=[46,26,16],
                                 failed_actual_edge=[16,15]),
                general_cycle_catalogue_complete=False)
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
