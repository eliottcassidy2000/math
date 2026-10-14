"""Exact signed cofactor defects and complete one-prime completions.

Run: python -B -X utf8 04-computation/experiments/giuga_korselt_defects_20261014f.py
No orbit discovery, probabilistic primality, floating point, or assert gates.
"""
from dataclasses import dataclass, asdict, replace
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import gcd, lcm, prod
import json


def integer(x, minimum=None):
    if type(x) is not int or (minimum is not None and x < minimum):
        raise ValueError("exact integer outside domain")
    return x


def sheet(epsilon):
    integer(epsilon)
    if epsilon not in (-1, 1):
        raise ValueError("sign must be -1 or +1")
    return epsilon


@lru_cache(None)
def _factor(n):
    answer = []
    p = 2
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            answer.append((p, e))
        p = 3 if p == 2 else p+2
    if n > 1:
        answer.append((n, 1))
    return tuple(answer)


def factor(n):
    return _factor(integer(n, 1))


def prime(n):
    integer(n)
    return n >= 2 and factor(n) == ((n, 1),)


def derivative(n):
    integer(n, 1)
    return sum(e*(n//p) for p, e in factor(n))


@dataclass(frozen=True)
class Core:
    source: int
    primes: tuple
    cofactor_sum: int
    period: int


def core(n):
    integer(n, 2)
    fs = factor(n)
    if any(e != 1 for _, e in fs):
        raise ValueError("core must be squarefree")
    ps = tuple(p for p, _ in fs)
    return Core(n, ps, sum(n//p for p in ps), lcm(*(p-1 for p in ps)))


def audit_core(c):
    if type(c) is not Core or type(c.primes) is not tuple:
        raise ValueError("exact Core required")
    integer(c.source, 2)
    integer(c.cofactor_sum, 1)
    integer(c.period, 1)
    for p in c.primes:
        integer(p, 2)
    if c != core(c.source):
        raise ValueError("source/factor/carry mismatch")
    return c


@dataclass(frozen=True)
class Signed:
    core: Core
    sign: int
    index: int


def signed(c, epsilon):
    audit_core(c)
    sheet(epsilon)
    q, r = divmod(c.cofactor_sum-epsilon, c.source)
    return None if r else Signed(c, epsilon, q)


def audit_signed(s):
    if type(s) is not Signed:
        raise ValueError("exact Signed required")
    audit_core(s.core)
    sheet(s.sign)
    integer(s.index, 0)
    if signed(s.core, s.sign) != s:
        raise ValueError("wrong signed defect/index")
    return s


def korselt(c):
    audit_core(c)
    return (c.source-1) % c.period == 0


def classical_giuga_sum_criterion(n):
    """The classical necessary-and-sufficient *congruence* predicate, not primality."""
    integer(n, 2)
    if any(e != 1 for _, e in factor(n)):
        return False
    c = core(n)
    return signed(c, 1) is not None and korselt(c)


def defect_vector(c, epsilon):
    audit_core(c)
    sheet(epsilon)
    return tuple((p, (c.source//p-epsilon) % p,
                  (c.source-1) % (p-1)) for p in c.primes)


@dataclass(frozen=True)
class OldJointGuard:
    source: Core
    residue: int
    modulus: int


def old_joint_guard(c):
    """Only old-factor obligations; this is not a completed Giuga counterexample."""
    audit_core(c)
    n, lam = c.source, c.period
    if gcd(n, lam) != 1:
        return None
    a = pow(c.cofactor_sum, -1, n)
    b = pow(n, -1, lam) if lam > 1 else 0
    t = ((b-a)*pow(n, -1, lam)) % lam if lam > 1 else 0
    return OldJointGuard(c, a+n*t, n*lam)


def audit_old_guard(g):
    if type(g) is not OldJointGuard:
        raise ValueError("exact OldJointGuard required")
    audit_core(g.source)
    integer(g.residue, 0)
    integer(g.modulus, 1)
    if old_joint_guard(g.source) != g:
        raise ValueError("old-prime guard/source mismatch")
    return g


def old_joint_accepts(g, block_product):
    audit_old_guard(g)
    integer(block_product, 1)
    return (block_product-g.residue) % g.modulus == 0


@dataclass(frozen=True)
class Completion:
    source: Core
    sign: int
    prime: int
    residue: int
    target: Signed


def complete_one(c, epsilon):
    """Complete over ALL new primes, using at most one candidate, no cutoff."""
    audit_core(c)
    sheet(epsilon)
    n = c.source
    residue = (epsilon*pow(c.cofactor_sum, -1, n)) % n
    q = n+1 if epsilon == -1 and residue == 1 else residue
    if not prime(q) or n % q == 0 or (n-epsilon) % q:
        return None
    target = signed(core(n*q), epsilon)
    if target is None:
        raise RuntimeError("derived completion lost its local guard")
    return Completion(c, epsilon, q, residue, target)


def audit_completion(r):
    if type(r) is not Completion:
        raise ValueError("exact Completion required")
    audit_core(r.source)
    sheet(r.sign)
    integer(r.prime, 2)
    integer(r.residue, 1)
    audit_signed(r.target)
    if complete_one(r.source, r.sign) != r:
        raise ValueError("completion belongs to a different source/sign")
    return r


def signed_step(s, epsilon):
    """Independent closed classification for already-signed sources."""
    audit_signed(s)
    sheet(epsilon)
    n, alpha = s.core.source, s.sign
    q = (n+1 if epsilon == -1 else n-1) if alpha == -1 else (
        5-n if n in (2, 3) and epsilon == -1 else None)
    return q if q is not None and prime(q) and n % q else None


@dataclass(frozen=True)
class Block:
    source: Signed
    added_primes: tuple
    sign: int
    target: Signed


def complete_block(s, new_primes, epsilon):
    audit_signed(s)
    sheet(epsilon)
    if type(new_primes) is not tuple or not new_primes:
        raise ValueError("nonempty ordered tuple of new primes required")
    for p in new_primes:
        integer(p, 2)
        if not prime(p):
            raise ValueError("nonprime block factor")
    if tuple(sorted(set(new_primes))) != new_primes:
        raise ValueError("canonical distinct block factors required")
    n, b = s.core.source, prod(new_primes)
    if gcd(n, b) != 1:
        raise ValueError("block repeats source factor")
    if (b-s.sign*epsilon) % n or any((n*(b//q)-epsilon) % q for q in new_primes):
        return None
    t = signed(core(n*b), epsilon)
    if t is None:
        raise RuntimeError("block CRT implication failed")
    return Block(s, new_primes, epsilon, t)


def signed_parents(s):
    """All parents obtained by removing one prime and retaining either sign."""
    audit_signed(s)
    answer = []
    for q in s.core.primes:
        n = s.core.source//q
        if n <= 1:
            continue
        c = core(n)
        for alpha in (-1, 1):
            parent = signed(c, alpha)
            if parent is not None:
                r = complete_one(c, s.sign)
                if r is not None and r.prime == q:
                    answer.append((parent, r))
    return tuple(answer)


def primary_fork(s):
    """Two authenticated completions, when both neighbouring primes exist."""
    audit_signed(s)
    if s.sign != -1 or s.index != 1:
        raise ValueError("index-one minus source required")
    lower, upper = complete_one(s.core, 1), complete_one(s.core, -1)
    return None if lower is None or upper is None else (lower, upper)


def power_sum_residue(c):
    """Exact CRT evaluation of sum_{a=1}^{n-1} a^(n-1)+1 mod n."""
    audit_core(c)
    n = c.source
    total = 0
    for p in c.primes:
        m = n//p
        value = (1-m) % p if (n-1) % (p-1) == 0 else 1
        total += value*m*pow(m, -1, p)
    return total % n


CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def rejects(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        return True
    return False


def main():
    minus, plus, carmichael = [], [], []
    for n in range(2, 10001):
        fs = factor(n)
        sf = all(e == 1 for _, e in fs)
        for epsilon in (-1, 1):
            local = all((n//p-epsilon) % p == 0 for p, _ in fs)
            check(((derivative(n)-epsilon) % n == 0) == local,
                  "multiplicity-aware derivative iff local")
            if not sf:
                check(not local, "repeated prime obstruction")
        if not sf:
            continue
        c = core(n)
        check(gcd(c.cofactor_sum, n) == 1, "squarefree derivative is a unit")
        if korselt(c) and len(c.primes) > 1:
            carmichael.append(n)
            check(gcd(n, c.period) == 1, "Carmichael normal support")
        for epsilon in (-1, 1):
            s = signed(c, epsilon)
            if s:
                (minus if epsilon == -1 else plus).append(n)
                check(Fraction(c.cofactor_sum-epsilon, n) == s.index,
                      "exact signed index")
                if epsilon == 1 and len(c.primes) > 1:
                    check(s.index >= 1, "Giuga index positive, not assumed one")
            r = complete_one(c, epsilon)
            # Independent finite list: every possible q divides n-epsilon.
            candidates = [q for q, _ in factor(n-epsilon)
                          if n % q and all((n*q//p-epsilon) % p == 0
                                           for p in c.primes+(q,))]
            check(candidates == ([] if r is None else [r.prime]),
                  "complete one-prime classifier versus divisor path")
            if r is not None:
                check(audit_completion(r) == r, "source-bound receipt")
            if s is not None:
                for target_sign in (-1, 1):
                    rr = complete_one(c, target_sign)
                    q = signed_step(s, target_sign)
                    check(q == (None if rr is None else rr.prime), "four-case classification")
                    if rr is not None:
                        delta = Fraction(s.sign*q+n-target_sign, n*q)
                        check(rr.target.index == s.index+delta, "exact index transition")
                        check(delta == (1 if n in (2, 3) and (s.sign, target_sign) == (1, -1)
                                        else 0), "only exceptional index jump")
    check(minus == [2, 6, 42, 1806], "declared minus census")
    giuga = [n for n in plus if not prime(n)]
    check(giuga == [30, 858, 1722], "declared Giuga census")

    literal_count = 0
    for n in range(2, 401):
        actual = (1+sum(pow(a, n-1, n) for a in range(1, n))) % n
        check((actual == 0) == classical_giuga_sum_criterion(n), "literal full power sum")
        if all(e == 1 for _, e in factor(n)):
            check(power_sum_residue(core(n)) == actual, "full CRT power-sum readout")
        literal_count += 1
    for n in range(2, 201):
        actual = all(pow(a, n, n) == a for a in range(n))
        expected = all(e == 1 for _, e in factor(n)) and korselt(core(n))
        check(actual == expected, "literal full-ring Korselt function")

    old_guard_controls = 0
    for n in range(2, 301):
        if any(e != 1 for _, e in factor(n)):
            continue
        c = core(n)
        g = old_joint_guard(c)
        for b in range(1, 101):
            direct = all((b*(n//p)-1) % p == 0 and (n*b-1) % (p-1) == 0
                         for p in c.primes)
            predicted = g is not None and old_joint_accepts(g, b)
            check(direct == predicted, "partial old-prime CRT guard iff")
            old_guard_controls += 1
        if g is not None:
            for t in (0, 1, 17):
                b = g.residue+g.modulus*t
                check(all((b*(n//p)-1) % p == 0 and (n*b-1) % (p-1) == 0
                          for p in c.primes), "positive old-prime progression")

    # All factors retained: Leibniz is valid even on overlapping prime supports.
    for n, m in product(range(1, 101), repeat=2):
        check(derivative(n*m) == m*derivative(n)+n*derivative(m), "Leibniz carrier")
        if gcd(n, m) == 1 and n > 1 and m > 1:
            if all(e == 1 for _, e in factor(n*m)):
                c, d = core(n), core(m)
                for alpha, beta, epsilon in product((-1, 1), repeat=3):
                    left = derivative(n*m)-epsilon
                    right = m*(c.cofactor_sum-alpha)+n*(d.cofactor_sum-beta)+m*alpha+n*beta-epsilon
                    check(left == right, "signed block defect identity")

    seeds = [2, 3, 5, 6, 7, 30, 42, 858, 1722, 1806, 47058, 66198]
    signed_seeds = [s for n in seeds for epsilon in (-1, 1)
                    if (s := signed(core(n), epsilon)) is not None]
    block_controls = 0
    for s in signed_seeds:
        for q, r in [(5, 7), (11, 13), (17, 19), (23, 31), (41, 43)]:
            if gcd(s.core.source, q*r) != 1:
                continue
            for epsilon in (-1, 1):
                b = complete_block(s, (q, r), epsilon)
                actual = signed(core(s.core.source*q*r), epsilon)
                check((b is None) == (actual is None), "old/new block gates iff")
                block_controls += 1
    for a, b in product(signed_seeds, repeat=2):
        n, m = a.core.source, b.core.source
        if gcd(n, m) != 1:
            continue
        for epsilon in (-1, 1):
            actual = signed(core(n*m), epsilon) is not None
            small, large = (a, b) if n < m else (b, a)
            expected = (large.core.source == small.core.source+1 and
                        large.sign == -small.sign and epsilon == small.sign)
            expected = expected or (small.core.source == 2 and large.core.source == 3
                                    and epsilon == -1)
            check(actual == expected, "two signed-factor adjacency classification")

    recursion = []
    for n in (2, 6, 42, 47058):
        for epsilon in (-1, 1):
            r = complete_one(core(n), epsilon)
            if r is not None:
                check(r.target.index == 1, "primary index preservation")
                check(any(rr == r for _, rr in signed_parents(r.target)), "inverse parent recovery")
                recursion.append(asdict(r))
    check(signed_parents(signed(core(858), 1)) == (), "858 missing from one-prime signed grammar")
    check(complete_one(core(15), 1).prime == 2, "unsuccessful source may complete")
    check(complete_one(core(30), 1) is None and complete_one(core(30), -1) is None,
          "Giuga no single-prime exit")
    check(signed(core(30), 1) is not None and not korselt(core(30)), "Giuga not Carmichael")
    check(korselt(core(561)) and signed(core(561), 1) is None, "Carmichael not Giuga")
    check(gcd(15, core(15).period) == 1 and not korselt(core(15)), "normal support not enough")
    check(prime(6-1) and prime(6+1) and not prime(1806-1) and not prime(1806+1),
          "both twin-branch outcomes occur")
    check(1807 == 13*139 and 1805 == 5*19*19, "unconditional termination witness for one seed")
    goldbach_witnesses = []
    for n in (6, 42, 47058):
        lower, upper = primary_fork(signed(core(n), -1))
        pair = (lower.prime, upper.prime)
        triple = pair+(3,)
        check(pair == (n-1, n+1) and sum(pair) == 2*n and all(prime(p) for p in pair),
              "source-specific binary Goldbach witness")
        check(sum(triple) == 2*n+3 and all(prime(p) for p in triple),
              "source-specific ternary Goldbach witness")
        goldbach_witnesses.append({"source": n, "binary": pair, "binary_sum": 2*n,
                                   "ternary": triple, "ternary_sum": 2*n+3})
    check(primary_fork(signed(core(1806), -1)) is None, "fork failure is not Goldbach failure")
    check(old_joint_guard(core(30)) is None, "Giuga30 cannot retain its factors in any Carmichael multiple")
    g561 = old_joint_guard(core(561))
    check(g561 is not None and old_joint_accepts(g561, g561.residue), "561 partial joint class")
    check(not classical_giuga_sum_criterion(561*g561.residue), "old guard not new-prime completion")

    valid = complete_one(core(6), 1)
    hostiles = [
        (core, (True,)), (core, (6.0,)), (core, (1,)), (core, (12,)),
        (complete_one, (core(6), True)), (complete_one, (core(6), 0)),
        (audit_core, (replace(core(6), source=6.0),)),
        (audit_core, (replace(core(6), primes=(2, 3.0)),)),
        (audit_core, (replace(core(6), cofactor_sum=5.0),)),
        (audit_core, (replace(core(6), period=2.0),)),
        (audit_signed, (replace(signed(core(6), -1), index=1.0),)),
        (audit_completion, (replace(valid, prime=5.0),)),
        (audit_completion, (replace(valid, residue=5.0),)),
        (audit_completion, (replace(valid, source=core(10)),)),
        (audit_completion, (replace(valid, sign=-1),)),
        (audit_completion, (replace(valid, target=signed(core(42), -1)),)),
        (complete_block, (signed(core(6), -1), (5, 5), 1)),
        (complete_block, (signed(core(6), -1), (7, 5), 1)),
        (complete_block, (signed(core(6), -1), (2, 5), 1)),
        (complete_block, (signed(core(6), -1), (5, 9), 1)),
        (audit_old_guard, (replace(g561, residue=float(g561.residue)),)),
        (audit_old_guard, (replace(g561, modulus=float(g561.modulus)),)),
        (audit_old_guard, (replace(g561, source=core(15)),)),
        (old_joint_accepts, (g561, True)),
        (primary_fork, (signed(core(30), 1),)),
    ]
    for fn, args in hostiles:
        check(rejects(fn, *args), "malformed/source-switch rejection")

    data = {"status": "PROVED scoped algebra; FINITE-EXACT declared census; classical Giuga conjecture OPEN",
            "sign_convention": "A(N)=index*N+sign; + is Giuga-or-prime; - is weak primary pseudoperfect",
            "census_bound": 10000, "minus_sources": minus, "composite_plus_sources": giuga,
            "carmichael_sources": carmichael, "literal_power_sum_sources": literal_count,
            "block_controls": block_controls, "recursions": recursion,
            "old_joint_guard_controls": old_guard_controls,
            "old_joint_561": asdict(g561),
            "goldbach_witnesses": goldbach_witnesses,
            "hostile_858": {"signed_parents": [], "factorization": factor(858)},
            "hostile_defect_vectors": {str(n): defect_vector(core(n), 1) for n in (30, 561)},
            "malformed_controls": len(hostiles), "checks": CHECKS}
    print(json.dumps(data, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
