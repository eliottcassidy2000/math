"""Source-marked unimodular lattice and fixed-count prime-block capacity.

Exact arithmetic only. The finite 56-prime obligation is not an extension.
Run with python -B [-O] -X utf8; stdout is the deterministic JSON receipt.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod
import json

import giuga_korselt_defects_20261014f as base
import adaptive_giuga_completion_20261014g as adaptive


def integer(x, minimum=None):
    return base.integer(x, minimum)


@dataclass(frozen=True)
class Lattice:
    core: base.Core
    mu: int
    index: int
    vector: tuple
    norm: int
    signature: int


def lattice(n, mu):
    c = base.core(integer(n, 2))
    base.sheet(mu)
    k = F(c.cofactor_sum + mu, n)
    if k.denominator != 1:
        return None
    t = sum((F(1, p) for p in c.primes if p % 2), F())
    v0 = mu*n*(k-t)
    v = (v0,) + tuple((F(p % 2)-v0)/p for p in c.primes)
    if any(x.denominator != 1 for x in v):
        raise ArithmeticError("unimodular characteristic vector was not integral")
    v = tuple(int(x) for x in v)
    norm = int(k)*v[0]**2 + 2*v[0]*sum(v[1:])
    norm += sum(p*x*x for p, x in zip(c.primes, v[1:]))
    return Lattice(c, mu, int(k), v, norm, len(c.primes)+mu)


def audit_lattice(s):
    if type(s) is not Lattice or type(s.vector) is not tuple:
        raise ValueError("exact Lattice required")
    base.audit_core(s.core)
    base.sheet(s.mu)
    integer(s.index, 0)
    integer(s.norm)
    integer(s.signature)
    for v in s.vector:
        integer(v)
    if s != lattice(s.core.source, s.mu):
        raise ValueError("source, index or characteristic marker mismatch")
    return s


def gram(s):
    audit_lattice(s)
    ps = s.core.primes
    return ((s.index,) + (1,)*len(ps),) + tuple(
        (1,) + tuple(p if i == j else 0 for j in range(len(ps)))
        for i, p in enumerate(ps))


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    out = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F()
        if j != i:
            a[i], a[j] = a[j], a[i]
            out = -out
        pivot = a[i][i]
        out *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i]/pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio*a[i][k]
    return out


@dataclass(frozen=True)
class BlockGuard:
    source: Lattice
    count: int
    jump: int
    allowed_count: bool
    product_modulus: int
    product_residue: int
    half_leaf_sum: int


def block_guard(s, count, jump):
    """Necessary guard for same-sign, same-parity squarefree extensions."""
    audit_lattice(s)
    integer(count, 0)
    integer(jump)
    if s.core.source % 2:
        return BlockGuard(s, count, jump, (count-jump) % 8 == 0, 1, 0, 0)
    n = s.core.source//2
    r = (1 + 2*s.mu*pow(n, -1, 16)*(count-jump)) % 16
    delta = (s.index+len(s.core.primes)-1) % 2
    delta2 = (delta+jump+count) % 2
    return BlockGuard(s, count, jump, True, 16, r, (delta-delta2) % 4)


def audit_guard(g):
    if type(g) is not BlockGuard or type(g.allowed_count) is not bool:
        raise ValueError("exact block guard required")
    audit_lattice(g.source)
    integer(g.count, 0)
    integer(g.jump)
    integer(g.product_modulus, 1)
    integer(g.product_residue, 0)
    integer(g.half_leaf_sum, 0)
    if g != block_guard(g.source, g.count, g.jump):
        raise ValueError("wrong block guard")
    return g


def accepts(g, added):
    audit_guard(g)
    if type(added) is not tuple or len(added) != g.count:
        raise ValueError("exact added-prime tuple required")
    for p in added:
        integer(p, 3)
        if p % 2 == 0 or not base.prime(p) or g.source.core.source % p == 0:
            raise ValueError("new distinct odd primes required")
    if len(set(added)) != len(added):
        raise ValueError("repeated new prime")
    if not g.allowed_count:
        return False
    if g.source.core.source % 2:
        return True
    return (prod(added) % 16 == g.product_residue and
            sum((p-1)//2 for p in added) % 4 == g.half_leaf_sum)


def first_primes(n, count, forbidden=()):
    integer(n, 2)
    integer(count, 0)
    if type(forbidden) is not tuple:
        raise ValueError("forbidden tuple required")
    for p in forbidden:
        integer(p, 2)
    ps, p = [], 2
    while len(ps) < count:
        if n % p and p not in forbidden and base.prime(p):
            ps.append(p)
        p += 1
    return tuple(ps)


def reciprocal(ps):
    return sum((F(1, p) for p in ps), F())


def capacity_upper(n, count, forbidden=()):
    return F(1, n) + reciprocal(first_primes(n, count, forbidden))


def largest_prime_bound(n, count, jump):
    """If m=jump is feasible, m < 1/n + sum 1/q forces max q < reciprocal gap."""
    integer(n, 2)
    integer(count, 1)
    integer(jump, 1)
    gap = F(jump) - capacity_upper(n, count-1)
    if gap <= 0:
        return None
    bound = 1/gap
    return (bound.numerator-1)//bound.denominator


def mandatory_primes(n, count, jump):
    integer(jump, 1)
    return tuple(p for p in first_primes(n, count)
                 if capacity_upper(n, count, (p,)) <= jump)


@dataclass(frozen=True)
class Capacity:
    pool: tuple
    count: int
    modulus: int
    residue: int
    selected: tuple
    value: F


def bounded_capacity(pool, count, modulus, residue):
    """Exact maximum over one FINITE supplied prime universe; no guard sufficiency."""
    if type(pool) is not tuple:
        raise ValueError("exact pool tuple required")
    integer(count, 0)
    integer(modulus, 1)
    integer(residue, 0)
    if residue >= modulus:
        raise ValueError("noncanonical residue")
    for p in pool:
        integer(p, 2)
        if not base.prime(p):
            raise ValueError("pool must consist of primes")
    if tuple(sorted(set(pool))) != pool:
        raise ValueError("ordered distinct pool required")
    groups = {}
    for p in pool:
        groups.setdefault(p % modulus, []).append(p)
    dp = {(0, 1 % modulus): (F(), ())}
    for r, ps in sorted(groups.items()):
        pref = [F()]
        for p in ps:
            pref.append(pref[-1]+F(1, p))
        nxt = {}
        for (used, rem), (value, chosen) in dp.items():
            power = rem
            for j in range(min(count-used, len(ps))+1):
                key = used+j, power
                item = value+pref[j], tuple(sorted(chosen+tuple(ps[:j])))
                old = nxt.get(key)
                if old is None or item[0] > old[0] or (item[0] == old[0] and item[1] < old[1]):
                    nxt[key] = item
                power = power*r % modulus
        dp = nxt
    item = dp.get((count, residue))
    return None if item is None else Capacity(pool, count, modulus, residue, item[1], item[0])


def audit_capacity(c):
    if type(c) is not Capacity or type(c.selected) is not tuple or type(c.value) is not F:
        raise ValueError("exact capacity receipt required")
    for p in c.selected:
        integer(p, 2)
    if c != bounded_capacity(c.pool, c.count, c.modulus, c.residue):
        raise ValueError("capacity not maximal in its marked universe")
    return c


def thirty_obligation():
    s = lattice(30, -1)
    mandatory = mandatory_primes(30, 56, 1)
    cap = largest_prime_bound(30, 56, 1)
    pool = tuple(p for p in range(7, cap+1) if base.prime(p))
    best = bounded_capacity(pool, 56, 240, 31)
    return s, mandatory, cap, best


def forced_core_obstruction():
    """Exact all-old-prime CRT exclusion of a 56-prime extension of 30."""
    forced = mandatory_primes(30, 56, 1)
    cap = largest_prime_bound(30, 56, 1)
    factors = (2, 3, 5) + forced
    m = prod(factors)
    carry = sum(m//p for p in factors)
    residue = pow(carry, -1, m)
    available = tuple(p for p in range(7, cap+1) if p not in forced and base.prime(p))
    remaining = 56-len(forced)
    maximum = prod(available[-remaining:])
    return {"factors": factors, "core": m, "cofactor_sum": carry,
            "residue": residue, "remaining": remaining,
            "prime_cap": cap, "largest_allowed": available[-remaining:],
            "maximum": maximum, "excluded": maximum < residue < m}


def next57_obligation():
    forced = mandatory_primes(30,57,1)
    factors = (2,3,5)+forced
    m = prod(factors)
    carry = sum(m//p for p in factors)
    deficit = 2*m-carry
    remaining = 57-len(forced)
    return {"original_source":30,"added_count":57,"target_index":2,
            "mandatory":forced,"core":m,"deficit":deficit,
            "remaining":remaining,"first_prime_upper":(remaining*m-1)//deficit,
            "pending_residue":pow(carry,-1,m),"status":"OPEN: no completion tree run"}


def tests():
    checks = 0
    def check(ok, label):
        nonlocal checks
        if not ok:
            raise ArithmeticError(label)
        checks += 1
    def rejects(fn):
        try:
            fn()
        except (ValueError, TypeError):
            check(True, "rejection")
        else:
            check(False, "malformed packet accepted")

    receipts = []
    for n in range(2, 2501):
        if any(e != 1 for _, e in base.factor(n)):
            continue
        for mu in (-1, 1):
            s = lattice(n, mu)
            if s is None:
                check((base.derivative(n)+mu) % n != 0, "nonintegral index")
                continue
            receipts.append(s)
            b = gram(s)
            check(determinant(b) == mu, "determinant")
            check(all((sum(row[j]*s.vector[j] for j in range(len(row)))-row[i]) % 2 == 0
                      for i, row in enumerate(b)), "characteristic")
            check((s.norm-s.signature) % 8 == 0, "van der Blij control")
            formula = F(s.index) if n % 2 else s.index+F(mu*n+2, 4)
            check(formula == s.norm, "characteristic norm")
    for source, target, mu in ((2,6,1),(6,42,1),(42,1806,1),(3,15,-1)):
        s, t = lattice(source, mu), lattice(target, mu)
        if t is None:
            check(source == 3 and target == 15, "bad premise control")
            continue
        added = tuple(p for p in t.core.primes if p not in s.core.primes)
        g = block_guard(s, len(added), t.index-s.index)
        check(accepts(g, added), "actual same-sign block")
    for n in (30,858,1722,66198,2214408306):
        s = lattice(n,-1)
        check(s.index == 1, "known Giuga source")
        check(s.norm % 8 == s.signature % 8, "large norm")
    check(lattice(15,-1) is None, "Mersenne15 is outside lattice domain")
    check(lattice(31,-1).index == 0 and lattice(161,-1) is None, "Collatz domain not invariant")
    check((3*107+1)//2 == 161, "literal domain exit")
    rejects(lambda: lattice(63,-1))
    s = lattice(30,-1)
    check(s.vector == (-14,7,5,3) and s.norm == -6 and s.signature == 2, "even marked vector")
    check(block_guard(s,56,1).product_residue == 15, "even product guard")
    check(block_guard(s,56,1).half_leaf_sum == 1, "even leaf guard")
    check(not block_guard(lattice(3,-1),56,1).allowed_count, "odd count guard")
    for pool in ((2,3,5,7),(3,5,7,11,13),(7,11,13,17,19,23)):
        for count in range(5):
            for modulus in (1,4,8,15,16):
                for r in range(modulus):
                    row = bounded_capacity(pool,count,modulus,r)
                    candidates = [ps for ps in combinations(pool,count) if prod(ps)%modulus == r]
                    expected = max(map(reciprocal,candidates), default=None)
                    check((None if row is None else row.value) == expected, "independent subset capacity")
                    if row:
                        check(prod(row.selected)%modulus == r and len(row.selected) == count, "retained argmax")
    s, forced, cap, best = thirty_obligation()
    baseline = first_primes(30,56)
    check(capacity_upper(30,55) < 1 < capacity_upper(30,56) < 2, "56 count/index boundary")
    check(cap == 793 and not base.prime(793), "strict prime cap")
    check(forced == tuple(p for p in range(7,168) if base.prime(p)) and len(forced) == 36, "forced core")
    check(capacity_upper(30,56,(173,)) > 1, "first optional prime")
    check(best.selected == tuple(sorted(set(baseline)-{271}|{317})), "joint residue maximizer")
    check(F(1,30)+best.value > 1, "capacity is not an impossibility proof")
    check(prod(best.selected)%240 == 31, "joint product residue")
    check(not accepts(block_guard(s,56,1), best.selected), "product guard alone loses second characteristic constraint")
    target = 30*prod(best.selected)
    cofactor_sum = sum(target//p for p in (2,3,5)+best.selected)
    check((cofactor_sum-1)%target != 0, "maximizer not a Giuga extension")
    obstruction = forced_core_obstruction()
    check(obstruction["excluded"], "all-old-prime CRT excludes56")
    check((obstruction["cofactor_sum"]*obstruction["residue"])%obstruction["core"] == 1,
          "CRT inverse")
    for p in obstruction["factors"]:
        check((obstruction["core"]//p*obstruction["residue"]-1)%p == 0, "every old prime")
    pending = adaptive.pending(base.core(obstruction["core"]), -1)
    check(pending.residue == obstruction["residue"] and pending.modulus == obstruction["core"],
          "independent adaptive pending consumer")
    check(obstruction["remaining"] == 20, "remaining labels")
    check(obstruction["core"] == 962947420735983927056946215901134429196419130606213075415963491270,
          "forced core independent reproduction")
    check(obstruction["residue"] == 241154391290894464663467115463743362761848622423921060203644724609,
          "CRT independent reproduction")
    check(obstruction["maximum"] == 1294399569600457862304602788432784367142015736802901139771,
          "product bound independent reproduction")
    check(largest_prime_bound(30,57,1) is None, "declared57 stopping boundary")
    next_row = next57_obligation()
    check(capacity_upper(30,57) < 2, "57 target index2")
    check(next_row["mandatory"] == tuple(p for p in range(7,104) if base.prime(p)) and
          next_row["remaining"] == 33 and next_row["first_prime_upper"] == 185,
          "57 first adaptive branch")
    pending57 = adaptive.pending(base.core(next_row["core"]),-1)
    check(pending57.residue == next_row["pending_residue"], "57 whole remaining product")
    malformed = (True,1.0,-1)
    for x in malformed:
        rejects(lambda x=x: block_guard(s,x,1))
        rejects(lambda x=x: bounded_capacity((3,5),x,8,1))
    rejects(lambda: audit_lattice(replace(s,index=True)))
    rejects(lambda: accepts(block_guard(s,1,0),(3,)))
    rejects(lambda: bounded_capacity((3,3),1,8,1))
    rejects(lambda: audit_capacity(replace(best,value=1.0)))
    # Mersenne prime sources are index-zero lattices, not composite Giuga receipts.
    return {
        "checks": checks,
        "signed_squarefree_sources_through_2500": len(receipts),
        "thirty": {
            "count_lower_bound": 57, "target_index_if_56": 2,
            "new_prime_upper_bound": cap, "mandatory": list(forced),
            "forced_core": str(30*prod(forced)), "remaining_count": 20,
            "product_residue":31, "product_modulus":240,
            "maximizing_block":list(best.selected),
            "capacity_upper": str(F(1,30)+best.value),
            "baseline_upper":str(capacity_upper(30,56)),
            "no_extension_claim":True,
        },
        "crt_exclusion": {k: ([str(x) for x in v] if type(v) is tuple else
                              str(v) if type(v) is int else v)
                          for k,v in obstruction.items()},
        "next57_obligation": {k: (list(v) if type(v) is tuple else str(v) if type(v) is int else v)
                              for k,v in next_row.items()},
        "lattice_30": {"vector":list(s.vector),"norm":s.norm,"signature":s.signature},
    }


if __name__ == "__main__":
    print(json.dumps(tests(), indent=2, sort_keys=True))
