"""Exact oriented wheel sums, collision correction, and the Apéry form.

All odd finite fields: elementary moment/collision identities.
Odd primes: Ahlgren--Ono Theorem 6 supplies the modular evaluation.
Nonprime field orders: direct exact tests of the Frobenius extension only.
"""
from dataclasses import dataclass, asdict, replace
from functools import lru_cache
from itertools import product
from math import comb, isqrt
import json


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer outside domain")
    return n


def prime(p):
    integer(p)
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def odd_prime(p):
    integer(p, 3)
    if p % 2 == 0 or not prime(p):
        raise ValueError("odd prime characteristic required")
    return p


def remainder(a, b, p):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    while len(a) >= len(b) and any(a):
        k, c = len(a) - len(b), a[-1] * pow(b[-1], -1, p) % p
        for j, coefficient in enumerate(b):
            a[k + j] = (a[k + j] - c * coefficient) % p
        while len(a) > 1 and a[-1] == 0:
            a.pop()
    return tuple(a)


@dataclass(frozen=True)
class Field:
    p: int
    modulus: tuple  # ascending coefficients, monic and irreducible

    def __post_init__(self):
        odd_prime(self.p)
        f = self.modulus
        if type(f) is not tuple or len(f) < 2 or any(
                type(c) is not int or not 0 <= c < self.p for c in f) or f[-1] != 1:
            raise ValueError("canonical monic polynomial required")
        degree = len(f) - 1
        for d in range(1, degree // 2 + 1):
            for coefficients in product(range(self.p), repeat=d):
                if remainder(f, coefficients + (1,), self.p) == (0,):
                    raise ValueError("reducible field modulus")

    @property
    def degree(self):
        return len(self.modulus) - 1

    @property
    def q(self):
        return self.p ** self.degree


def source(field):
    if type(field) is not Field:
        raise ValueError("a marked finite field is required")
    return field


@lru_cache(maxsize=None)
def _tables(field):
    p, degree, q = field.p, field.degree, field.q
    digits = [tuple((x // p ** i) % p for i in range(degree)) for x in range(q)]

    def encode(coefficients):
        return sum(c * p ** i for i, c in enumerate(coefficients))

    add, sub, mul = [], [], []
    for x in digits:
        ar, sr, mr = [], [], []
        for y in digits:
            ar.append(encode(tuple((a + b) % p for a, b in zip(x, y))))
            sr.append(encode(tuple((a - b) % p for a, b in zip(x, y))))
            z = [0] * (2 * degree - 1)
            for i, a in enumerate(x):
                for j, b in enumerate(y):
                    z[i + j] = (z[i + j] + a * b) % p
            mr.append(encode(remainder(z, field.modulus, p)))
        add.append(tuple(ar)); sub.append(tuple(sr)); mul.append(tuple(mr))
    chi = [-1] * q
    chi[0] = 0
    for x in range(1, q):
        chi[mul[x][x]] = 1
    inverse = (0,) + tuple(next(y for y in range(1, q) if mul[x][y] == 1)
                           for x in range(1, q))
    return tuple(add), tuple(sub), tuple(mul), tuple(chi), inverse


def tables(field):
    return _tables(source(field))


def moment(field):
    """Sum_y chi(y) (sum_x chi(x(x-1)(x-y)))^2, including singular y."""
    _, sub, mul, chi, _ = tables(field)
    traces = tuple(sum(chi[mul[mul[x][sub[x][1]]][sub[x][y]]]
                       for x in range(field.q)) for y in range(field.q))
    return sum(chi[y] * t * t for y, t in enumerate(traces))


def wheel(field, injective=True, free_hub=False, reversed_edges=0):
    """Fixed labelled W4 orientation; reverse exactly the supplied edge count."""
    source(field)
    if type(injective) is not bool or type(free_hub) is not bool:
        raise ValueError("exact Boolean mode required")
    integer(reversed_edges)
    if reversed_edges > 8:
        raise ValueError("W4 has eight edges")
    _, sub, _, chi, _ = tables(field)
    value = (field.q - 1) * (moment(field) + (field.q - 4 if injective else 0))
    return value * (field.q if free_hub else 1) * chi[sub[0][1]] ** reversed_edges


def direct_normalized(field, injective=True):
    """Independent literal triple sum with hub0 and first rim vertex1."""
    source(field)
    if type(injective) is not bool:
        raise ValueError("exact Boolean mode required")
    _, sub, mul, chi, _ = tables(field)
    value = 0
    for x, y, z in product(range(field.q), repeat=3):
        if injective and len({0, 1, x, y, z}) != 5:
            continue
        v = 1
        for a in (x, y, z, sub[1][x], sub[x][y], sub[y][z], sub[z][1]):
            v = mul[v][a]
        value += chi[v]
    return value


def collision_sums(field):
    """The two opposite-rim diagonals and their intersection, normalized."""
    source(field)
    _, sub, mul, chi, _ = tables(field)
    first = second = both = 0
    for x, z in product(range(field.q), repeat=2):
        if x not in (0, 1) and z not in (0, 1):
            first += chi[mul[x][z]]
    for x, y in product(range(field.q), repeat=2):
        if x not in (0, 1) and y not in (0, x):
            second += chi[y]
    for x in range(field.q):
        if x not in (0, 1):
            both += 1
    return first, second, both


def calabi_count(field):
    """Literal count of four nonzero x with sum(x+1/x)=0, via pair convolution."""
    add, sub, _, _, inverse = tables(field)
    q = field.q
    frequency = [0] * q
    for x in range(1, q):
        frequency[add[x][inverse[x]]] += 1
    pairs = [0] * q
    for a, b in product(range(q), repeat=2):
        pairs[add[a][b]] += frequency[a] * frequency[b]
    return sum(pairs[a] * pairs[sub[0][a]] for a in range(q))


@lru_cache(maxsize=None)
def _eta_coefficients(limit):
    coefficients = [0] * (limit + 1)
    if limit:
        coefficients[1] = 1
    for step in (2, 4):
        for exponent in range(step, limit + 1, step):
            for _ in range(4):
                for n in range(limit, exponent - 1, -1):
                    coefficients[n] -= coefficients[n - exponent]
    return tuple(coefficients)


def eta_coefficients(limit):
    integer(limit)
    return _eta_coefficients(limit)


def four_square_coefficient(n):
    """Jacobi-product formula; all four coordinates positive and odd."""
    integer(n, 1)
    target = 4 * n
    pairs = {}
    for a in range(1, isqrt(target) + 1, 2):
        for b in range(1, isqrt(target - a * a) + 1, 2):
            key = a * a + b * b
            pairs[key] = pairs.get(key, 0) + (-1) ** ((a - 1) // 2) * a
    return sum(value * pairs.get(target - key, 0) for key, value in pairs.items())


def prime_wheel(p):
    """CITED Ahlgren--Ono plus the elementary injective collision correction."""
    odd_prime(p)
    return -(p - 1) * (eta_coefficients(p)[p] + 4)


def frobenius_trace(p, degree):
    """Trace of the modular local factor, not a proved general wheel evaluator."""
    odd_prime(p); integer(degree, 1)
    a = eta_coefficients(p)[p]
    before, current = 2, a
    for _ in range(1, degree):
        before, current = current, a * current - p ** 3 * before
    return current


def apery(n):
    integer(n)
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))


@dataclass(frozen=True)
class WheelReceipt:
    field: Field
    normalized_all: int
    anchored_all: int
    anchored_injective: int
    collision: int


def compile_receipt(field):
    source(field)
    q, s = field.q, moment(field)
    return WheelReceipt(field, s, (q - 1) * s,
                        (q - 1) * (s + q - 4), (q - 1) * (4 - q))


def audit_receipt(receipt, field):
    source(field)
    if type(receipt) is not WheelReceipt or type(receipt.field) is not Field:
        raise ValueError("typed receipt required")
    for value in (receipt.normalized_all, receipt.anchored_all,
                  receipt.anchored_injective, receipt.collision):
        if type(value) is not int:
            raise ValueError("exact signed counts required")
    if receipt.field != field or receipt != compile_receipt(field):
        raise ValueError("wrong field or altered wheel count")
    return receipt


def main():
    checks = 0

    def check(condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(f"control {checks} failed")

    def reject(call):
        try:
            call()
        except (ValueError, TypeError):
            check(True)
        else:
            check(False)

    primes = [p for p in range(3, 200, 2) if prime(p)]
    coefficients = eta_coefficients(199)
    for n in range(1, 200):
        check(coefficients[n] == four_square_coefficient(n))
    rows = []
    for p in primes:
        f = Field(p, (0, 1)); s = moment(f)
        r = compile_receipt(f)
        check(audit_receipt(r, f) == r)
        check(s == -coefficients[p] - p)
        check(r.anchored_injective == prime_wheel(p))
        check(collision_sums(f) == (1, 1, p - 2))
        check(calabi_count(f) == p ** 3 - 2 * p ** 2 - 7 - coefficients[p])
        check((coefficients[p] - apery((p - 1) // 2)) % (p * p) == 0)
        check((r.anchored_injective // (p - 1) + p ** 3 + 5) % 32 == 0)
        check(coefficients[p] ** 2 <= 4 * p ** 3)
        check(wheel(f, free_hub=True) == p * r.anchored_injective)
        check(wheel(f, reversed_edges=1) == (-1) ** ((p - 1) // 2) * prime_wheel(p))
        if p <= 31:
            check(direct_normalized(f, False) == s)
            check(direct_normalized(f, True) * (p - 1) == prime_wheel(p))
        rows.append({"p": p, "a_p": coefficients[p], "all_map": r.anchored_all,
                     "injective": r.anchored_injective, "normalized_collision": 4 - p})

    extensions = [(3, (1, 0, 1)), (3, (1, 2, 0, 1)), (3, (2, 1, 0, 0, 1)),
                  (5, (2, 0, 1)), (5, (1, 1, 0, 1)), (7, (1, 0, 1)),
                  (11, (1, 0, 1)), (13, (2, 0, 1))]
    extension_rows = []
    for p, polynomial in extensions:
        f = Field(p, polynomial); q = f.q
        r = compile_receipt(f); t = frobenius_trace(p, f.degree)
        check(audit_receipt(r, f) == r)
        check(collision_sums(f) == (1, 1, q - 2))
        check(r.anchored_all - r.anchored_injective == (q - 1) * (4 - q))
        check(calabi_count(f) == q ** 3 - 2 * q ** 2 + q - 7 + r.normalized_all)
        check(r.normalized_all == -q - t)  # finite-exact extension control, not proof
        if q <= 49:
            check(direct_normalized(f, False) == r.normalized_all)
            check(direct_normalized(f, True) * (q - 1) == r.anchored_injective)
        add, sub, mul, chi, inverse = tables(f)
        for x in range(q):
            check(add[x][0] == x and sub[x][x] == 0 and mul[x][1] == x)
            check(x == 0 or mul[x][inverse[x]] == 1)
            for y in range(q):
                check(chi[mul[x][y]] == chi[x] * chi[y])
        extension_rows.append({"p": p, "degree": f.degree, "q": q,
                               "modulus": list(polynomial), "moment": r.normalized_all,
                               "injective": r.anchored_injective,
                               "frobenius_trace": t, "fourier_a_q": eta_coefficients(q)[q]})

    # Literal field model change: two irreducible F9 presentations agree.
    f9 = Field(3, (1, 0, 1)); f9b = Field(3, (2, 1, 1))
    check(wheel(f9) == wheel(f9b) == 272)
    check(eta_coefficients(9)[9] == -11 and frobenius_trace(3, 2) == -38)
    check(-(9 - 1) * (eta_coefficients(9)[9] + 4) == 56 != wheel(f9))
    check(prime_wheel(7) == -168 and wheel(Field(7, (0, 1)), False) == -186)
    check(prime_wheel(11) == 400)
    # Full anchored enumeration, independent of scaling normalization.
    for p in (3, 5, 7):
        f = Field(p, (0, 1)); _, sub, mul, chi, _ = tables(f)
        for hub in range(p):
            full = inj = 0
            for rim in product(range(p), repeat=4):
                v = 1
                for i in range(4):
                    v = mul[v][sub[rim[i]][hub]]
                    v = mul[v][sub[rim[i]][rim[(i + 1) % 4]]]
                full += chi[v]
                if len(set((hub,) + rim)) == 5:
                    inj += chi[v]
            check(full == wheel(f, False) and inj == wheel(f))

    f7 = Field(7, (0, 1)); receipt = compile_receipt(f7)
    hostile = [lambda: Field(True, (0, 1)), lambda: Field(9, (0, 1)),
               lambda: Field(2, (0, 1)), lambda: Field(3, (True, 0, 1)),
               lambda: Field(3, (0, 0, 1)), lambda: Field(3, [0, 1]),
               lambda: prime_wheel(7.0), lambda: prime_wheel(9),
               lambda: eta_coefficients(True), lambda: frobenius_trace(3, 0),
               lambda: frobenius_trace(3, True), lambda: wheel(f7, 1),
               lambda: wheel(f7, reversed_edges=True), lambda: wheel(f7, reversed_edges=9),
               lambda: audit_receipt(replace(receipt, normalized_all=float(receipt.normalized_all)), f7),
               lambda: audit_receipt(replace(receipt, anchored_injective=0), f7),
               lambda: audit_receipt(receipt, Field(11, (0, 1))),
               lambda: audit_receipt(compile_receipt(f9), f9b)]
    for call in hostile:
        reject(call)
    output = {"checks": checks,
              "scope": {"elementary": "all odd finite fields; anchored oriented W4",
                        "modular_theorem": "all odd primes, Ahlgren--Ono Theorem 6",
                        "prime_power_trace_identity": "finite-exact tests only in this packet",
                        "prime_limit": 199, "direct_prime_limit": 31,
                        "extension_fields": len(extensions), "hostiles": len(hostile)},
              "prime_rows": rows, "extension_rows": extension_rows,
              "hostile_F9": {"actual_injective": 272, "naive_fourier_formula": 56,
                              "frobenius_trace": -38, "fourier_a9": -11}}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
