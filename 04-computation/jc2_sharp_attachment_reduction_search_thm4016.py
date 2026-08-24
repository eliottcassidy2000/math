#!/usr/bin/env python3
"""Exact finite-field reduction search for the THM-4007/4008 attachment point."""

from hashlib import sha256
from math import isqrt


def require(condition: bool, detail) -> None:
    if not condition:
        raise RuntimeError(detail)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def inv(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def add(P, Q, p: int):
    """Addition on y^2=x^3+1 over F_p; None is infinity."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = 3 * x1 * x1 * inv(2 * y1, p) % p
    else:
        lam = (y2 - y1) * inv(x2 - x1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return x3, y3


def mul(n: int, P, p: int):
    R = None
    Q = P
    while n:
        if n & 1:
            R = add(R, Q, p)
        Q = add(Q, Q, p)
        n >>= 1
    return R


def point_order(P, p: int) -> int:
    # Hasse gives a strict stopping bound; direct addition is an independent
    # order computation and does not assume a group-order formula.
    bound = p + 1 + 2 * isqrt(p) + 2
    R = None
    for n in range(1, bound + 1):
        R = add(R, P, p)
        if R is None:
            return n
    raise RuntimeError((p, P, bound))


def vp(n: int, ell: int) -> int:
    e = 0
    while n % ell == 0:
        n //= ell
        e += 1
    return e


def compatible_orders(p: int, mp: int, q: int, mq: int) -> bool:
    """Necessary compatibility if both are reductions of one torsion point.

    At good residue characteristic p, reduction preserves the prime-to-p
    order and can only lose p-power order.
    """
    primes = {r for r in range(2, max(mp, mq) + 1) if is_prime(r)}
    for ell in primes:
        ep = vp(mp, ell)
        eq = vp(mq, ell)
        if ell != p and ell != q and ep != eq:
            return False
        if ell == p and ell != q and ep > eq:
            return False
        if ell == q and ell != p and eq > ep:
            return False
    return True


def main() -> None:
    records = []
    for p in range(5, 500):
        if not is_prime(p) or p in (7,):
            continue
        a = 43 * inv(84, p) % p
        b = 127 * inv(84, p) % p
        xs = [x for x in range(p) if x**3 % p == a]
        ys = [y for y in range(p) if y*y % p == b]
        if not xs or not ys:
            continue
        orders = sorted({point_order((x, y), p) for x in xs for y in ys})
        for x in xs:
            for y in ys:
                require((y*y - x**3 - 1) % p == 0, (p, x, y))
        records.append((p, tuple(xs), tuple(ys), tuple(orders)))

    record_text = "\n".join(repr(rec) for rec in records)
    print(f"DEGREE_ONE_RECORDS {len(records)}")
    print(f"RECORDS_SHA256 {sha256(record_text.encode('ascii')).hexdigest()}")
    print("FIRST_INCOMPATIBLE_UNIQUE_CUBE_ROOT_PAIR")
    unique = [(p, xs[0], ys[0], orders[0]) for p, xs, ys, orders in records
              if len(xs) == 1 and len(orders) == 1]
    for i, (p, xp, yp, mp) in enumerate(unique):
        for q, xq, yq, mq in unique[i + 1:]:
            if not compatible_orders(p, mp, q, mq):
                print((p, xp, yp, mp), (q, xq, yq, mq))
                # Exact witness multiples and curve checks.
                Pp = (xp, yp)
                Pq = (xq, yq)
                require(mul(mp, Pp, p) is None, (p, mp))
                require(all(mul(mp // ell, Pp, p) is not None
                            for ell in range(2, mp + 1)
                            if is_prime(ell) and mp % ell == 0), (p, mp))
                require(mul(mq, Pq, q) is None, (q, mq))
                require(all(mul(mq // ell, Pq, q) is not None
                            for ell in range(2, mq + 1)
                            if is_prime(ell) and mq % ell == 0), (q, mq))
                print("INDEPENDENT REDUCTION SEARCH PASSED")
                return
    raise RuntimeError("no incompatible pair found")


if __name__ == "__main__":
    main()
