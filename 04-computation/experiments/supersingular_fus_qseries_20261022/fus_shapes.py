#!/usr/bin/env python3
"""F, U, S at the shapes p^2, pq, p^2 q (opus-2026-10-22-ssj).

Definitions (arithmetic_braids_20260917_divisors.md, DB1), for N >= 2:
    D(N) = {d : d | N, 1 < d < N}
    F(N) = |D(N)|
    S(N) = #{d in D(N) : d squarefree}
    U(N) = #{d in D(N) : d prime}
    defect delta(N) = F - U - S.

Parts
  1. F, U, S at p^2, pq, p^2 q by trial division over all primes below PMAX
     (all ordered pairs of distinct primes), against the closed forms
     (1,1,1), (2,2,2), (4,2,3); controls p, p^3, p^2 q r; hostiles
     (the collisions q = p, the literal exponent reading p^q r).
  2. Exhaustive level sets of delta for 2 <= N <= NMAX by a divisor sieve,
     typed by exponent profile: delta = -1 iff p^2 or p^2 q; delta = -2 iff pq;
     delta = 0 iff p, p^3, p^2 q r; delta < 0 iff p^2, p^2 q, or squarefree
     composite.
  3. The prime-adjunction recursion: for a prime r not dividing N,
        delta(N r) = 2 delta(N) + U(N) - [N squarefree] - [N prime],
     checked for all N <= RMAX and the three smallest primes r not dividing N.
  4. The Klein action (C, Q) on D(p^2 q r): orbit representatives in the
     closed quadrant f >= 0, g >= 0 of the (f, g) = (a-b-c, b-c) chart of
     THM-4733 are exactly {p, p^2, pq, p^2 q}.
Run: python3 -B fus_shapes.py > fus_shapes.out
"""
from __future__ import annotations

import sys
from collections import Counter

PMAX = 400
NMAX = 2_000_000
RMAX = 20_000


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise SystemExit("FAILED: " + msg)


def primes_below(n: int) -> list[int]:
    sieve = bytearray([1]) * n
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, n, i)))
    return [i for i in range(n) if sieve[i]]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def is_squarefree(n: int) -> bool:
    i = 2
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True


def fus_trial(N: int) -> tuple[int, int, int]:
    """F, U, S by listing every divisor 1 < d < N with trial division."""
    F = U = S = 0
    for d in range(2, N):
        if N % d:
            continue
        F += 1
        if is_squarefree(d):
            S += 1
        if is_prime(d):
            U += 1
    return F, U, S


def fus_divisor_pairs(N: int) -> tuple[int, int, int]:
    """Same counts, enumerating divisors in pairs (d, N/d) up to sqrt(N)."""
    divs = set()
    i = 1
    while i * i <= N:
        if N % i == 0:
            divs.add(i)
            divs.add(N // i)
        i += 1
    D = [d for d in divs if 1 < d < N]
    return len(D), sum(1 for d in D if is_prime(d)), sum(1 for d in D if is_squarefree(d))


def part1() -> None:
    print("== Part 1: F, U, S at p^2, pq, p^2 q (trial division, all primes < %d) ==" % PMAX)
    P = primes_below(PMAX)
    closed = {"p^2": (1, 1, 1), "pq": (2, 2, 2), "p^2q": (4, 2, 3)}
    counts = Counter()
    for p in P:
        got = fus_divisor_pairs(p * p)
        require(got == closed["p^2"], f"p^2 at p={p}: {got}")
        counts["p^2"] += 1
    small = [p for p in P if p < 60]
    for p in P:
        for q in P:
            if q == p:
                continue
            if p * q < 10 ** 9:
                got = fus_divisor_pairs(p * q)
                require(got == closed["pq"], f"pq at {p},{q}: {got}")
                counts["pq"] += 1
            if p * p * q < 10 ** 9:
                got = fus_divisor_pairs(p * p * q)
                require(got == closed["p^2q"], f"p^2q at {p},{q}: {got}")
                counts["p^2q"] += 1
    # independent slow path (plain trial division over every d < N) on small primes
    slow = 0
    for p in small:
        require(fus_trial(p * p) == closed["p^2"], f"slow p^2 {p}")
        slow += 1
        for q in small:
            if q != p:
                require(fus_trial(p * q) == closed["pq"], f"slow pq {p},{q}")
                if p * p * q < 200_000:
                    require(fus_trial(p * p * q) == closed["p^2q"], f"slow p^2q {p},{q}")
                    slow += 1
                slow += 1
    for k in ("p^2", "pq", "p^2q"):
        F, U, S = closed[k]
        print(f"  {k:5s}: (F, U, S) = ({F}, {U}, {S}),  F - U - S = {F - U - S:+d};"
              f"  confirmed on {counts[k]} prime choices")
    print(f"  slow independent trial-division path agrees on {slow} cases (primes < 60)")
    # controls and hostiles
    print("  controls (balanced shapes):")
    for name, N in (("p = 7", 7), ("p^3 = 27", 27), ("p^3 = 8", 8), ("p^2qr = 60", 60),
                    ("p^2qr = 2^2*7*11", 4 * 7 * 11)):
        F, U, S = fus_trial(N)
        require(F == U + S, f"control {name}")
        print(f"    {name:18s} (F,U,S) = ({F},{U},{S})  balanced")
    print("  hostiles:")
    for name, N, want in (("q = p in pq -> p^2 (p=5)", 25, (1, 1, 1)),
                          ("q = p in p^2q -> p^3 (p=5)", 125, (2, 1, 1)),
                          ("literal p^q r, p=2,q=3,r=5 -> 2^3*5", 40, (6, 2, 3)),
                          ("pqr = 30", 30, (6, 3, 6))):
        got = fus_trial(N)
        require(got == want, f"hostile {name}: {got}")
        F, U, S = got
        print(f"    {name:36s} (F,U,S) = ({F},{U},{S}), defect {F - U - S:+d}")
    print("  (the collision q = p sends p^2 q to the balanced shape p^3: the defect jumps -1 -> 0)")


def part2() -> None:
    print("\n== Part 2: exhaustive defect level sets for 2 <= N <= %d ==" % NMAX)
    n = NMAX + 1
    # smallest prime factor sieve
    spf = list(range(n))
    i = 2
    while i * i < n:
        if spf[i] == i:
            for k in range(i * i, n, i):
                if spf[k] == k:
                    spf[k] = i
        i += 1
    by_def = {}
    shape_examples = {}
    for N in range(2, n):
        m = N
        exps = []
        while m > 1:
            p = spf[m]
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            exps.append(e)
        r = len(exps)
        tau = 1
        for e in exps:
            tau *= e + 1
        sqfree = all(e == 1 for e in exps)
        prime = (r == 1 and exps[0] == 1)
        F = tau - 2
        S = 2 ** r - 1 - (1 if sqfree else 0)
        U = r - (1 if prime else 0)
        d = F - U - S
        prof = tuple(sorted(exps, reverse=True))
        if d <= 0:
            by_def.setdefault(d, Counter())[prof] += 1
            shape_examples.setdefault((d, prof), N)
    for d in sorted(by_def):
        profs = by_def[d]
        txt = ", ".join(f"{prof}:{cnt}" for prof, cnt in sorted(profs.items()))
        print(f"  defect {d:+d}: profiles(count) {txt}")
    want = {
        -1: {(2,), (2, 1)},
        0: {(1,), (3,), (2, 1, 1)},
        -2: {(1, 1)},
        -3: {(1, 1, 1)},
        -4: {(1, 1, 1, 1)},
        -5: {(1, 1, 1, 1, 1)},
        -6: {(1, 1, 1, 1, 1, 1)},
        -7: {(1, 1, 1, 1, 1, 1, 1)},
    }
    for d, shapes in want.items():
        if d in by_def:
            require(set(by_def[d]) == shapes, f"level set {d}: {set(by_def[d])}")
    require(set(by_def) == set(want), f"negative/zero levels present: {sorted(by_def)}")
    print("  verdict: delta = -1 exactly at p^2 and p^2 q; delta = -2 exactly at pq;")
    print("           delta = 0 exactly at p, p^3, p^2 q r; every other delta < 0 is a squarefree")
    print("           composite with delta = -omega (DB2 of arithmetic_braids, re-run exhaustively).")


def fus_formula(exps: list[int]) -> tuple[int, int, int]:
    r = len(exps)
    tau = 1
    for e in exps:
        tau *= e + 1
    sqfree = all(e == 1 for e in exps)
    prime = (r == 1 and exps[0] == 1)
    return tau - 2, r - (1 if prime else 0), 2 ** r - 1 - (1 if sqfree else 0)


def factor_exps(N: int) -> dict[int, int]:
    out = {}
    p = 2
    while p * p <= N:
        while N % p == 0:
            out[p] = out.get(p, 0) + 1
            N //= p
        p += 1
    if N > 1:
        out[N] = out.get(N, 0) + 1
    return out


def part3() -> None:
    print("\n== Part 3: prime-adjunction recursion delta(N r) = 2 delta(N) + U(N) - [N sqfree] - [N prime] ==")
    checks = 0
    for N in range(2, RMAX + 1):
        fe = factor_exps(N)
        F, U, S = fus_formula(list(fe.values()))
        d = F - U - S
        sqfree = all(e == 1 for e in fe.values())
        prime = (len(fe) == 1 and list(fe.values())[0] == 1)
        rs = []
        r = 2
        while len(rs) < 3:
            if N % r and is_prime(r):
                rs.append(r)
            r += 1
        for r in rs:
            F2, U2, S2 = fus_formula(list(factor_exps(N * r).values()))
            d2 = F2 - U2 - S2
            pred = 2 * d + U - (1 if sqfree else 0) - (1 if prime else 0)
            require(d2 == pred, f"recursion at N={N}, r={r}")
            checks += 1
    print(f"  recursion verified on {checks} pairs (N <= {RMAX}, three primes r per N)")
    print("  chains:  p  --(x q)-->  pq        :  0  -> -2   (away from balance)")
    print("           p^2 --(x q)--> p^2 q     : -1  -> -1   (2(-1) + U(p^2) = -2 + 1)")
    print("           p^2 q --(x r)--> p^2 q r : -1  ->  0   (2(-1) + U(p^2 q) = -2 + 2)")
    print("           p^2 --(x p)--> p^3       : -1  ->  0   (exponent raise: tau 3 -> 4)")


def part4() -> None:
    print("\n== Part 4: Klein orbits on D(p^2 q r) and the quadrant representatives ==")
    pts = [(a, b, c) for a in range(3) for b in range(2) for c in range(2)
           if (a, b, c) not in ((0, 0, 0), (2, 1, 1))]
    C = lambda v: (2 - v[0], 1 - v[1], 1 - v[2])
    Q = lambda v: (v[0], v[2], v[1])
    name = {}
    for a, b, c in pts:
        s = ("p" if a == 1 else ("p^2" if a == 2 else "")) + ("q" if b else "") + ("r" if c else "")
        name[(a, b, c)] = s
    seen = set()
    reps = []
    for v in pts:
        if v in seen:
            continue
        orb = {v, C(v), Q(v), C(Q(v))}
        seen |= orb
        quad = [w for w in orb if (w[0] - w[1] - w[2]) >= 0 and (w[1] - w[2]) >= 0]
        require(len(quad) == 1, f"orbit {orb} quadrant reps {quad}")
        reps.append(quad[0])
        print(f"  orbit {sorted(name[w] for w in orb)}  ->  quadrant representative {name[quad[0]]}")
    require(sorted(name[w] for w in reps) == sorted(["p", "p^2", "pq", "p^2q"]), "reps")
    print("  the four Klein orbits have quadrant representatives {p, p^2, pq, p^2 q} = p * Div(pq):")
    print("  the three asked shapes plus p, i.e. exactly the divisors of p^2 q divisible by p.")


if __name__ == "__main__":
    print(f"python {sys.version.split()[0]}")
    part1()
    part2()
    part3()
    part4()
    print("\nALL CHECKS PASSED")
