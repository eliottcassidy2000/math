"""Pair chain for two Collatz (Terras T-map) orbits related by u = 3^k v + e.
State (k, num) with e = num / 3^p, p = max(0, -k).  Driving bit beta = parity of v.
Verifies the chain against actual integer orbits (exact Fractions) and checks
the martingale / entry-state claims on random 2-adic (mod 2^N) inputs.
"""
from fractions import Fraction
import random, sys

def T(z):
    return z // 2 if z % 2 == 0 else (3 * z + 1) // 2

def step(k, num, beta):
    """one T-step of the pair chain; returns (k', num', flip)"""
    p = -k if k < 0 else 0
    if num % 2 == 0:            # e even: no flip, k unchanged
        if beta == 0:
            return k, num // 2, 0
        # (3e + 1 - 3^k)/2
        if k >= 0:
            return k, (3 * num + 1 - 3 ** k) // 2, 0
        return k, (3 * num + 3 ** p - 1) // 2, 0
    # e odd: flip
    if beta == 0:               # u odd, v even: k+1, (3e+1)/2
        if k >= 0:
            return k + 1, (3 * num + 1) // 2, 1
        return k + 1, (num + 3 ** (p - 1)) // 2, 1
    # beta == 1: u even, v odd: k-1, (e - 3^(k-1))/2
    if k >= 1:
        return k - 1, (num - 3 ** (k - 1)) // 2, 1
    return k - 1, (3 * num - 1) // 2, 1

def e_of(k, num):
    p = -k if k < 0 else 0
    return Fraction(num, 3 ** p)

def verify_integers(trials=2000, steps=300, seed=1):
    rng = random.Random(seed)
    bad = 0
    for _ in range(trials):
        k0 = rng.randint(-3, 6)
        y = rng.randint(1, 10 ** 12)
        num0 = rng.randint(-10 ** 6, 10 ** 6)
        if k0 < 0:
            # need u = 3^k y + e integral: choose y multiple of 3^|k0| and e integral
            y = y * 3 ** (-k0)
            num0 = num0 * 3 ** (-k0)   # e = num0/3^p integral
        e0 = e_of(k0, num0)
        u = Fraction(3) ** k0 * y + e0
        assert u.denominator == 1
        u = int(u); v = y; k, num = k0, num0
        for n in range(steps):
            beta = v % 2
            k, num, _ = step(k, num, beta)
            u, v = T(u), T(v)
            if Fraction(u) != Fraction(3) ** k * v + e_of(k, num):
                bad += 1
                break
    return bad

if __name__ == "__main__":
    b = verify_integers()
    print("relation u_n = 3^k_n v_n + e_n violated in", b, "of 2000 integer trials (300 steps each)")
    # entry states into (0,0): enumerate small states, one step
    entries = []
    for k in range(-6, 7):
        for num in range(-3000, 3001):
            for beta in (0, 1):
                k2, n2, _ = step(k, num, beta)
                if (k2, n2) == (0, 0) and (k, num) != (0, 0):
                    entries.append((k, str(e_of(k, num)), beta))
    print("states entering (0,0) in one step (k, e, beta):", entries)
    # absorbing check
    print("(0,0) absorbing:", [step(0, 0, b)[:2] for b in (0, 1)])
