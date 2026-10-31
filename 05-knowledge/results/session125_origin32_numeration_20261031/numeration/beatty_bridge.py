"""Exact Beatty division supplement; no floating point, old files unchanged."""

from math import isqrt


def natural(n):
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError("argument must be an integer")
    if n < 0:
        raise ValueError("argument must be nonnegative")


def shift(n):
    """ceil(beta*n), beta=(3+sqrt(5))/2; canonical CG append-zero."""
    natural(n)
    return 0 if n == 0 else (3 * n + isqrt(5 * n * n)) // 2 + 1


def parent(n):
    """floor(n/beta), computed exactly by an integer square root."""
    natural(n)
    return 0 if n == 0 else (3 * n - isqrt(5 * n * n) - 1) // 2


def gap(n):
    natural(n)
    return shift(n + 1) - shift(n)


def divide(n):
    """Return unique (p,d) with n=shift(p)+d and 0<=d<gap(p)."""
    p = parent(n)
    return p, n - shift(p)


def encode(n):
    """Canonical CG word from repeated Beatty division, not weight lookup."""
    natural(n)
    if n == 0:
        return "0"
    digits = []
    while n:
        n, digit = divide(n)
        digits.append(str(digit))
    return "".join(reversed(digits))


def decode(word):
    """Read legal digits using only exact Beatty multiplication."""
    if not isinstance(word, str):
        raise TypeError("word must be a string")
    n, locked = 0, False
    for char in word:
        if char not in "012" or (locked and char == "2"):
            raise ValueError("word is not a legal Chung--Graham representation")
        d = int(char)
        n = shift(n) + d
        if d == 0:
            locked = False
        elif d == 2:
            locked = True
    return n


def free_state(n):
    """Legality state of the canonical prefix, without computing its digits."""
    return gap(n) == 3


def free_count(start, length):
    """Number of Free prefixes in [start,start+length)."""
    natural(start)
    natural(length)
    return shift(start + length) - shift(start) - 2 * length


def orbit_norm(n):
    """Positive integer D(C) parameter for n>0; zero at zero."""
    p = shift(n)
    return n * n - 3 * n * p + p * p


def iterate_shift(n, steps):
    natural(n)
    natural(steps)
    for _ in range(steps):
        n = shift(n)
    return n


def orbit_core(n):
    """Unique (primitive positive core, shift depth); convention (0,0) at 0."""
    natural(n)
    if n == 0:
        return 0, 0
    depth = 0
    while True:
        p, digit = divide(n)
        if digit:
            return n, depth
        n, depth = p, depth + 1
