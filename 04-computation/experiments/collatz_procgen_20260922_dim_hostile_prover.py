"""Exact hostility prover / descent finder for negative 2-adic rationals in the forward E-game.

E-moves (2-adic): odd v -> 3v+1 (forced); even v -> v/2 or 3v+1.  x is HOSTILE (x in Bad_inf(E)) iff no
E-path from x ever reaches running multiplier R = 3^a/2^b < 1 at a halving.

Values are kept exactly as v = N/3^e (3 !| N when e > 0); v is 2-adically even iff N is even.

Certificate logic (proved in collatz_procgen_20260922_exceptional_dimension.md, Lemma S):
  SAFE  : at an integer state v <= -1 with running multiplier R, every future value is an integer <= -1,
          so every future multiplier satisfies R' >= 1/|v|, and R' >= 1/(|v|-1/3) if v is odd.
          Hence R > |v| (or v odd and R > |v|-1/3) certifies that the branch never descends.
  ESCAPE: a value v >= 0 leads to a descent (followed explicitly: Collatz on the positive value).
  DESCENT: a prefix with 3^a < 2^b.
prove(x) returns ('hostile', nodes) iff every branch was certified SAFE (a complete proof), or
('descends', witness), or ('undecided', nodes) at the node limit.  For |x| < 3/2 the unsafe tree is
finite (Lemma S(c)), so the search always terminates.
route(x) searches the excursion route M^a then Collatz-with-few-excursions, used for |x| > 3/2.
"""
from fractions import Fraction as F
import sys

def positive_descent(N, e, a, b, cap=200000):
    """follow Collatz from v=N/3^e >= 0 until 3^a < 2^b; returns (a,b)."""
    for _ in range(cap):
        if N == 0:
            while 3**a >= 2**b: b += 1
            return (a, b)
        if N % 2 == 0:
            N //= 2; b += 1
            if 3**a < 2**b: return (a, b)
        else:
            if e > 0: N = N + 3**(e-1); e -= 1
            else: N = 3*N + 1
            a += 1
    return None

def prove(x, nodelimit=5_000_000):
    x = F(x); assert x < 0 and x.denominator % 2 == 1
    N0, q = x.numerator, x.denominator
    e0 = 0
    while q % 3 == 0: q //= 3; e0 += 1
    assert q == 1, "denominator must be a power of 3"
    stack = [(N0, e0, 0, 0)]
    nodes = 0; maxa = 0
    while stack:
        N, e, a, b = stack.pop()
        nodes += 1
        if nodes > nodelimit: return ('undecided', nodes)
        if b >= 1 and 3**a < 2**b: return ('descends', ('prefix', a, b, str(F(N, 3**e))))
        if N >= 0: return ('descends', ('escape', a, b, str(F(N, 3**e)), positive_descent(N, e, a, b)))
        if e == 0:
            if 3**a > (-N) * 2**b: continue
            if N % 2 != 0 and 3**(a+1) > (-3*N - 1) * 2**b: continue
        maxa = max(maxa, a)
        # x3 move
        if e > 0: stack.append((N + 3**(e-1), e-1, a+1, b))
        else: stack.append((3*N + 1, 0, a+1, b))
        if N % 2 == 0: stack.append((N//2, e, a, b+1))   # halve (explored first)
    return ('hostile', nodes)

def route(x, amax=80, extra=2, steps=400):
    """descent search: M^a (always legal) then Collatz on negatives with <= extra optional x3 moves."""
    x = F(x)
    for a0 in range(1, amax+1):
        v0 = 3**a0 * x + F(3**a0 - 1, 2)
        def dfs(v, a, b, k, t):
            if b >= 1 and 3**a < 2**b: return (a, b)
            if t > steps: return None
            if v >= 0:
                if v.denominator != 1: return None
                return positive_descent(v.numerator, 0, a, b)
            if v.numerator % 2 == 0:
                r = dfs(v/2, a, b+1, k, t+1)
                if r: return r
                if k > 0: return dfs(3*v+1, a+1, b, k-1, t+1)
                return None
            return dfs(3*v+1, a+1, b, k, t+1)
        r = dfs(v0, a0, 0, extra, 0)
        if r: return ('descends', ('route', a0, r))
    return ('no-route-found', amax)

if __name__ == '__main__':
    for s in sys.argv[1:]:
        x = F(s)
        res = prove(x) if abs(x) < F(3, 2) else route(x)
        print(f"{s:>16} = {float(x):.6f}: {res}", flush=True)
