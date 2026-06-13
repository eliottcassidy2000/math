#!/usr/bin/env python3
"""THM-491 / HYP-2436/2437 (claudebox-2026-06-11-S6): LRC n=14 vs n=19 through the recursive lens.
- the shell ramification tower 27->9->3 and its Pisano signature pi(p^k)=p^{k-1}pi(p)
- the freed-clock formula (exact Q) for AP perturbations
- n=19 cyclotomic transversal core; the Goldbach-comet singular series local-global contrast
Pure python, exact. M(S)=max_t min_i ||v_i t||; canon normalization, LRC: M>=1/n."""
from fractions import Fraction as F
from math import gcd

def norm(x):
    r = x - int(x)
    if r < 0: r += 1
    return min(r, 1 - r)
def candidate_denoms(S):
    D = set()
    for v in S:
        D.add(abs(v))
        for w in S:
            if v != w: D.add(abs(v - w)); D.add(abs(v + w))
    D.discard(0); return D
def M_exact(S):
    best, bt = F(0), F(0)
    for d in candidate_denoms(S):
        for k in range(1, d):
            t = F(k, d); m = min(norm(v * t) for v in S)
            if m > best: best, bt = m, t
    return best, bt
def pisano(m):
    if m == 1: return 1
    a, b = 0, 1
    for i in range(6 * m * m):
        a, b = b, (a + b) % m
        if a == 0 and b == 1: return i + 1
def factor(m):
    f, d, fs = m, 2, {}
    while d * d <= f:
        while f % d == 0: fs[d] = fs.get(d, 0) + 1; f //= d
        d += 1
    if f > 1: fs[f] = fs.get(f, 0) + 1
    return fs

def part1_tower():
    print("== Part 1: the ramification tower + Pisano signature ==")
    print("  27 ->(/3) 9 ->(/3) 3 : non-unit (p|v) core of shell p^k = shell p^{k-1}")
    assert sorted(set((r // 3) % 9 for r in range(27) if r % 3 == 0 and r)) == list(range(1, 9))
    assert pisano(9) == 3 * pisano(3) and pisano(27) == 3 * pisano(9)
    print(f"  Pisano tower pi(3,9,27) = {pisano(3)},{pisano(9)},{pisano(27)} = 8,24,72 (p^{{k-1}}*pi(p))")
    ram = [(n, 2*n-1, *list(factor(2*n-1).items())[0]) for n in range(4, 201)
           if len(factor(2*n-1)) == 1 and list(factor(2*n-1).values())[0] >= 2]
    print(f"  ramified n<=200 (depth=v_p(2n-1)): " +
          ", ".join(f"{n}({p}^{e})" for n, m, p, e in ram))
    assert (14, 27, 3, 3) in ram and min(n for n,m,p,e in ram if e>=3) == 14
    print("  => n=14 (27=3^3) is the SMALLEST depth-3 shell")

def part2_freedclock():
    print("\n== Part 2: the freed-clock formula (exact Q) ==")
    for n in (14, 19):
        AP = list(range(1, n)); freed = {}
        for j in range(1, n):
            S = [v for v in AP if v != j] + [n]
            M, t = M_exact(S); freed[j] = M
        clean = [j for j in range(1, n) if freed[j] == F(1, j)]
        blocked = [j for j in range(2, n-1) if n % j == 0 and freed[j] != F(1, j)]  # proper divisors >1
        print(f"  n={n}: drop j frees to 1/j for j in {clean}; divisor-blocked j (j|n): {blocked or 'none (n prime)'}")
        assert freed[n-1] == F(1, n-1)
        if n == 14: assert blocked == [2, 7]
        if n == 19: assert blocked == []  # 19 prime: no proper-divisor clock blocked

def part3_n19_goldbach():
    print("\n== Part 3: n=19 cyclotomic core + Goldbach-comet local-global ==")
    m = 37
    def transversal(S):
        seen = set()
        for s in S:
            r = s % m
            if gcd(r, m) != 1: return False
            k = min(r, (m - r) % m)
            if k in seen: return False
            seen.add(k)
        return len(seen) == (m - 1) // 2
    import random
    random.seed(1)
    c = sum(transversal(random.sample(range(1, m), 18)) for _ in range(4000))
    print(f"  n=19: transversals among 4000 random 18-subsets mod 37 = {c} (THM-420: non-transversal=>loose)")
    assert M_exact(list(range(1, 19)))[0] == F(1, 19)
    assert M_exact(list(range(1, 18)) + [19])[0] == F(1, 18)
    print("  AP {1..18} tight 1/19; skip-18 frees to 1/18 (exact)")
    def rad(n):
        s, d = set(), 2
        while d * d <= n:
            while n % d == 0: s.add(d); n //= d
            d += 1
        if n > 1: s.add(n)
        return s
    sing = lambda n: F(1).__class__(1) if False else __import__('functools').reduce(
        lambda a, p: a * F(p-1, p-2) if p > 2 else a, rad(n), F(1))
    print(f"  Goldbach singular factor: n=14 -> {sing(14)} (7-wing 6/5); n=19 -> {sing(19)} (18/17);")
    print(f"    n=27=3^3 -> {sing(27)} = SAME as n=3 ({sing(3)}): comet reads RADICAL, blind to the tower")
    assert sing(27) == sing(3) == sing(9)

if __name__ == '__main__':
    part1_tower(); part2_freedclock(); part3_n19_goldbach()
    print("\nall checks passed")
