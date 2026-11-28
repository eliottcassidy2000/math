"""Exact witnesses, two residual average intervals, and four-rail cube test."""
from fractions import Fraction
from math import isqrt, lcm
from pathlib import Path
import json

N = 200
witnesses = json.loads(Path(__file__).with_name("witnesses.json").read_text())
witnesses = {int(h): (word, flip) for h, (word, flip) in witnesses.items()}

def tournament(word, flip):
    arcs, root = [0], 0
    for i, move in enumerate(word, 1):
        assert move in "PKM"
        arcs.append(0 if move == "P" else 1 << root)
        for x in range(i):
            if x != root or move == "P":
                arcs[x] |= 1 << i
        if move == "M":
            root = i
    if flip:
        i, j = flip
        assert 0 <= i < j < len(arcs)
        arcs[i] ^= 1 << j
        arcs[j] ^= 1 << i
    return arcs

def count(arcs):
    n = len(arcs)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for i, v in enumerate(dp[mask]):
            if not v:
                continue
            rest = arcs[i] & ~mask
            while rest:
                bit = rest & -rest
                rest ^= bit
                dp[mask | bit][bit.bit_length() - 1] += v
    return sum(dp[-1])

def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))

def in_join_closure(h, known):
    if h in known or h == 1:
        return True
    return any(h % d == 0 and in_join_closure(h // d, known)
               for d in known if 1 < d < h)

def filtration(values):
    L = 1
    M = Fraction(1)
    for t in range(1, N + 1):
        if t in values:
            L = lcm(L, t)
        M += Fraction(1, L)
    return L, M

def common_decimal(lo, hi):
    for k in range(1, 1000):
        scale = 10**k
        a, b = lo * scale // 1, hi * scale // 1
        if a != b:
            break
        result = str(a)
    return result[0] + "." + result[1:]

def cube_rails():
    slopes = (495, 567)
    offsets = (303, 305)
    ratios = []
    for a in slopes:
        for b in slopes:
            L = b
            while L <= a**3:
                L *= 2
            ratios.append(Fraction(L, a**3))
    assert min(ratios) > (1 + Fraction(306, 495 * 2**5))**3
    for a in slopes:
        for b in slopes:
            assert a % 2 and b % 2 and a**3 != b
    for k in range(1, 5):
        for a in slopes:
            for u in offsets:
                h = a * 2**k + u
                j = 1
                while 495 * 2**j + 303 < (h + 1)**3:
                    assert all(not h**3 <= b * 2**j + v < (h + 1)**3
                               for b in slopes for v in offsets)
                    j += 1
    print("No cubic shell contains a second count from the four rails.")

if __name__ == "__main__":
    needed = ({p for p in range(3, N+1, 2) if prime(p) and p != 7}
              | {7*p for p in range(5, N+1, 2) if prime(p) and p != 7}
              | {5, 9, 35, 49, 63, 343})
    assert set(witnesses) == needed
    for h, (word, flip) in witnesses.items():
        assert count(tournament(word, flip)) == h
    print(len(witnesses), "literal tournament witnesses checked")
    for i in range(2):
        target = {d for d in range(1, N+1, 2)
                  if (d not in (7,21) if i == 0 else d not in (1,3) and d % 7)}
        witnessed = {d for d in target if in_join_closure(d * (7 if i else 1), needed)}
        L, M = filtration(target)
        assert filtration(witnessed) == (L, M)
        result = common_decimal(M, M + Fraction(15 * N, L))
        print(f"D{i}: {result}... ({len(result)-2} certified places)")
    cube_rails()
