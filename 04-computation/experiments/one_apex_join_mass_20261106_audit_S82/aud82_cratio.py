# aud82_cratio.py (independent audit of S82): C(B)/H(B) for all tournaments of order <= 6, with
# C(B) = sum over A subset V of H(B[A]) * H(B[V - A]), H(empty) = 1 (labelled two-bin covers, Session197 B_B(2)).
# Also identifies strong tournaments (via reachability) and reports the non-power-of-2 ratios by order and strongness.
import sys
from fractions import Fraction
from itertools import combinations

def ham_counts(n, out):
    end = [[0]*n for _ in range(1 << n)]
    for v in range(n): end[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = end[S][v]
            if not c: continue
            nb = out[v] & ~S
            while nb:
                w = (nb & -nb).bit_length() - 1; nb &= nb - 1
                end[S | (1 << w)][w] += c
    return [1 if S == 0 else sum(end[S]) for S in range(1 << n)]

def strong(n, out):
    for s in range(n):
        seen = 1 << s; stack = [s]
        while stack:
            v = stack.pop(); nb = out[v] & ~seen
            while nb:
                w = (nb & -nb).bit_length() - 1; nb &= nb - 1; seen |= 1 << w; stack.append(w)
        if seen != (1 << n) - 1: return False
    return True

def is_pow2(fr):
    if fr.denominator != 1: return False
    x = fr.numerator
    return x > 0 and x & (x - 1) == 0

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for n in range(1, NMAX + 1):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ratios = {}
    for bits in range(1 << len(pairs)):
        out = [0]*n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        hs = ham_counts(n, out); full = (1 << n) - 1
        H = hs[full]; C = sum(hs[A] * hs[full & ~A] for A in range(1 << n))
        r = Fraction(C, H); st = strong(n, out)
        ratios.setdefault((r, H, C), set()).add(st)
    allr = sorted(set(k[0] for k in ratios))
    nonp = sorted(set(k for k in ratios if not is_pow2(k[0])))
    print(f"order {n}: {len(allr)} distinct C/H ratios; non-power-of-2 ratios (C/H, H, C, strong?):",
          [(str(k[0]), k[1], k[2], sorted(ratios[k])) for k in nonp][:40])
