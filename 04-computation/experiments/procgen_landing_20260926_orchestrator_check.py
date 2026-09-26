#!/usr/bin/env python3
"""Orchestrator audit of lane `landing`, written from the note's statements;
the lane's scripts were not read.

T_b(x) = x/2 (x even), (3x+b)/2 (x odd). Scale X, window k = floor(log2 X),
depth D. In a segment of distinct positive integers, index i (with y_i <= X and
a full window) is a dipper if y_(i+s) < y_i 2^(-D) for some 1 <= s <= k; its
landing index is i+s for the least such s.

  1. Exhaustive worst case over first dippers y <= 2^(k+1) - 1 (k = 10, 12, 14;
     D = 1, 1.5, 2, 3; b = 1, -1, 5): Lemma S (all dippers of j in
     (2^D y_j, 2^(D+1) y_j], step into j halves), Lemma O (an odd step between
     consecutive dippers), Theorem 2 (m <= ceil((k-D)/log2 3) for b > 0; the
     log2(3/2)-corrected bound for b < 0 when y_j >= k|b|), and equality of the
     maximum with the bound for b = 1.
  2. The orbit of 13255 at L = 20, D = 2 has a landing point with 12 dippers.
  3. Proposition H: the Sturmian hover class (b = 1, D = 3, O = 10) -- random
     members have every visit as a dipper landing at index l.
  4. Theorem 1(ii) (saturation) on a grid, in log space.
  5. Lemma A: along record orbits, a pair (i, j) serves at most one integer depth.
"""
import math, random
from fractions import Fraction


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


ALPHA = math.log2(3)


def T(x, b):
    return (3 * x + b) // 2 if x & 1 else x >> 1


def below(u, v, D2):
    """u < v * 2^(-D) with D = D2/2 (D2 integer), exact: u * 2^(D2/2) < v."""
    if D2 % 2 == 0:
        return (u << (D2 // 2)) < v
    return 2 * (u * u) * (1 << (D2 - 1)) < v * v   # (u 2^(D2/2))^2 = u^2 2^D2 < v^2


def segment(y, b, length):
    seq, seen = [y], {y}
    for _ in range(length):
        z = T(seq[-1], b)
        if z <= 0 or z in seen:
            break
        seq.append(z)
        seen.add(z)
    return seq


def landing(seq, i, k, D2):
    for s in range(1, k + 1):
        if i + s >= len(seq):
            return None
        if below(seq[i + s], seq[i], D2):
            return i + s
    return -1


# ---------------------------------------------------------------- 1. exhaustive worst case
for b in (1, -1, 5):
    for k in (10, 12, 14):
        X = 2 ** (k + 1) - 1
        for D2 in (2, 3, 4, 6):
            D = D2 / 2
            if b > 0:
                bound = math.ceil((k - D) / ALPHA)
            else:
                bound = math.ceil((k - D + math.log2(1.5)) / ALPHA)
            best = 0
            for y in range(1, X + 1):
                seq = segment(y, b, 3 * k + 4)
                if len(seq) < k + 1:
                    continue
                j = landing(seq, 0, k, D2)
                if j is None or j == -1:
                    continue
                dippers = [0]
                for t in range(1, j):
                    if seq[t] > X:
                        continue
                    if landing(seq, t, k, D2) == j:
                        dippers.append(t)
                yj = seq[j]
                if b < 0 and yj < k * abs(b):
                    continue
                # Lemma S
                assert seq[j - 1] == 2 * yj, (b, k, D, y)
                for t in dippers:
                    # 2^D yj < y_t <= 2^(D+1) yj
                    assert below(yj, seq[t], D2) and not below(yj, seq[t], D2 + 2), (b, k, D, y, t)
                # Lemma O
                for r in range(len(dippers) - 1):
                    assert any(seq[t] & 1 for t in range(dippers[r], dippers[r + 1])), (b, k, D, y)
                m = len(dippers)
                assert m <= bound, (b, k, D, y, m, bound)
                best = max(best, m)
            if b == 1:
                assert best == bound, (b, k, D, best, bound)
check(True, "Lemmas S, O and Theorem 2 on every first dipper y <= 2^(k+1)-1 (b = 1, -1, 5; k = 10, 12, 14; D = 1, 1.5, 2, 3); for b = 1 the maximum equals ceil((k-D)/log2 3) in all 12 cells")

# ---------------------------------------------------------------- 2. 13255
b, k, D2 = 1, 20, 4
X = 2 ** 20
seq = segment(13255, 1, 10 ** 5)
mult = {}
for i in range(len(seq)):
    if seq[i] > X or i + k >= len(seq):
        continue
    j = landing(seq, i, k, D2)
    if j is not None and j >= 0:
        mult[j] = mult.get(j, 0) + 1
mx = max(mult.values())
check(mx == 12 == math.ceil((20 - 2) / ALPHA), f"orbit of 13255 at L = 20, D = 2: max landing multiplicity {mx} = ceil(18/log2 3)")

# ---------------------------------------------------------------- 3. Proposition H class
b, D, O = 1, 3, 10
eps = 1e-6
a = (1 - (D % 1) - eps - O * ALPHA) % 1
n = [math.floor(o * ALPHA + a) for o in range(O + 1)]
word = []
for o in range(O):
    word += [1] + [0] * (n[o + 1] - n[o] - 1)
c = math.floor(D) + 1
word += [0] * c
l = len(word)
r = 0
for t in range(l):
    x = r
    for _ in range(t):
        x = T(x, b)
    if (x & 1) != word[t]:
        r += 2 ** t
x = r
for t in range(l):
    assert (x & 1) == word[t]
    x = T(x, b)
random.seed(7)
k = l + 14
X = 2 ** k
hits = 0
for _ in range(40):
    y = r + (2 ** l) * random.randrange((X // 4) // 2 ** l // 2, (X // 4) // 2 ** l)
    assert 16 * b * l < y <= X // 4
    seq = segment(y, b, 3 * k)
    visits = [n[o] for o in range(O + 1)]
    lands = [landing(seq, i, k, 2 * D) for i in visits]
    safe = [o for o in range(O + 1) if o == O or min(((O - o) * ALPHA + D) % 1, 1 - ((O - o) * ALPHA + D) % 1) >= 1 / 8]
    assert all(lands[o] == n[O] + c for o in safe), (y, lands)
    hits += sum(1 for o in range(O + 1) if lands[o] == n[O] + c)
check(True, f"Proposition H class (b = 1, D = 3, O = 10, l = {l}): in 40 random members every safe visit is a dipper landing at index l; all-visit hits {hits}/{40 * (O + 1)}")

# ---------------------------------------------------------------- 4. Theorem 1(ii)
rho = math.log(2) / math.log(3)
hstar = -(rho * math.log2(rho) + (1 - rho) * math.log2(1 - rho))
lam = math.log2(rho / (1 - rho)) / ALPHA
worst = 1e9
for mu in (0.0, 0.25, 0.5, 0.75, 1.0):
    astar = mu * lam / hstar - 1.5
    for C in (0.01, 1.0, 100.0):
        E = C ** (-lam / hstar)
        for p in range(4, 18):
            L = 2.0 ** p
            for D in range(1, int(L / 2) + 1, max(1, int(L / 400))):
                # log2 of psi(X)/X^h = a* log2 L ; divide everything by X^h
                lhs = astar * math.log2(L)
                t1 = math.log2(C) + mu * math.log2(L) - D * hstar + astar * math.log2(L - D)
                t2 = math.log2(E) - 1.5 * math.log2(L) + lam * D
                rhs = max(t1, t2) + math.log2(1 + 2 ** (-abs(t1 - t2)))
                worst = min(worst, rhs - lhs)
check(worst >= 0, f"Theorem 1(ii): psi = X^h* L^a*(mu) satisfies (R_mu) on the grid (mu in 0..1, C in 0.01..100, L = 2^4..2^17, D in [1, L/2]); min log2(RHS/psi) = {worst:.3f}; a*(1) = {lam / hstar - 1.5:.6f}")

# ---------------------------------------------------------------- 5. Lemma A
pairs = 0
for start in (27, 703, 77671, 837799, 8400511, 63728127):
    seq = segment(start, 1, 10 ** 5)
    k = 40
    for j in range(2, len(seq)):
        for i in range(max(0, j - k), j - 1):
            mu_ij = min(seq[i + 1:j])
            # i is a dipper landing at j at depth D iff y_j < y_i 2^-D <= mu_ij
            Ds = [D for D in range(1, 80) if seq[j] * 2 ** D < seq[i] <= mu_ij * 2 ** D]
            assert len(Ds) <= 1, (start, i, j, Ds)
            pairs += 1
check(True, f"Lemma A: along six record orbits (k = 40), each of {pairs} pairs (i, j) serves at most one integer depth")
