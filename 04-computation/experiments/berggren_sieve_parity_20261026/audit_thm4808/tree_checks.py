#!/usr/bin/env python3
"""Independent audit of THM-4808 (a)-(f) on the integer tree (own code, written from the theorem text).

Word convention: w = X_1 X_2 ... X_k applied root first, (m,n) = X_k ... X_1 (2,1), M_w = X_k ... X_1.
Descent: m > 3n -> D^-1, 2n < m < 3n -> A^-1, n < m < 2n -> U^-1 (run-length compressed).
Checks
 (a)  chi = alpha n - gamma m with (alpha, gamma) = M_w (1,0); the box 0 <= alpha <= m, 0 <= gamma <= n, and the
      stronger CONE LEMMA (proof in REPORT.md): alpha >= gamma >= 0 and m - alpha >= n - gamma >= 0.
 (b)  the THM's continued-fraction digit automaton equals #A; and the audit's closed forms
      chi = (-1)^(L + [X_1 = U]), L + 1 = number of regular CF partial quotients of m/n (last digit >= 2), and
      chi = (-1)^(m + c_0 + ... + c_L)  (two-line proof in REPORT.md; uses THM-4057 Theorem 4.1).
 (c)  (16M+2, 4M+1): descent word A D^(M-1) U D, residues, triple == (3,4,5) mod M, chi = -1, M_w == R mod M.
 (d)  Q-recursions (symbolic) and the trailing-A-run sign law; A-ray Pell sign.
 (e)  level identity (brute force over all words, k <= 7, symbolic).
 (f)  node partition (root / A-child / interior of exactly one ray), the ray pairing, interleaving, the exact
      identity S(X) = 1 + N_U + N_D + E(X) with E from the local terms, delta in {0,1}.
"""
import sys, itertools
from math import gcd, isqrt
import sympy as sp

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
XS = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]


def descent(m, n):
    """run-length word, root first: list of (letter, run).  Also returns M_w as (a, b, c, d)."""
    runs = []
    P = (1, 0, 0, 1)
    while (m, n) != (2, 1):
        if m > 3 * n:
            # largest r with m - 2 r n > 2 n (stay > 2n so that the next region is decided) but at least 1
            r = max(1, (m - 2 * n - 1) // (2 * n))
            while m - 2 * r * n < n:
                r -= 1
            m -= 2 * r * n
            runs.append(('D', r)); X = (1, 2 * r, 0, 1)
        elif m > 2 * n:
            m, n = n, m - 2 * n
            runs.append(('A', 1)); X = (2, 1, 1, 0)
        else:
            d = m - n
            r = max(1, -(-(2 * n - m) // d))     # ceil((2n - m)/d)
            r = min(r, n - 1) if d == 1 else r
            m, n = m - r * d, n - r * d
            runs.append(('U', r)); X = (1 + r, -r, r, 1 - r)
        a, b, c, d_ = P; x1, x2, x3, x4 = X
        P = (a * x1 + b * x3, a * x2 + b * x4, c * x1 + d_ * x3, c * x2 + d_ * x4)
    runs.reverse()
    # merge adjacent equal letters
    merged = []
    for L, r in runs:
        if merged and merged[-1][0] == L and L != 'A':
            merged[-1] = (L, merged[-1][1] + r)
        else:
            merged.append((L, r))
    return merged, P


def apply_word(runs, v):
    x, y = v
    for L, r in runs:
        for _ in range(r):
            if L == 'U': x, y = 2 * x - y, x
            elif L == 'A': x, y = 2 * x + y, x
            else: x, y = x + 2 * y, y
    return x, y


def cf(m, n):
    d = []
    while n:
        d.append(m // n); m, n = n, m % n
    return d


def automaton(m, n):
    """THM-4808 (b) exactly as stated: returns the number of counted even heads."""
    d = cf(m, n)
    count = 0
    while True:
        h, rest = d[0], d[1:]
        if not rest:                       # integer slope: must be even; stop
            assert h % 2 == 0
            return count
        if h % 2 == 0:                     # even head with digits remaining
            count += 1; d = rest
        else:
            if len(rest) == 1:             # [odd; c1]: stop
                return count
            d = [rest[1] + 1] + rest[2:]   # odd head, >= 2 further digits: next head c2 + 1


def main():
    tot = ok_a = ok_box = ok_cone = ok_b = ok_cf = ok_d = ok_sum = ok_par = 0
    for m in range(2, MMAX + 1):
        for n in range(1, m):
            if (m - n) % 2 == 0 or gcd(m, n) != 1:
                continue
            tot += 1
            runs, P = descent(m, n)
            nA = sum(r for L, r in runs if L == 'A')
            chi = -1 if nA % 2 else 1
            a_, b_, c_, d_ = P
            assert (2 * a_ + b_, 2 * c_ + d_) == (m, n)
            al, ga = a_, c_
            ok_a += (al * n - ga * m == chi) and (a_ * d_ - b_ * c_ == chi)
            ok_box += (0 <= al <= m and 0 <= ga <= n)
            ok_cone += (al >= ga >= 0 and m - al >= n - ga >= 0)
            ok_b += automaton(m, n) == nA
            digits = cf(m, n)
            L = len(digits) - 1
            first_U = bool(runs) and runs[0][0] == 'U'
            ok_cf += chi == (-1) ** (L + first_U)
            # closed form (audit): #A = m + sum of partial quotients (mod 2), from
            #   m = #U + #A (mod 2)  [U, A swap the parities of (m, n); D keeps them; root m = 2]
            #   sum c_i = 2 + #U + 2#A + 2#D  [Stern-Brocot depth cocycle, THM-4057 Thm 4.1 in THM-4808 letters]
            ok_sum += chi == (-1) ** (m + sum(digits))
            nU = sum(r for Lr, r in runs if Lr == 'U'); nD = sum(r for Lr, r in runs if Lr == 'D')
            ok_par += (m - nU - nA) % 2 == 0 and sum(digits) == 2 + nU + 2 * nA + 2 * nD
            # (d) trailing A-run sign law
            r = 0; k = len(runs) - 1
            while k >= 0 and runs[k][0] == 'A':
                r += 1; k -= 1
            base = 1 if (k >= 0 and runs[k][0] == 'D') else -1
            A, B = m * m - n * n, 2 * m * n
            ok_d += (1 if A > B else -1) == (-1) ** r * base
    print(f"Euclid pairs with m <= {MMAX}: {tot}")
    print(f"(a) chi = alpha n - gamma m = det M_w: {ok_a}/{tot}; box 0<=alpha<=m, 0<=gamma<=n: {ok_box}/{tot}; "
          f"cone lemma alpha>=gamma>=0, m-alpha>=n-gamma>=0: {ok_cone}/{tot}")
    print(f"(b) THM automaton = #A: {ok_b}/{tot}; closed form chi = (-1)^(L + [first letter U]): {ok_cf}/{tot}")
    print(f"(b'') CLOSED FORM chi = (-1)^(m + c_0 + c_1 + ... + c_L): {ok_sum}/{tot}; "
          f"lemmas m = #U + #A (mod 2) and sum c_i = 2 + #U + 2#A + 2#D: {ok_par}/{tot}")
    print(f"(d) sign(a-b) = (-1)^(trailing A-run) * (+1 after D, -1 after U or empty): {ok_d}/{tot}")
    # (d) symbolic recursions
    m_, n_ = sp.symbols('m n')
    Q = lambda x, y: sp.expand(x ** 2 - y ** 2 - 2 * x * y)
    P_ = lambda x, y: sp.expand(x ** 2 - y ** 2 + 2 * x * y)
    rec = (sp.expand(Q(2 * m_ + n_, m_) + Q(m_, n_)) == 0, sp.expand(Q(m_ + 2 * n_, n_) - P_(m_, n_)) == 0,
           sp.expand(Q(2 * m_ - n_, m_) + P_(m_, n_)) == 0, Q(2, 1) == -1,
           sp.expand(Q(m_, n_) - ((m_ - n_) ** 2 - 2 * n_ ** 2)) == 0)
    print(f"(d) Q(A.) = -Q, Q(D.) = P, Q(U.) = -P, Q(2,1) = -1, Q = (m-n)^2 - 2n^2: {rec}")
    x = (2, 1); pell = True
    for k in range(40):
        A, B = x[0] ** 2 - x[1] ** 2, 2 * x[0] * x[1]
        pell &= abs(A - B) == 1 and (-1) ** k == -(A - B)
        x = (2 * x[0] + x[1], x[0])
    print(f"(d) A-ray k = 0..39: |a-b| = 1 and chi = -(a-b): {pell}")
    # (c)
    badc = 0
    Rm = sp.Matrix([[-1, 4], [0, 1]])
    for M in range(1, 5001):
        m, n = 16 * M + 2, 4 * M + 1
        runs, P = descent(m, n)
        word = ''.join(L * r for L, r in runs)
        good = word == 'A' + 'D' * (M - 1) + 'UD'
        good &= (m % M, n % M) == (2 % M, 1 % M)
        a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
        good &= (a % M, b % M, c % M) == (3 % M, 4 % M, 5 % M)
        good &= P[0] * P[3] - P[1] * P[2] == -1
        good &= all((P[i] - (-1, 4, 0, 1)[i]) % M == 0 for i in range(4))
        good &= apply_word(runs, (2, 1)) == (m, n) if M <= 300 else True
        badc += not good
    Dm = sp.Matrix([[1, 2], [0, 1]]); Um = sp.Matrix([[2, -1], [1, 0]]); Am = sp.Matrix([[2, 1], [1, 0]])
    Rcomp = Dm * Um * Dm.inv() * Am
    print(f"(c) M = 1..5000 failures: {badc}; R = D U D^-1 A = {Rcomp.tolist()}, R(2,1) = {list(Rcomp * sp.Matrix([2, 1]))}, "
          f"det R = {Rcomp.det()}")
    # (e)
    X_, Y_ = sp.symbols('x y')
    e_ok = True
    for k in range(0, 8):
        tot_poly = 0
        for w in itertools.product('UAD', repeat=k):
            tot_poly += (-1) ** w.count('A') * X_ ** w.count('D') * Y_ ** w.count('U')
        e_ok &= sp.expand(tot_poly - (X_ + Y_ - 1) ** k) == 0
    print(f"(e) sum chi x^#D y^#U over |w| = k equals (x+y-1)^k for k = 0..7: {e_ok}")
    # (f)
    for X in XS:
        S = 0; E = 0; bad_delta = 0; bad_inter = 0; nodes = 0
        cls = {'root': 0, 'A-child': 0, 'Dray-interior': 0, 'Uray-interior': 0}
        stack = [((2, 1), 1, '')]
        while stack:
            (m, n), s, last = stack.pop()
            if m * m + n * n > X:
                continue
            nodes += 1; S += s
            cls['root' if last == '' else {'A': 'A-child', 'D': 'Dray-interior', 'U': 'Uray-interior'}[last]] += 1
            wU, wA, wD = (2 * m - n, m), (2 * m + n, m), (m + 2 * n, n)
            # D-rays of wU and wA (j >= 1), U-rays of wD and wA (j >= 1)
            def cnt(v, step):
                j = 0
                while (v[0] + (j + 1) * step[0]) ** 2 + (v[1] + (j + 1) * step[1]) ** 2 <= X:
                    j += 1
                return j
            dD = cnt(wU, (2 * m, 0)) - cnt(wA, (2 * m, 0))
            dU = cnt(wD, (m + n, m + n)) - cnt(wA, (m + n, m + n))
            bad_delta += dD not in (0, 1) or dU not in (0, 1)
            if X == XS[0]:
                for j in range(0, 6):
                    hU = (2 * m - n + 2 * j * m) ** 2 + m * m; hA = (2 * m + n + 2 * j * m) ** 2 + m * m
                    hU1 = (2 * m - n + 2 * (j + 1) * m) ** 2 + m * m
                    gD = (m + 2 * n + j * (m + n)) ** 2 + (n + j * (m + n)) ** 2
                    gA = (2 * m + n + j * (m + n)) ** 2 + (m + j * (m + n)) ** 2
                    gD1 = (m + 2 * n + (j + 1) * (m + n)) ** 2 + (n + (j + 1) * (m + n)) ** 2
                    bad_inter += not (hU < hA < hU1 and gD < gA < gD1)
            E += s * (dD + dU - (1 if wA[0] ** 2 + wA[1] ** 2 <= X else 0))
            stack += [(wU, s, 'U'), (wA, -s, 'A'), (wD, s, 'D')]
        NU = sum(1 for k in range(1, isqrt(X) + 2) if (k + 2) ** 2 + (k + 1) ** 2 <= X)
        ND = sum(1 for k in range(1, isqrt(X) + 2) if (2 * k + 2) ** 2 + 1 <= X)
        # root-ray interiors must equal N_U + N_D: count nodes with n = 1 (D ray) or n = m - 1 (U ray), minus root
        print(f"(f) X = {X}: nodes {nodes}, S = {S}, 1 + N_U + N_D + E = {1 + NU + ND + E} (N_U = {NU}, N_D = {ND}, "
              f"E = {E}), identity {S == 1 + NU + ND + E}; delta outside {{0,1}}: {bad_delta}; "
              f"interleaving violations (j<6, X=1e3 nodes): {bad_inter}; classes {cls}")
    worst = 0.0
    from math import sqrt
    for X in [10 ** k for k in range(2, 25)] + [5 * 10 ** k for k in range(2, 24)]:
        NU = max(0, (isqrt(2 * X - 1) - 3) // 2)
        while 2 * (NU + 1) ** 2 + 6 * (NU + 1) + 5 <= X: NU += 1
        while NU > 0 and 2 * NU ** 2 + 6 * NU + 5 > X: NU -= 1
        ND = max(0, (isqrt(X - 1) - 2) // 2)
        while (2 * (ND + 1) + 2) ** 2 + 1 <= X: ND += 1
        while ND > 0 and (2 * ND + 2) ** 2 + 1 > X: ND -= 1
        worst = max(worst, abs(NU + ND - (0.5 + 1 / sqrt(2)) * sqrt(X)))
    print(f"(f) max |N_U + N_D - (1/2 + 1/sqrt2) sqrt X| over X = 10^2..10^24 (1-5 grid): {worst:.3f}")


if __name__ == '__main__':
    main()
