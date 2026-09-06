"""Exact controls for the rational-handoff obstruction; no producer imports.

The infinite theorem is analytic/algebraic, not inferred from this finite bank.
All arithmetic is integral or rational; all checks remain active under -O.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a or [0]


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    return trim(c)


def scale(a, c):
    return trim([x*c for x in a])


def mul(a, b):
    c = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def power(a, n):
    z = [1]
    for _ in range(n):
        z = mul(z, a)
    return z


def reflect(a):
    z = [0]
    for k, c in enumerate(a):
        z = add(z, scale(power([1, -1], k), c))
    return z


def affine(a, slope, intercept):
    z = [0]
    for k, c in enumerate(a):
        z = add(z, scale(power([intercept, slope], k), c))
    return z


def value(a, x):
    z = 0
    for c in reversed(a):
        z = z*x+c
    return z


def row_power(box):
    d = len(box)-1
    z = [0]
    for k, c in enumerate(box):
        z = add(z, [0]*(d-k)+scale(power([1, -1], k), c))
    return z


def modrank(a, prime):
    a = [[x % prime for x in r] for r in a]
    rank = 0
    for j in range(len(a[0])):
        piv = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if piv is None:
            continue
        a[rank], a[piv] = a[piv], a[rank]
        inv = pow(a[rank][j], -1, prime)
        a[rank] = [x*inv % prime for x in a[rank]]
        for i in range(rank+1, len(a)):
            u = a[i][j]
            a[i] = [(x-u*y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def annulus(m):
    """Reconstruct THM-3343 literal lexicographic donor by Hamming counts."""
    M = 2*m
    donor0, donor1 = [0]*m, [0]*m
    for w in range(1, M):
        def C(n, k):
            return comb(n, k) if 0 <= k <= n else 0
        A0, A1 = C(m-1, w-1), C(m-1, w-m)
        defect = sum((1 if (n-m) % 2 else -1) *
                     (C(M-n-1, w-1)+C(M-n-1, w-n))
                     for n in range(m+1, M))
        need(0 <= defect <= A0+A1, 'rotation capacity')
        need((A0+A1-defect) % 2 == 0, 'donor parity')
        take = (A0+A1-defect)//2
        if 0 <= w-1 < m:
            donor0[w-1] = min(take, A0)
        if 0 <= M-1-w < m:
            donor1[M-1-w] = max(0, take-A0)
    rows = {m: (row_power(donor0), row_power(donor1))}
    for n in range(m+1, M):
        rows[n] = ([int((n-m) % 2 == 1)],)*2
    for n, (W, V) in rows.items():
        d = m-1 if n == m else 0
        need(len(W) <= d+1 and len(V) <= d+1, 'deadline degree')
    return rows


def orientation_sum(rows, stop):
    F, G = [0], [0]
    for n in range(1, stop):
        W, V = rows[n]
        F = add(F, [0]*n+mul([1, -1], W))
        G = add(G, [0]*n+mul([1, -1], V))
    return F, G


def main():
    # Exact Bernstein box -> integral power coefficients and parity for all d<=4.
    boxes = 0
    for d in range(5):
        for b in product(*(range(comb(d, k)+1) for k in range(d+1))):
            W = row_power(b)
            Delta = add(scale(W, 2), [-1])
            need(all(isinstance(c, int) for c in W), 'integer power germ')
            need(Delta[0] % 2 == 1 and all(c % 2 == 0 for c in Delta[1:]),
                 'signed parity invariant')
            boxes += 1
    need(boxes == 782, 'complete box count')

    # Denominators with constant one cannot vanish after reflection mod ell.
    reflection_cases = 0
    for prime, degree in ((2, 4), (3, 2), (5, 1)):
        denominators = [[1]+list(t) for t in product(range(prime), repeat=degree)]
        for B in denominators:
            for D in denominators:
                R = mul(B, reflect(D))
                need(any(c % prime for c in R), 'nonzero denominator product')
                reflection_cases += 1
    # Content, rather than constant-term normalization, is sufficient.
    need(any(c % 2 for c in mul([2, 1], reflect([2, 1]))),
         'primitive even-constant denominators survive reflection')
    for prime in (2, 3, 5):
        for mask in product(range(prime), repeat=3):
            if not any(mask):
                continue
            for slope in range(1, prime):
                for intercept in range(prime):
                    need(any(c % prime for c in affine(mask, slope, intercept)),
                         'unit affine substitution preserves nonzero denominator')
    # F_d(t)=t/(1-dt) has integral Taylor coefficients. Nonunit slope 4 permits
    # F_3(4x+1)+F_6(-x)=-1/2, showing the affine prime exclusion is necessary.
    A, B = [1, 4], [-2, -12]
    C, D = [0, -1], [1, 6]
    need(add(scale(add(mul(A, D), mul(C, B)), 2), mul(B, D)) == [0],
         'nonunit affine exact half-constant hostile')

    # A compressed signed boundary can pass every tail/parity test while
    # its first actual row is impossible. K_N=(1-2p)^(2N)-p^N-q^N, N>=2.
    # This does NOT refute the orientation-resolved rational-state theorem.
    for N in range(2, 10):
        K = add(add(power([1, -2], 2*N), scale(power([0, 1], N), -1)),
                scale(power([1, -1], N), -1))
        parity = add(add([1], power([0, 1], N)), power([1, -1], N))
        need(all(c % 2 == 0 for c in add(K, parity)), 'compressed residual parity')
        need(value(K, 0) == value(K, 1) == 0, 'compressed residual endpoints')
        for j in range(1, 16):
            x = Q(j, 16)
            need(abs(value(K, x)) <= x**N+(1-x)**N, 'compressed residual tail capacity')
    K2 = add(add(power([1, -2], 4), scale(power([0, 1], 2), -1)),
             scale(power([1, -1], 2), -1))
    need(value(K2, Q(1, 4)) / (Q(1, 4)*Q(3, 4)) == -3,
         'compressed residual violates first-row signed capacity two')

    # A genuine stationary, admissible slope-two rule: parity of m-1 fresh bits.
    # W_(m+1)=(2p-1)W_m+(1-p), W_1=1; V_m=1-W_m.
    W = [1]
    for m in range(1, 13):
        d = m-1
        raw = [comb(d, k) if k % 2 == 0 else 0 for k in range(d+1)]
        need(W == row_power(raw), 'literal parity-row Bernstein box')
        W = add(mul([-1, 2], W), [1, -1])
    # F=p(1+p)/(1+2p), G=u^2/(1+2u); exact failure 16/35 at p=1/3.
    need(Q(1, 3)*(1+Q(1, 3))/(1+2*Q(1, 3)) +
         Q(2, 3)**2/(1+2*Q(2, 3)) == Q(16, 35), 'stationary fair-bias hostile')
    # Cleared exact-fairness numerator; mod 2 equals B(p)B(1-p).
    B = [1, 2]
    numerator = add(add(scale(mul([0, 1, 1], reflect(B)), 2),
                        scale(mul(reflect([0, 0, 1]), B), 2)),
                    scale(mul(B, reflect(B)), -1))
    need(numerator != [0] and any(x % 2 for x in numerator), 'matrix parity witness')

    # von Neumann is exactly fair and rational, but not critical-deadline bounded.
    Fvn, Gvn = [0, 1, Q(-1, 2)], [0, 0, Q(1, 2)]
    need(add(Fvn, reflect(Gvn)) == [Q(1, 2)], 'von Neumann rational fairness')
    need(any(isinstance(c, Q) and c.denominator == 2 for c in Fvn+Gvn),
         'von Neumann fails integer-germ premise')
    for k in range(12):
        word = '0011' + '00'*k + '01'
        pairs = [word[j:j+2] for j in range(0, len(word), 2)]
        need(next(j for j in range(1, len(word)) if word[j] != word[0]) == 2,
             'fixed critical value two')
        need(all(v in ('00', '11') for v in pairs[:-1]) and pairs[-1] == '01',
             'unbounded post-disagreement wait')

    # Reconstruct inherited shell-spill control, without its implementation.
    rows = {
        1: ([1], [0]), 2: ([0], [1]), 3: ([0, 1], [1, -1]),
        4: ([1], [0]), 5: ([1, -1], [0, 1]), 6: ([0], [1]), 7: ([1], [0])}
    F, G = orientation_sum(rows, 8)
    target = scale(add(add([1], scale(power([0, 1], 8), -1)),
                       scale(power([1, -1], 8), -1)), Q(1, 2))
    need(add(F, reflect(G)) == target, 'cross-shell exact fairness')
    shell_vectors = {2: [0]*9, 4: [0]*9}
    for digits in product((0, 1), repeat=8):
        if len(set(digits)) == 1:
            continue
        n = next(j for j in range(1, 8) if digits[j] != digits[0])
        head = (digits[0] == 0 if n in (1, 4, 7) else
                digits[0] == 1 if n in (2, 6) else
                digits[4] == 0 if n == 3 else digits[6] == 1)
        if n >= 2:
            shell_vectors[2 if n < 4 else 4][sum(digits)] += 2*int(head)-1
    need(shell_vectors[2] == [0, 0, -2, -4, 0, 4, 2, 0, 0], 'inherited donor residual')
    need(shell_vectors[4] == [-c for c in shell_vectors[2]], 'exact residual handoff')
    need(rows[3][0] == [0, 1], 'interior row three spills beyond time four')

    # Positive infinite architecture, inherited THM-3343: reconstruct five annuli.
    rows = {}
    for m in (1, 2, 4, 8, 16, 32):
        rows.update(annulus(m))
        F, G = orientation_sum(rows, 2*m)
        target = scale(add(add([1], scale(power([0, 1], 2*m), -1)),
                           scale(power([1, -1], 2*m), -1)), Q(1, 2))
        need(add(F, reflect(G)) == target, 'complete dyadic prefix bisection')
    F, G = orientation_sum(rows, 64)
    ranks = []
    for size in (2, 4, 8, 12, 16, 24, 30):
        ranks.append([size] + [modrank([[a[i+j+1] if i+j+1 < len(a) else 0
                                        for j in range(size)] for i in range(size)], 101)
                               for a in (F, G)])
    need(ranks[-1][1] >= 20 and ranks[-1][2] >= 20, 'positive nonstationary rank control')
    print('Rational handoff obstruction: PROOF CANDIDATE; exact companion PASS')
    print('Complete Bernstein boxes:', boxes)
    print('Prime-reflected denominator controls:', reflection_cases)
    print('Affine sharpness: F_3(4x+1)+F_6(-x)=-1/2; bad prime 2 divides slope 4')
    print('Aggregate-state hostile: tail/parity valid, first-row signed load -3 < -2')
    print('Stationary admissible slope-two rule: H(1/3)=16/35 != 1/2')
    print('von Neumann: rational exact fairness, half-integral germs, T(2) unbounded')
    print('Inherited THM-3337 opposite shell vectors:', shell_vectors[2], shell_vectors[4])
    print('Inherited THM-3343 prefix checks through horizon 64; Hankel ranks mod 101:')
    print(json.dumps(ranks, separators=(',', ':')))
    print('No uniform slope below two is proved; fixed polynomial-state architecture is excluded.')
    print('PASS', GATES, 'active exact gates')


if __name__ == '__main__':
    main()
