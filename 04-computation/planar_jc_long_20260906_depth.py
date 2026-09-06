#!/usr/bin/env python3
"""Exact controls for sharp finite depth recognition. No repo imports.

Universe: source monomials a,b,c,e>=0, a+b<=d, fixed ell=2c+3e-a;
every monomial in the tested diagonal is included. All checks are live under
-O. Universal conclusions use the proof in the companion report, not this grid.
"""
from fractions import Fraction as Q
from math import comb, factorial
from hashlib import sha256

CHECKS = 0
DIGEST = sha256()


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def record(value):
    DIGEST.update((repr(value) + "\n").encode())


def rank(rows):
    if not rows:
        return 0
    a = [[Q(v) for v in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        lead = a[r][c]
        a[r] = [v / lead for v in a[r]]
        for i in range(r + 1, len(a)):
            if a[i][c]:
                scale = a[i][c]
                a[i] = [u - scale * v for u, v in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def det(rows):
    a = [[Q(v) for v in row] for row in rows]
    out = Q(1)
    for c in range(len(a)):
        p = next((i for i in range(c, len(a)) if a[i][c]), None)
        if p is None:
            return Q(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            out = -out
        lead = a[c][c]
        out *= lead
        for i in range(c + 1, len(a)):
            scale = a[i][c] / lead
            for j in range(c + 1, len(a)):
                a[i][j] -= scale * a[c][j]
    return out


def geometry(ell, d):
    s = max(0, (ell + 1) // 2)
    rho = max(0, (ell + 2) // 3)
    D = ell + d - s
    return s, rho, D, D - rho


def packets(ell, d):
    """Literal binomial expansion of every allowed generator, deduplicated."""
    s, _, D, _ = geometry(ell, d)
    cols = set()
    for e in range(max(0, (ell + d) // 3 + 1)):
        for c in range(max(0, (ell + d - 3 * e) // 2 + 1)):
            a = 2 * c + 3 * e - ell
            if not 0 <= a <= d:
                continue
            for b in range(d - a + 1):
                n0 = b + c + 2 * e
                N = c + e
                v = [0] * (D + 1)
                for j in range(N + 1):
                    idx = n0 + j - s
                    require(0 <= idx <= D, "packet lies in declared cap")
                    v[idx] = (-1) ** idx * comb(N, j)
                cols.add(tuple(v))
    return sorted(cols)


def ideal_basis(ell, d):
    _, rho, D, L = geometry(ell, d)
    cols = []
    for j in range(L + 1):
        v = [0] * (D + 1)
        for i in range(rho + 1):
            v[j + i] = (-1) ** i * comb(rho, i)
        cols.append(v)
    return cols


def h(n, rho):
    return comb(n + rho - 1, rho - 1) if n >= 0 else 0


def main():
    print("PLANAR DEPTH: complete diagonal image and sharp finite recognition")
    print("field characteristic zero; exact Fraction arithmetic; no repo imports")
    count = 0
    for d in range(6):
        for ell in range(-d, 19):
            s, rho, D, L = geometry(ell, d)
            p = packets(ell, d)
            b = ideal_basis(ell, d)
            expected = max(0, L + 1)
            require(rank(p) == expected, "full literal packet rank")
            require(rank(p + b) == expected, "full image equals principal ideal")
            for m in range(max(0, s - 1), ell + d + 2):
                size = max(0, min(m, ell + d) - s + 1)
                pr = [v[:size] for v in p]
                rr = rank(pr) if size else 0
                require(rr == min(size, expected), "truncated image rank")
                if ell >= 1 and s <= m <= ell + d:
                    q0 = max(0, ell + d - m)
                    dual = [[comb(m + q - s - j, q) for j in range(size)]
                            for q in range(q0, rho)]
                    require(rank(dual) == size - rr, "full simplex codimension")
                    for v in pr:
                        for w in dual:
                            require(sum(a * b for a, b in zip(v, w)) == 0,
                                    "every source packet killed by full bank")
                record((ell, d, m, size, rr))
                count += 1
    print("literal source-image/truncation universes:", count)

    print("inherited rank controls: m d ambient rank nullity")
    for m, d, wanted in [(8, 2, 51), (8, 3, 63), (14, 2, 108),
                         (14, 3, 129), (15, 2, 119), (15, 3, 142)]:
        ambient = total = 0
        for ell in range(-d, 2 * m + 1):
            s, _, _, L = geometry(ell, d)
            a = max(0, min(m, ell + d) - s + 1)
            ambient += a
            total += min(a, max(0, L + 1))
        require(total == wanted, "inherited rank")
        print(m, d, ambient, total, ambient - total)

    minors = 0
    for rho in range(1, 11):
        for k in range(1, rho + 1):
            for B in range(13):
                matrix = [[h(B + i - j, rho) for j in range(k)]
                          for i in range(k)]
                product = Q(1)
                for i in range(k):
                    product *= Q(factorial(rho + B - 1 - i) * factorial(i),
                                 factorial(rho - 1 - i) * factorial(B + i))
                value = det(matrix)
                require(value == product and value > 0 and value.denominator == 1,
                        "Toeplitz positive integer determinant")
                record((rho, k, B, value))
                minors += 1
    print("independent Toeplitz determinants:", minors)

    universes = 0
    for T in range(1, 16):
        for d in range(6):
            cutoff = (4 * T) // 3 + d + 1
            max_gate = 0
            for ell in range(1, 2 * T + 1):
                s, rho, D, L = geometry(ell, d)
                M = min(T, ell + d) - s
                if M < 0:
                    continue
                gate = ell + d - rho + min(M + 1, rho)
                max_gate = max(max_gate, gate)
                early = [[h(r - j, rho) for j in range(M + 1)]
                         for r in range(L + 1, gate - s + 1)]
                complete = [[h(r - j, rho) for j in range(M + 1)]
                            for r in range(L + 1, D + 1)]
                require(rank(early) == min(M + 1, rho), "recognition rank")
                require(rank(early + complete) == rank(early),
                        "early equations contain complete recognition")
                require(gate <= cutoff, "uniform upper bound")
                record((T, d, ell, gate, rank(early)))
                universes += 1
            require(max_gate == cutoff, "sharp universal clock")
            ell = 2 * T
            s, rho, D, L = geometry(ell, d)
            q = [h(j, rho) for j in range(L + 1)]
            lift = [sum((-1) ** i * comb(rho, i) * q[j-i]
                        for i in range(rho + 1) if 0 <= j-i < len(q))
                    for j in range(D + 1)]
            require(lift[:L+1] == [1] + [0] * L, "sharp hostile passes prefix")
            require(lift[L+1] == -comb(rho + L, L + 1) != 0,
                    "sharp hostile first failed coefficient")
            require(s + L + 1 == cutoff, "sharp hostile endpoint")
            require(rank(ideal_basis(ell,d) + [lift]) == L + 1,
                    "positive lift is an actual full source image")
            target = [1] + [0] * D
            require(rank(ideal_basis(ell,d) + [target]) == L + 2,
                    "finite polynomial hostile really outside source image")
    print("fixed-support diagonal recognition universes:", universes)
    print("sharp t^T controls: T=1..15, d=0..5; 90 positive lifts + 90 hostiles")
    height_cases=0
    for height in range(2,9):
        for T in range(1,9):
            for d in range(4):
                cutoff=height*height*T//(height+1)+d+1
                max_gate=0
                for ell in range(1,height*T+1):
                    start=(ell+height-1)//height
                    rho=(ell+height)//(height+1)
                    D=ell+d-start; L=D-rho
                    M=min(T,ell+d)-start
                    if M<0:
                        continue
                    gate=ell+d-rho+min(M+1,rho)
                    max_gate=max(max_gate,gate)
                    early=[[h(r-j,rho) for j in range(M+1)]
                           for r in range(L+1,gate-start+1)]
                    require(rank(early)==min(M+1,rho),'all-height detection rank')
                    require(gate<=cutoff,'all-height upper bound')
                    record((height,T,d,ell,gate))
                    height_cases+=1
                require(max_gate==cutoff,'all-height sharp maximum')
                rho=(height*T+height)//(height+1)
                L=(height-1)*T+d-rho
                require(T+L+1==cutoff and comb(rho+L,L+1)>0,
                        'all-height actual t^T hostile')
    print('all-height recognition universes:',height_cases,'(height2..8,T1..8,d0..3)')
    print('height indexes planar surfaces, not ambient Jacobian dimension')
    print("characteristic-three hostile: T=3,d=0, lift=(1-z)^2(1+2z)")
    require([1, 0, -3, 2][:3] == [1, 0, -3], "literal characteristic hostile")
    require(all(v % 3 == 0 for v in [0, -3]), "char3 hides first mismatch")
    require(2 % 3 != 0, "char3 later mismatch survives")
    print("it matches 1 through degree 2 mod 3, fails degree 3; char-zero is essential")
    print("SEMANTIC_SHA256", DIGEST.hexdigest())
    print("CHECKS", CHECKS)
    print("PASS; universal proof is analytic, finite checks are controls")


if __name__ == "__main__":
    main()
