#!/usr/bin/env python3
"""Independent audit of creation_decoder C3-C5, C9 (degree theorem) and the
transport reading q_-(n) = q_+(-n).  Written from the note's statements only;
does not import or read codex scripts.

q_{s,d}(n) = sum_{j<d} (T_s^j(n) mod 2) 2^j, T_s(n) = n/2 (even), (3n+s)/2 (odd).
ANF degree via the binary Moebius transform.
"""
from itertools import product


def T(n, s):
    return n // 2 if n % 2 == 0 else (3 * n + s) // 2


def q(n, s, d):
    out = 0
    x = n
    for j in range(d):
        out |= (x & 1) << j
        x = T(x, s)
    return out


def anf_degree(truth):
    """truth: list of 0/1 of length 2^m indexed by input integer (bit i = x_i).
    returns algebraic degree of the multilinear F2 polynomial."""
    a = list(truth)
    m = (len(a)).bit_length() - 1
    for i in range(m):
        step = 1 << i
        for x in range(len(a)):
            if x & step:
                a[x] ^= a[x ^ step]
    deg = -1
    for mono, c in enumerate(a):
        if c:
            deg = max(deg, bin(mono).count("1"))
    return deg


def perm_sign(p):
    n = len(p)
    seen = [False] * n
    sgn = 1
    for i in range(n):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                sgn = -sgn
    return sgn


def main():
    # depth-3 tables
    for s in (1, -1):
        print("q_%+d,3 table:" % s, [q(n, s, 3) for n in range(8)])
    qp = [q(n, 1, 3) for n in range(8)]
    qm = [q(n, -1, 3) for n in range(8)]
    lin_m = all(qm[x ^ y] == qm[x] ^ qm[y] for x in range(8) for y in range(8))
    lin_p = all(qp[x ^ y] == qp[x] ^ qp[y] for x in range(8) for y in range(8))
    print("q_- linear over F2^3:", lin_m, "; q_+ linear:", lin_p)
    polar_ok = all((qp[x] ^ qp[y] ^ qp[x ^ y]) == 4 * (((x & 1) & ((y >> 1) & 1)) ^ (((x >> 1) & 1) & (y & 1)))
                   for x in range(8) for y in range(8))
    print("polar defect q+(x)^q+(y)^q+(x^y) = 4(x0y1+x1y0):", polar_ok)
    # Fano lines
    lines = set()
    for x in range(1, 8):
        for y in range(x + 1, 8):
            lines.add(frozenset((x, y, x ^ y)))
    img = {L: frozenset(qp[v] for v in L) for L in lines}
    fixed = sorted(sorted(L) for L in lines if img[L] == L)
    moved = sorted((sorted(L), sorted(img[L]), img[L] in lines) for L in lines if img[L] != L)
    print("Fano lines fixed by q+:", fixed)
    print("Fano lines moved by q+ (line, image, image is a line?):", moved)

    # transport: q_-(n) = q_+(-n mod 2^d)
    print("\nd  deg(top bit q+)  deg(top bit q-)  all-bit degs q+ | q-   sign(q+) sign(q-)  transport ok")
    for d in range(1, 14):
        N = 1 << d
        QP = [q(n, 1, d) for n in range(N)]
        QM = [q(n, -1, d) for n in range(N)]
        transport = all(QM[n] == QP[(-n) % N] for n in range(N))
        bij = len(set(QP)) == N and len(set(QM)) == N
        degp = [anf_degree([(QP[n] >> b) & 1 for n in range(N)]) for b in range(d)]
        degm = [anf_degree([(QM[n] >> b) & 1 for n in range(N)]) for b in range(d)]
        print(d, degp[-1], degm[-1], degp, degm, perm_sign(QP), perm_sign(QM), transport and bij)

    # the negation permutation itself: top-bit degree
    print("\nnegation mod 2^d top-bit ANF degree:",
          [anf_degree([(((-n) % (1 << d)) >> (d - 1)) & 1 for n in range(1 << d)]) for d in range(1, 13)])


if __name__ == "__main__":
    main()
