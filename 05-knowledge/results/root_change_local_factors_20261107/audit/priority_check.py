"""priority_check.py - auditor's exact checks that the THM-4916 identities coincide with, or follow in one line from,
statements already in the repository (see AUDIT.md section 7).

(P1) (1b) is Session223 cutoff/REPORT.md eq. (2) (commit 1c524e45a9): for a nonempty O/Z block C = c_1..c_n,
     F(s.JC) = alpha D + beta E + gamma U with h = a_C + 2^n, r = c_C + q_C + e_C, alpha = 1 - i_C + h - r,
     beta = i_C - 1, gamma = 1 + h. With m = tail integer of the word J C Z (bits c_1..c_n then a top 1), Session223's
     (alpha, beta, gamma) equal (1b)'s coefficients: h = m, r = Q(m).
(P2) (1c), (1e) from KNOWN identities: s.JZ = s.OJ = (D+2U, D, D+2U, D) (Session205; THM-4888 (D)1; Session199
     hidden-structure Theorem 2: balanced sources are exactly W.OJ), J Z^r = O^r J (Session223 bridges section 1), a
     balanced state (h,e,h,e) followed by O^k has F = (2^k+1) h (the mechanism of Session223's rectangle
     F(t.Z^a O^b) = (2^a+1)(2^b+1)), and K((h,e,h,e),(h',e',h',e')) = h h' (Session199 physical-adjoint section 4).
(P3) (1d): the defect a = D-E-U+d is minus Session199's positive coordinate f = U-D+E-d, whose moves are
     O: f -> a199 + 2f + e + 2p, Z: f -> f, J: f -> 0 (Session199 hidden-structure dependency BALANCED_SOURCES.md).
"""
import itertools


def mv(s, c):
    D, E, U, d = s
    if c == 'O':
        return (D, 0, 2 * U, d)
    if c == 'Z':
        return (D - E + U, d, 2 * U - d, d)
    return (D - E + U, D - E, D + U, D)


def act(s, w):
    for c in w:
        s = mv(s, c)
    return s


def F(s):
    return s[0] - s[1] + s[2]


def Qrec(m):
    if m == 0:
        return 0
    if m % 2 == 0:
        return Qrec(m // 2)
    k = m // 2
    return Qrec(k) + k + (k % 2)


ok1 = True; n1 = 0
for n in range(1, 11):
    for bits in itertools.product((0, 1), repeat=n):
        b = list(bits)                                   # b[0] = c_1 (first letter), 1 = Z
        q = sum(b[i] * 2 ** (n - 1 - i) for i in range(n))
        a = sum(b[i] * 2 ** i for i in range(n))
        c = sum(b[i] * b[j] * 2 ** (i - 1 - j) for i in range(n) for j in range(i)) + sum(b[i] * b[i - 1] for i in range(1, n))
        iC, eC = b[0], b[-1]
        h = a + 2 ** n; r = c + q + eC
        alpha, beta, gamma = 1 - iC + h - r, iC - 1, 1 + h
        m = a + 2 ** n                                    # tail integer of J C Z
        Qm = Qrec(m)
        mine = ((m - Qm), 0, (m + 1)) if m % 2 else ((m + 1 - Qm), -1, (m + 1))
        ok1 &= (alpha, beta, gamma) == mine and h == m and r == Qm
        # and literally on a random state
        s = (17, 5, 40, 7)
        ok1 &= F(act(s, 'J' + ''.join('Z' if x else 'O' for x in b))) == alpha * s[0] + beta * s[1] + gamma * s[2]
        n1 += 1
print('(P1) Session223 eq. (2) coefficients == THM-4916 (1b) coefficients (h = m, r = Q(m)), all blocks of length 1..10:', ok1, n1)

ok2 = True; n2 = 0
states = {act((1, 0, 2, 1), ''.join(w)) for L in range(0, 8) for w in itertools.product('OZJ', repeat=L)}
for s in states:
    D, E, U, d = s
    ok2 &= act(s, 'JZ') == act(s, 'OJ') == (D + 2 * U, D, D + 2 * U, D)
    for r in range(0, 6):
        ok2 &= act(s, 'J' + 'Z' * r) == act(s, 'O' * r + 'J')
    h, e = D + 2 * U, D
    for k in range(1, 8):
        ok2 &= F(act((h, e, h, e), 'O' * k)) == (2 ** k + 1) * h          # balanced then O^k
        ok2 &= F(act(s, 'JZ' + 'O' * k)) == (2 ** k + 1) * (D + 2 * U)      # = (1c) with j - 1 = k
    n2 += 1
QM = [[0, 0, 0, 1], [0, 1, 0, -1], [0, 0, 1, -1], [1, -1, -1, 1]]
K = lambda v, w: sum(v[i] * QM[i][j] * w[j] for i in range(4) for j in range(4))
for h in range(3, 40, 2):
    for e in range(1, h):
        for h2 in range(3, 30, 4):
            ok2 &= K((h, e, h, e), (h2, 1, h2, 1)) == h * h2
jfree = [s for s in {act((1, 0, 2, 1), ''.join(w)) for L in range(0, 7) for w in itertools.product('OZ', repeat=L)}]
for s in jfree:
    x, eps, y, one = s
    for p in range(1, 6):
        v = act(s, 'J' + 'Z' * p)
        ok2 &= v == (x + 2 ** p * y, x, x + 2 ** p * y, x)              # w J Z^p is balanced, h = x + 2^p y
print('(P2) JZ = OJ, JZ^r = O^r J, balanced O^k product, K = hh\' on balanced pairs, wJZ^p balanced with h = x + 2^p y:', ok2, n2)

ok3 = True
for s in states:
    D, E, U, d = s
    a199, f, e, p = D - E - d, U - D + E - d, E, d - E
    aS15 = D - E - U + d
    ok3 &= aS15 == -f
    O_, Z_, J_ = act(s, 'O'), act(s, 'Z'), act(s, 'J')
    fO = O_[2] - O_[0] + O_[1] - O_[3]; fZ = Z_[2] - Z_[0] + Z_[1] - Z_[3]; fJ = J_[2] - J_[0] + J_[1] - J_[3]
    ok3 &= fO == a199 + 2 * f + e + 2 * p and fZ == f and fJ == 0
print('(P3) a = D-E-U+d = -f (Session199 positive coordinate) with Session199\'s O/Z/J laws for f:', ok3, len(states))
