"""identities_audit.py - independent auditor's exact-integer checks of THM-4916 (1a)-(1e).

Everything is recomputed from the move definitions by literal letter-by-letter evaluation; no closed form from the
packet is used to produce the left-hand sides.

  O(D,E,U,d) = (D,0,2U,d),  Z = (D-E+U,d,2U-d,d),  J = (D-E+U,D-E,D+U,D),  t = (1,0,2,1)
  value(word) = final D;  F(s) = D - E + U;  tail m: binary digits of m from the lowest bit are the letters (0=O, 1=Z)
  Q(0) = 0, Q(2m) = Q(m), Q(2m+1) = Q(m) + m + (m mod 2)
  K(v,w) = v^T Q w with Session199's matrix [[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]]

usage: python3 identities_audit.py   (deterministic; seed fixed)
"""
import itertools
import random
from math import gcd

random.seed(20261107)


def O(s):
    D, E, U, d = s
    return (D, 0, 2 * U, d)


def Z(s):
    D, E, U, d = s
    return (D - E + U, d, 2 * U - d, d)


def J(s):
    D, E, U, d = s
    return (D - E + U, D - E, D + U, D)


MOVE = {'O': O, 'Z': Z, 'J': J}
T = (1, 0, 2, 1)


def act(s, word):
    for c in word:
        s = MOVE[c](s)
    return s


def value(s, word):
    return act(s, word)[0]


def F(s):
    D, E, U, d = s
    return D - E + U


def Qrec(m):
    if m == 0:
        return 0
    if m % 2 == 0:
        return Qrec(m // 2)
    k = m // 2
    return Qrec(k) + k + (k % 2)


def tailword(m):
    w = ''
    while m:
        w += 'Z' if m & 1 else 'O'
        m >>= 1
    return w


def dagger(w):
    sw = {'O': 'Z', 'Z': 'O', 'J': 'J'}
    return ''.join(sw[c] for c in reversed(w))


QM = [[0, 0, 0, 1], [0, 1, 0, -1], [0, 0, 1, -1], [1, -1, -1, 1]]


def K(v, w):
    return sum(v[i] * QM[i][j] * w[j] for i in range(4) for j in range(4))


def words_upto(n, alphabet='OZJ', maxJ=None):
    for L in range(n + 1):
        for w in itertools.product(alphabet, repeat=L):
            w = ''.join(w)
            if maxJ is None or w.count('J') <= maxJ:
                yield w


def report(name, ok, n):
    print('%-78s %s (%d cases)' % (name, 'OK' if ok else 'FAIL', n))
    assert ok, name


# reachable states with at most two J, words up to length 9, plus random integer 4-tuples (the identities are formal)
reach = {}
for w in words_upto(9, maxJ=2):
    reach.setdefault(act(T, w), w)
states = list(reach.items())
rnd = [(random.randint(-50, 50), random.randint(-50, 50), random.randint(-50, 50), random.randint(-50, 50)) for _ in range(400)]
print('reachable states (<= 2 J, length <= 9): %d distinct; random integer tuples: %d' % (len(states), len(rnd)))

# Q: the incremental rule used in census_lf_audit.c, Q(m + 2^L) = Q(m) + b_{L-1} + rev_L(m), against the recursion
ok = True; n = 0
for m in range(1, 1 << 14):
    L = m.bit_length() - 1
    base = m - (1 << L)
    last = (base >> (L - 1)) & 1 if L >= 1 else 0
    rev = int(format(base, '0%db' % L)[::-1], 2) if L >= 1 else 0
    ok &= Qrec(m) == Qrec(base) + last + rev; n += 1
report('Q(m + 2^L) = Q(m) + b_{L-1} + rev_L(m) for 1 <= m < 2^14', ok, n)
ok = all(Qrec((1 << j) + 1) == (1 << (j - 1)) for j in range(2, 40)) and Qrec(3) == 2
report('Q(2^j + 1) = 2^(j-1) for 2 <= j < 40 (and Q(3) = 2, so j = 1 is excluded)', ok, 38)
ok = all(Qrec((1 << j) - 1) == (1 << j) - 2 for j in range(1, 40))
report('Q(2^j - 1) = 2^j - 2 (so the m = 2^j - 1 tail gives D + 2^j U)', ok, 39)

# tail formula H_s(m) = h + R m - d Q(m) - e (m mod 2) for the seed (h,R,e,d) = (D,U,E,d) of a state s
ok = True; n = 0
for s, _ in states[:3000] + [(r, None) for r in rnd]:
    D, E, U, d = s
    for m in range(0, 1 << 9):
        lhs = value(s, tailword(m))
        rhs = D + U * m - d * Qrec(m) - E * (m % 2)
        ok &= lhs == rhs; n += 1
report('tail formula value(s . tail m) = h + Rm - dQ(m) - e(m mod 2), (h,R,e,d)=(D,U,E,d)', ok, n)

# (1a) seed confinement, and gcd(R,d) = gcd(D,U) of the pre-J state for a root-change seed
ok = True; n = 0; okg = True
for s, _ in states[:3000] + [(r, None) for r in rnd]:
    D, E, U, d = s
    g = gcd(U, d)
    for m in range(0, 1 << 8):
        if g:
            ok &= (value(s, tailword(m)) - (D - E * (m % 2))) % g == 0
        n += 1
    sj = J(s)
    okg &= gcd(sj[2], sj[3]) == gcd(D, U)
report('(1a) H_s(m) = h - e (m mod 2) mod gcd(R,d)', ok, n)
report('(1a) for s = J(p): gcd(R,d) = gcd(D_p, U_p)', okg, len(states[:3000]) + len(rnd))
# the packet's phrase "odd tails land in U (mod gcd(D,U))": U = 0 mod gcd(D,U), i.e. odd tails are DIVISIBLE by gcd(D,U)
ok = True; n = 0
for p, _ in states[:3000]:
    g = gcd(p[0], p[2])
    for m in range(1, 1 << 8, 2):
        ok &= value(J(p), tailword(m)) % g == 0; n += 1
report('(1a) odd tails of J(p) are divisible by gcd(D_p,U_p) (U = 0 mod gcd)', ok, n)

# (1b) root-change tail in pre-change coordinates, and tail confinement
ok = True; okc = True; n = 0
for p, _ in states[:3000] + [(r, None) for r in rnd]:
    D, E, U, d = p
    for m in range(0, 1 << 9):
        N = value(p, 'J' + tailword(m))
        q = Qrec(m)
        if m % 2:
            ok &= N == (m - q) * D + (m + 1) * U
            okc &= N % gcd(q + 1, m + 1) == 0 and gcd(m - q, m + 1) == gcd(q + 1, m + 1)
        else:
            ok &= N == (m + 1 - q) * D + (m + 1) * U - E
        n += 1
report('(1b) value(p J tail m) = (m-Q)D + (m+1)U (m odd), (m+1-Q)D + (m+1)U - E (m even)', ok, n)
report('(1b) gcd(Q(m)+1, m+1) = gcd(m-Q(m), m+1) divides every odd-tail value', okc, n)

# (1c) Fermat-type multipliers, for ALL states (0, 1 and 2 J) and random tuples, j = 2..14; j = 1 fails
ok = True; n = 0; j1fail = 0
for p, _ in states + [(r, None) for r in rnd]:
    D, E, U, d = p
    for j in range(2, 15):
        ok &= value(p, 'JZ' + 'O' * (j - 1) + 'Z') == ((1 << (j - 1)) + 1) * (D + 2 * U); n += 1
    j1fail += value(p, 'JZZ') != 2 * (D + 2 * U)
report('(1c) value(p J Z O^(j-1) Z) = (2^(j-1)+1)(D+2U), j = 2..14, every state', ok, n)
print('     (j = 1, i.e. p J Z Z, differs from 2(D+2U) for %d of %d states: the j >= 2 hypothesis is needed)' % (j1fail, len(states) + len(rnd)))
# membership: for p with <= 1 J the word p J Z O^(j-1) Z has <= 2 J and ends in Z, so its value lies in F_<=2
one_J = [(p, w) for p, w in states if w.count('J') <= 1]
ok = all(('J' + 'Z' + 'O' * (j - 1) + 'Z').count('J') + w.count('J') <= 2 for p, w in one_J for j in range(2, 6))
report('(1c) realizing word has <= 2 J when the prefix word of p has <= 1 J', ok, len(one_J) * 4)

# (1d) defect dynamics and the K expansion
ok = True; n = 0
for s, _ in states[:5000] + [(r, None) for r in rnd]:
    D, E, U, d = s
    a = D - E - U + d
    aZ = (lambda x: x[0] - x[1] - x[2] + x[3])(Z(s)); aJ = (lambda x: x[0] - x[1] - x[2] + x[3])(J(s)); aO = (lambda x: x[0] - x[1] - x[2] + x[3])(O(s))
    ok &= aZ == a and aJ == 0 and aO == a + E - U; n += 1
report('(1d) a = D-E-U+d: Z fixes a, J resets a to 0, O sends a to a + E - U', ok, n)
ok = True; n = 0
pool = [s for s, _ in states[:600]] + rnd[:200]
for v in pool:
    for w in random.sample(pool, 40):
        av = v[0] - v[1] - v[2] + v[3]; aw = w[0] - w[1] - w[2] + w[3]
        ok &= K(v, w) == v[2] * w[2] + v[1] * w[1] - v[3] * w[3] + av * w[3] + aw * v[3]; n += 1
report('(1d) K(v,w) = U_vU_w + E_vE_w - d_vd_w + a_v d_w + a_w d_v', ok, n)
ok = True; n = 0
W7 = list(words_upto(6))
for al in W7:
    for be in random.sample(W7, 40):
        ok &= K(act(T, al), act(T, be)) == F(act(T, be + dagger(al))); n += 1
report("(1d) Session199: K(t.alpha, t.beta) = F(t.(beta alpha^dagger)) (all alpha <= 6 letters)", ok, n)
mn = min(K(act(T, al), act(T, be)) for al in W7[:400] for be in W7[:400])
print('     smallest K value on these pairs: %d (K >= 3 always: K >= U_vU_w/2 + d_vd_w >= 3), so 1 in F_<=2 is not a K-value' % mn)

# (1e) product family, built literally as a word and evaluated letter by letter
ok = True; okJ = True; n = 0
jfree = [(s, w) for s, w in states if 'J' not in w]
for _ in range(3000):
    (s1, w1), (s2, w2) = random.choice(jfree), random.choice(jfree)
    p, q = random.randint(1, 7), random.randint(1, 7)
    x, eps, y, one = s1; x2, eps2, y2, one2 = s2
    assert one == 1 and one2 == 1
    alpha = w1 + 'J' + 'Z' * p
    beta = w2 + 'J' + 'Z' * q
    word = beta + dagger(alpha)            # the prefix; the counted word is word + 'Z'
    val = value(T, word + 'Z')
    ok &= val == (x + (1 << p) * y) * (x2 + (1 << q) * y2) == K(act(T, alpha), act(T, beta))
    okJ &= (word + 'Z').count('J') == 2 and word.count('J') == 2
    n += 1
report('(1e) value(beta alpha^dagger Z) = (x+2^p y)(x\'+2^q y\'), alpha = w J Z^p, beta = w\' J Z^q', ok, n)
report('(1e) the realizing word has exactly two J', okJ, n)
# p = 0 is genuinely excluded
bad = 0
for _ in range(500):
    (s1, w1), (s2, w2) = random.choice(jfree), random.choice(jfree)
    x, eps, y, _o = s1; x2, eps2, y2, _o2 = s2
    if value(T, w2 + 'J' + 'Z' + dagger(w1 + 'J') + 'Z') != (x + y) * (x2 + 2 * y2):
        bad += 1
print('     (p = 0: the product formula fails for %d of 500 random pairs, so p, q >= 1 is needed)' % bad)
print('ALL CHECKS PASSED')
