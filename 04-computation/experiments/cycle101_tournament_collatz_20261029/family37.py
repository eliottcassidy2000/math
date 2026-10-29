"""The (k,H) = (3,7) rational 3x+1 cycles (denominator 2^7 - 3^3 = 101) as weight-3 binary necklaces of length 7.

A valuation word w = (v0, v1, v2), sum 7, is the U_3 cycle x = C(w)/D with C = 9*2^s0 + 3*2^s1 + 2^s2,
s = (0, v0, v0+v1), D = 2^7 - 27 = 101 (collatz_fermat_reframes_20261029, C1).  Its T-map parity word is the
characteristic vector of S = {s0, s1, s2} in Z_7 (positions of odd steps), so rotation classes of words are
rotation classes of 3-subsets of Z_7.  The decimation group Z_7^* (S -> mS, the multiplier group of circulant
tournaments) and reversal (S -> -S, the converse) act on the five classes.  Also: the one-letter future
spectrum (append a letter v to each of the 15 marked words; new D = 2^(7+v) - 81) and the integer forcing
cycles of 3x+101.
"""
import itertools, math, json
from fractions import Fraction

a, k, H = 3, 3, 7
D = 2 ** H - a ** k


def C_of(word, a=3):
    s, out, kk = 0, 0, len(word)
    for j, v in enumerate(word):
        out += a ** (kk - 1 - j) * 2 ** s
        s += v
    return out


def subset_of(word, H):
    s, S = 0, []
    for v in word:
        S.append(s % H)
        s += v
    return frozenset(S)


def word_of_subset(S, H):
    """Rotate so 0 in S; return the gap word starting at 0."""
    pts = sorted(S)
    return tuple((pts[(i + 1) % len(pts)] - pts[i]) % H or H for i in range(len(pts)))


def canon_subset(S, H):
    return min(tuple(sorted((x - t) % H for x in S)) for t in range(H))


words = [w for w in itertools.product(range(1, H + 1), repeat=k) if sum(w) == H]
classes = {}
for w in words:
    classes.setdefault(canon_subset(subset_of(w, H), H), []).append(w)
out = {'D': D, 'classes': []}
print('(k,H) = (3,7): D = 2^7 - 3^3 =', D, '; words:', len(words), '; necklaces:', len(classes))
for key, ws in sorted(classes.items()):
    rows = []
    for w in ws:
        C = C_of(w)
        x = Fraction(C, D)
        # verify the cycle step: U(x) = (3x+1)/2^v0 is the next rotation
        nxt = Fraction(3 * C + D, 2 ** w[0] * D)
        assert nxt == Fraction(C_of(w[1:] + w[:1]), D)
        rows.append((w, C))
    rev = canon_subset(frozenset((-x) % H for x in key), H)
    dec = {m: canon_subset(frozenset((m * x) % H for x in key), H) for m in range(1, H)}
    print('  S =', key, 'words', [r[0] for r in rows], 'numerators', [r[1] for r in rows],
          '| reversal ->', rev, '| x2 ->', dec[2], '| x3 ->', dec[3])
    out['classes'].append({'S': list(key), 'words': [list(r[0]) for r in rows], 'C': [r[1] for r in rows],
                           'reversal': list(rev), 'dec': {m: list(v) for m, v in dec.items()}})
# decimation orbits
orb = {}
for key in classes:
    o = frozenset(canon_subset(frozenset((m * x) % H for x in key), H) for m in range(1, H))
    orb[o] = 1
print('Z_7^* decimation orbits on the five necklaces:', [sorted(o) for o in orb])
QR = {1, 2, 4}
for key in classes:
    trans = [frozenset((x + t) % H for x in key) for t in range(H)]
    if frozenset(QR) in trans:
        print('  necklace', key, 'contains a translate of QR = {1,2,4}: Paley P7 connection set')
    if frozenset({3, 5, 6}) in trans:
        print('  necklace', key, 'contains a translate of NQR = {3,5,6}: the converse P7^op')
# one-letter future spectrum
print('\nOne-letter futures: append v to each marked word w; C(wv) = 3 C(w) + 2^7, D(wv) = 2^(7+v) - 81')
fut = {}
for v in range(1, 7):
    Dn = 2 ** (H + v) - a ** (k + 1)
    spec = {}
    for w in words:
        C = C_of(w)
        Cn = C_of(w + (v,))
        assert Cn == 3 * C + 2 ** H
        g = math.gcd(Cn, Dn)
        spec.setdefault(Dn // g, []).append((w, C))
    fut[v] = {'D': Dn, 'factor': None, 'spectrum': {str(dd): [list(x[0]) for x in ws] for dd, ws in spec.items()}}
    print('  v=%d: D = %d; reduced denominators -> numerators C(w):' % (v, Dn),
          {dd: sorted(c for _, c in ws) for dd, ws in sorted(spec.items())})
out['futures'] = fut
# 3x+101 integer cycles from the family, and the fixed point 101 (scaled trivial cycle)
print('\n3x+101 integer cycles of the family (numerators):', sorted(sorted(c for _, c in [(w, C_of(w)) for w in ws]) for ws in classes.values()))
json.dump(out, open('out/family37.json', 'w'), indent=1, default=str)
