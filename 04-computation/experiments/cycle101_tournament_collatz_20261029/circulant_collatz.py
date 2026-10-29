"""Circulant tournaments as rational Collatz cycles.

For an odd prime q and a connection set S of Z_q (|S| = (q-1)/2, S and -S partition Z_q^*), read the indicator of S
as a cyclic T-parity word (1 = odd step).  The ones sit at the elements of S; the gap word w (k = (q-1)/2 letters,
sum H = q) defines the rational 3x+1 cycle x = C(w)/D with D = 2^q - 3^((q-1)/2) (classical necklace <-> cycle
correspondence).  Translation S -> S + t rotates the word (same cycle); multipliers S -> mS act by decimation;
the converse S -> -S reverses the word.  For each q this lists: D and its factorization, every connection set,
its necklace class, the cycle (reduced denominator, integer points of 3x + d), the tournament isomorphism class
(multiplier orbit), the cycle sum, and whether converse pairs share sums.  Also the 'append a letter 2' map from
Z_q-words to Z_(q+2)-words."""
import itertools, math, json
from fractions import Fraction
from sympy import factorint

def C_of(w, a=3):
    s, out, kk = 0, 0, len(w)
    for j, v in enumerate(w):
        out += a ** (kk - 1 - j) * 2 ** s
        s += v
    return out

def gap_word(S, q):
    pts = sorted(S)
    return tuple((pts[(i + 1) % len(pts)] - pts[i]) % q or q for i in range(len(pts)))

def canon(S, q):
    return min(tuple(sorted((x - t) % q for x in S)) for t in range(q))

def cycle_info(S, q):
    w = gap_word(S, q)
    H, k = sum(w), len(w)
    D = 2 ** H - 3 ** k
    pts = [Fraction(C_of(w[i:] + w[:i]), D) for i in range(k)]
    d = pts[0].denominator
    assert all(p.denominator == d for p in pts)
    ints = sorted(int(p * d) for p in pts)
    return w, D, d, ints, sum(pts)

out = {}
for q in (3, 5, 7, 11, 13):
    D = 2 ** q - 3 ** ((q - 1) // 2)
    pairs = [(x, q - x) for x in range(1, (q + 1) // 2)]
    sets = [frozenset(c) for c in itertools.product(*pairs)]
    classes = {}
    for S in sets:
        classes.setdefault(canon(S, q), []).append(S)
    # tournament isomorphism classes = multiplier orbits of connection sets (Turner/Elspas at prime order)
    iso = {}
    for S in sets:
        key = min(tuple(sorted((m * x) % q for x in S)) for m in range(1, q))
        iso.setdefault(key, []).append(S)
    print('q=%d: D = 2^%d - 3^%d = %d = %s; %d connection sets, %d necklace classes, %d tournament classes'
          % (q, q, (q - 1) // 2, D, factorint(D), len(sets), len(classes), len(iso)))
    rows = []
    for key, Ss in sorted(classes.items()):
        S = Ss[0]
        w, Dw, d, ints, sm = cycle_info(S, q)
        tclass = [i for i, (k2, v) in enumerate(sorted(iso.items())) if S in v][0]
        rev = canon(frozenset((-x) % q for x in S), q)
        rows.append({'class': key, 'sets': [sorted(x) for x in Ss], 'word': w, 'reduced_den': d, 'cycle_of_3x+d': ints,
                     'sum': str(sm), 'tournament_class': tclass, 'reversal_class': rev})
        print('   class %s: connection sets %s | word %s | reduced denominator %d | 3x+%d cycle %s | sum %s | tournament #%d | reversal -> %s'
              % (key, [sorted(x) for x in Ss], w, d, d, ints, sm, tclass, 'self' if rev == key else str(rev)))
    out[q] = rows
# 'append 2': Z_5 words to Z_7 words
print('\nappend a letter 2 to the Z_5 circulant cycles (all markings):')
for S in [frozenset({1, 2}), frozenset({1, 3})]:
    w = gap_word(S, 5)
    for i in range(len(w)):
        ww = w[i:] + w[:i] + (2,)
        pos = [0]
        for v in ww[:-1]:
            pos.append(pos[-1] + v)
        print('   %s marked %s + (2) -> %s, Z_7 class %s' % (sorted(S), w[i:] + w[:i], ww, canon(frozenset(pos), 7)))
json.dump({str(k): v for k, v in out.items()}, open('out/circulant_collatz.json', 'w'), indent=1, default=str)
