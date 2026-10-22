"""Audit S50 - a5: claims of THM-4744 (3) and of the results note sections 1-3, checked on my own inventory
(inventory_375.json from a3_ear_census.py; complete for every strong class with H <= 375).
  (I1) strong polynomials P_d for every odd d < 375 compared term by term with the session's printed lists
       (h_patterns.out for d < 125, strong_ears_375.out for 125 <= d < 375; the files are read as data);
  (I2) the sequence law with my strong inventory: total class counts A000568 for n <= 8 (validation), the note's
       N(n,h) table for n = 10, 11, 12, and the note's 'small polynomials' (all n <= 30);
  (I3) converse pairing at orders 10 and 11; H = 1 (mod 4) for self-converse even-order classes; |Aut| divides H;
  (I4) second floor by order; the claim 'at odd m >= 7 the second floor is 3^((m-1)/2)' tested at m = 13 and 15;
  (I5) late seven: no strong class of order <= 6 has 7 | H; multiples of 7 at n = 7 with multiplicities; uniquely
       attained values at n = 7 and the regular 7-tournaments; holes of n = 7 in (45, 150)."""
import os, re, json, time
from math import comb
from collections import defaultdict, Counter
import alib as L

T0 = time.time()
FAIL = []
HERE = os.path.dirname(os.path.abspath(__file__))
SESS = os.path.dirname(HERE)


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


inv = json.load(open(os.path.join(HERE, 'inventory_375.json')))
orders = {int(k): [(int(c), h, sc, aut) for c, h, sc, aut in v] for k, v in inv['orders'].items()}
full8 = [(int(c), h) for c, h in inv['uncapped8']]

# ------------------------------------------------------------------ (I1)
P = defaultdict(Counter)
for m, rows in orders.items():
    for c, h, sc, aut in rows:
        if h < 375:
            P[h][m] += 1
zeros = [d for d in range(1, 375, 2) if not P.get(d)]
check(zeros == [1, 7, 21], f'(I1) odd d < 375 with P_d = 0: {zeros}')


def parse(path, lo, hi):
    out = {}
    for line in open(path, encoding='utf-8'):
        mm = re.match(r'\s*P_(\d+) = (.*)$', line)
        if mm:
            d = int(mm.group(1))
            if lo <= d < hi:
                terms = Counter()
                for t in mm.group(2).split('+'):
                    t = t.strip()
                    tm = re.match(r'(\d*)z\^(\d+)', t)
                    terms[int(tm.group(2))] += int(tm.group(1)) if tm.group(1) else 1
                out[d] = terms
    return out


sess = parse(os.path.join(SESS, 'h_patterns.out'), 25, 125)
sess.update(parse(os.path.join(SESS, 'strong_ears_375.out'), 125, 375))
small = {3: {3: 1}, 5: {4: 1}, 9: {5: 2}, 11: {5: 1}, 13: {5: 1}, 15: {5: 2, 6: 2}, 17: {6: 2}, 19: {6: 2}, 23: {6: 4}}
for d, v in small.items():
    sess[d] = Counter(v)
mism = [d for d in range(3, 375, 2) if d not in (7, 21) and Counter(P[d]) != sess.get(d)]
check(not mism and len(sess) == 184, f'(I1) my P_d = the session\'s P_d for all {len(sess)} odd d in [3, 373] other than 7, 21' + (f'; MISMATCH {mism[:8]}' if mism else ''))
note_sel = {25: {6: 2, 7: 1}, 27: {6: 2, 7: 3}, 35: {7: 2}, 45: {6: 2, 7: 5, 8: 5}, 49: {7: 7, 8: 1}, 63: {8: 2},
            75: {7: 7, 8: 2, 9: 3}, 125: {7: 1, 8: 10, 9: 2, 10: 1}, 135: {7: 8, 8: 14, 9: 14, 10: 10},
            147: {7: 1, 8: 10, 9: 10, 10: 2}}
check(all(Counter(P[d]) == Counter(v) for d, v in note_sel.items()), '(I1) the note\'s table of selected P_d (section 1) is correct')

# ------------------------------------------------------------------ (I2)
strong_by_order = defaultdict(Counter)       # complete for H < 375
for m, rows in orders.items():
    for c, h, sc, aut in rows:
        if h < 375:
            strong_by_order[m][h] += 1


def tuples(strong, hmax, nmax):
    """B[(k, t, h)] = number of ordered k-tuples of nontrivial strong classes, total size t, product h <= hmax."""
    B = defaultdict(int)
    B[(0, 0, 1)] = 1
    frontier = {(0, 0, 1): 1}
    while frontier:
        nxt = defaultdict(int)
        for (k, t, h), c in frontier.items():
            for m, dist in strong.items():
                if t + m > nmax:
                    continue
                for h1, c1 in dist.items():
                    if h * h1 <= hmax:
                        nxt[(k + 1, t + m, h * h1)] += c * c1
        for key, c in nxt.items():
            B[key] += c
        frontier = nxt
    return B


def N_of(B, n):
    res = defaultdict(int)
    for (k, t, h), c in B.items():
        if n - t >= 0:
            res[h] += c * comb(n - t + k, k)
    return res


# validation: full strong inventories to order 8 give A000568
full = defaultdict(Counter)
for m in range(3, 8):
    for c, h, sc, aut in orders[m]:
        full[m][h] += 1
for c, h in full8:
    full[8][h] += 1
Bfull = tuples(full, 10 ** 9, 8)
tot = [sum(N_of(Bfull, n).values()) for n in range(1, 9)]
check(tot == [1, 1, 2, 4, 12, 56, 456, 6880], f'(I2) sequence law on my complete strong inventories (orders <= 8) gives the class counts {tot} (A000568)')
B = tuples(strong_by_order, 374, 40)
table = {63: (6, 8, 10), 75: (70, 121, 192), 105: (59, 85, 115), 125: (43, 67, 98), 135: (146, 278, 490),
         147: (70, 123, 192), 189: (102, 148, 198), 195: (126, 238, 396), 223: (82, 117, 152)}
mine = {h: tuple(N_of(B, n).get(h, 0) for n in (10, 11, 12)) for h in table}
check(mine == table, f'(I2) the note\'s N(n,h) table (n = 10, 11, 12) reproduced: {mine}')
polys = {3: lambda n: n - 2, 9: lambda n: (n - 4) * (n - 1) // 2, 27: lambda n: n * (n - 5) * (n - 4) // 6,
         33: lambda n: n * n - 5 * n - 1, 45: lambda n: 6 * n * n - 74 * n + 233, 49: lambda n: 8 * n - 49,
         63: lambda n: 2 * (n - 7), 75: lambda n: 4 * (n - 6) ** 2}
onset = {3: 2, 9: 4, 27: 6, 33: 6, 45: 7, 49: 7, 63: 7, 75: 8}
rows = []
for h, f in polys.items():
    bad = [n for n in range(onset[h], 31) if N_of(B, n).get(h, 0) != f(n)]
    rows.append((h, bad[:3]))
good = [h for h, bad in rows if not bad]
badp = {h: bad for h, bad in rows if bad}
print(f'     small polynomials of the note: correct for all n <= 30: {good}; WRONG (first n): {badp}')
true45 = lambda n: 12 * n - 75 + 4 * comb(max(n - 6, 0), 2) + 8 * comb(max(n - 7, 0), 2) + 3 * comb(max(n - 7, 0), 3)
true75 = lambda n: 12 * n - 80 + 8 * comb(max(n - 7, 0), 2) + 6 * comb(max(n - 8, 0), 2) + 3 * comb(max(n - 8, 0), 3)
check(all(N_of(B, n).get(45, 0) == true45(n) for n in range(7, 31)) and all(N_of(B, n).get(75, 0) == true75(n) for n in range(8, 31)),
      '(I2) the true laws: N(n,45) = 12n-75+4C(n-6,2)+8C(n-7,2)+3C(n-7,3) (n >= 7), N(n,75) = 12n-80+8C(n-7,2)+6C(n-8,2)+3C(n-8,3) (n >= 8); '
      f'e.g. N(10,45) = {N_of(B, 10)[45]} (note polynomial gives {polys[45](10)}), N(10,75) = {N_of(B, 10)[75]} (note polynomial gives {polys[75](10)})')

# ------------------------------------------------------------------ (I3)
ok_pair = True; ok_mod4 = True; ok_aut = True
for m in (10, 11):
    byh = defaultdict(lambda: [0, 0])
    for c, h, sc, aut in orders[m]:
        byh[h][0] += 1
        byh[h][1] += 1 if sc else 0
        ok_aut &= h % aut == 0
        if m % 2 == 0 and sc:
            ok_mod4 &= h % 4 == 1
    ok_pair &= all((a - s) % 2 == 0 for a, s in byh.values())
nsc10 = sum(1 for c, h, sc, aut in orders[10] if sc)
check(ok_pair, '(I3) converse pairing P_d = SC_d (mod 2) at orders 10 and 11 (H <= 375)')
check(ok_mod4, f'(I3) every self-converse strong 10-class with H <= 375 ({nsc10} classes) has H = 1 (mod 4)')
check(ok_aut, '(I3) |Aut| divides H for every class at orders 10, 11')
ok8 = True
for c, h in full8:
    adj = L.decode(8, c)
    if L.canon(L.converse(adj)) == c:
        ok8 &= h % 4 == 1
check(ok8, '(I3) every self-converse strong 8-class (complete census) has H = 1 (mod 4)')

# ------------------------------------------------------------------ (I4)
sec = {m: sorted(set(h for c, h, sc, aut in orders[m]))[1] for m in range(5, 12)}
check([sec[m] for m in range(5, 12)] == [11, 17, 27, 49, 81, 135, 243], f'(I4) second-lowest strong count, orders 5..11: {[sec[m] for m in range(5, 12)]}')
mult = {m: sum(1 for c, h, sc, aut in orders[m] if h == sec[m]) for m in (7, 9, 11)}
allt = {m: L.canon(L.chain((2,) * ((m - 1) // 2))) in {c for c, h, sc, aut in orders[m] if h == sec[m]} for m in (7, 9, 11)}
print(f'     classes attaining the second floor at m = 7, 9, 11: {mult}; the all-triangle chain among them: {allt}')
h13 = L.H(L.chain((2, 2, 2, 3, 3))); h15 = L.H(L.chain((2, 2, 2, 2, 3, 3)))
check(h13 == 675 and 625 < h13 < 3 ** 6 and 1875 < h15 < 3 ** 7,
      f'(I4) the note\'s "at odd m >= 7 the second floor is the all-triangle chain 3^((m-1)/2)" FAILS: H(IC(2,2,2,3,3)) = {h13} lies in '
      f'(f(13), 3^6) = (625, 729); H(IC(2,2,2,2,3,3)) = {h15} lies in (f(15), 3^7) = (1875, 2187)')

# ------------------------------------------------------------------ (I5)
check(all(h % 7 for m in range(3, 7) for c, h, sc, aut in orders[m]), '(I5) no strong class of order <= 6 has 7 | H (complete inventories)')
m7 = Counter(h for c, h, sc, aut in orders[7] if h % 7 == 0)
check(dict(sorted(m7.items())) == {35: 2, 49: 7, 77: 10, 91: 15, 105: 4, 133: 6, 147: 1, 175: 1, 189: 1},
      f'(I5) multiples of 7 among strong 7-classes: {dict(sorted(m7.items()))}')
N7 = N_of(tuples(full, 10 ** 9, 7), 7)
nonstrong_max = max(h for h in N_of(tuples({m: full[m] for m in range(3, 7)}, 10 ** 9, 7), 7))
check(nonstrong_max == 45 and sum(N7.values()) == 456, f'(I5) every non-strong 7-tournament has H <= {nonstrong_max}; 456 classes at n = 7 (so every 7-tournament with 7 | H is strong)')
uniq = sorted(h for h, c in N7.items() if c == 1)
check(uniq == [1, 61, 85, 87, 125, 127, 147, 155, 157, 171, 175, 189], f'(I5) uniquely attained values at n = 7: {uniq}')
reg = sorted((h, aut) for c, h, sc, aut in orders[7] if all(bin(a).count('1') == 3 for a in L.decode(7, c)))
check(reg == [(171, 3), (175, 7), (189, 21)], f'(I5) regular 7-tournaments (H, |Aut|): {reg}')
holes = [h for h in range(47, 150, 2) if N7.get(h, 0) == 0]
primes_attained = [h for h in range(47, 150, 2) if all(h % q for q in range(2, int(h ** 0.5) + 1)) and N7.get(h, 0) > 0]
check(holes == [63, 107, 119, 149], f'(I5) odd holes of n = 7 in (45, 150): {holes}')
print(f'     primes in (45, 150) that ARE attained at n = 7 (none can borrow a product factor there either): {primes_attained}')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
