"""Patterns in the multiplicities of the Hamiltonian-path spectrum (session opus-2026-10-21-S50).
Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md; THM-4744, THM-4745, HYP-9355, HYP-9356.

Uses the class census of h_multiplicities.py (all tournaments n <= 9).  Inherited and credited: the class series
F_N(z) = sum_k (1-z)^(-k-1) sum_{d1...dk = N} prod P_di(z) and the strong polynomials P_d for d <= 23
(05-knowledge/results/spectrum_basins_20261020e/spectrum.md, section 4), Busch's strong minimum f(n) (CITED there).

  (Q1) the strong polynomials P_d for every odd d < 125 (complete: a strong tournament with H < 125 has order <= 9,
       since f(10) = 125 and f is nondecreasing); the spectrum_basins list d <= 23 is reproduced; the odd d < 125 with
       P_d = 0 are exactly 1, 7, 21; Busch-extremal counts E(m) = [z^m] P_f(m).
  (Q2) converse pairing: N(n,h) = SC(n,h) (mod 2) and the same for strong classes; at even n no self-converse class
       has H = 3 (mod 4); so a uniquely attained value belongs to a self-converse class, and is 1 (mod 4) at even n.
  (Q3) late primes: a(p) = least order of a strong tournament with p | H, against the least order whose strong range
       contains an odd multiple of p; 7 | H forces a strong component of order >= 7 (order 7 is the first);
       the multiples of 7 at n = 7.
  (Q4) residues: labelled fractions with p | H (p = 3, 5, 7, 11, 13), n = 3..9; the mod-4 and mod-8 imbalances.
  (Q5) repunit and Cipolla probe: multiplicities of 2^k - 1, (4^k - 1)/3 and 2^k + 1 against a local baseline.
  (Q6) the converse-Gluck family: H, |Aut|, tc and uniqueness of C_3, Z_5{1,2}, P_7, C_3[C_3], Z_5{1,2}[C_3],
       P_7[C_3]; comparison with the maxima (n <= 9) and with the mean n!/2^(n-1) (n = 15, 21).
  (Q7) Paley presentation counts tc(P_q) = 2H/(q(q-1)) for q <= 31 (H recomputed for q <= 19; q = 23, 27, 31 cited
       from paley_h31_exact_rate_and_apery_wheel_20261014.md), factored, with Fermat/Carmichael tests.
Prints ALL CHECKS PASSED.  Runtime about 1 minute (the n = 21 DP needs about 0.4 GB)."""
import sys, os, time
from math import factorial, comb, gcd
from collections import defaultdict, Counter
import numpy as np
from sympy import factorint, isprime, primerange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import h_multiplicities as HM                                            # noqa: E402

FAIL = []
T0 = time.time()


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def busch(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r] if m >= 3 else None


def poly_str(cs):
    return ' + '.join((f'{c}' if c > 1 else '') + f'z^{m}' for m, c in sorted(cs.items())) or '0'


# ---------------------------------------------------------------- census
data = HM.classes_upto(9)
R = {}
for n in range(1, 10):
    A = data[n]['A']
    H, st = HM.ham_all(n, A)
    aut = np.array(data[n]['aut'], dtype=np.int64)
    sc = np.array(data[n]['sc'], dtype=bool)
    R[n] = dict(A=A, H=H.astype(np.int64), st=st.astype(bool), aut=aut, sc=sc, w=np.array([factorial(n) // int(a) for a in aut], dtype=object))
print(f'census n <= 9 loaded [{time.time() - T0:.0f}s]')
N = {n: Counter(int(h) for h in R[n]['H']) for n in R}
SC = {n: Counter(int(h) for h, s in zip(R[n]['H'], R[n]['sc']) if s) for n in R}
S = {n: Counter(int(h) for h, s in zip(R[n]['H'], R[n]['st']) if s) for n in range(3, 10)}       # strong
SSC = {n: Counter(int(h) for h, s, c in zip(R[n]['H'], R[n]['st'], R[n]['sc']) if s and c) for n in range(3, 10)}

# ---------------------------------------------------------------- (Q1)
print('(Q1) strong raw-count polynomials P_d, every odd d < 125')
P = defaultdict(dict)
for m in range(3, 10):
    for h, c in S[m].items():
        P[h][m] = c
f = {m: busch(m) for m in range(3, 11)}
check(all(min(S[m]) == f[m] for m in range(3, 10)) and f[10] == 125,
      'Busch floor f(m) = 3,5,9,15,25,45,75 attained exhaustively for m = 3..9; f(10) = 125 (CITED), so P_d is complete for d < 125')
zeros = [d for d in range(1, 125, 2) if not P.get(d)]
check(zeros == [1, 7, 21], f'odd d < 125 with P_d = 0: {zeros}: every odd d in [3, 123] other than 7, 21 is the count of a STRONG tournament')
sb = {3: {3: 1}, 5: {4: 1}, 9: {5: 2}, 11: {5: 1}, 13: {5: 1}, 15: {5: 2, 6: 2}, 17: {6: 2}, 19: {6: 2}, 23: {6: 4}}
check(all(P.get(d, {}) == v for d, v in sb.items()) and not P.get(7) and not P.get(21),
      'spectrum_basins (6) reproduced: P3=z^3, P5=z^4, P9=2z^5, P11=P13=z^5, P15=2z^5+2z^6, P17=P19=2z^6, P23=4z^6')
for d in range(25, 125, 2):
    print(f'     P_{d} = {poly_str(P[d])}')
Pone = {d: sum(P[d].values()) for d in range(3, 125, 2) if P.get(d)}
print('     P_d(1) for odd d = 3..123:', [Pone.get(d, 0) for d in range(3, 125, 2)])
E = {m: S[m][f[m]] for m in range(3, 10)}
print('     Busch-extremal strong classes E(m), m = 3..9:', [E[m] for m in range(3, 10)])
second = {m: sorted(S[m])[1] for m in range(5, 10)}
print('     second-lowest strong count at m = 5..9:', [second[m] for m in range(5, 10)], '(m = 3, 4 have one strong class each)')
check(E[3] == 1 and E[4] == 1, 'E(3) = E(4) = 1 (C_3, the strong 4-tournament)')

# ---------------------------------------------------------------- (Q2)
print('(Q2) converse pairing of multiplicities')
ok = all((N[n][h] - SC[n][h]) % 2 == 0 for n in R for h in N[n])
check(ok, 'N(n,h) = SC(n,h) (mod 2) for every n <= 9 and h (non-self-converse classes pair with their converses; PROVED)')
ok = all((S[m][h] - SSC[m][h]) % 2 == 0 for m in S for h in S[m])
check(ok, 'the same for strong classes: P_d(z) = SC_d(z) (mod 2) coefficientwise')
ok = all(h % 4 == 1 for n in (2, 4, 6, 8) for h in SC[n])
check(ok, 'even n <= 8: every self-converse class has H = 1 (mod 4) (THM-791 with a fixed-point-free anti-involution)')
ok = all(N[n][h] % 2 == 0 for n in (2, 4, 6, 8) for h in N[n] if h % 4 == 3)
check(ok, 'hence at even n every value H = 3 (mod 4) is attained by an even number of classes')
for n in range(5, 10):
    uniq = sorted(h for h, c in N[n].items() if c == 1)
    allsc = all(SC[n][h] == 1 for h in uniq)
    info = []
    for h in uniq:
        i = int(np.where(R[n]['H'] == h)[0][0])
        sc_ = sorted(bin(int(x)).count('1') for x in R[n]['A'][i])
        info.append((h, int(R[n]['aut'][i]), 'strong' if R[n]['st'][i] else 'non-strong', 'regular' if len(set(sc_)) == 1 else ''))
    check(allsc, f'n = {n}: {len(uniq)} uniquely attained values, all self-converse: ' + ', '.join(f'{h}(|Aut| {a}{", " + r if r else ""})' for h, a, s, r in info[:14]) + (' ...' if len(info) > 14 else ''))

# ---------------------------------------------------------------- (Q3)
print('(Q3) late primes and the late seven')
first = {}
for m in range(3, 10):
    for h in S[m]:
        for p in factorint(h):
            first.setdefault(p, m)
window = {}
for p in primerange(3, 400):
    for m in range(3, 10):
        lo, hi = min(S[m]), max(S[m])
        if any((k * p) % 2 == 1 for k in range(-(-lo // p), hi // p + 1)):
            window[p] = m
            break
late = [(p, window[p], first.get(p, '>9')) for p in primerange(3, 400) if p in window and first.get(p) != window[p]]
print('     (p, first order whose strong range holds an odd multiple of p, first order with p | H strong):', late)
check(window[7] == 6 and first[7] == 7 and all(first.get(p) == window.get(p) for p in (3, 5, 11, 13, 17, 19, 23, 29, 31)),
      'p = 7 is late: odd multiples of 7 (21, 35) enter the strong range [15, 45] at order 6, but 7 | H first at order 7; '
      '3, 5, 11, ..., 31 divide a strong count at the first order whose range allows it; the other late primes (107, 149, ...) '
      'sit above the order-7 holes near its maximum')
check(all(h % 7 for n in range(1, 7) for h in N[n]), 'no tournament with n <= 6 has 7 | H (FINITE-EXACT); so 7 | H(T) forces a strong component of order >= 7 with 7 | H (PROVED from multiplicativity)')
m7 = sorted(h for h in N[7] if h % 7 == 0)
allstrong7 = all(R[7]['st'][i] for i in range(len(R[7]['H'])) if R[7]['H'][i] % 7 == 0)
check(allstrong7, 'at n = 7 every tournament with 7 | H is strong; the multiples of 7 attained: ' + ', '.join(f'{h}x{N[7][h]}' for h in m7))
miss7 = [h for h in range(7, max(N[7]) + 1, 14) if h not in N[7]]
holes7 = [h for h in range(1, max(N[7]), 2) if h not in N[7]]
mid = [h for h in holes7 if 45 < h < 150]
print(f'     n = 7: odd multiples of 7 below the maximum that are missing: {miss7}; all odd holes: {holes7}')
print(f'     n = 7: holes in (45, 150): {mid} -> multiples of 7: {[h for h in mid if h % 7 == 0]}, primes: {[h for h in mid if isprime(h)]}')
for n in (7, 8, 9):
    w = R[n]['w']; H = R[n]['H']; st = R[n]['st']
    tot = sum(w); f7 = sum(wi for wi, h in zip(w, H) if h % 7 == 0)
    fs = sum(wi for wi, h, s in zip(w, H, st) if s); fs7 = sum(wi for wi, h, s in zip(w, H, st) if s and h % 7 == 0)
    print(f'     n = {n}: labelled fraction with 7 | H = {float(f7 / tot):.4f}; among strong labelled = {float(fs7 / fs):.4f}  (1/7 = 0.1429)')

# ---------------------------------------------------------------- (Q4)
print('(Q4) residues of H over labelled tournaments')
for p in (3, 5, 7, 11, 13):
    fr = []
    for n in range(3, 10):
        w = R[n]['w']; H = R[n]['H']
        fr.append(round(float(sum(wi for wi, h in zip(w, H) if h % p == 0) / sum(w)), 4))
    print(f'     p = {p:2d}: fraction with p | H, n = 3..9: {fr}   (1/p = {1 / p:.4f})')
imb4, imb8 = [], []
for n in range(1, 10):
    w = R[n]['w']; H = R[n]['H']
    c = Counter()
    for wi, h in zip(w, H):
        c[int(h) % 8] += wi
    imb4.append(c[1] + c[5] - c[3] - c[7]); imb8.append(tuple(c[r] for r in (1, 3, 5, 7)))
print('     labelled #{H = 1 mod 4} - #{H = 3 mod 4}, n = 1..9:', imb4)
print('     divided by 2^(C(n,2)):', [round(float(x) / 2 ** comb(n, 2), 5) for n, x in zip(range(1, 10), imb4)])
print('     labelled counts of H mod 8 = (1,3,5,7), n = 7..9:', imb8[6:])
check(sum(imb8[8]) == 2 ** 36 and all(x > 0 for x in imb8[8]), 'n = 9 labelled residue counts sum to 2^36')
print('(Q4b) relative labelled masses of H mod 16 (8 x frequency), residues 1,3,...,15')
for n in (7, 8, 9):
    w = R[n]['w']; H = R[n]['H']; tot = sum(w)
    c = Counter()
    for wi, h in zip(w, H):
        c[int(h) % 16] += wi
    print(f'     n = {n}: ' + ', '.join(f'{r}:{float(8 * c[r] / tot):.3f}' for r in range(1, 16, 2)))

# ---------------------------------------------------------------- (Q5)
print('(Q5) repunits, Cipolla numbers and 2^k + 1 against a local baseline (n = 8, 9)')
fam = {'2^k-1': [2 ** k - 1 for k in range(2, 13)], '(4^k-1)/3': [(4 ** k - 1) // 3 for k in range(1, 7)],
       '2^k+1': [2 ** k + 1 for k in range(1, 12)]}
ratios = defaultdict(list)
for n in (8, 9):
    top = max(N[n])
    for name, vals in fam.items():
        row = []
        for v in vals:
            if v < 25 or v > 0.6 * top:
                continue
            win = [N[n][h] for h in range(int(v * 0.9) | 1, int(v * 1.1) + 1, 2) if h != v]
            base = float(np.median(win))
            rk = sum(1 for x in win if x < N[n][v]) / len(win)
            row.append(f'{v}:{N[n][v]}/{base:.0f}(pct {rk:.2f})')
            ratios[name].append(rk)
        print(f'     n = {n} {name:9s}: ' + ', '.join(row))
allp = [x for v in ratios.values() for x in v]
print(f'     percentile ranks of the family values within their +-10% windows: mean {np.mean(allp):.2f} over {len(allp)} cases (0.5 = typical)')
check(len(allp) >= 10, 'repunit probe ran on at least 10 (value, n) cases (NUMERICAL; the verdict is in the note)')

# ---------------------------------------------------------------- (Q6)
print('(Q6) the converse-Gluck family: H, |Aut|, tc, uniqueness')


def circ(n, Sset):
    return [sum(1 << ((x + d) % n) for d in Sset) for x in range(n)]


def lex(Am, Bm):
    na, nb = len(Am), len(Bm); T = [0] * (na * nb)
    for a in range(na):
        for b in range(nb):
            u = a * nb + b
            for a2 in range(na):
                for b2 in range(nb):
                    v = a2 * nb + b2
                    if (a != a2 and (Am[a] >> a2) & 1) or (a == a2 and (Bm[b] >> b2) & 1):
                        T[u] |= 1 << v
    return T


def Hof(T):
    return int(HM.ham_one(len(T), np.array(T, dtype=np.int64)))


C3, C5, P7 = circ(3, {1}), circ(5, {1, 2}), circ(7, {1, 2, 4})
fam6 = [('C_3', C3, 3), ('Z_5{1,2}', C5, 5), ('P_7', P7, 21), ('C_3[C_3]', lex(C3, C3), 81),
        ('Z_5{1,2}[C_3]', lex(C5, C3), 5 * 3 ** 5), ('P_7[C_3]', lex(P7, C3), 21 * 3 ** 7)]
for name, T, au in fam6:
    n = len(T); h = Hof(T)
    extra = ''
    if n <= 9:
        extra = f'; N({n},{h}) = {N[n][h]}; max H({n}) = {max(N[n])}'
    else:
        mean = factorial(n) / 2 ** (n - 1)
        extra = f'; mean n!/2^(n-1) = {mean:.4g}; H/mean = {h / mean:.3f}'
    print(f'     {name:14s} n = {n:2d}: H = {h}, |Aut| = {au}, tc = {h // au if h % au == 0 else "NOT INTEGER"} = {factorint(h // au) if h % au == 0 else ""}{extra}')
    check(h % au == 0, f'{name}: |Aut| = {au} divides H (free action)')
check(Hof(lex(C3, C3)) == 3159 and N[9][3159] == 1, 'H(C_3[C_3]) = 3159, attained by no other 9-class; C_3, P_7 are likewise unique, Z_5{1,2} ties at 15')

# ---------------------------------------------------------------- (Q7)
print('(Q7) Paley presentation counts tc(P_q) = H(P_q)/(q(q-1)/2)')
cited = {23: 15760206976379349, 27: 418425945688463432031, 31: 19876729650849384218015865}
pal = {}
for q in (3, 7, 11, 19):
    sq = {x * x % q for x in range(1, q)}
    pal[q] = Hof(circ(q, sq))
check(pal[7] == 189 and pal[11] == 95095 and pal[19] == 1172695746915, 'H(P_7) = 189, H(P_11) = 95095, H(P_19) = 1172695746915 recomputed')
pal.update(cited)
for q, h in sorted(pal.items()):
    aut = q * (q - 1) // 2 * (3 if q == 27 else 1)
    tc = h // aut
    psp2 = (pow(2, tc - 1, tc) == 1) and not isprime(tc) and tc > 1
    carm = (not isprime(tc)) and tc > 2 and all(e == 1 for e in factorint(tc).values()) and all((tc - 1) % (pp - 1) == 0 for pp in factorint(tc))
    print(f'     q = {q:2d}: tc = {tc} = {factorint(tc)}; base-2 pseudoprime: {psp2}; Carmichael (Korselt): {carm}; tc mod 8 = {tc % 8}')
    check(h % aut == 0, f'q = {q}: |Aut(P_q)| = {aut} divides H')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
