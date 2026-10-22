"""Audit S50 - a6: numerical claims of the results note sections 3 and 5 and of HYP-9355's evidence.
  (M1) Paley: H(P_q), q = 3, 7, 11, 19 by my DP; q = 23, 27, 31 from paley_h31_exact_rate_and_apery_wheel_20261014.md;
       tc under the affine normalisation 2H/(q(q-1)) and under |Aut| (q = 27: Frobenius), factorisations, base-2
       Fermat pseudoprimality, Korselt, residues mod 8;
  (M2) converse-Gluck lifts: H(Z_5{1,2}[C_3]), H(P_7[C_3]) by my DP, |Aut| divisibility, ratios to n!/2^(n-1);
       H(Z_9{1,2,3,5}), H(C_3[C_3]); order-9 maximum and uniqueness of 3159 from the mac-mini gentourng census (data);
  (M3) HYP-9355 evidence: labelled fractions with p | H (p = 3,5,7,11,13) and mod-8 / mod-16 masses:
       n = 7 by enumerating all 2^21 labelled tournaments with my DP; n = 8, 9 from the mac-mini census (n!/|Aut| weights).
       Test of the HYP-9355 title's monotone 2-adic claim at n = 7, 8, 9 (mod 16)."""
import os, time, math
from math import factorial
from collections import Counter
import numpy as np
import numba as nb
from sympy import factorint, isprime
import alib as L

T0 = time.time()
FAIL = []
HERE = os.path.dirname(os.path.abspath(__file__))
MM = os.path.join(os.path.dirname(os.path.dirname(HERE)), 'tournament_pseudoprimes_20261021')


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def circ(n, S):
    return [sum(1 << ((x + d) % n) for d in S) for x in range(n)]


def lex(A, B):
    na, nb_ = len(A), len(B)
    T = [0] * (na * nb_)
    for a in range(na):
        for b in range(nb_):
            for a2 in range(na):
                for b2 in range(nb_):
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (B[b] >> b2) & 1):
                        T[a * nb_ + b] |= 1 << (a2 * nb_ + b2)
    return T


def psp2(t):
    return t > 1 and not isprime(t) and pow(2, t - 1, t) == 1


def korselt(t):
    if t < 3 or isprime(t):
        return False
    f = factorint(t)
    return all(e == 1 for e in f.values()) and all((t - 1) % (p - 1) == 0 for p in f)


# ------------------------------------------------------------------ (M1)
Hp = {}
for q in (3, 7, 11, 19):
    sq = {x * x % q for x in range(1, q)}
    Hp[q] = L.H(circ(q, sq))
check(Hp == {3: 3, 7: 189, 11: 95095, 19: 1172695746915}, f'(M1) H(P_q) by my DP: {Hp}')
Hp.update({23: 15760206976379349, 27: 418425945688463432031, 31: 19876729650849384218015865})
note_tc = {3: 1, 7: 9, 11: 1729, 19: 6857869865, 23: 62293308207033, 31: 42745655163116955307561}
rows = []; ok = True; any_psp = []
for q, h in sorted(Hp.items()):
    aff = q * (q - 1) // 2
    assert h % aff == 0
    t_aff = h // aff
    t_aut = h // (aff * (3 if q == 27 else 1))
    rows.append(f'q={q}: 2H/(q(q-1)) = {t_aff} (mod 8: {t_aff % 8}, psp2 {psp2(t_aff)}, Korselt {korselt(t_aff)})'
                + (f'; H/|Aut| = {t_aut} (mod 8: {t_aut % 8}, psp2 {psp2(t_aut)})' if q == 27 else ''))
    if q in note_tc:
        ok &= note_tc[q] == t_aff
    ok &= t_aff % 8 == 1
    if psp2(t_aff) or (q == 27 and psp2(t_aut)):
        any_psp.append(q)
for r in rows:
    print('     ' + r)
check(ok and Hp[27] // 1053 == 397365570454381227, '(M1) the note\'s tc values (q = 3,7,11,19,23,31) and the q = 27 Frobenius value 397365570454381227 are correct; '
      'every AFFINE-normalised value is 1 (mod 8) (the printed q = 27 value is 3 (mod 8))')
check(any_psp == [11], f'(M1) base-2 pseudoprimes among the presentation counts (both normalisations at 27): q in {any_psp}; 1729 = 7*13*19 is Carmichael')

# ------------------------------------------------------------------ (M2)
C3, C5, P7 = circ(3, {1}), circ(5, {1, 2}), circ(7, {1, 2, 4})
h15 = L.H(lex(C5, C3))
check(h15 == 197030475 and h15 % 1215 == 0, f'(M2) H(Z_5{{1,2}}[C_3]) = {h15}, divisible by 5*3^5 = 1215')
code, aut15 = L.canon(lex(C5, C3), want_aut=True)
check(aut15 == 1215, f'(M2) |Aut(Z_5{{1,2}}[C_3])| = {aut15} (canon leaf count)')
h21 = L.H(lex(P7, C3))
check(h21 == 123156134945235 and h21 % 45927 == 0, f'(M2) H(P_7[C_3]) = {h21}, divisible by 21*3^7 = 45927')
r15 = h15 / (factorial(15) / 2 ** 14); r21 = h21 / (factorial(21) / 2 ** 20)
r19 = Hp[19] / (factorial(19) / 2 ** 18); r23 = Hp[23] / (factorial(23) / 2 ** 22)
check(round(r15, 2) == 2.47 and round(r21, 2) == 2.53 and round(r19, 2) == 2.53 and round(r23, 2) == 2.56,
      f'(M2) H/(n!/2^(n-1)): {r15:.4f} (n=15), {r21:.4f} (n=21); Paley {r19:.4f} (19), {r23:.4f} (23)')
hz9 = L.H(circ(9, {1, 2, 3, 5})); hcc = L.H(lex(C3, C3))
check(hz9 == 3357 and hcc == 3159, f'(M2) H(Z_9{{1,2,3,5}}) = {hz9}, H(C_3[C_3]) = {hcc}')


def read_mm(n):
    rows = []
    import gzip
    path = os.path.join(MM, f'census_n{n}.txt')
    fh = open(path) if os.path.exists(path) else gzip.open(path + '.gz', 'rt')
    for ln in fh:
        if ln.strip():
            H, C, s, A, AA, c3 = map(int, ln.split())
            rows.append((H, s, A, AA))
    return rows


mm9 = read_mm(9)
top = max(r[0] for r in mm9)
check(top == 3357 and [r[2] for r in mm9 if r[0] == top] == [9] and sum(1 for r in mm9 if r[0] == 3159) == 1,
      f'(M2) gentourng census n = 9 (mac-mini data): max H = {top}, attained once with |Aut| 9; H = 3159 attained by exactly one class')
mm7 = read_mm(7); mm5 = read_mm(5)
check(sum(1 for r in mm7 if r[0] == 189) == 1 and sum(1 for r in mm5 if r[0] == 15) == 2, '(M2) max H(7) = 189 unique; max H(5) = 15 attained twice')


# ------------------------------------------------------------------ (M3)
@nb.njit(parallel=True, cache=True)
def all_labelled(n, out):
    E = n * (n - 1) // 2
    tot = 1 << E
    for code in nb.prange(tot):
        adj = np.zeros(n, dtype=np.int64)
        e = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> e) & 1:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                e += 1
        out[code] = L.ham_dp(n, adj)


out7 = np.zeros(1 << 21, dtype=np.int64)
all_labelled(7, out7)
dist = {7: Counter(out7.tolist())}
for n in (8, 9):
    c = Counter()
    for H, s, A, AA in read_mm(n):
        c[H] += factorial(n) // A
    dist[n] = c
check(sum(dist[8].values()) == 2 ** 28 and sum(dist[9].values()) == 2 ** 36 and sum(dist[7].values()) == 2 ** 21,
      '(M3) labelled totals 2^21, 2^28, 2^36 (n = 7 enumerated; n = 8, 9 from the census with n!/|Aut| weights)')
note = {3: (0.3147, 0.3297, 0.3388), 5: (0.1681, 0.1806, 0.2005), 7: (0.1086, 0.1377, 0.1420), 11: (0.0884, 0.0882, 0.0919),
        13: (0.0985, 0.0756, 0.0764)}
ok = True
for p, vals in note.items():
    mine = tuple(round(sum(w for h, w in dist[n].items() if h % p == 0) / 2 ** (n * (n - 1) // 2), 4) for n in (7, 8, 9))
    ok &= mine == vals
check(ok, '(M3) HYP-9355 table of labelled fractions with p | H (n = 7, 8, 9) reproduced to 4 decimals')
m8 = {n: tuple(round(4 * sum(w for h, w in dist[n].items() if h % 8 == r) / 2 ** (n * (n - 1) // 2), 3) for r in (1, 3, 5, 7)) for n in (7, 8, 9)}
check(m8 == {7: (1.179, 0.985, 1.029, 0.807), 8: (1.064, 0.979, 1.019, 0.938), 9: (1.015, 0.992, 1.0, 0.993)},
      f'(M3) HYP-9355 mod-8 relative masses reproduced: {m8}')


def dist2(r):
    d = (r - 1) % 16
    if d == 0:
        return 0.0
    v = 0
    while d % 2 == 0:
        d //= 2; v += 1
    return 2.0 ** (-v)


for n in (7, 8, 9):
    tot = 2 ** (n * (n - 1) // 2)
    m16 = {r: 8 * sum(w for h, w in dist[n].items() if h % 16 == r) / tot for r in range(1, 16, 2)}
    groups = {}
    for r, v in m16.items():
        groups.setdefault(dist2(r), []).append(round(v, 3))
    keys = sorted(groups)
    mono = all(min(groups[keys[i]]) > max(groups[keys[i + 1]]) for i in range(len(keys) - 1))
    print(f'     n = {n}: mod-16 masses grouped by 2-adic distance |r-1|_2: ' +
          '; '.join(f'{k}: {groups[k]}' for k in keys) + f'  -> strictly decreasing with distance: {mono}')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
