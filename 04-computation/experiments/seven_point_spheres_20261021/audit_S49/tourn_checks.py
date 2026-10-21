"""THM-4743 audit: (a) Paley q <= 19 with groups computed from the tournament; (b) circulant census odd n <= 15;
(c) primality of C3, T5, P7; (d) cheap HYP-9353 counterexample hunt (Z3^2 Cayley, unequal substitutions, stacks)."""
import sys, time, itertools
import numpy as np
from tlib import *

t0 = time.time()
# ---------------- (a) Paley
print('(a) Paley P_q, groups by backtracking (no affine assumption)')
for q in (3, 7, 11, 19):
    sq = {x * x % q for x in range(1, q)}
    T = circ(q, sq)
    aut, ca, anti, cn = group_pm(T, q, 100000)
    regA = count_regular(q, nonid(aut), ca - 1)
    els = np.vstack([nonid(aut), anti]) if cn else nonid(aut)
    regPM = count_regular(q, els, els.shape[0])
    print(f'   q={q}: |Aut|={ca} |anti|={cn}; regular subsets Aut {regA} (orbits {regA // ca}), Aut+- {regPM} (orbits {regPM // (ca + cn)})', flush=True)

# ---------------- (c) primality
def is_prime_tourn(T):
    n = len(T)
    for M in range(1, 1 << n):
        k = bin(M).count('1')
        if k < 2 or k == n:
            continue
        ok = True
        for x in range(n):
            if (M >> x) & 1:
                continue
            outs = int(T[x]) & M
            if outs != 0 and outs != M:
                ok = False; break
        if ok:
            return False
    return True
print('(c) prime:', {nm: is_prime_tourn(T) for nm, T in (('C3', circ(3, {1})), ('T5', circ(5, {1, 2})), ('P7', circ(7, {1, 2, 4})))})

# ---------------- (b) circulants
print('(b) circulant tournaments, odd n <= 15 (exhaustive subsets, full Aut+-)')
C3 = circ(3, {1}); T5 = circ(5, {1, 2}); P7 = circ(7, {1, 2, 4})
fam = {3: [C3], 5: [T5], 7: [P7], 9: [lexprod(C3, C3)], 15: [lexprod(T5, C3)]}
for n in range(3, 16, 2):
    h = (n - 1) // 2
    exc = []; reps = []
    for bits in range(1 << h):
        S = {d if (bits >> (d - 1)) & 1 else n - d for d in range(1, h + 1)}
        T = circ(n, S)
        aut, ca, anti, cn = group_pm(T, n, 200000)
        els = np.vstack([nonid(aut), anti])
        if count_regular(n, els, els.shape[0]) == 0:
            exc.append(sorted(S))
            if not any(iso_exists(T, R) for R in reps):
                reps.append(T)
    infam = [any(iso_exists(R, F) for F in fam.get(n, [])) for R in reps]
    print(f'   n={n}: {1 << h} connection sets, exceptional {len(exc)}, iso classes {len(reps)}, in family {infam}  [{time.time() - t0:.0f}s]', flush=True)

# ---------------- (d) counterexample hunt
print('(d) HYP-9353 hunt')
rng = np.random.default_rng(9353)
cands = {}
# Cayley tournaments on Z3^2
Z32 = [(a, b) for a in range(3) for b in range(3)]
idx = {g: i for i, g in enumerate(Z32)}
pairs = [((1, 0), (2, 0)), ((0, 1), (0, 2)), ((1, 1), (2, 2)), ((1, 2), (2, 1))]
for bits in range(16):
    Sset = [pairs[i][(bits >> i) & 1] for i in range(4)]
    rows = []
    for g in Z32:
        r = 0
        for s in Sset:
            r |= 1 << idx[((g[0] + s[0]) % 3, (g[1] + s[1]) % 3)]
        rows.append(r)
    cands[f'Cay(Z3^2,{bits})'] = from_rows(rows)
K1 = from_rows([0]); TT2 = from_rows([2, 0])
C33 = lexprod(C3, C3)
cands['C3[C3,C3,T5] (n=11)'] = subst(C3, [C3, C3, T5])
cands['C3[C3,T5,T5] (n=13)'] = subst(C3, [C3, T5, T5])
cands['C3[C3,C3,C3[C3]] (n=15)'] = subst(C3, [C3, C3, C33])
cands['C3[C3[C3],C3[C3],C3] (n=21)'] = subst(C3, [C33, C33, C3])
cands['C3[K1,K1,C3] (n=5)'] = subst(C3, [K1, K1, C3])
cands['T5[C3,C3,C3,C3,K1] (n=13)'] = subst(T5, [C3, C3, C3, C3, K1])
cands['T5[C3,K1,K1,K1,K1] (n=7)'] = subst(T5, [C3, K1, K1, K1, K1])
cands['P7[C3,K1x6] (n=9)'] = subst(P7, [C3] + [K1] * 6)
cands['P7[C3,C3,K1,C3,K1,K1,K1] (n=13)'] = subst(P7, [C3, C3, K1, C3, K1, K1, K1])     # C3 on the Fano line {0,1,3}
cands['TT2[C3] stack (n=6)'] = lexprod(TT2, C3)
cands['C3=>C3=>C3 (n=9)'] = subst(from_rows([6, 4, 0]), [C3, C3, C3])
cands['C3=>T5=>C3 (n=11)'] = subst(from_rows([6, 4, 0]), [C3, T5, C3])
cands['P7=>P7conv? (n=14)'] = subst(TT2, [P7, circ(7, {3, 5, 6})])
cands['C3[T5] (n=15)'] = lexprod(C3, T5)
cands['T5[C3] (n=15)'] = lexprod(T5, C3)
cands['C3[C3[C3]] (n=27)'] = lexprod(C3, C33)
for nm, T in cands.items():
    assert is_tournament(T), nm
    n = len(T)
    S, t = random_regular_search(T, 3000, rng)
    print(f'   {nm}: n={n}; regular Aut+- subset found: {S is not None} (after {t} random tries)', flush=True)
print(f'done {time.time() - t0:.0f}s')
