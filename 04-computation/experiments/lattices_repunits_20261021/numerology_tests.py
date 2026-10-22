"""Numerology tests and known-value checks around 1728/1729, 24, Heegner j-values, exotic spheres and repunits
(session opus-2026-10-21-S52; note 05-knowledge/results/lattices_modular_repunits_20261021.md).

(J) j = E4^3/Delta = q^-1 + 744 + 196884 q + ...; 744 = 720 + 24 = 3 * 248 (theta_{E8^3}/eta^24: roots + rank =
    dim E8^3); j(i) = 1728 = 12^3; 1729 = 1728 + 1 = 12^3 + 1^3 = 10^3 + 9^3; 1/Delta = 1/eta^24 counts the 24
    transverse oscillators of the bosonic string (26 = 24 + 2); the 24 of q^(1/24) is |Mp_2(Z)^ab| = 24 (eta's
    multiplier eta(tau + 1) = e^(2 pi i/24) eta(tau) generates its character group), against |SL_2(Z)^ab| = 12.
    Heegner values j((1 + sqrt(-d))/2) for d = 3, 7, 11, 19, 43, 67, 163 are (negative) cubes; e^(pi sqrt 163) is
    within 7.5e-13 of 640320^3 + 744.  All KNOWN.
(P) H(P_q)/|Aut(P_q)| for the Paley tournaments with H known exactly (q = 3, 7, 11, 19, 23, 27, 31): 1, 9, 1729 at
    q = 3, 7, 11, so value - 1 = 8 m^3 (m = 0, 1, 6).  Test: does it continue?  (It fails at q = 19: NUMEROLOGY, the
    same q where the prior split-prime reading failed.)  Also H(P_q)/|Aff(P_q)| = 1 (mod 8) (HYP-9273) at every known q
    (|Aff| = q(q-1)/2; |Aut(P_27)| = 3 |Aff| adds the Frobenius).
(K) Kervaire-Milnor: |bP_{4m}| = 2^(2m-2) (2^(2m-1) - 1) num(4 |B_2m| / m) (the framing of MISTAKE-563), whose odd
    Mersenne factor 2^(2m-1) - 1 runs over the odd-exponent repunits 7, 31, 127, 511, 2047, 8191, ... ; the Euclid
    factor 2^(2m-2)(2^(2m-1) - 1) is a perfect number exactly when 2^(2m-1) - 1 is a Mersenne prime (m = 2, 3, 4, 7,
    9, 10 below 16).  3 never divides 2^(2m-1) - 1 (the arithmetic of THM-4749 (2): Phi_2(2) = 3 divides 2^odd + 1).
    Adams: im J in stem 4s - 1 is cyclic of order denom(B_2s/4s) = 24, 240, 504, 480, 264, 65520, which is the
    numerator of the q-coefficient -4s/B_2s of E_2s; the 691 of 65520/691 (the Leech/Niemeier theta series) is
    num(B_12/24) and divides |bP_24| (exotic 23-spheres).  KNOWN.
(H) Perfect binary 1-error-correcting codes exist exactly at the repunit lengths n = 2^r - 1 (Hamming bound
    (1 + n) | 2^n  <=>  n + 1 = 2^r); the Golay code: 1 + 23 + 253 + 1771 = 2^11.  KNOWN.
(Q) The projective-plane sizes q^2 + q + 1 = 111 in base q (q a prime power) against the strong H-spectrum below
    2100 (complete strong inventory as in S50/S51, 723,544 classes): only 7 = |PG(2,2)| and 21 = |PG(2,4)| are
    missing; every other size up to 1893 is attained within a factor 2 of its local median (31 = |PG(2,5)|: 4 classes
    against a median of 6).  The band has little power (it holds 98.6% of all odd d; audit_S52 regresses log S(d) on
    log d and 3|d, 9|d, 5|d and finds no projective-plane outlier).  KNOWN design side: a projective plane of order q
    extends to a 3-design only if (q + 2) | 12, i.e. q in {2, 4, 10}, and order 10 does not exist (Lam 1989), so
    {7, 21} = {111_2, 111_4} are exactly the point sets of the extendable planes, which extend to S(3,4,8) (8 points,
    extended Hamming -> E8) and S(5,8,24) (24 points, Golay -> Leech).  As an H statement this is NUMEROLOGY.
HYP-9273 (H(P_q)/|Aff| = 1 mod 8) is the S33/S34 hypothesis with its pre-registered q = 31 pass; it is re-checked here.
Prints ALL CHECKS PASSED."""
import os, sys, math, time
from fractions import Fraction
from collections import Counter
from multiprocessing import Pool
import mpmath
from sympy import bernoulli, factorint, isprime, Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h_spectrum_multiplicities_20261021', 'audit_S50'))

FAIL = []
T0 = time.time()
mpmath.mp.dps = 50


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def main():
    print('(J) j, 1728, 1729, Delta^-1, Heegner')
    N = 12

    def sig(k, n):
        return sum(d ** k for d in range(1, n + 1) if n % d == 0)
    E4 = [1] + [240 * sig(3, n) for n in range(1, N + 2)]
    pr = [1] + [0] * (N + 1)
    for n in range(1, N + 2):
        for _ in range(24):
            pr = [pr[i] - (pr[i - n] if i >= n else 0) for i in range(N + 2)]
    # 1/prod(1-q^n)^24 coefficients (q * (1/Delta)): p_24(n)
    inv = [1] + [0] * (N + 1)
    for n in range(1, N + 2):
        for _ in range(24):
            for i in range(n, N + 2):
                inv[i] += inv[i - n]
    E43 = [sum(E4[i] * sum(E4[j] * E4[k - i - j] for j in range(k - i + 1)) for i in range(k + 1)) for k in range(N + 2)]
    jq = [sum(E43[i] * inv[k - i] for i in range(k + 1)) for k in range(N + 2)]   # q * j(q)
    check(jq[:4] == [1, 744, 196884, 21493760] and jq[1] == E43[1] + inv[1] == 720 + 24 == 3 * 248,
          f'q j(q) = {jq[:4]} ... (744 = 720 roots of E8^3 + 24 = 3 * 248; 196884 = 196883 + 1)')
    snf = lambda rel: [int(x) for x in smith_normal_form(Matrix(rel), domain=ZZ).diagonal()]
    # SL_2(Z) = <S, T | S^4, S^2 = (ST)^3>; Mp_2(Z) = <S, T | S^8, S^2 = (ST)^3>; abelianise (s, t)
    check(snf([[1, 3], [4, 0]]) == [1, 12] and snf([[1, 3], [8, 0]]) == [1, 24],
          'SL_2(Z)^ab = Z/12, Mp_2(Z)^ab = Z/24 (Smith normal forms): the 24 of q^(1/24) is the order of eta\'s multiplier')
    check(inv[:6] == [1, 24, 324, 3200, 25650, 176256], f'1/eta^24 = q^-1 (1 + 24 q + 324 q^2 + 3200 q^3 + ...): {inv[:6]}')
    ji = mpmath.kleinj(mpmath.mpc(0, 1)) * 1728
    check(abs(ji - 1728) < 1e-30 and 1729 == 12 ** 3 + 1 == 10 ** 3 + 9 ** 3, 'j(i) = 1728 = 12^3; 1729 = 12^3 + 1^3 = 10^3 + 9^3')
    heeg = {3: 0, 7: -15, 11: -32, 19: -96, 43: -960, 67: -5280, 163: -640320}
    ok = True
    for d, c in heeg.items():
        tau = (1 + mpmath.sqrt(-d)) / 2
        jv = mpmath.kleinj(tau) * 1728
        ok &= abs(jv - c ** 3) < mpmath.mpf(10) ** -15 * max(1, abs(c) ** 3)
    ram = mpmath.e ** (mpmath.pi * mpmath.sqrt(163)) - (640320 ** 3 + 744)
    check(ok and abs(ram) < 1e-12, f'j((1+sqrt(-d))/2) = {[(d, c) for d, c in heeg.items()]} cubed; e^(pi sqrt 163) - 640320^3 - 744 = {mpmath.nstr(ram, 5)}')

    print('(P) H(P_q)/|Aut(P_q)| and H(P_q)/|Aff(P_q)|')
    H = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349,
         27: 418425945688463432031, 31: 19876729650849384218015865}
    aff = {q: q * (q - 1) // 2 for q in H}
    aut = dict(aff)
    aut[27] = 3 * aff[27]                      # x -> a x^sigma + b, sigma in Gal(F_27/F_3)
    rows = []
    cube_ok = []

    def icbrt(x):
        r = round(x ** (1 / 3)) if x > 0 else 0
        while r ** 3 > x:
            r -= 1
        while (r + 1) ** 3 <= x:
            r += 1
        return r
    for q in sorted(H):
        assert H[q] % aut[q] == 0
        t, ta = H[q] // aut[q], H[q] // aff[q]
        is_cube = (t - 1) % 8 == 0 and icbrt((t - 1) // 8) ** 3 == (t - 1) // 8
        cube_ok.append(is_cube)
        rows.append((q, t, ta))
        print(f'     q = {q:2d}: H/|Aut| = {t}; (H/|Aut| - 1)/8 a cube: {is_cube}; H/|Aff| mod 8 = {ta % 8}')
    check(cube_ok[:3] == [True, True, True] and not any(cube_ok[3:]), 'H/|Aut| - 1 = 8 m^3 at q = 3, 7, 11 (m = 0, 1, 6) and at no later q: NUMEROLOGY')
    check(all(ta % 8 == 1 for q, t, ta in rows), 'H(P_q)/|Aff(P_q)| = 1 (mod 8) at every known q, including q = 27 (HYP-9273)')

    print('(K) Kervaire-Milnor bP_{4m}')
    ref = {2: 28, 3: 992, 4: 8128, 5: 261632, 7: 67100672}
    perfect_m = []
    ok = True
    for m in range(2, 16):
        B = abs(Fraction(bernoulli(2 * m)))
        num = (4 * B / m).numerator
        bp = 2 ** (2 * m - 2) * (2 ** (2 * m - 1) - 1) * num
        mer = 2 ** (2 * m - 1) - 1
        if m in ref:
            ok &= bp == ref[m]
        euclid = 2 ** (2 * m - 2) * mer
        if isprime(mer):
            perfect_m.append(m)
        ok &= mer % 3 != 0
        print(f'     m = {m:2d}: |bP_{4 * m}| = {bp} = 2^{2 * m - 2} * {mer} * {num}; Mersenne factor {mer} = {dict(factorint(mer))}{" (prime: Euclid factor %d perfect)" % euclid if isprime(mer) else ""}')
    check(ok and perfect_m == [2, 3, 4, 7, 9, 10], f'reference orders 28, 992, 8128, 261632, 67100672 reproduced; perfect Euclid factors at m = {perfect_m} (2m - 1 = 3, 5, 7, 13, 17, 19); 3 never divides 2^(2m-1) - 1')
    imj = {1: 24, 2: 240, 3: 504, 4: 480, 5: 264, 6: 65520}
    ok = True
    for s_ in range(1, 7):
        B = Fraction(bernoulli(2 * s_))
        coef = -Fraction(4 * s_) / B                      # E_2s = 1 + coef sum sigma_(2s-1)(n) q^n
        ok &= (B / (4 * s_)).denominator == imj[s_] == abs(coef.numerator)
        print(f'     stem {4 * s_ - 1:2d}: |im J| = denom(B_{2 * s_}/{4 * s_}) = {(B / (4 * s_)).denominator}; E_{2 * s_} q-coefficient = {coef}')
    check(ok and (Fraction(bernoulli(12)) / 24).numerator == -691 and (2 ** 10 * 2047 * 691) % 691 == 0,
          'im J orders 24, 240, 504, 480, 264, 65520 = numerators of the Eisenstein coefficients -24, 240, -504, 480, -264, 65520/691; 691 = |num(B_12/24)| divides |bP_24| (Adams; KNOWN)')

    print('(H) perfect binary codes at repunit lengths')
    perfect_len = [n for n in range(1, 4100) if (2 ** n) % (1 + n) == 0]
    check(perfect_len == [2 ** r - 1 for r in range(1, 13)], f'(1 + n) | 2^n for n < 4100 exactly at n = {perfect_len} = 2^r - 1 (Hamming)')
    check(1 + 23 + math.comb(23, 2) + math.comb(23, 3) == 2 ** 11, 'Golay: 1 + 23 + 253 + 1771 = 2^11 (perfect, radius 3)')

    print('(Q) projective-plane sizes in the strong H-spectrum below 2100')
    import alib as L
    from b3_ear_census_625 import grow
    LIM = 2100
    inv_ = {3: {L.canon(L.chain((2,))): 3}}
    with Pool(12) as pool:
        n = 4
        while True:
            par = [c for c, h in inv_[n - 1].items() if h <= LIM]
            got, _ = grow(pool, n, par, cap=LIM) if par else ({}, None)
            got = {c: h for c, h in got.items() if h <= LIM}
            if not got:
                break
            inv_[n] = got
            n += 1
    S = Counter(h for dd in inv_.values() for h in dd.values())
    print(f'     strong inventory below {LIM}: {sum(S.values())} classes [{time.time() - T0:.0f}s]')
    check(sum(S.values()) == 723544, 'strong inventory below 2100 = 723,544 classes (as in S51 fan_uniqueness.out / repunit_degrees.out)')
    pp = [q for q in range(2, 46) if len(factorint(q)) == 1]
    sizes = [(q, q * q + q + 1) for q in pp if q * q + q + 1 < LIM]
    # local relative count: S(d) / median of S over odd d within +-10%
    rel = {}
    for q, d in sizes:
        win = sorted(S[x] for x in range(int(d * 0.9) | 1, int(d * 1.1) + 1, 2) if x != d)
        med = win[len(win) // 2] if win else 0
        rel[q] = (S[d], med)
        print(f'     q = {q:2d}: |PG(2,q)| = {d:4d}: S = {S[d]}, local median {med}')
    missing = [d for q, d in sizes if S[d] == 0]
    check(missing == [7, 21], f'missing projective-plane sizes: {missing} (q = 2, 4: the two H-spectrum gaps)')
    off = [q for q, (s_, med) in rel.items() if q not in (2, 4) and not (med and 0.5 * med <= s_ <= 2 * med)]
    check(off == [] and rel[5] == (4, 6), f'every other size within a factor 2 of its local median (31 = |PG(2,5)|: {rel[5][0]} vs {rel[5][1]}); the projective-plane reading stops at q = 2, 4: NUMEROLOGY')
    band = []
    for d in range(9, 1900, 2):
        if d == 21:
            continue
        win = sorted(S[x] for x in range(int(d * 0.9) | 1, int(d * 1.1) + 1, 2) if x != d)
        med = win[len(win) // 2] if win else 0
        band.append(bool(med) and 0.5 * med <= S[d] <= 2 * med)
    print(f'     base rate: {sum(band) / len(band):.3f} of all odd d in [9, 1899] (except 21) lie in the same factor-2 band: the band test has little power')
    ext = [q for q in range(2, 2000) if Fraction((q * q + q + 2) * (q * q + q + 1), q + 2).denominator == 1]
    check(ext == [2, 4, 10] and [q * q + q + 1 for q in ext] == [int('111', q) for q in ext] == [7, 21, 111]
          and Fraction(math.comb(24, 5), math.comb(8, 5)) == 759 and Fraction(math.comb(8, 3), math.comb(4, 3)) == 14,
          'extendable projective planes: (q + 2) | 12 forces q in {2, 4, 10} (q = 10 does not exist, Lam 1989); point sets 111_2 = 7, 111_4 = 21; extensions S(3,4,8) (14 blocks) and S(5,8,24) (759 octads): KNOWN; as H statement NUMEROLOGY')

    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')


if __name__ == '__main__':
    main()
