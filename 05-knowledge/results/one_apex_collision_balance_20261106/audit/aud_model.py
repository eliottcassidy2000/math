"""Independent audit of the HYP-9482 Poisson-pair model constants (auditor's own derivation and code).

Model: on Y_n, g(v)/2^n -> Z = XY, X, Y i.i.d. uniform on [1/2, 1).  Converse classes ~ Poisson with mean
mu(t) = f(t)/4 per odd value at t = N/2^n, f the density of Z.  Odd values per unit t: 2^(n-1) = 2|Y|.

Checks (mpmath, 30 digits):
  (1) f(t) = 4 ln(4t) on [1/4,1/2], 4 ln(1/t) on [1/2,1] by direct convolution; int f = 1; int f^2 = 24 - 32 ln 2.
  (2) E^Y/|Y| -> 2 + 8 int mu^2 = 14 - 16 ln 2; non-converse part 12 - 16 ln 2; smoothed slope 12 - 16 ln 2.
  (3) values/|Y| -> 2 int (1 - e^-mu) = 3/4 - (ln 2)/2.
  (4) values with exactly q classes per |Y|: closed form vs direct numerical integration of 2 int e^-mu mu^q/q!.
  (5) cross layers: consecutive layers j, j-1 overlap on t in [1/4,1/2] (t = N/2^j); cross part per |X_n| ->
      4 * 2 * int_{1/4}^{1/2} 4 mu(t) mu(2t) dt ... = 12 ln 2 - 8, and E/|X| -> 6 - 4 ln 2.
  (6) tail at |Y| = 2^24: expected numbers of values of multiplicity 12, 14, 16, 18, 20 and
      P(max multiplicity >= 18) = 1 - exp(-sum_{q>=9} N_q) (values independent under the model).
  (7) Monte Carlo of the model itself at |Y| = 2^24 (n = 26): replicates of the multiplicity histogram,
      giving the distribution of the maximum and of the counts at 12..18 (100 replicates).
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 30
L2 = mp.log(2)


def f_conv(t):
    # density of XY, X,Y ~ U[1/2,1): f(t) = int 2 * 2 * [t/x in [1/2,1)] / x dx over x in [1/2,1)
    lo = max(mp.mpf(1) / 2, t)          # t/x < 1  <=> x > t
    hi = min(mp.mpf(1), 2 * t)          # t/x >= 1/2 <=> x <= 2t
    if hi <= lo:
        return mp.mpf(0)
    return 4 * mp.log(hi / lo)


def f_closed(t):
    if t < mp.mpf(1) / 4 or t > 1:
        return mp.mpf(0)
    return 4 * mp.log(4 * t) if t <= mp.mpf(1) / 2 else 4 * mp.log(1 / t)


def main():
    out = []
    md = max(abs(f_conv(t) - f_closed(t)) for t in [mp.mpf(k) / 997 for k in range(250, 997)])
    out.append('(1) max |convolution - closed form| on a grid of [1/4,1): %s' % mp.nstr(md, 5))
    I1 = mp.quad(f_closed, [0.25, 0.5, 1])
    I2 = mp.quad(lambda t: f_closed(t) ** 2, [0.25, 0.5, 1])
    out.append('(1) int f = %s ; int f^2 = %s ; 24 - 32 ln 2 = %s' % (mp.nstr(I1, 20), mp.nstr(I2, 20), mp.nstr(24 - 32 * L2, 20)))
    mu = lambda t: f_closed(t) / 4
    Imu2 = mp.quad(lambda t: mu(t) ** 2, [0.25, 0.5, 1])
    EY = 2 + 8 * Imu2
    out.append('(2) E^Y/|Y| = 2 + 8 int mu^2 = %s ; 14 - 16 ln 2 = %s ; non-converse 8 int mu^2 = %s ; 12 - 16 ln 2 = %s'
               % (mp.nstr(EY, 15), mp.nstr(14 - 16 * L2, 15), mp.nstr(8 * Imu2, 15), mp.nstr(12 - 16 * L2, 15)))
    out.append('(2) smoothed slope: each nonzero even offset d contributes 2|Y| * int 4 mu(t)^2 dt per |Y| (offsets +-2..+-K: K of them) -> '
               '8 int mu^2 = %s' % mp.nstr(8 * Imu2, 15))
    vals = 2 * mp.quad(lambda t: 1 - mp.e ** (-mu(t)), [0.25, 0.5, 1])
    out.append('(3) values/|Y| = %s ; 3/4 - ln2/2 = %s' % (mp.nstr(vals, 15), mp.nstr(mp.mpf(3) / 4 - L2 / 2, 15)))
    tot = 0
    fr = {}
    for q in range(1, 13):
        closed = 2 * (L2 ** (q + 1) / (4 * mp.factorial(q + 1)) + mp.quad(lambda u: mp.e ** (-2 * u) * u ** q / mp.factorial(q), [0, L2]))
        direct = 2 * mp.quad(lambda t: mp.e ** (-mu(t)) * mu(t) ** q / mp.factorial(q), [0.25, 0.5, 1])
        fr[q] = closed
        tot += closed * (2 * q) ** 2
        out.append('(4) q=%2d (mult %2d): closed %s  direct %s' % (q, 2 * q, mp.nstr(closed, 12), mp.nstr(direct, 12)))
    out.append('(4) energy from the histogram sum_q (2q)^2 frac_q = %s (should be 14 - 16 ln 2 = %s)'
               % (mp.nstr(tot, 12), mp.nstr(14 - 16 * L2, 12)))
    # (5) cross layers
    cross_pair = mp.quad(lambda t: 4 * mu(t) * mu(2 * t), [0.25, 0.5])   # E[r_j r_{j-1}] per odd value, integrated in t
    # ordered cross pairs between layers j and j-1: 2 * (2^(j-1) odd values per unit t) * cross_pair = 2^j cross_pair
    # sum over j <= n: ~ 2^(n+1) cross_pair; per |X_n| = 2^(n-1): 4 cross_pair
    out.append('(5) cross-layer per |X_n| = 4 * int_{1/4}^{1/2} 4 mu(t) mu(2t) dt = %s ; 12 ln 2 - 8 = %s'
               % (mp.nstr(4 * cross_pair, 15), mp.nstr(12 * L2 - 8, 15)))
    out.append('(5) E/|X| = (14 - 16 ln 2) + (12 ln 2 - 8) = %s ; 6 - 4 ln 2 = %s'
               % (mp.nstr(14 - 16 * L2 + 4 * cross_pair, 15), mp.nstr(6 - 4 * L2, 15)))
    out.append('(5) note: 4 - 4 ln 2 - (12 - 16 ln 2) = %s' % mp.nstr(4 - 4 * L2 - (12 - 16 * L2), 15))
    # (6) tail at |Y| = 2^24
    Y = mp.mpf(2) ** 24
    for q in (6, 7, 8, 9, 10):
        out.append('(6) |Y| = 2^24: expected number of values with multiplicity %d = %s' % (2 * q, mp.nstr(fr[q] * Y, 6)))
    lam18 = sum(fr[q] for q in range(9, 13)) * Y
    out.append('(6) expected number with multiplicity >= 18 = %s ; P(max >= 18) = 1 - exp(-that) = %s'
               % (mp.nstr(lam18, 6), mp.nstr(1 - mp.e ** (-lam18), 6)))
    lam16 = sum(fr[q] for q in range(8, 13)) * Y
    out.append('(6) expected number with multiplicity >= 16 = %s ; P(max >= 16) = %s' % (mp.nstr(lam16, 6), mp.nstr(1 - mp.e ** (-lam16), 6)))
    for (nn, YY) in ((22, 2 ** 20), (24, 2 ** 22)):
        l18 = sum(fr[q] for q in range(9, 13)) * YY
        out.append('(6) n=%d (|Y| = 2^%d): P(max multiplicity >= 18) = %s' % (nn, int(math.log2(YY)), mp.nstr(1 - mp.e ** (-l18), 4)))
    out.append('(6) observed at n = 26: 356, 30, 3, 1 values of multiplicity 12, 14, 16, 18; ratios to model: %s'
               % ', '.join('%.2f' % (o / float(fr[q] * Y)) for o, q in ((356, 6), (30, 7), (3, 8), (1, 9))))
    Y28 = mp.mpf(2) ** 26
    out.append('(6) |Y| = 2^26 (n = 28): model counts at multiplicity 12, 14, 16, 18 = %s; >= 18: %s; >= 20: %s'
               % (', '.join(mp.nstr(fr[q] * Y28, 5) for q in (6, 7, 8, 9)),
                  mp.nstr(sum(fr[q] for q in range(9, 13)) * Y28, 4), mp.nstr(sum(fr[q] for q in range(10, 13)) * Y28, 4)))
    out.append('(6) observed at n = 28 (audit run): 1306, 130, 10, 0 values of multiplicity 12, 14, 16, 18; ratios to model: %s'
               % ', '.join('%.2f' % (o / float(fr[q] * Y28)) for o, q in ((1306, 6), (130, 7), (10, 8), (0, 9))))
    # (7) Monte Carlo of the model at n = 26
    rng = np.random.default_rng(20261106)
    n = 26
    odd = np.arange(2 ** (n - 2) + 1, 2 ** n, 2, dtype=np.float64)   # odd values in (2^(n-2), 2^n)
    t = odd / 2.0 ** n
    muv = np.where(t <= 0.5, np.log(4 * t), np.log(1 / t))
    reps = 100
    maxes = []
    cnt = np.zeros((reps, 10))
    for r in range(reps):
        K = rng.poisson(muv)
        mult = 2 * K
        maxes.append(int(mult.max()))
        for i, mm in enumerate((12, 14, 16, 18)):
            cnt[r, i] = (mult == mm).sum()
    maxes = np.array(maxes)
    out.append('(7) model Monte Carlo n = 26, %d replicates: P(max >= 18) = %.3f, P(max >= 20) = %.3f; mean counts at 12,14,16,18 = %s'
               % (reps, (maxes >= 18).mean(), (maxes >= 20).mean(), ', '.join('%.2f' % x for x in cnt[:, :4].mean(axis=0))))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
