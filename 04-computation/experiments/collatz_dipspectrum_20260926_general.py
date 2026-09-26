#!/usr/bin/env python3
"""collatz_dipspectrum_20260926_general.py -- control for the general dip-spectrum theorem
(session collatz-exponent-atlas-20260926, opus; Theorem 4 of the dip-spectrum note, THM-4487).

For a Conway/Matthews-Watts map g(x) = (p_i x + q_i)/m on x = i mod m (p_i coprime to m), Theorem 4
says that the exponent of D_g(X, gamma) = #{n <= X : g^j(n) >= n^gamma for 0 <= j <= floor(log_m n)} is

    E_g(gamma) = max { H_m(pi) : pi a law on Z/m, sum_i pi_i log_m (p_i/m) >= gamma - 1 },

the constrained maximum entropy of the multiplier law (base-m entropy), for every 0 < gamma <= 1 when
max p_i > m (the carries are handled multiplicatively: no additive carry condition is needed). The maximiser is the tilted law pi_i ~ (p_i/m)^lambda when the constraint
is active, and uniform (E = 1) when the uniform law already satisfies it; that boundary,
sum_i (1/m) log_m(p_i/m) = gamma - 1, is Korec's threshold for 3n+1 (gamma = log_4 3).

Three checks, for the m = 3 map g = x/3, (2x+1)/3, (4x+1)/3 (multipliers 1/3, 2/3, 4/3; gamma_0 = log_3 2):
 (1) E_g(gamma) by a 1-parameter search over the tilt lambda (and, for 3n+1, agreement with
     h(max(1/2, gamma/log_2 3)) to five decimals);
 (2) exact orbit counts D_g(3^T - 1, gamma) for T <= tmax, next to the carry-free word model
     Dw = #{n : prod_(l<j) a_(i_l) * n >= n^gamma for all j <= t} (same residue words, carries dropped):
     D/Dw ~ 1 is the carry claim of the proof;
 (3) the exact number W_t(gamma) of residue words of length t with every prefix product >= 3^(t(gamma-1)),
     by dynamic programming on (position, #1 + 2 #2) for t up to 400: log_3 W_t / t -> E_g(gamma) slowly
     (finite-size effect of the multinomial), which is why the counts at X = 3^12 sit far below X^E_g.
Usage: python3 collatz_dipspectrum_20260926_general.py [tmax=12]
"""
import math, sys

LOG3_2 = math.log(2) / math.log(3)


def H_m(pi, m):
    return -sum(p * math.log(p) for p in pi if p > 0) / math.log(m)


def E_g(m, ps, gamma, grid=40000):
    """max entropy over tilted laws pi_i ~ a_i^lambda meeting sum pi_i log_m a_i >= gamma - 1."""
    la = [math.log(p / m) for p in ps]
    if sum(la) / m / math.log(m) >= gamma - 1:
        return 1.0
    best = 0.0
    for k in range(grid + 1):
        lam = 20.0 * k / grid
        w = [math.exp(lam * x) for x in la]
        Z = sum(w)
        pi = [x / Z for x in w]
        if sum(pi[i] * la[i] for i in range(m)) / math.log(m) >= gamma - 1:
            best = max(best, H_m(pi, m))
    return best


P3 = {0: (1, 0), 1: (2, 1), 2: (4, 1)}
LA3 = {0: math.log(1 / 3), 1: math.log(2 / 3), 2: math.log(4 / 3)}


def orbit_and_model_exponents(n, t):
    """(min_j log g^j(n)/log n, min_j log(M_j n)/log n) over 1 <= j <= t."""
    x = n
    mn = n
    logM = 0.0
    logMmin = 0.0
    for _ in range(t):
        r = x % 3
        p, q = P3[r]
        x = (p * x + q) // 3
        if x < mn:
            mn = x
        logM += LA3[r]
        if logM < logMmin:
            logMmin = logM
    ln = math.log(n)
    return math.log(mn) / ln, 1.0 + logMmin / ln


def words_dp(t, gamma):
    """# words in {0,1,2}^t with s_j log_3 2 - j >= t (gamma - 1) for all 1 <= j <= t, s_j = #1 + 2 #2 in prefix j."""
    thr = t * (gamma - 1)
    cur = {0: 1}
    for j in range(1, t + 1):
        nxt = {}
        for s, c in cur.items():
            for step in (0, 1, 2):
                s2 = s + step
                if s2 * LOG3_2 - j >= thr - 1e-12:
                    nxt[s2] = nxt.get(s2, 0) + c
        cur = nxt
    return sum(cur.values())


def main():
    tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("== (1) 3n+1 check: E_g(gamma) for m = 2, p = (1, 3) versus h(max(1/2, gamma/log_2 3)) ==")
    h = lambda p: 0.0 if p <= 0 or p >= 1 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)
    for gamma in (0.75, 0.7925, 0.85, 0.91, 0.97, 1.0):
        print("   gamma=%.4f: E_g = %.5f   h(rho) = %.5f" % (gamma, E_g(2, (1, 3), gamma), h(max(0.5, gamma / math.log2(3)))))
    ps = (1, 2, 4)
    thr = 1 + sum(math.log(p / 3) for p in ps) / 3 / math.log(3)
    print("== m = 3 map g = x/3, (2x+1)/3, (4x+1)/3: additive-carry threshold log_3(4/3) = %.5f (not needed); uniform-law threshold gamma_0 = %.5f = log_3 2 ==" % (math.log(4 / 3) / math.log(3), thr))
    gammas = [0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.97, 1.0]
    Eg = [E_g(3, ps, g) for g in gammas]
    print("   gamma:      " + "  ".join("%7.3f" % g for g in gammas))
    print("   E_g(gamma): " + "  ".join("%7.4f" % e for e in Eg))
    print("== (3) exact word counts W_t(gamma) (carry-free, DP): log_3 W_t / t ==")
    for t in (12, 14, 20, 30, 50, 100, 200, 400):
        print("   t=%3d:      " % t + "  ".join("%7.4f" % (math.log(words_dp(t, g)) / (t * math.log(3)) if words_dp(t, g) > 0 else float('nan')) for g in gammas))
    print("== (2) exact orbit counts D(3^T - 1, gamma) [top row] and carry-free model Dw [second row], with log D / log X ==")
    blk = {t: [0] * len(gammas) for t in range(4, tmax + 1)}
    blkw = {t: [0] * len(gammas) for t in range(4, tmax + 1)}
    for t in range(4, tmax + 1):
        for n in range(3 ** t, 3 ** (t + 1)):
            e, ew = orbit_and_model_exponents(n, t)
            for gi, g in enumerate(gammas):
                if e >= g - 1e-12:
                    blk[t][gi] += 1
                if ew >= g - 1e-12:
                    blkw[t][gi] += 1
    cum = [0] * len(gammas)
    cumw = [0] * len(gammas)
    rows = {}
    for t in range(4, tmax + 1):
        for gi in range(len(gammas)):
            cum[gi] += blk[t][gi]
            cumw[gi] += blkw[t][gi]
        rows[t + 1] = (list(cum), list(cumw))
    for T in range(6, tmax + 2):
        c, cw = rows[T]
        print("   X=3^%2d D : " % T + "  ".join("%9d" % v for v in c))
        print("          Dw: " + "  ".join("%9d" % v for v in cw))
        print("     logD/logX:" + "  ".join("%7.4f" % (math.log(v) / (T * math.log(3)) if v > 0 else float('nan')) for v in c))
    for T0 in range(6, tmax - 2):
        T1 = T0 + 4
        if T1 in rows:
            c0, c1 = rows[T0][0], rows[T1][0]
            print("   slope 3^%2d->3^%2d: " % (T0, T1) + "  ".join("%7.4f" % (math.log(c1[gi] / c0[gi]) / (4 * math.log(3)) if c0[gi] > 0 and c1[gi] > 0 else float('nan')) for gi in range(len(gammas))))


if __name__ == '__main__':
    main()
