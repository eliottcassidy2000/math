#!/usr/bin/env python3
"""Sturmian test for the transversality foundry (collatz-procgen-20260922).

Bernstein's 2-adic number of a parity word, for Sturmian (mechanical) words,
for the maps 3x+1, 3x-1, 5x+1 and Mahler's ceil(3a/2), with three checks:

  (1) residue Phi(w) mod 2^N (N >= 10^4 bits) and the least positive
      representative (a positive integer below 2^N would have to equal it);
  (2) rational reconstruction: no rational a/b with |a|, b <= 2^((N-1)/2)-1
      is congruent to Phi(w) mod 2^N (Wang's bound);
  (3) repetition certificates (Theorem R of the note): for a prefix u v^e of w
      with lambda = |common prefix of w and u v^inf|, the rational
      Phi(u v^inf) = P/Q gives, IF Phi(w) = a/b, the bound
      max(|a|,b) >= 2^(lambda - log2(|P|+|Q|)) =: 2^gain.
      Theorem S (proved in the note) says gain -> infinity along the
      convergent denominators for every Sturmian word when mu < golden ratio.

Every quantity is an exact integer; the only floating point is in printing
log2 of exact integers.  Mechanical words are generated from a 520-bit
rational approximation of the slope with an explicit margin check, so each
bit printed is the exact bit of the real mechanical word.

Affine 2-adic shift maps are given as (m0, r0, m1, r1):
  T(x) = (m_e x + r_e)/2 for x = e mod 2, with m_e odd and r_e = e mod 2.
  Collatz 3x+1: (1,0,3,1);  3x-1: (1,0,3,-1);  5x+1: (1,0,5,1);
  Mahler ceil(3a/2) (THM-2228): (3,0,3,1).
"""
import math
import random
from fractions import Fraction

import mpmath

MAPS = {
    "3x+1": (1, 0, 3, 1),
    "3x-1": (1, 0, 3, -1),
    "5x+1": (1, 0, 5, 1),
    "mahler": (3, 0, 3, 1),
}

PREC_BITS = 520
mpmath.mp.prec = PREC_BITS + 40


# ------------------------------------------------------------------ slopes
def cf_value(prefix, period):
    """Value of [0; prefix, period, period, ...] (period nonempty -> quadratic)."""
    # solve the periodic tail t = [period..., t] exactly with mpmath
    def ev(terms, tail):
        v = tail
        for a in reversed(terms):
            v = a + 1 / v
        return v
    if period:
        # t = ev(period, t)  ->  fixed point iteration (contracting for positive CF)
        t = mpmath.mpf(period[0]) + 1
        for _ in range(4000):
            t_new = ev(period, t)
            if abs(t_new - t) < mpmath.mpf(2) ** (-(PREC_BITS + 20)):
                t = t_new
                break
            t = t_new
        tail = t
    else:
        tail = mpmath.inf
    v = tail
    for a in reversed(prefix):
        v = a + 1 / v
    return 1 / v


def slope_library():
    lib = {}
    lib["log_3 2 (critical q=3)"] = mpmath.log(2) / mpmath.log(3)
    lib["[0;1,1,1,(2)] quadratic, below crit."] = cf_value([1, 1, 1], [2])
    lib["[0;1,1,1,2,2,(3)] quadratic"] = cf_value([1, 1, 1, 2, 2], [3])
    lib["[0;1,1,1,2,2,3,1,5,2,23,2,(3)] quadr., just above"] = cf_value([1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2], [3])
    lib["[0;1,1,1,(3)] quadratic, above crit."] = cf_value([1, 1, 1], [3])
    lib["1/phi=[0;(1)] quadratic"] = (mpmath.sqrt(5) - 1) / 2
    lib["sqrt2/2=[0;1,(2)] quadratic"] = mpmath.sqrt(2) / 2
    lib["[0;1,(1,2)] quadratic, near 0.73"] = cf_value([1], [1, 2])
    lib["[0;1,(12)] quadratic, near 0.92"] = cf_value([1], [12])
    lib["log_5 2 (critical q=5)"] = mpmath.log(2) / mpmath.log(5)
    lib["[0;2,3,(9)] quadratic, near log_5 2"] = cf_value([2, 3], [9])
    lib["[0;2,(3)] quadratic, near 0.434"] = cf_value([2], [3])
    lib["sqrt2-1=[0;(2)] quadratic, 0.414"] = mpmath.sqrt(2) - 1
    lib["[0;3,(3)] quadratic, 0.303 (Mahler-safe range)"] = cf_value([3], [3])
    return lib


def cf_terms(x, n):
    out = []
    for _ in range(n):
        a = int(mpmath.floor(x))
        out.append(a)
        f = x - a
        if f == 0:
            break
        x = 1 / f
    return out


def convergents(alpha, n):
    """Denominators q_k of the convergents of alpha (alpha in (0,1))."""
    terms = cf_terms(alpha, n)
    qs = []
    q_prev, q = 0, 1  # q_{-1}=0, q_0=1 for [a0;...]
    for a in terms[1:]:
        q_prev, q = q, a * q + q_prev
        qs.append(q)
    return terms, qs


def to_fraction(x):
    """Rational approximation of an mpf with |error| < 2^-PREC_BITS."""
    return Fraction(int(mpmath.nint(x * (mpmath.mpf(2) ** PREC_BITS))), 1 << PREC_BITS)


# ------------------------------------------------------------------ words
def mechanical_word(alpha, rho, n, upper=False):
    """Exact bits s(i) = floor((i+1)a+r) - floor(i a + r) (ceil if upper), i < n.

    alpha, rho are mpf values; they are replaced by 520-bit rationals A, R and
    every i*A+R (i <= n) is checked to be farther than 2^-400 from an integer,
    except at i = 0 when rho itself is an exact integer (then exact anyway).
    The true slope differs from A by < 2^-520, so for i <= 10^6 the floor
    values of the true word and of the rational word coincide.
    """
    A = to_fraction(alpha)
    R = to_fraction(rho)
    num_a, den = A.numerator * (R.denominator), A.denominator * R.denominator
    num_r = R.numerator * A.denominator
    margin = den >> 400  # 2^-400 in units of 1/den
    fl = []
    val = num_r
    rho_is_int = (R.denominator == 1)
    for i in range(n + 1):
        q, rem = divmod(val, den)
        exact_zero = (i == 0 and rho_is_int)  # i*alpha+rho = rho is an exact integer
        if not exact_zero and not (margin < rem < den - margin):
            raise ValueError(f"margin check failed at i={i}")
        if upper:
            fl.append(q + (1 if rem else 0))
        else:
            fl.append(q)
        val += num_a
    return [fl[i + 1] - fl[i] for i in range(n)]


def critical_halving_word(rho, n_letters):
    """In-house discrepancy note words: K_j = floor(j*log2 3 + rho) - floor(rho),
    k_j = K_j - K_(j-1) in {1,2}; T-parity vector = concat of 1 0^(k_j - 1)."""
    alpha = mpmath.log(3) / mpmath.log(2)
    ks = mechanical_word(alpha - 1, rho, n_letters)  # k_j - 1 in {0,1}
    w = []
    for b in ks:
        w.append(1)
        if b:
            w.append(0)
    return w


def thue_morse(n):
    return [bin(i).count("1") & 1 for i in range(n)]


def random_word(n, density, seed):
    rnd = random.Random(seed)
    return [1 if rnd.random() < density else 0 for _ in range(n)]


def eventually_periodic(u, v, n):
    w = list(u)
    while len(w) < n:
        w.extend(v)
    return w[:n]


def balanced_random_critical(n_letters, seed, width=3.0):
    """Critical-slope bounded-discrepancy halving word that is NOT Sturmian:
    k_j in {1,2} chosen at random subject to 0 <= K_j - j*log2(3) < width."""
    rnd = random.Random(seed)
    a = math.log2(3)
    K, w = 0, []
    for j in range(1, n_letters + 1):
        opts = [k for k in (1, 2) if 0 <= K + k - j * a < width]
        k = rnd.choice(opts)
        K += k
        w.append(1)
        if k == 2:
            w.append(0)
    return w


# ------------------------------------------------------------------ 2-adic arithmetic
def affine_prefix(word, amap):
    """(M, R) with T^s(x) = (M x + R)/2^s along the word (s = len(word))."""
    m0, r0, m1, r1 = amap
    M, R = 1, 0
    p2 = 1
    for e in word:
        if e:
            M, R = m1 * M, m1 * R + r1 * p2
        else:
            M, R = m0 * M, m0 * R + r0 * p2
        p2 <<= 1
    return M, R


def phi_residue(word, amap, N):
    """Phi(w) mod 2^N from the first N letters."""
    M, R = affine_prefix(word[:N], amap)
    mod = 1 << N
    return (-R * pow(M, -1, mod)) % mod


def approximant(word, t, p, amap):
    """Exact P/Q = Phi(u v^inf) for u = word[:t], v = word[t:t+p]."""
    Mu, Ru = affine_prefix(word[:t], amap)
    Mv, Rv = affine_prefix(word[t:t + p], amap)
    D = (1 << p) - Mv
    P = (1 << t) * Rv - Ru * D
    Q = Mu * D
    # Phi(u v^inf) = P/Q; sanity: Q odd and nonzero
    assert Q % 2 != 0
    return P, Q


def common_prefix_len(word, t, p):
    """lambda = length of the common prefix of word and word[:t] (word[t:t+p])^inf."""
    n = len(word)
    k = t + p
    while k < n and word[k] == word[k - p]:
        k += 1
    return k  # may equal n (truncated)


def rat_recon(x, N):
    """Wang rational reconstruction mod 2^N with |a|, b <= floor(sqrt(2^(N-1)))."""
    M = 1 << N
    bound = math.isqrt(M >> 1)
    r0, r1 = M, x % M
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    a, b = r1, s1
    if b == 0 or abs(b) > bound:
        return None, bound
    if b < 0:
        a, b = -a, -b
    if b % 2 == 0 or (a - b * x) % M != 0:
        return None, bound
    return Fraction(a, b), bound


def v2(n):
    return (n & -n).bit_length() - 1 if n else math.inf


# ------------------------------------------------------------------ certificates
def repetition_certificates(word, amap, periods, n_clusters=4):
    """For each period p: candidates u = empty, and u = prefix ending just after each
    of the first defect clusters (maximal runs of positions i with w[i+p] != w[i]).
    These are options A and B of the proof of Theorem S.  Returns, per p, the best
    (p, gain, t, lambda, P, Q) among candidates whose repetition ends inside the window."""
    out = []
    n = len(word)
    for p in periods:
        if 3 * p > n:
            break
        starts = [0]
        i, found = 0, 0
        while i + p < n and found < n_clusters:
            if word[i] != word[i + p]:
                j = i
                while j + 1 + p < n and word[j + 1] != word[j + 1 + p]:
                    j += 1
                starts.append(j + 1)
                found += 1
                i = j + 1
                continue
            i += 1
        best = None
        for t in starts:
            if t + 2 * p >= n:
                continue
            lam = common_prefix_len(word, t, p)
            if lam >= n:
                continue  # repetition runs to the end of the window: lambda unknown
            P, Q = approximant(word, t, p, amap)
            gain = lam - math.log2(abs(P) + abs(Q))
            if best is None or gain > best[0]:
                best = (gain, t, lam, P, Q)
        if best:
            out.append((p,) + best)
    return out


def check_approximant_valuation(word, amap, t, p, lam, P, Q, N):
    """v_2(Phi(w) - P/Q) must equal lambda exactly (isometry), when lambda < N."""
    x = phi_residue(word, amap, N)
    mod = 1 << N
    diff = (x - P * pow(Q, -1, mod)) % mod
    return v2(diff) == lam


def mu_of(alpha, amap):
    m0, _, m1, _ = amap
    rate = (1 - alpha) * math.log2(m0) + alpha * math.log2(m1)
    return max(1.0, rate)


def parity_agreement(cand, word, amap):
    """Length of the common prefix of the parity vector of the rational cand = a/b
    (b odd) under amap and the given word (exact integer iteration of numerators)."""
    m0, r0, m1, r1 = amap
    a, b = cand.numerator, cand.denominator
    for i, e in enumerate(word):
        if (a & 1) != e:
            return i
        if e:
            a = (m1 * a + r1 * b) // 2
        else:
            a = (m0 * a + r0 * b) // 2
    return len(word)


def rational_certificate(word_fn, amap, N, check_len):
    """Wang reconstruction of Phi(w) mod 2^N (unique candidate of height <= 2^h,
    h = floor((N-1)/2)), then the candidate's parity vector is compared with w
    up to check_len letters (word_fn(L) returns the first L letters of w).
      ('NONE', h, lam)  the unique candidate (if any) differs from w at letter lam,
                        so NO rational a/b with |a|, b <= 2^h equals Phi(w);
      ('AGREES', h, L)  the candidate agrees with w on check_len letters (the word
                        contains a repetition running past the window)."""
    w = word_fn(N)
    x = phi_residue(w, amap, N)
    cand, bound = rat_recon(x, N)
    h = bound.bit_length() - 1
    if cand is None:
        return "NONE", h, None, None
    w_long = word_fn(check_len)
    lam = parity_agreement(cand, w_long, amap)
    hb = max(abs(cand.numerator), cand.denominator).bit_length()
    if lam < check_len:
        return "NONE", h, lam, hb
    return "AGREES", h, check_len, hb


def theorem_s_bound(qn, qn1, mu):
    """Proved lower bound (up to the O(log) constants of the note) for the best
    gain at the convergent pair (q_n, q_{n+1}): q_{n+1}/mu - (mu-1) q_n."""
    return qn1 / mu - (mu - 1) * qn


def sturmian_tests(lib):
    rng = random.Random(20260922)
    rand_rho = Fraction(rng.getrandbits(60), 1 << 60)
    rhos = [("0", mpmath.mpf(0), False), ("0 upper", mpmath.mpf(0), True),
            ("1/2", mpmath.mpf(1) / 2, False),
            ("rnd", mpmath.mpf(rand_rho.numerator) / mpmath.mpf(2) ** 60, False)]
    tests = []
    q3 = ["log_3 2 (critical q=3)", "[0;1,1,1,(2)] quadratic, below crit.", "[0;1,1,1,2,2,(3)] quadratic",
          "[0;1,1,1,2,2,3,1,5,2,23,2,(3)] quadr., just above", "[0;1,1,1,(3)] quadratic, above crit.",
          "1/phi=[0;(1)] quadratic", "sqrt2/2=[0;1,(2)] quadratic",
          "[0;1,(1,2)] quadratic, near 0.73", "[0;1,(12)] quadratic, near 0.92"]
    for name in q3:
        for rn, rho, up in rhos:
            tests.append(("3x+1", name, rn, rho, up))
    for name in ["log_3 2 (critical q=3)", "sqrt2/2=[0;1,(2)] quadratic"]:
        tests.append(("3x-1", name, "0", mpmath.mpf(0), False))
    for name in ["log_5 2 (critical q=5)", "[0;2,3,(9)] quadratic, near log_5 2",
                 "[0;2,(3)] quadratic, near 0.434", "sqrt2-1=[0;(2)] quadratic, 0.414",
                 "1/phi=[0;(1)] quadratic", "sqrt2/2=[0;1,(2)] quadratic"]:
        for rn, rho, up in rhos[:3]:
            tests.append(("5x+1", name, rn, rho, up))
    tests.append(("mahler", "[0;3,(3)] quadratic, 0.303 (Mahler-safe range)", "0", mpmath.mpf(0), False))
    tests.append(("mahler", "sqrt2-1=[0;(2)] quadratic, 0.414", "0", mpmath.mpf(0), False))
    return tests, rand_rho


def analyse_word(word_fn, amap, N, periods, check_len):
    """Certificates for one word; word_fn(L) returns its first L letters.
    Repetition certificates use the long window (check_len letters), so that
    repetitions longer than N (large partial quotients) are measured exactly."""
    w = word_fn(N)
    x = phi_residue(w, amap, N)
    kind, h, lam, hb = rational_certificate(word_fn, amap, N, check_len)
    w_long = word_fn(check_len)
    certs = repetition_certificates(w_long, amap, periods)
    ok = True
    for (p, gain, t, lam_, P, Q) in certs[-2:]:
        prec = min(lam_ + 64, len(w_long))
        ok &= check_approximant_valuation(w_long[:prec], amap, t, p, lam_, P, Q, prec)
    return x, (kind, h, lam, hb), certs, ok


def sturmian_report(N=20000, check_factor=8, printer=print):
    """Run the Sturmian test; return summary rows (for the main script)."""
    lib = slope_library()
    tests, rand_rho = sturmian_tests(lib)
    rows = []
    check_len = check_factor * N
    printer(f"S4.1 Sturmian words: Phi(w) mod 2^{N} (exact).  Rational certificate: Wang reconstruction at N bits gives the"
            f" unique candidate of height <= 2^{(N - 1) // 2}; its parity vector is compared with w up to {check_len} letters.")
    printer(f"     'NONE@lam(hb)': the candidate (height 2^hb) leaves w at letter lam, so no rational of height <= 2^{(N - 1) // 2} equals Phi(w).")
    printer(f"     intercept 'rnd' = {rand_rho.numerator}/2^60; 'U' = upper mechanical word.")
    printer("     cols: map | slope | intercept | 1-density | mu | 1-mu(mu-1) | bitlen(Phi mod 2^N) | rational certificate |"
            " best gain [p,t,lambda] | gain at q_n (proof bound q_{n+1}/mu-(mu-1)q_n) for the last 4 q_n")
    for mapname, sname, rn, rho, up in tests:
        amap = MAPS[mapname]
        alpha = lib[sname]
        cache = {}

        def word_fn(L, alpha=alpha, rho=rho, up=up, cache=cache):
            if cache.get("n", 0) < L:
                cache["w"] = mechanical_word(alpha, rho, L, upper=up)
                cache["n"] = L
            return cache["w"][:L]
        terms, qs = convergents(alpha, 60)
        periods = [qq for qq in qs if 3 * qq <= check_len]
        x, (kind, h, lam, hb), certs, ok = analyse_word(word_fn, amap, N, periods, check_len)
        a = float(alpha)
        mu = mu_of(a, amap)
        best = max(certs, key=lambda c: c[1]) if certs else None
        qpos = {qq: i for i, qq in enumerate(qs)}
        tail = []
        for c in certs[-4:]:
            i = qpos.get(c[0])
            b = theorem_s_bound(c[0], qs[i + 1], mu) if i is not None and i + 1 < len(qs) else float('nan')
            tail.append(f"{c[0]}:{c[1]:.0f}({b:.0f})")
        if kind == "NONE":
            cert_txt = f"NONE@{lam}({hb})" if lam is not None else "NONE(no cand.)"
        else:
            cert_txt = f"AGREES>={lam}({hb})"
        m0_, _, m1_, _ = amap
        crit = 1 / math.log2(m1_) if m0_ == 1 else None
        regime = "-" if crit is None else ("crit" if abs(a - crit) < 1e-15 else ("super" if a > crit else "sub"))
        row = dict(map=mapname, slope=sname, alpha=a, rho=rn, upper=up, density=sum(word_fn(N)) / N, mu=mu, regime=regime,
                   margin=1 - mu * (mu - 1), bitlen=x.bit_length(), kind=kind, h=h, lam=lam, hb=hb,
                   best_gain=best[1] if best else None, best_p=best[0] if best else None,
                   best_t=best[2] if best else None, best_lam=best[3] if best else None,
                   gains=[(c[0], c[1]) for c in certs], iso_ok=ok,
                   bounds=[(c[0], theorem_s_bound(c[0], qs[qpos[c[0]] + 1], mu)) for c in certs if c[0] in qpos and qpos[c[0]] + 1 < len(qs)])
        rows.append(row)
        printer(f"  {mapname:6s} | {sname[:36]:36s} | {regime:5s} | {rn + (' U' if up else ''):7s} | {row['density']:.5f} | {mu:.4f} |"
                f" {row['margin']:+.3f} | {x.bit_length():6d} | {cert_txt:22s} |"
                + (f" {best[1]:7.1f} [{best[0]},{best[2]},{best[3]}] | {' '.join(tail)}" if best else " (no certificate in window)")
                + f"{'' if ok else '  ISOMETRY CHECK FAILED'}")
    return rows


def control_report(N=20000, check_factor=8, printer=print):
    """Non-Sturmian controls: eventually periodic (must reconstruct), Thue-Morse,
    random words, in-house critical halving words, random bounded-discrepancy words."""
    amap = MAPS["3x+1"]
    rows = []
    L = check_factor * N

    def pv27(n):
        x0, w = 27, []
        while len(w) < n:
            w.append(x0 % 2)
            x0 = x0 // 2 if x0 % 2 == 0 else (3 * x0 + 1) // 2
        return w
    ctrls = [
        ("eventually periodic u=1101, v=10011", eventually_periodic([1, 1, 0, 1], [1, 0, 0, 1, 1], L)),
        ("eventually periodic: parity vector of 27", pv27(L)),
        ("Thue-Morse (density 1/2)", thue_morse(L)),
        ("random, density 1/2 (seed 1)", random_word(L, 0.5, 1)),
        ("random, density 0.7 (seed 2)", random_word(L, 0.7, 2)),
        ("critical halving word rho=0 (in-house sec. 3)", critical_halving_word(mpmath.mpf(0), L)[:L]),
        ("critical halving word rho=1/2 (in-house sec. 3)", critical_halving_word(mpmath.mpf(1) / 2, L)[:L]),
        ("random bounded critical discrepancy, width 3 (seed 3)", balanced_random_critical(L, 3)[:L]),
    ]
    printer("S4.2 Controls (3x+1): cols: word | 1-density | rational certificate | best repetition gain, periods 1..400 (p,t,lambda)")
    for name, wl in ctrls:
        word_fn = (lambda n, wl=wl: wl[:n])
        kind, h, lam, hb = rational_certificate(word_fn, amap, N, L)
        certs = repetition_certificates(wl[:N], amap, list(range(1, 401)), n_clusters=2)
        best = max(certs, key=lambda c: c[1]) if certs else None
        if kind == "NONE":
            cert_txt = f"NONE@{lam}({hb})" if lam is not None else "NONE(no cand.)"
        else:
            cand = rat_recon(phi_residue(wl[:N], amap, N), N)[0]
            cert_txt = f"AGREES>={lam}: {cand}" if hb < 64 else f"AGREES>={lam}({hb})"
        rows.append(dict(name=name, density=sum(wl[:N]) / N, kind=kind, lam=lam, hb=hb, best=best[:4] if best else None))
        printer(f"  {name:52s} | {sum(wl[:N]) / N:.4f} | {cert_txt:28s} | "
                + (f"{best[1]:7.1f} (p={best[0]}, t={best[2]}, lam={best[3]})" if best else "none"))
    return rows


def deep_report(N=100000, check_factor=3, printer=print):
    """One deep row: the critical slope log_3 2, intercept 0, 3x+1, at N = 10^5 bits."""
    lib = slope_library()
    alpha = lib["log_3 2 (critical q=3)"]
    amap = MAPS["3x+1"]
    cache = {}

    def word_fn(L):
        if cache.get("n", 0) < L:
            cache["w"] = mechanical_word(alpha, mpmath.mpf(0), L)
            cache["n"] = L
        return cache["w"][:L]
    check_len = check_factor * N
    terms, qs = convergents(alpha, 60)
    periods = [qq for qq in qs if 3 * qq <= check_len]
    x, (kind, h, lam, hb), certs, ok = analyse_word(word_fn, amap, N, periods, check_len)
    best = max(certs, key=lambda c: c[1])
    cert_txt = (f"NONE@{lam}({hb})" if lam is not None else "NONE(no cand.)") if kind == "NONE" else f"AGREES>={lam}({hb})"
    printer(f"S4.3 Deep row: 3x+1, slope log_3 2, intercept 0, N = {N} bits, words to {check_len} letters:"
            f" bitlen(Phi mod 2^N) = {x.bit_length()}; rational certificate {cert_txt} (no rational of height <= 2^{h});"
            f" best repetition gain {best[1]:.1f} at p={best[0]}, t={best[2]}, lambda={best[3]}; gains "
            + " ".join(f"{c[0]}:{c[1]:.0f}" for c in certs) + ("" if ok else "  ISOMETRY CHECK FAILED"))
    return dict(kind=kind, h=h, lam=lam, hb=hb, best=best[:4], certs=[(c[0], c[1]) for c in certs])


if __name__ == "__main__":
    sturmian_report()
    control_report()
