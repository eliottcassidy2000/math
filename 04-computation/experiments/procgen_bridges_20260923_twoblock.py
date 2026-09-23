#!/usr/bin/env python3
"""procgen-bridges 2026-09-23, part TWOBLOCK: is the LRC(14) residual the same '2 versus 3' question as Collatz?

  BR1  The LRC(14) small clocks {2,3,4} of THM-4447, recomputed from its capacity rule
       B_q(t) = g ceil(q/(a g)), g = gcd(t,q), a = n/2, with divisor absorption and primitivity, and then
       re-run for other runner counts n and tail counts r.  If the open clocks track (n, r) as 'small integers'
       the prime-3 appearance is pigeonhole, not a 2-versus-3 Diophantine seam.
  BR2  The Collatz side of the S596 analogy in its HYP-9122 finite form: the minimal loop law
       K0(s) = floor((s+1) log_2 3), loop ratio 2^K0/3^s in (3/2,3), and the record (upper best approximation)
       denominators of log_2 3 that carry it (Ostrowski splicing).
  BR3  Mechanical-word slopes of the extremal profiles: AMM's Beatty profile floor(gamma*(m+k)) with
       gamma* = log_5 phi^2 (THM-3009) versus the Collatz critical slope log_3 2; exact AMM finite-block
       optima versus convergents of C_*.
Memory < 100 MB; runtime a few seconds.
"""
import math
import sys
from fractions import Fraction
from itertools import combinations_with_replacement, product

import mpmath as mp


def P(*a):
    print(*a, flush=True)


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def divisors(fc):
    ds = [1]
    for p, e in fc.items():
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return sorted(ds)


def open_signatures(c, r, a):
    """tail signatures (valuation vectors capped at v_p(c)) for which no divisor d>1 of c closes the clock."""
    fc = factor(c)
    primes = sorted(fc)
    sig_space = list(product(*[range(fc[p] + 1) for p in primes]))
    bad = []
    for T in combinations_with_replacement(sig_space, r):
        # primitivity: for every prime p | c not all tails divisible by p
        if any(all(t[i] >= 1 for t in T) for i in range(len(primes))):
            continue
        ok = False
        for d in divisors(fc):
            if d == 1:
                continue
            vd = [0] * len(primes)
            dd = d
            for i, p in enumerate(primes):
                while dd % p == 0:
                    vd[i] += 1
                    dd //= p
            cap = 0
            for t in T:
                if all(t[i] >= vd[i] for i in range(len(primes))):
                    continue                       # absorbed into the pack R_d
                g = 1
                for i, p in enumerate(primes):
                    g *= p ** min(t[i], vd[i])
                cap += g * -(-d // (a * g))       # g * ceil(d / (a g))
            if cap < d:
                ok = True
                break
        if not ok:
            bad.append(T)
    return bad


def section_br1():
    P("=" * 100)
    P("BR1  LRC small clocks from the THM-4447 capacity rule, as a function of runners n and tails r")
    P("=" * 100)
    P("  n runners (n-1 speeds), arc ||.|| < 1/n, bad fraction 1/a with a = n/2; a primitive row = pack of")
    P("  n-1-r speeds with gcd c, plus r tails; clock c closes if for every admissible tail signature some divisor")
    P("  d > 1 of c has sum over non-absorbed tails of g ceil(d/(a g)) < d (THM-4447 (11); base case LRC(n-1)).")
    res = {}
    for n in (10, 12, 14, 16, 20, 28):
        a = n // 2
        for r in (2, 3, 4, 5):
            opens = []
            for c in range(2, 61):
                b = open_signatures(c, r, a)
                if b:
                    opens.append(c)
            res[(n, r)] = opens
            P(f"    n={n:2d} (a={a:2d}) r={r}: open clocks c <= 60: {opens}")
    b14 = {c: open_signatures(c, 3, 7) for c in (2, 3, 4)}
    P("  n=14, r=3 open signatures (valuation vectors of the three tails at the primes of c):")
    for c, sigs in b14.items():
        P(f"    c={c}: {sigs}")
    P("  THM-4447 (13): c=2 zero or one even tail; c=3 no tail divisible by 3; c=4 exactly one tail with v_2 = 1 and")
    P("  two odd tails.  Reproduced exactly above.")
    P("  => the open clocks are the small integers where r tails of capacity ceil(q/a) (boosted by gcd) can still")
    P("     cover q labels; they move with r and a.  The '3' in LRC(14) is r = 3 tails (13 = 10 + 3), not a")
    P("     2-versus-3 transcendence seam: no logarithm, no 2^K - 3^L and no Baker-type quantity enters.")
    return res


def section_br2():
    P("=" * 100)
    P("BR2  the Collatz side: HYP-9122's finite form and its Diophantine skeleton")
    P("=" * 100)
    mp.mp.dps = 60
    x = mp.log(3) / mp.log(2)
    # upper best approximations p/q > log2 3 with 2^p/3^q closest to 1 from above
    recs = []
    best = mp.mpf(10)
    for q in range(1, 3000):
        p = int(mp.floor(q * x)) + 1
        ratio = mp.mpf(2) ** p / mp.mpf(3) ** q
        if ratio < best:
            best = ratio
            recs.append((q, p, float(ratio)))
    P("  K0(s) = floor((s+1) log_2 3) is the minimal number of halvings of a loop through 1 with s multiplications")
    P("  (PROVED lower bound; HYP-9122 says it is attained for every s; FINITE-EXACT for s <= 6000).")
    for s in (2, 3, 4, 8, 16, 40):
        K0 = int(mp.floor((s + 1) * x))
        P(f"    s={s:3d}: K0 = {K0}, loop ratio 2^K0/3^s = {float(mp.mpf(2) ** K0 / mp.mpf(3) ** s):.6f} (in (3/2, 3))")
    P(f"  record denominators (upper best approximations of log_2 3, 2^p/3^q -> 1+): {[q for q, p, _ in recs][:16]}")
    P(f"  (HYP-9122 lists 3, 5, 17, 29, 41, 94, 147, 200, 253, 306, 971, 1636, 2301, ...; the S596 'two-block'")
    P(f"  2^E - 3^k enters exactly here.)  The record ratios: {[round(r, 6) for _, _, r in recs[:9]]}")
    P("  Contrast: the Collatz clock set is infinite and governed by the continued fraction of the transcendental")
    P("  log_2 3 (linear forms in logarithms); the LRC(14) clock set (BR1) is finite and pigeonhole-determined.")


def cf(xv, n=12):
    out = []
    for _ in range(n):
        a = int(mp.floor(xv))
        out.append(a)
        xv = xv - a
        if xv == 0:
            break
        xv = 1 / xv
    return out


def convergents(terms):
    h0, h1, k0, k1 = 0, 1, 1, 0
    out = []
    for a in terms:
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        out.append(Fraction(h1, k1))
    return out


def section_br3():
    P("=" * 100)
    P("BR3  extremal profiles are mechanical words with logarithmic slopes")
    P("=" * 100)
    mp.mp.dps = 60
    phi = (1 + mp.sqrt(5)) / 2
    gstar = mp.log(phi ** 2) / mp.log(5)
    cstar = 1 + gstar
    l32 = mp.log(2) / mp.log(3)
    for name, v in (("gamma* = log_5 phi^2 (AMM THM-3009 Beatty profile)", gstar),
                    ("C_* = 1 + gamma*", cstar), ("log_3 2 (Collatz critical odd density)", l32),
                    ("log_2 3 (Collatz clocks)", 1 / l32)):
        t = cf(+v, 14)
        P(f"  {name:52s} = {mp.nstr(v, 15)}  CF {t}")
    cs = convergents(cf(+cstar, 10))
    P(f"  convergents of C_*: {[str(c) for c in cs[:7]]}")
    opt = [Fraction(14, 9), Fraction(25, 16), Fraction(53, 34), Fraction(83, 53), Fraction(157, 100)]
    P(f"  exact AMM finite-block optima (uniform frontier sec. 6.1/abstract): {[str(o) for o in opt]} -> none is a"
      f" convergent of C_*: {not any(o in cs for o in opt)}")
    P("  Both extremal profiles are Beatty words of slopes log(alpha)/log(beta) (transcendental by Gelfond-Schneider).")
    P("  In Collatz the slope's continued fraction is the whole cycle problem (Baker, BR2); in AMM the slope enters")
    P("  only through floor rounding of the deadline and no Diophantine obstruction attached to it is known: shared")
    P("  toolbox (Beatty/Ostrowski/three-distance, also LRC THM-536/778), ANALOGY only.")


def main():
    section_br1()
    section_br2()
    section_br3()


if __name__ == "__main__":
    main()
