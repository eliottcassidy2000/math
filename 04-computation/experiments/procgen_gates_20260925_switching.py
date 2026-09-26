#!/usr/bin/env python3
"""procgen_gates_20260925_switching.py -- PART D of the cycle-gate equidistribution lane.

Switching bound (a Weyl-type differencing step on the product structure of the DP).
Pair the positions (0,1), (2,3), ...; a pair is MIXED in w if it reads 01 or 10.  Swapping the two
letters of a mixed pair keeps the number of ones before it, so the index i_k of its one is fixed and
c_w changes by delta_k = q^(a-1-i_k) 2^(2k).  Grouping words by their skeleton gives, for every
modulus m and every h (each factor is <= cos(pi/m) when gcd(m, 2q) = 1 and h != 0 mod m),

   |sum_{w in W(p,a)} e(h c_w / m)|  <=  sum_w prod_{k mixed in w} |cos(pi h delta_k(w) / m)|  =:  C(p,a) B_m(h).

D1  exact residue counts of c_w modulo small primes l (dividing the gate, or not) by a DP over
    residues: equidistribution with the proven rate E cos(pi/l)^K for gcd(l, 2q) = 1; the coarse moduli
    2 and q (c_w = w_0 mod 2, c_w = 2^(s_{a-1}) != 0 mod q) are the exceptions;
D2  the certificate B_M(h) against |S(h)|/C at the modulus M = |G| itself: B is large exactly on the
    structured frequencies and exponentially small elsewhere, but its size e^(-c p) is far above the
    square-root scale 2^(-h(theta) p/2) that a count of integral points would need;
D3  the dense divisor regime against published data: exact primitive-cycle counts of the (3m+d)-maps for
    the eleven d of Belaga-Mignotte 2006 table (20), per clock, against the equidistribution prediction.
stdout = results, stderr = timing/memory.  Checks raise on failure.
"""
import math
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_gates_20260925_core import (S_dp, carry, factor_small, gate, lean_malloc, report_mem,  # noqa: E402
                                         words)

RNG = random.Random(20260925)


def residue_counts(p, a, q, ell):
    """exact counts n_r = #{w in W(p,a) : c_w = r mod ell} (Python ints via object arrays)."""
    v = [np.zeros(ell, dtype=object) for _ in range(a + 1)]
    v[0][0] = 1
    for t in range(p):
        lo = max(0, a - (p - t))
        for i in range(min(t, a - 1), lo - 1, -1):
            term = (pow(q, a - 1 - i, ell) * pow(2, t, ell)) % ell
            v[i + 1] = v[i + 1] + np.roll(v[i], term)
        for i in range(lo):
            v[i] = np.zeros(ell, dtype=object)
    return v[a]


def switching_B(p, a, q, m, hs):
    """B_m(h) = (1/C) sum_w prod_{mixed pairs} |cos(pi h delta / m)| by a DP over pairs (float)."""
    C = math.comb(p, a)
    out = []
    for h in hs:
        v = np.zeros(a + 1)
        v[0] = 1.0
        for k in range(p // 2):
            t = 2 * k
            nv = v.copy()                       # pair 00
            nv[2:] += v[:-2]                    # pair 11
            for i in range(0, a):               # mixed pair: one with index i at position t or t+1
                if v[i] == 0.0:
                    continue
                d = (h * pow(q, a - 1 - i, m) * pow(2, t, m)) % m
                nv[i + 1] += 2.0 * abs(math.cos(math.pi * d / m)) * v[i]
            v = nv
        if p % 2:
            nv = v.copy()
            nv[1:] += v[:-1]
            v = nv
        out.append(v[a] / C)
    return np.array(out)


def expected_cos_power(p, a, ell):
    """E cos(pi/ell)^K for a uniform word, K = number of mixed pairs (exact, by counting pair types)."""
    c = math.cos(math.pi / ell)
    npairs = p // 2
    odd = p % 2
    tot = 0.0
    C = math.comb(p, a)
    # choose K mixed pairs, J pairs '11', rest '00' (plus a free last letter when p is odd)
    for K in range(0, npairs + 1):
        for J in range(0, npairs - K + 1):
            for last in range(0, odd + 1):
                if K + 2 * J + last != a:
                    continue
                ways = math.comb(npairs, K) * math.comb(npairs - K, J) * 2 ** K
                tot += ways * c ** K
    return tot / C


def d1_divisors():
    print("=" * 110)
    print("D1. Carries modulo small primes l (exact residue counts): n_0 l / C and max_r |n_r l / C - 1|, with the")
    print("    proven switching bound max_{h != 0} |S_l(h)|/C <= E cos(pi/l)^K (K = number of mixed pairs)")
    print("=" * 110)
    print(f"  {'q':>2} {'(p,a)':>9} {'l':>5} {'l|G':>4} {'n_0 l/C':>9} {'max dev':>10} {'max|S_l|/C':>11} {'bound':>10}")
    rows = []
    for q, clocks in ((3, [(19, 12), (27, 17), (38, 24), (46, 29), (65, 41), (84, 53), (149, 94), (40, 20), (60, 30)]),
                      (5, [(28, 12), (35, 15), (37, 16), (40, 17), (60, 26)])):
        for (p, a) in clocks:
            G = gate(p, a, q)
            fac, _ = factor_small(G, 2000)
            ells = sorted({f for f, _ in fac if f > 3 and f <= 1500})
            ctrl = [ell for ell in (7, 11, 13, 101) if G % ell and ell not in ells][:2]
            for ell in ells + ctrl:
                n = residue_counts(p, a, q, ell)
                C = math.comb(p, a)
                assert sum(n) == C
                dev = max(abs(int(x) * ell - C) for x in n) / C
                # Fourier coefficients mod ell from the counts, centred first (sum_r zeta^(hr) = 0 for h != 0),
                # so the float error is relative to the deviation, not to C
                phases = np.exp(2j * np.pi * np.arange(ell) / ell)
                cen = np.array([float(int(x) * ell - C) for x in n]) / (ell * float(C))
                mS = max(abs(np.sum(cen * phases ** hh)) for hh in range(1, ell))
                bnd = expected_cos_power(p, a, ell)
                assert mS <= bnd * (1 + 1e-9) + 1e-15, (q, p, a, ell, mS, bnd)
                rows.append((q, p, a, ell))
                print(f"  {q:>2} {str((p, a)):>9} {ell:>5} {'yes' if G % ell == 0 else 'no':>4} "
                      f"{int(n[0]) * ell / C:9.6f} {dev:10.3e} {mS:11.3e} {bnd:10.3e}")
    # coarse moduli 2 and 3: exact statements checked on all words of small clocks
    bad = 0
    for q in (3, 5):
        for p in range(1, 15):
            for a in range(1, p + 1):
                for pos in words(p, a):
                    c = carry(pos, a, q)
                    if c % 2 != (1 if pos[0] == 0 else 0):
                        bad += 1
                    if c % q != pow(2, pos[-1], q):
                        bad += 1
    assert bad == 0
    print("  coarse moduli (all words, p <= 14): c_w = w_0 mod 2, and c_w = 2^(s_{a-1}) != 0 mod q; so the carries are")
    print("  never equidistributed mod 2 or mod q (the Tao-type coarse irregularity), while every modulus coprime to 2q")
    print("  is reached at the proven exponential rate.")


def d2_certificate():
    print()
    print("=" * 110)
    print("D2. Switching certificate at the gate modulus M = |G|: |S(h)|/C <= B_M(h) (checked), on structured h =")
    print("    u 2^-j 3^k and on random h.  typical(B) = median over 200 random h; sqrt-scale = C^(-1/2)")
    print("=" * 110)
    for (p, a) in [(16, 10), (19, 12), (24, 15), (27, 17), (46, 29), (65, 41), (84, 53)]:
        M = abs(gate(p, a))
        C = math.comb(p, a)
        inv2 = pow(2, -1, M)
        struct = []
        for j in (0, 1, 2, p // 2, p):
            for k in (0, 1, a // 2, a - 1):
                for u in (1, -1):
                    struct.append(u * pow(inv2, j, M) * pow(3, k, M) % M)
        struct = sorted(set(h for h in struct if h))
        rnd = [RNG.randrange(1, M) for _ in range(200)]
        Ss = np.abs(S_dp(p, a, struct, 3, normalize=True))
        Bs = switching_B(p, a, 3, M, struct)
        Sr = np.abs(S_dp(p, a, rnd, 3, normalize=True))
        Br = switching_B(p, a, 3, M, rnd)
        assert np.all(Ss <= Bs * (1 + 1e-9) + 1e-12) and np.all(Sr <= Br * (1 + 1e-9) + 1e-12)
        lC = math.log2(C)
        k = int(np.argmax(Ss))
        print(f"  ({p},{a}) M = {M}: structured max |S|/C = {Ss.max():.3e} (B = {Bs[k]:.3e});  random: "
              f"median |S|/C = {np.median(Sr):.2e}, median B = {np.median(Br):.2e}, max B = {Br.max():.2e};  "
              f"sqrt-scale 2^{-lC / 2:.1f} = {2 ** (-lC / 2):.2e}")


# Belaga-Mignotte 2006 (DMTCS proc. AG, 249-260), table (20): numbers omega(d) of primitive cycles of the
# (3m+d)-maps for the eleven d <= 19999 with omega(d) > 160 (primary text read, section 6).
BM_TABLE = {7463: 162, 18359: 164, 7727: 198, 15655: 207, 10289: 214, 9823: 241, 17021: 258, 14197: 329,
            13085: 335, 6487: 534, 14303: 944}


def td_cycles_on_clock(d, p, a, chunk=1 << 19):
    """primitive T_d-cycles (T_d(y) = y/2, (3y + d)/2 on positive integers) with p steps, a odd steps and
    primitive period p: enumerate the odd least points y in the perigee window
    d 3^(a-1)/G <= y <= d/(2^(p/a) - 3)  (the product identity prod(3 + d/y_i) = 2^p over the odd points)."""
    G = 2 ** p - 3 ** a
    assert G > 0 and G % d == 0
    lo = max(1, -(-(d * 3 ** (a - 1)) // G))
    r = 2 ** (p / a) - 3
    hi = int(d / r) + 2
    found = []
    killed = 0
    LIM = 1 << 61
    y0 = lo if lo % 2 else lo + 1
    while y0 <= hi:
        y1 = min(hi, y0 + 2 * chunk)
        Y0 = np.arange(y0, y1 + 1, 2, dtype=np.int64)
        Y = Y0.copy()
        odd = np.zeros_like(Y)
        mn = Y0.copy()
        alive = np.ones(len(Y), dtype=bool)
        first = np.zeros(len(Y), dtype=np.int32)
        for s in range(1, p + 1):
            b = Y & 1
            odd += b
            Y = np.where(b == 1, (3 * Y + d) >> 1, Y >> 1)
            big = Y > LIM
            if big.any():
                killed += int((big & alive).sum())
                alive &= ~big
                Y = np.where(big, 1, Y)
            np.minimum(mn, Y, out=mn)
            ret = (Y == Y0) & (first == 0)
            first[ret] = s
        ok = alive & (first == p) & (odd == a) & (mn == Y0)
        found.extend(int(v) for v in Y0[ok] if math.gcd(int(v), d) == 1)
        y0 = y1 + 1 if (y1 + 1) % 2 else y1 + 2
    return found, killed, lo, hi


def phi_ratio(d):
    fac, cof = factor_small(d, 10 ** 6)
    ps = [f for f, _ in fac] + ([cof] if cof > 1 else [])
    r = 1.0
    for f in ps:
        r *= 1 - 1 / f
    return r


def lyndon(p, a):
    from procgen_gates_20260925_core import lyndon_count
    return lyndon_count(p, a)


def d3_belaga_mignotte(pmax=250, wmax=40_000_000):
    print()
    print("=" * 110)
    print("D3. The dense divisor regime against published data: primitive cycles of T_d(y) = y/2, (3y+d)/2 for the")
    print("    eleven d of Belaga-Mignotte 2006 table (20).  A primitive T_d-cycle with clock (p,a) is a word with")
    print("    (G/d) | c_w and gcd(c_w/(G/d), d) = 1, so it lives on a clock with d | G.  Per clock: the exact count")
    print(f"    (perigee-window enumeration, all clocks p <= {pmax} with d | G > 0) and the equidistribution")
    print("    prediction L(p,a) (d/G) prod_(l|d)(1-1/l); a/p measures the distance to the critical line log_3 2.")
    print("=" * 110)
    tot_ok = 0
    for d, om in sorted(BM_TABLE.items(), key=lambda t: t[1]):
        clocks = []
        for p in range(2, pmax + 1):
            for a in range(1, p):
                G = 2 ** p - 3 ** a
                if G > 0 and G % d == 0:
                    clocks.append((p, a))
        total = 0
        pred_tot = 0.0
        lines = []
        skipped = []
        for (p, a) in clocks:
            G = 2 ** p - 3 ** a
            r = 2 ** (p / a) - 3
            hi = d / r
            pred = lyndon(p, a) * d / G * phi_ratio(d)
            if hi - d * 3 ** (a - 1) / G > wmax:
                skipped.append((p, a, pred))
                continue
            f, k, lo, hi_ = td_cycles_on_clock(d, p, a)
            assert k == 0, (d, p, a, k)
            total += len(f)
            pred_tot += pred
            if len(f) or pred >= 0.5:
                lines.append(f"({p},{a}) a/p={a / p:.3f}: {len(f)} vs {pred:.1f}")
        # extension along the main family k (p0, a0) to p <= 600 (the only clocks near the critical line
        # with d | G beyond p = 250; the other lattice points there have a/p < 0.2 and no eligible cycle)
        p0, a0 = max(((p, a) for (p, a) in clocks), key=lambda c: lyndon(*c) * d / (2 ** c[0] - 3 ** c[1]))
        ext = 0
        k = pmax // p0 + 1
        while k * p0 <= 600:
            f, kk, _, _ = td_cycles_on_clock(d, k * p0, k * a0)
            assert kk == 0
            ext += len(f)
            k += 1
        total += ext
        status = "MATCH" if total == om else f"differs by {om - total:+d}"
        tot_ok += total == om
        print(f"  d = {d:>5}: exact total {total:>4} | Belaga-Mignotte omega(d) = {om:>4} [{status}] | equidistribution "
              f"prediction {pred_tot:7.1f} ({len(clocks)} clocks p <= {pmax}, {len(skipped)} skipped windows; "
              f"family {p0},{a0} extended to p <= 600: +{ext})")
        print("      " + ";  ".join(lines[:6]))
        if skipped:
            print("      skipped (window > %.0e): " % wmax + ", ".join(f"({p},{a}) pred {pr:.2g}" for p, a, pr in skipped))
    print(f"  {tot_ok} of {len(BM_TABLE)} totals reproduce the published omega(d).")


def main():
    lean_malloc()
    t = time.time()
    d1_divisors()
    print(f"[D1 {time.time() - t:.1f}s]", file=sys.stderr, flush=True)
    d2_certificate()
    print(f"[D2 {time.time() - t:.1f}s]", file=sys.stderr, flush=True)
    d3_belaga_mignotte()
    report_mem("end")
    print(f"[D total {time.time() - t:.1f}s]", file=sys.stderr)


if __name__ == "__main__":
    main()
