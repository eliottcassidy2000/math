#!/usr/bin/env python3
"""procgen-bridges 2026-09-23, part DEFECTS: the owner's inequality principle made concrete.

Every Cauchy-Schwarz / AM-GM / Jensen / triangle / union bound has an exact defect identity.  This script
computes the defect on the extremal objects of the three programs.

  DF1  AMM, THM-4467 box majorant |W_m(p)| <= S(p)^(d_m) (a triangle inequality in the Bernstein basis)
       (a) per level: over box-admissible parity-correct coefficient vectors the bound is attained up to a
           factor in [2/pi, 1] (zonotope support function), so there is no exponential per-level defect;
       (b) for the golden (separately balanced, lacunary) class the block defect rate is
           delta(p) = gamma log S(p) - log|1-p| = -log|p(1-p)| on the box boundary: zero exactly at the
           golden point p = -1/phi when gamma = gamma*, and nowhere else.
  DF2  Collatz, the drift is an AM-GM defect: for T = (x/2, (qx+1)/2) the arithmetic mean of the size
       factor is (q+1)/4 and the geometric mean sqrt(q)/2; q = 3 is the unique odd q with AM = 1.
  DF3  LRC(14), union/Bonferroni bounds for 13 speeds with arcs ||v t|| < 1/14, computed exactly (rational
       sweep): the spectrum of the multiplicity M(t), BONF_j, the exact defect identity
       P(M=0) - BONF_j = E[C(M-1, j) 1{M>=1}] (j odd), the integrality floor of Var(M), and the location
       of the defect mass near rationals a/q, and the excised Bonferroni bound (the LRC lever test).
  DF4  LRC(14), THM-4449's owner identity E(T) = mu(F_T) + mu(Omega_T) recomputed independently (exact sweep)
       and generalized: pair energy - union mass = E[(|K_1||K_2| - 1) 1_F] for any odd tail set.
  DF5  dyadic shadows: GF(2) Hankel/linear-complexity profiles of the lacunary skeletons
       sum w^(2^j) (AMM, Lemma P), sum x^(k^2) and sum x^(k^3) (Collatz theta and cube-swap numbers).
Memory < 200 MB; runtime about 2 minutes.
"""
import math
import sys
import time
from fractions import Fraction
from itertools import combinations
from math import comb

import numpy as np

T0 = time.time()
PHI = (1 + 5 ** 0.5) / 2
GSTAR = math.log(PHI) / math.log(5 ** 0.5)
GAMMA1 = 0.3775246


def P(*a):
    print(*a, flush=True)


def tick(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------------------------- DF1
def r_boundary(theta, g):
    """t > 0 with t (t + |1 - t e^{i theta}|)^g = 1 : the boundary of Omega_0(g) in direction theta"""
    lo, hi = 0.0, 1.0
    e = complex(math.cos(theta), math.sin(theta))
    f = lambda t: t * (t + abs(1 - t * e)) ** g - 1
    while f(hi) < 0:
        hi *= 2
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def section_df1():
    P("=" * 100)
    P("DF1  AMM: the Bernstein-box triangle inequality |W_m(p)| <= S(p)^(d_m)")
    P("=" * 100)
    P("  (a) per-level sharpness: max over e_k in [-C(d,k), C(d,k)], e_k == C(d,k) mod 2 (extreme points +-C(d,k)")
    P("      satisfy the parity) of |sum e_k p^(d-k) q^k| equals max_theta sum_k C(d,k)|p|^(d-k)|q|^k |cos(theta - arg z_k)|.")
    pts = [("golden point -1/phi", complex(-1 / PHI, 0)), ("boundary, theta=2.2 (gamma1)", None),
           ("boundary, theta=1.3 (gamma1)", None), ("interior 0.3+0.6i", complex(0.3, 0.6))]
    out = []
    for name, p in pts:
        if p is None:
            th = float(name.split("theta=")[1].split()[0])
            p = r_boundary(th, GAMMA1) * complex(math.cos(th), math.sin(th))
        q = 1 - p
        S = abs(p) + abs(q)
        row = []
        for d in (50, 200, 800, 2000):
            k = np.arange(d + 1)
            logw = (np.array([math.lgamma(d + 1) - math.lgamma(i + 1) - math.lgamma(d - i + 1) for i in range(d + 1)])
                    + (d - k) * math.log(abs(p)) + k * math.log(abs(q)) - d * math.log(S))
            wts = np.exp(logw)
            ang = (d - k) * np.angle(p) + k * np.angle(q)
            best = max(float(np.sum(wts * np.abs(np.cos(t - ang)))) for t in np.linspace(0, math.pi, 721))
            row.append(best)
        out.append((name, p, row))
        P(f"      {name:30s} p={p.real:+.4f}{p.imag:+.4f}i  ratio max|E_m|/S^d at d=50,200,800,2000: "
          + ", ".join(f"{x:.4f}" for x in row))
    P(f"      all ratios lie in [2/pi, 1] = [{2 / math.pi:.4f}, 1]: the per-level bound is exponentially sharp; any")
    P("      exponential gain must come from coupling ACROSS levels (the fairness identity), not from one level.")
    P("  (b) golden class (separately balanced blocks, deviation eps_N w^N): per-unit-N defect of the box bound on")
    P("      the boundary |p| S(p)^gamma = 1 of Omega_0(gamma) is delta(p) = gamma log S - log|1-p| = -log|p(1-p)|.")
    for g, nm in ((GSTAR, "gamma* (golden)"), (0.5, "gamma = 1/2"), (GAMMA1, "gamma_1 = 0.3775")):
        ths = np.linspace(0.02, math.pi, 4000)
        vals = []
        for th in ths:
            r = r_boundary(th, g)
            p = r * complex(math.cos(th), math.sin(th))
            vals.append(-math.log(abs(p * (1 - p))))
        vals = np.array(vals)
        i = int(np.argmin(vals))
        frac_small = float(np.mean(vals < 0.05))
        th_min = ths[i]
        P(f"      {nm:20s}: min delta = {vals[i]:+.6f} at arg p = {th_min:.4f} (pi = {math.pi:.4f});"
          f" boundary fraction with delta < 0.05: {frac_small:.3f}")
    P("      => at gamma* the box inequality is TIGHT for the lacunary class exactly at the golden point and loose")
    P("      everywhere else; below gamma* it is violated there (Theorem B).  Every sub-golden construction must")
    P("      manufacture a defect at w = -1: the golden-zero handoff states S_N(-1) = 0 (uniform frontier, sec. 6)")
    P("      are exactly such a structured defect.")


# ---------------------------------------------------------------------------------------------- DF2
def section_df2():
    P("=" * 100)
    P("DF2  Collatz: the drift is exactly an AM-GM defect")
    P("=" * 100)
    P("  T(x) = x/2 (even), (qx+1)/2 (odd); for Haar-random x in Z_2 the parity bits are i.i.d. fair, so the size")
    P("  factor per step is 1/2 or q/2 with probability 1/2 each (ignoring the +1).")
    for q in (1, 3, 5, 7):
        am = (1 + q) / 4
        gm = math.sqrt(q) / 2
        P(f"    q={q}: AM = {am:.4f}, GM = {gm:.4f}, log AM = {math.log(am):+.5f}, log GM = {math.log(gm):+.5f},"
          f" AM-GM defect log(AM/GM) = {math.log(am / gm):.5f}")
    P("  q = 3 is the unique odd q with AM = 1: the mean size is conserved and the whole typical descent")
    P(f"  log(2/sqrt3) = {math.log(2 / math.sqrt(3)):.5f} nats/step is the AM-GM gap (second-order estimate")
    P(f"  Var(log factor)/2 = (log 3)^2/8 = {math.log(3) ** 2 / 8:.5f}).")
    beta = math.log(2) / math.log(3)
    h = -beta * math.log2(beta) - (1 - beta) * math.log2(1 - beta)
    P(f"  The structure this discards: the words whose empirical GM is >= 1, i.e. odd-step density >= log_3 2 =")
    P(f"  {beta:.5f}; their entropy h(log_3 2) = {h:.5f} bits is the PROVED dimension of the exceptional set.  The same")
    P("  AM and GM hold for 3x-1, where cycles exist: the defect is SHEET-blind (foundry typology).  So the Jensen")
    P("  defect is not a lever for Collatz; it is the descent mechanism itself, and its equality cases (GM >= 1,")
    P("  near-cycles at convergent clocks of log_2 3) are exactly where the problem lives.")


# ---------------------------------------------------------------------------------------------- DF3
def multiplicity_sweep(V, half=Fraction(1, 14), extra=()):
    """exact sweep of M(t) = #{v : ||v t|| < half} on [0,1).  Returns intervals (a, b, M) with Fraction ends.
    `extra` are additional neutral breakpoints (used to make clock neighbourhoods exact)."""
    ev = []
    for v in V:
        for j in range(v):
            ev.append((Fraction(j) / v + half / v, -1))            # leave arc around j/v
            ev.append((Fraction(j + 1) / v - half / v, +1))        # enter arc around (j+1)/v
    for x in extra:
        if 0 < x < 1:
            ev.append((x, 0))
    ev.sort()
    M = len(V)
    cur = Fraction(0)
    out = []
    for x, s in ev:
        if x > cur:
            out.append((cur, x, M))
            cur = x
        M += s
    if cur < 1:
        out.append((cur, Fraction(1), M))
    return out


def clock_nbhd_points(Qs, c):
    pts = []
    for q in Qs:
        r = Fraction(c) / (14 * q)
        for a in range(0, q + 1):
            pts += [Fraction(a, q) - r, Fraction(a, q) + r]
    return pts


def in_clock_nbhd(t, Qs, c):
    for q in Qs:
        r = Fraction(c) / (14 * q)
        a = round(t * q)
        if abs(t - Fraction(a, q)) < r:
            return True
    return False


def clock_of(t, Qmax):
    """smallest q <= Qmax with |t - a/q| < 1/(28 q) for some a (circle), else 0"""
    for q in range(1, Qmax + 1):
        a = round(t * q)
        if abs(t - a / q) < 1.0 / (28 * q):
            return q
    return 0


def section_df3():
    P("=" * 100)
    P("DF3  LRC(14): union/Bonferroni/second-moment defects for 13 speeds (exact rational sweep)")
    P("=" * 100)
    rng = np.random.default_rng(2026)
    sets = [("AP13 = {1..13} (tight)", list(range(1, 14))),
            ("{1..12, 5460} (frontier control)", list(range(1, 13)) + [5460]),
            ("26*{1..12} u {339} (frontier control)", [26 * i for i in range(1, 13)] + [339])]
    for s in range(3):
        sets.append((f"random 13-set from [1,80] #{s}", sorted(int(x) for x in rng.choice(np.arange(1, 81), 13, replace=False))))
    for name, V in sets:
        iv = multiplicity_sweep(V)
        tick(f"DF3 sweep {name} ({len(iv)} intervals)")
        dist = {}
        for a, b, M in iv:
            dist[M] = dist.get(M, Fraction(0)) + (b - a)
        pM0 = dist.get(0, Fraction(0))
        EM = sum(M * w for M, w in dist.items())
        EM2 = sum(M * M * w for M, w in dist.items())
        S = [sum(comb(M, k) * w for M, w in dist.items()) for k in range(0, 14)]
        bonf = [sum((-1) ** k * S[k] for k in range(j + 1)) for j in range(0, 8)]
        var = EM2 - EM * EM
        frac = EM - math.floor(EM)
        floor_var = frac * (1 - frac)
        P(f"  {name}: {len(iv)} intervals; P(M=0) = {float(pM0):.6f}; E[M] = {EM} ; Var(M) = {float(var):.5f}"
          f" (integrality floor {{m}}(1-{{m}}) = {float(floor_var):.5f})")
        spec = sorted(dist.items())
        P("      spectrum P(M=k): " + " ".join(f"{k}:{float(w):.4f}" for k, w in spec))
        P("      BONF_j (j=1..7): " + " ".join(f"{float(bonf[j]):+.4f}" for j in range(1, 8)))
        # exact defect identity at j = 5 and j = 3
        for j in (3, 5):
            dfe = sum(comb(M - 1, j) * w for M, w in dist.items() if M >= 1)
            lhs = pM0 - bonf[j]
            P(f"      j={j}: P(M=0) - BONF_j = {float(lhs):.6f} = E[C(M-1,{j}) 1(M>=1)] = {float(dfe):.6f}"
              f" (exact equality: {lhs == dfe})")
        # localization of the j=5 defect by clocks
        loc = {}
        tot = Fraction(0)
        for a, b, M in iv:
            if M >= 6:
                wgt = comb(M - 1, 5) * (b - a)
                q = clock_of(float((a + b) / 2), 30)
                loc[q] = loc.get(q, Fraction(0)) + wgt
                tot += wgt
        if tot > 0:
            parts = sorted(loc.items())
            small = sum(w for q, w in loc.items() if 1 <= q <= 4)
            P("      where the BONF_5 defect lives (least q <= 30 with |t - a/q| < 1/(28q) at the midpoint; 0 = none): "
              + " ".join(f"q={q}:{float(w / tot):.3f}" for q, w in parts)
              + f"; clocks 1..4 carry {float(small / tot):.3f}")
        # the lever test: Bonferroni restricted to the complement of small-clock neighbourhoods (pointwise valid)
        best = None
        cells = []
        for Qs in ((1,), (1, 2), (1, 2, 3, 4), (1, 2, 13, 26)):
            for c in (Fraction(1), Fraction(1, 2), Fraction(1, 4), Fraction(1, 8)):
                iv2 = multiplicity_sweep(V, extra=clock_nbhd_points(Qs, c))
                keep = [(a, b, M) for a, b, M in iv2 if not in_clock_nbhd((a + b) / 2, Qs, c)]
                p0 = sum((b - a) for a, b, M in keep if M == 0)
                b5 = sum((b - a) * sum((-1) ** k * comb(M, k) for k in range(6)) for a, b, M in keep)
                cells.append((Qs, c, p0, b5))
                if best is None or b5 > best[3]:
                    best = (Qs, c, p0, b5)
        P("      excised Bonferroni BONF5 on the complement of {|t - a/q| < c/(14q), q in Qs} (a valid lower bound"
          " for P(M=0)); Qs, then c = 1, 1/2, 1/4, 1/8:")
        for Qs0 in ((1,), (1, 2), (1, 2, 3, 4), (1, 2, 13, 26)):
            P(f"        Qs={Qs0}: " + " ".join(f"{float(b5):+.4f}" for Qs, c, p0, b5 in cells if Qs == Qs0))
        P(f"        best: Qs={best[0]}, c={best[1]}: BONF5 = {float(best[3]):+.5f} <= P(M=0 off nbhd) = {float(best[2]):.5f}"
          f" <= P(M=0) = {float(pM0):.5f}  (raw BONF5 = {float(bonf[5]):+.4f})")
        tick(f"DF3 excision {name}")
        pc = 1 - pM0
        pz = EM * EM / EM2
        cond_cv2 = (EM2 / pc) / ((EM / pc) ** 2) - 1 if pc > 0 else 0
        P(f"      [wrong direction for LRC, shown for the ledger] Cauchy-Schwarz/PZ: P(M>=1) = {float(pc):.6f} >="
          f" E[M]^2/E[M^2] = {float(pz):.6f}; defect P(M>=1) CV^2/(1+CV^2), conditional CV^2 = {float(cond_cv2):.5f}")


# ---------------------------------------------------------------------------------------------- DF4
def two_lift_sweep(T, half=Fraction(1, 14)):
    """quotient phase y in [0,1): lifts y/2 and (y+1)/2; tail t kills lift x iff ||t x|| < half.
    Returns intervals (a, b, K1, K2) with K1, K2 frozensets of killing tails."""
    pts = {Fraction(0), Fraction(1)}
    for t in T:
        for lift in (0, 1):
            # ||t (y+lift)/2|| < half  <=>  t (y+lift)/2 in (j - half, j + half)  <=> y in ((2(j-half))/t - lift, ...)
            for j in range(0, t + 2):
                for e in (j - half, j + half):
                    y = 2 * e / t - lift
                    if 0 < y < 1:
                        pts.add(y)
    pts = sorted(pts)
    out = []
    for a, b in zip(pts[:-1], pts[1:]):
        m = (a + b) / 2
        K1 = frozenset(t for t in T if abs(t * m / 2 - round(t * m / 2)) < half)
        K2 = frozenset(t for t in T if abs(t * (m + 1) / 2 - round(t * (m + 1) / 2)) < half)
        out.append((a, b, K1, K2))
    return out


def section_df4():
    P("=" * 100)
    P("DF4  LRC(14): THM-4449's pair energy vs physical union mass, recomputed and generalized")
    P("=" * 100)
    P("  For odd tails no tail kills both lifts, so y fails iff K1, K2 are both nonempty, and the number of pairs")
    P("  {a,b} killing both lifts is |K1||K2|.  Hence, exactly, pair energy - union mass = E[(|K1||K2| - 1) 1_F].")
    for T in ((1, 11, 121), (1, 9, 23), (1, 7, 11), (3, 5, 7, 11), (1, 9, 23, 37)):
        iv = two_lift_sweep(T)
        F = sum((b - a) for a, b, K1, K2 in iv if K1 and K2)
        E = sum((b - a) * len(K1) * len(K2) for a, b, K1, K2 in iv if K1 and K2)
        Om = sum((b - a) * (len(K1) * len(K2) - 1) for a, b, K1, K2 in iv if K1 and K2)
        both = any(K1 & K2 for a, b, K1, K2 in iv)
        P(f"    T={T}: mu(F_T) = {F} = {float(F):.6f}; E(T) = {E} = {float(E):.6f}; defect mu(Omega) = {Om}"
          f" ({float(Om / E) if E else 0:.3f} of E); any tail killing both lifts: {both}")
    P("  THM-4449 values: (1,11,121): E = 124/847 = 108/847 + 16/847; sharp caps mu(F) <= 214/1449 at (1,9,23) and")
    P("  72/539 at (1,7,11).  The defect is the (1,2)/(2,1) owner-count locus: an exactly calculable, structured part")
    P("  of the union bound, already harvested by THM-4449; the four-tail rows show it grows with the tail count.")


# ---------------------------------------------------------------------------------------------- DF5
def berlekamp_massey_profile(s):
    """linear complexity profile L_1..L_n of a 0/1 sequence over GF(2) (python ints as polynomials)"""
    n = len(s)
    C, B = 1, 1
    L, m = 0, 1
    prof = []
    for i in range(n):
        d = s[i]
        c = C >> 1
        k = 1
        while c:
            if c & 1:
                d ^= s[i - k]
            c >>= 1
            k += 1
        if d == 0:
            m += 1
        elif 2 * L <= i:
            Tt = C
            C ^= B << m
            L = i + 1 - L
            B = Tt
            m = 1
        else:
            C ^= B << m
            m += 1
        prof.append(L)
    return prof


def section_df5(n=6000):
    P("=" * 100)
    P("DF5  dyadic shadows: GF(2) linear-complexity (Hankel) profiles of the lacunary skeletons")
    P("=" * 100)
    P("  A 0/1 sequence has all Hankel determinants odd iff its profile is perfect (L_i = ceil(i/2)), iff all")
    P("  partial quotients of its GF(2) continued fraction have degree 1.  Wang-Massey (1986, CITED, not re-read):")
    P("  s_1 = 1 and s_(2i+1) = s_(2i) + s_i characterizes perfect profiles; the powers-of-2 indicator satisfies it.")
    seqs = {
        "powers of 2 (AMM skeleton, positions 1,2,4,...)": [1 if (j & (j - 1)) == 0 else 0 for j in range(1, n + 1)],
        "squares (Collatz theta, positions 1,4,9,...)": [1 if math.isqrt(j) ** 2 == j else 0 for j in range(1, n + 1)],
        "cubes (Collatz Y3, positions 1,8,27,...)": [1 if round(j ** (1 / 3)) ** 3 == j else 0 for j in range(1, n + 1)],
    }
    rng = np.random.default_rng(5)
    seqs["i.i.d. fair bits (control)"] = [int(x) for x in rng.integers(0, 2, n)]
    for name, s in seqs.items():
        prof = berlekamp_massey_profile(s)
        dev = max(abs(L - (i + 2) // 2) for i, L in enumerate(prof))
        jumps = [prof[i] - prof[i - 1] for i in range(1, n) if prof[i] != prof[i - 1]]
        degs = {}
        for jmp in jumps:
            degs[jmp] = degs.get(jmp, 0) + 1
        top = sorted(degs.items())[:6]
        perfect = all(L == (i + 2) // 2 for i, L in enumerate(prof))
        P(f"    {name:48s} n={n}: perfect profile: {perfect}; max |L_i - ceil(i/2)| = {dev};"
          f" jump sizes (count): {top}")
        tick(f"DF5 {name}")
    P("  => the AMM parity skeleton is maximally structured over GF(2) (automatic; Artin-Schreier X^2+X = w);")
    P("     the Collatz theta skeletons are GF(2)-random-like (non-automatic; cf. HC-TH6 for cubes).  The shared")
    P("     'dyadic' label hides opposite behaviour: finite-state parity in AMM, no parity structure in Collatz.")


def main():
    quick = "--quick" in sys.argv
    section_df1()
    tick("DF1")
    section_df2()
    section_df3()
    tick("DF3")
    section_df4()
    tick("DF4")
    section_df5(n=2000 if quick else 6000)
    tick("DF5")


if __name__ == "__main__":
    main()
