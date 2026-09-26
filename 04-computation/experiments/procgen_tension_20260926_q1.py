#!/usr/bin/env python3
"""
procgen_tension_20260926_q1.py -- Q1: periodic ranks = potential certificates = class (i);
the rank defect equals the maximum cycle mean; the flow (circulation) side; the Bernoulli and
finite-bank obstructions read as expanding cycles.
Run through procgen_tension_20260926_run.py.
"""
import random
from fractions import Fraction
from math import log

import numpy as np

from procgen_tension_20260926_lib import (
    all_minus, as_cycle, below_c, best_lower, check, claim, cycle_in_edges, density, karp, karp_batch,
    LOG2, LOG3, periodic_point, potential, potential_np, potential_ok, residue_of, rho_all,
    sigma_k_mask, signed_word, step, succ, tight_edges, v2)


def expanding_witness(k, mask, cyc, big=10 ** 12):
    """integer n > big following the cycle cyc of G_sigma for one period: n = x_gamma mod 2^(k+p)
    shifted up.  Returns (n, T^p(n))."""
    p = len(cyc)
    x = periodic_point(signed_word(k, mask, cyc))
    m = k + p
    r0 = residue_of(x, m)
    n = r0 + (1 << m) * (big // (1 << m) + 1)
    y = n
    for j in range(p):
        check(y % (1 << k) == cyc[j], "witness leaves the cycle at step %d" % j)
        y = step(y, k, mask)
    return n, y


def section_A(levels_exh=(2, 3, 4, 5), sample_levels=(6, 7, 8), nsample=60, seed=20260926):
    print("A. (a) class (i)  <=>  (c) integer potential below c  <=>  (b) periodic rank; witnesses")
    rng = random.Random(seed)
    tot_i = tot_non = 0
    for k in levels_exh:
        R = rho_all(k)
        NM = len(R)
        cls1 = [m for m in range(NM) if below_c(R[m])]
        # (a) => (c): potential at F = rho_max for every class-(i) strategy, checked on every edge
        for m in cls1:
            psi = potential(k, m, R[m])
            check(psi is not None and potential_ok(k, m, R[m], psi), "potential at rho_max k=%d m=%d" % (k, m))
        # not (a) => not (c): no potential at F_{2^k}, the largest possible cycle density below c
        FD = best_lower(1 << k)
        non = [m for m in range(NM) if not below_c(R[m])]
        test_non = non if k <= 4 else rng.sample(non, 3000)
        for m in test_non:
            check(potential(k, m, FD) is None, "no potential below c for non-(i) k=%d m=%d" % (k, m))
        # (c) => (b): explicit rank, margin eps = |mu|/2 above n_1 (exact inequality chain) and a
        # floating-point sanity check of R(Tn) - R(n) <= -eps on 2000 consecutive integers
        for m in (cls1 if k <= 4 else rng.sample(cls1, 200)):
            F = R[m]
            q, r = F.numerator, F.denominator
            psi = potential(k, m, F)
            mu = q / r * LOG3 - LOG2
            check(mu < 0, "mu < 0")
            eps = -mu / 2
            n1 = int(1 / (3 * eps)) + 2               # 1/(3n) <= |mu|/2 for n >= n1
            check(1 / (3 * n1) <= -mu / 2, "n_1 bound")
            h = [LOG3 / r * v for v in psi]
            worst = -1e9
            for n in range(n1, n1 + 2000):
                t = step(n, k, m)
                d = log(t) - log(n) + h[t % (1 << k)] - h[n % (1 << k)]
                worst = max(worst, d)
            check(worst <= -eps + 1e-9, "numerical rank descent k=%d m=%d worst=%g eps=%g" % (k, m, worst, eps))
        # (b) => (a): every non-(i) strategy has an expanding cycle with an integer witness
        # n > 10^12, T^p(n) = n (mod 2^k), T^p(n) > n: every periodic rank increases over the period
        for m in (non if k <= 4 else rng.sample(non, 3000)):
            F = R[m]
            psi = potential(k, m, F)
            check(psi is not None and potential_ok(k, m, F, psi), "potential at rho_max (non-(i))")
            cyc = cycle_in_edges(tight_edges(k, m, F, psi))
            check(cyc is not None, "tight subgraph has a cycle")
            cyc = as_cycle(k, m, cyc)
            check(density(cyc) == F and not below_c(density(cyc)), "tight cycle is an expanding max cycle")
            n, y = expanding_witness(k, m, cyc)
            check(y % (1 << k) == n % (1 << k) and y > n, "integer witness rises over a full period")
        tot_i += len(cls1)
        tot_non += len(non)
        claim(True, "level %d: %d class-(i) strategies carry an integer potential at F = rho_max (every edge checked); "
              "%s non-(i) strategies have no potential at F_{2^k} = %s and an integer witness n > 10^12 with "
              "T^p(n) = n mod 2^k and T^p(n) > n" % (k, len(cls1), ('all %d' % len(non)) if k <= 4 else '3000 sampled',
                                                     best_lower(1 << k)))
    # (c) => (b) for the near-critical approximants sigma_8, sigma_10, sigma_12 (rho_max = F_k = 5/8)
    for k in (8, 10, 12):
        m = sigma_k_mask(k)
        F = best_lower(k)
        psi = potential_np(k, m, F)
        check(psi is not None and potential_ok(k, m, F, psi), "potential for sigma_k")
        q, r = F.numerator, F.denominator
        mu = q / r * LOG3 - LOG2
        eps = -mu / 2
        n1 = int(1 / (3 * eps)) + 2
        h = [LOG3 / r * v for v in psi]
        worst = max(log(step(n, k, m)) - log(n) + h[step(n, k, m) % (1 << k)] - h[n % (1 << k)]
                    for n in range(n1, n1 + 2000))
        check(worst <= -eps + 1e-9, "numerical rank descent for sigma_%d" % k)
    claim(True, "sigma_8, sigma_10, sigma_12 (rho_max = 5/8): the rank log n + (log 3/8) psi(n mod 2^k) descends by "
          "|mu|/2 = %.5f at every step on 2000 consecutive integers above n_1 (floating point)" % eps)
    # sampled higher levels, through Karp on random strategies with sigma(1)=+, sigma(-1)=- (else trivially non-(i))
    for k in sample_levels:
        H = 1 << (k - 1)
        cnt_i = 0
        for _ in range(nsample):
            m = rng.getrandbits(H)
            m &= ~1                   # sigma(1) = +
            m |= 1 << (H - 1)         # sigma(-1) = -
            F = karp(k, m) if k <= 7 else None
            if F is None:
                rn, rd = karp_batch(k, [m])
                F = Fraction(int(rn[0]), int(rd[0]))
            psi = potential(k, m, F)
            check(psi is not None and potential_ok(k, m, F, psi), "potential at rho_max (sampled)")
            if below_c(F):
                cnt_i += 1
            else:
                cyc = as_cycle(k, m, cycle_in_edges(tight_edges(k, m, F, psi)))
                check(density(cyc) == F, "tight cycle density")
                n, y = expanding_witness(k, m, cyc)
                check(y % (1 << k) == n % (1 << k) and y > n, "sampled witness")
        claim(True, "level %d: %d random strategies with sigma(1)=+, sigma(-1)=-: potential at rho_max always exists; "
              "%d class (i); every other one has an expanding tight cycle with an integer witness" % (k, nsample, cnt_i))


def section_B():
    print("B. the rank defect: min over periodic h of the worst step = max cycle mean; Collatz = log(3/2) at -1")
    # B1: exact min-max for every strategy at levels 2..4 and all class-(i) at level 5:
    # the potential at F = rho_max makes every edge satisfy w_F(s) + psi(t) - psi(s) <= 0 with equality on a
    # cycle, i.e. max_edges (w + dh) = mu = rho log 3 - log 2 exactly, and no h does better (cycle sum).
    for k in (2, 3, 4):
        R = rho_all(k)
        for m in range(len(R)):
            psi = potential(k, m, R[m])
            check(psi is not None and potential_ok(k, m, R[m], psi), "B1 potential")
            cyc = cycle_in_edges(tight_edges(k, m, R[m], psi))
            check(cyc is not None and density(as_cycle(k, m, cyc)) == R[m], "B1 tight cycle")
    claim(True, "levels 2-4, all 276 strategies: the optimal periodic correction h = (log 3/r) psi attains "
          "max over edges of (w(s) + h(t) - h(s)) = rho_max log 3 - log 2, with equality on a cycle of density rho_max")
    # B2: Collatz at every level 2..16: the -1 self-loop, rho_max = 1, h = 0 optimal (psi = 0 at F = 1)
    for k in range(2, 17):
        M = 1 << k
        check((M - 1) in succ(k, 0, M - 1), "Collatz self-loop at -1")
        psi = potential(k, 0, Fraction(1))
        check(psi is not None and max(psi) == 0, "psi = 0 at F = 1 for Collatz")
        if k <= 8:
            check(karp(k, 0) == 1, "Collatz rho_max = 1")
    claim(True, "Collatz, k = 2..16: node 2^k - 1 carries the self-loop (density 1); the least potential at F = 1 "
          "is psi = 0, so h = 0 is an optimal periodic correction and the least achievable defect is exactly "
          "log(3/2) (rho_max = 1 re-derived by exact Karp for k <= 8)")
    # B3: the Bernoulli witness 2^H - 1 is the -1 self-loop: T^j(2^H - 1) = -1 mod 2^k for j <= H - k
    for k in range(2, 12):
        for H in range(k, k + 40):
            n = 2 ** H - 1
            for j in range(H - k + 1):
                check(n % (1 << k) == (1 << k) - 1, "Bernoulli witness on the -1 loop")
                n = step(n, k, 0)
    claim(True, "the Bernoulli-boundary witnesses n_H = 2^H - 1 stay on the node -1 (mod 2^k) for H - k + 1 "
          "consecutive steps (k = 2..11, H = k..k+39): they are the integer shadow of the expanding self-loop")


def section_B_banks():
    print("B'. finite banks of dyadic counters plus any periodic correction cannot pay for every expanding orbit")
    # bank: all centers used by the kuratowski/bernoulli notes and more (rationals, not positive integers)
    bank = [Fraction(-1), Fraction(-5), Fraction(-7), Fraction(-10), Fraction(-17), Fraction(-25), Fraction(-13, 9),
            Fraction(-5, 3), Fraction(0), Fraction(1, 3), Fraction(-1, 5), Fraction(1, 5), Fraction(-19, 11)]
    bank += [1 - 2 * Fraction(4, 3) ** d for d in range(1, 8)]
    K = 12     # the periodic correction h may depend on n mod 2^12
    rows = []
    for L in (5, 7, 9, 13, 23):
        word = [1] * (L - 1) + [0]            # necklace 1^(L-1) 0, density (L-1)/L > log_3 2 for L >= 3
        x0 = periodic_point(word)             # Collatz signs (all +)
        orbit = [x0]
        for b in word[:-1]:
            y = orbit[-1]
            orbit.append((3 * y + 1) / 2 if b else y / 2)
        back = orbit[-1] / 2                  # last letter is 0
        check(back == x0 and len(set(orbit)) == L, "orbit of 1^(L-1)0 closes up")
        check(all(z not in bank for z in orbit), "orbit avoids the bank")
        V = max(v2(z - beta) for z in orbit for beta in bank)
        N = K + L + V + 2
        n = residue_of(x0, N) + (1 << N) * 10 ** 6
        y = n
        for j in range(L):
            y = step(y, K, 0)
        same_h = (y % (1 << K)) == (n % (1 << K))
        same_counters = all(v2(y - beta) == v2(n - beta) for beta in bank)
        check(same_h and same_counters and y > n, "bank witness L=%d" % L)
        rows.append((L, V, n.bit_length()))
    claim(True, "for a bank of %d rational centers (incl. -1, -5, -17, -13/9, 1 - 2(4/3)^d, d <= 7) and any h mod 2^12: "
          "the Collatz orbits of 1^(L-1)0, L = 5, 7, 9, 13, 23, avoid the bank and have integer shadows n with "
          "T^L(n) > n and identical h-values and counters at n and T^L(n); so every rank "
          "a log n + h(n mod 2^12) + sum c_i v2(n - beta_i), a > 0, increases over one period" % len(bank))
    print("    (L, max counter along the orbit, bits of the witness) =", rows)


def section_C():
    print("C. the flow side: circulations, the Haar circulation, LP duality")
    from scipy.optimize import linprog
    # C1: Haar (uniform) edge flow is a circulation iff sigma is constant; then its density is 1/2
    for k in (2, 3, 4):
        NM = 1 << (1 << (k - 1))
        for m in range(NM):
            indeg = [0] * (1 << k)
            for s in range(1 << k):
                for t in succ(k, m, s):
                    indeg[t] += 1
            uniform_circ = all(d == 2 for d in indeg)
            check(uniform_circ == (m in (0, all_minus(k))), "uniform circulation <=> constant sigma, k=%d m=%d" % (k, m))
            check(min(indeg) >= 1 and max(indeg) <= 3, "in-degrees in {1,2,3}")
    for k in range(5, 15):
        indeg = [0] * (1 << k)
        for s in range(1 << k):
            for t in succ(k, 0, s):
                indeg[t] += 1
        check(all(d == 2 for d in indeg), "Collatz de Bruijn in-degree 2")
    claim(True, "levels 2-4 (all 276 strategies): the uniform edge flow is a circulation iff sigma is constant "
          "(Collatz or 3n-1); in-degrees always lie in {1,2,3}; Collatz is 2-in/2-out at k = 2..14, so its Haar "
          "circulation exists and has odd density exactly 1/2 (half of the nodes are odd)")
    # C2: LP duality, numerically: max over normalized circulations of the odd density = rho_max
    worst = 0.0
    cnt = 0
    for k in (2, 3, 4):
        R = rho_all(k)
        for m in range(len(R)):
            E = [(s, t) for s in range(1 << k) for t in succ(k, m, s)]
            nE = len(E)
            A = np.zeros(((1 << k) + 1, nE))
            for i, (s, t) in enumerate(E):
                A[s, i] -= 1
                A[t, i] += 1
                A[1 << k, i] = 1
            b = np.zeros((1 << k) + 1)
            b[-1] = 1
            c = -np.array([1.0 if s & 1 else 0.0 for s, t in E])
            res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)] * nE, method="highs")
            check(res.status == 0, "LP solved")
            worst = max(worst, abs(-res.fun - float(R[m])))
            cnt += 1
    check(worst < 1e-9, "LP duality numerically")
    claim(True, "levels 2-4 (%d strategies): the LP maximum of the odd density over normalized circulations equals "
          "the exact Karp rho_max (max deviation %.1e; HiGHS, floating point)" % (cnt, worst))


def run():
    section_A()
    section_B()
    section_B_banks()
    section_C()


if __name__ == "__main__":
    run()
