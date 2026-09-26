#!/usr/bin/env python3
"""
procgen_seven_20260926_potential.py -- exact corrected potentials for the negative-integer (top-lift)
adversary of THM-4486 (Theorem N) at every level; lane "seven" of session collatz-procgen-20260922,
2026-09-26.

U_k (level k, H = 2^(k-1)) has the nodes u = 1..H and the edges
    u even -> u/2,        u odd -> rho((q u + 1)/2), rho((q u - 1)/2),
rho(m) = the representative of m mod H in [1, H].  It is the Min graph of the arena when Max always takes
the lift with top bit 1 (THM-4486, Theorem N); a lower bound on the density of every cycle of U_k is a
lower bound on rho*(q,k) (Lemma G2 with the real potential log Phi).

Target density F = a0/p0 (the density of the adversary's critical cycle).  Put
    B = 2^(-beta),  A = 2^(beta (p0 - a0)/a0)   (an integer power of 2 for the chosen beta),
so that A^a0 B^(p0 - a0) = 1.  A positive function Phi on the positive integers with
    (E) Phi(u/2) <= B Phi(u)  (u even),     (O) Phi(v) <= A Phi(u)  (u odd, v a U_k-successor of u)
forces every cycle (a odd steps, p steps) to satisfy 1 <= A^a B^(p-a), i.e. a/p >= F.

Construction: Phi(u) = u^beta for u >= U0, and on 1 <= u < U0 the least solution of (E), (O) for the
unreduced successors (q u +- 1)/2 (exact rationals, Kleene iteration).  Hand proof (the note, section 6)
that (E), (O) then hold on U_k for every k with H >= (q U0 + 1)/2, given the finite checks below:
    C1 (E), (O) on 1 <= u < U0 with the unreduced successors;
    C2 Phi(w) <= w^beta for U0/2 <= w < U0      (even tail edges 2w -> w);
    C3 (q U0 + 1)^beta <= A 2^beta U0^beta      (odd tail edges, unwrapped; (q u + 1)/u decreases);
    C4 max_{u < U0} Phi(u) <= A U0^beta          (odd tail edges wrapped into the corrected range);
    C5 Phi > 0.
"""
from fractions import Fraction as Fr


def least_potential(q, a0, p0, beta, U0, max_passes=10000):
    eA = beta * (p0 - a0)
    if eA % a0:
        raise ValueError("beta (p0 - a0)/a0 must be an integer")
    A = Fr(2 ** (eA // a0))
    B = Fr(1, 2 ** beta)
    Phi = [Fr(0)] * U0          # index u, 1 <= u < U0

    def get(v):
        return Phi[v] if v < U0 else Fr(v ** beta)

    passes = 0
    while True:
        passes += 1
        if passes > max_passes:
            raise RuntimeError("no convergence (a cycle of density < F in the corrected range?)")
        changed = False
        for u in range(U0 - 1, 0, -1):
            if u % 2 == 0:
                need = get(u // 2) / B
            else:
                need = max(get((q * u + 1) // 2), get((q * u - 1) // 2)) / A
            if need > Phi[u]:
                Phi[u] = need
                changed = True
        if not changed:
            break
    return Phi, A, B, passes


def check_potential(q, a0, p0, beta, U0, Phi, A, B):
    """the finite checks C1-C5 (exact); returns a dict of booleans and data"""
    def get(v):
        return Phi[v] if v < U0 else Fr(v ** beta)
    c1 = True
    for u in range(1, U0):
        if u % 2 == 0:
            c1 &= get(u // 2) <= B * Phi[u]
        else:
            c1 &= get((q * u + 1) // 2) <= A * Phi[u] and get((q * u - 1) // 2) <= A * Phi[u]
    c2 = all(Phi[w] <= Fr(w ** beta) for w in range(U0 // 2, U0))
    c3 = (q * U0 + 1) ** beta <= A * 2 ** beta * U0 ** beta
    mx = max(Phi[1:])
    c4 = mx <= A * U0 ** beta
    c5 = all(Phi[u] > 0 for u in range(1, U0))
    cmax = max(Phi[u] / Fr(u ** beta) for u in range(1, U0))
    return {'C1': c1, 'C2': c2, 'C3': c3, 'C4': c4, 'C5': c5, 'cmax': cmax,
            'c1': Phi[1], 'Hmin': (q * U0 + 1 + 1) // 2}


def tight_cycle_check(q, cyc):
    """cyc = list of positive integers forming a closed walk of u -> u/2, (q u +- 1)/2; returns (a, p)"""
    p = len(cyc)
    a = 0
    for i, u in enumerate(cyc):
        v = cyc[(i + 1) % p]
        if u % 2 == 0:
            assert v == u // 2, (u, v)
        else:
            assert v in ((q * u + 1) // 2, (q * u - 1) // 2), (u, v)
            a += 1
    return a, p


def check_Uk_edges(q, a0, p0, beta, U0, Phi, A, B, k):
    """direct exact check of Phi(v) <= lambda Phi(u) on every edge of U_k (nodes 1..H, H = 2^(k-1));
    independent of the case analysis of the hand proof"""
    H = 1 << (k - 1)

    def P(v):
        return Phi[v] if v < U0 else Fr(v ** beta)

    def rho(m):
        r = m % H
        return H if r == 0 else r
    for u in range(1, H + 1):
        if u % 2 == 0:
            if P(u // 2) > B * P(u):
                return False, u
        else:
            for m in ((q * u + 1) // 2, (q * u - 1) // 2):
                if P(rho(m)) > A * P(u):
                    return False, u
    return True, None
