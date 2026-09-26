#!/usr/bin/env python3
"""Orchestrator audit of lane `seven`, Theorem Q (all-level floors 5/16 for q = 9 and
2/7 for q = 11), written from the note's statement; the lane's scripts were not read.
Potential: Phi(u) = u^beta for u >= U0; on u < U0 the least positive solution of
  (E) Phi(u/2) <= B Phi(u) (u even),  (O) Phi((q u +- 1)/2) <= A Phi(u) (u odd),
with B = 2^-beta, A = 2^(beta (p0 - a0)/a0), F = a0/p0. Checks C1..C5 exactly
(Fractions), and Proposition S on sample rational cycles.
"""
from fractions import Fraction

def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)

def audit(q, a0, p0, beta, U0):
    B = Fraction(1, 2 ** beta)
    Aexp = Fraction(beta * (p0 - a0), a0)
    assert Aexp.denominator == 1
    A = Fraction(2 ** int(Aexp))
    assert A ** a0 * B ** (p0 - a0) == 1
    def boundary(m):
        return Fraction(m) ** beta
    # least positive solution by upward iteration from 0 (constraints give lower bounds)
    Phi = [Fraction(0)] * U0
    def val(m):
        return boundary(m) if m >= U0 else Phi[m]
    changed = True
    it = 0
    while changed:
        changed = False
        it += 1
        for u in range(U0 - 1, 0, -1):
            if u % 2 == 0:
                lb = val(u // 2) / B
            else:
                lb = max(val((q * u + 1) // 2), val((q * u - 1) // 2)) / A
            if lb > Phi[u]:
                Phi[u] = lb
                changed = True
        assert it < 5000
    # C1: (E), (O) on u < U0 (unreduced successors)
    for u in range(1, U0):
        if u % 2 == 0:
            assert val(u // 2) <= B * Phi[u]
        else:
            for m in ((q * u + 1) // 2, (q * u - 1) // 2):
                assert val(m) <= A * Phi[u]
    # C2: Phi(w) <= w^beta on U0/2 <= w < U0
    assert all(Phi[w] <= Fraction(w) ** beta for w in range(U0 // 2, U0))
    # C3: (q U0 + 1)^beta <= A 2^beta U0^beta
    assert Fraction(q * U0 + 1) ** beta <= A * 2 ** beta * Fraction(U0) ** beta
    # C4: max Phi(u) <= A U0^beta
    assert max(Phi[1:]) <= A * Fraction(U0) ** beta
    # C5: Phi > 0
    assert all(Phi[u] > 0 for u in range(1, U0))
    ratio = max(Phi[u] / Fraction(u) ** beta for u in range(U0 // 2, U0))
    return it, float(ratio)

it, r = audit(11, 2, 7, 2, 1024)
check(True, f"Theorem Q, q = 11: the least fixed point on u < 1024 satisfies C1-C5 exactly ({it} sweeps; max Phi/u^2 on [512,1024) = {r:.4f}); floor 2/7 at every level")
it, r = audit(9, 5, 16, 5, 16384)
check(True, f"Theorem Q, q = 9: the least fixed point on u < 16384 satisfies C1-C5 exactly ({it} sweeps; max Phi/u^5 on [8192,16384) = {r:.4f}); floor 5/16 at every level")

# critical cycles
def cyc_of(q, start, signs_choice):
    pass
c11 = [1, 6, 3, 16, 8, 4, 2]
for i, u in enumerate(c11):
    v = c11[(i + 1) % len(c11)]
    assert (u % 2 == 0 and v == u // 2) or (u % 2 == 1 and v in ((11 * u + 1) // 2, (11 * u - 1) // 2))
c9 = [1, 5, 22, 11, 50, 25, 112, 56, 28, 14, 7, 32, 16, 8, 4, 2]
for i, u in enumerate(c9):
    v = c9[(i + 1) % len(c9)]
    assert (u % 2 == 0 and v == u // 2) or (u % 2 == 1 and v in ((9 * u + 1) // 2, (9 * u - 1) // 2))
check(sum(u % 2 for u in c11) * 7 == 2 * len(c11) and sum(u % 2 for u in c9) * 16 == 5 * len(c9),
      "the critical sign-choice cycles (1,6,3,16,8,4,2) of 11u-+1 (density 2/7) and the 16-cycle of 9u-+1 (density 5/16) are genuine")

# Proposition S: sgn rule breaks expanding rational cycles: product identity on sample cycles
import itertools
def orbit_sgn(x, q, n):
    out = []
    for _ in range(n):
        out.append(x)
        if x.numerator % 2 == 0:
            x = x / 2
        else:
            x = (q * x + (1 if x > 0 else -1)) / 2
    return out
cnt = 0
for q in (3, 5, 7, 9):
    for num in range(-40, 41):
        for den in (1, 3, 5, 7):
            x0 = Fraction(num, den)
            if x0 == 0:
                continue
            orb = orbit_sgn(x0, q, 400)
            # detect a cycle
            seen = {}
            for i, x in enumerate(orb):
                if x in seen:
                    cyc = orb[seen[x]:i]
                    a = sum(1 for y in cyc if y.numerator % 2)
                    p = len(cyc)
                    assert q ** a < 2 ** p, (q, cyc)   # contracting
                    cnt += 1
                    break
                seen[x] = i
check(True, f"Proposition S: every cycle reached under the real-sign rule (q = 3,5,7,9; {cnt} sample orbits) is contracting (q^a < 2^p)")
