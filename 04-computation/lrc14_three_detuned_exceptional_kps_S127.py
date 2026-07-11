# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-10-S127: the d=3 detuned EXCEPTIONAL SET.
#
# Context: LRCDetunedD3.lean proves ThreeDetunedClearing in the GENERIC regime
#   Sum_j badCount(delta_j, g) < g   <=>   Sum_j (floor(q_j/7)+1)/q_j < 1     (q_j = g/gcd(delta_j,g) >= 2)
# (the d_j = g/q_j cancel: d_j*(floor(q_j/7)+1) summed vs g = d_j*q_j gives the q-only condition).
# This script enumerates the EXCEPTIONAL triples where the generic bound FAILS (sum >= 1) --
# the (q1,q2,q3) that need the separate mod-lift residual, the d=3 analogue of opus's d=2 (2,2) pair.
from fractions import Fraction as F

def term(q):
    """The per-unit-g bad-branch density for one detuned coordinate with parameter q."""
    return F(q // 7 + 1, q)   # (floor(q/7)+1)/q

def main():
    QMAX = 60
    exc = []
    for q1 in range(2, QMAX):
        if term(q1) < F(1, 3):
            break                       # 3*term(q1) < 1 => no triple with min-coord q1 can fail
        for q2 in range(q1, QMAX):
            for q3 in range(q2, QMAX):
                s = term(q1) + term(q2) + term(q3)
                if s >= 1:
                    exc.append((q1, q2, q3, s))
    print(f"d=3 detuned exceptional set (q1<=q2<=q3, qj>=2, sum>=1), qj<{QMAX}: {len(exc)}")
    for q1, q2, q3, s in exc:
        print(f"  ({q1},{q2},{q3}): sum={s} = {float(s):.4f}")
    print()
    print("STRUCTURE:")
    print(f"  min-coord q1 values appearing: {sorted(set(e[0] for e in exc))}  (q1 >= 4 => generic bound ALWAYS holds)")
    print("  (2,2,*) : term(2)+term(2) = 1/2+1/2 = 1 already, so exceptional for ALL q3 -- an INFINITE")
    print("            family = the DOUBLE-HALF-HARMONIC, the d=3 analogue of opus's d=2 (2,2) residual;")
    print("            two detuned speeds at q=2 (half-integer vs g) need the mod-2g lift, not the count.")
    n22  = sum(1 for e in exc if e[0] == 2 and e[1] == 2)
    n23p = sum(1 for e in exc if e[0] == 2 and e[1] >= 3)
    n3   = sum(1 for e in exc if e[0] == 3)
    print(f"  (2,2,*): {n22} shown (INFINITE in truth);  (2,>=3,*): {n23p} (FINITE, q3<=42);  (3,*,*): {n3} (just (3,3,3)).")
    print("  => generic d=3 closes EVERYTHING except (2,2,*) [infinite, mod-2g] + a finite small-q set.")

if __name__ == "__main__":
    main()
