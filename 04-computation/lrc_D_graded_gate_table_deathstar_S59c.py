#!/usr/bin/env python3
"""
death-star-2026-07-19-S59c -- HYP-7900 part A: the D-GRADED GATE TABLE.

For the canonical-D single-far family F_D(N) = {1..N}\\{N-1} u {D(N-1)}, compute
exact M for D in {3,4,5,6}, N in 8..100.  Per (D,N): the value, the winning
modulus q* and multiplier a*, the branch label (b = q* - x when q* is an
x-pair-sum), whether M equals the slack-1 rung D/((N+1)D-1), and whether M lies
in the first-gap window (1/(N+1), 2/(2N+1)).

Purposes:
  (1) FIT the full D=4 gate G4 = {N : F_4(N) attains the rung in-window}.
      Conjecture (THM-1285 see-saw): G4 = {N == 1 mod 30, N != 1 mod 7}
      -> predicted G4 cap [8,100] = {31, 61, 91}.  N=61 already CONFIRMED
      inline (M = 4/247).
  (2) Extend the D=3 gate verification (HYP-4516, verified N<=37) to N=100.
  (3) Probe D=5 (binder 2D-1 = 9 = 3^2 COMPOSITE -- pattern may require the
      binder prime, predicting F_5 never attains its rung) and D=6 (binder 11;
      predicted first opening ~ N=211, outside this sweep -- separate probe).
Also tabulates the KILLING branch for every non-member (which competitor
overshoots), the raw mechanism data for the gate law.
"""
from fractions import Fraction as F
import sys, time
sys.path.insert(0, '04-computation')
from lrc_singlefar_absorption_atlas_deathstar_S59 import M_exact_wit

log = lambda s="": print(s, flush=True)

def branch_label(qstar, x, D):
    d = qstar - x
    if 0 < d <= 2*D + 5: return f"b={d}"
    if qstar == 2*x: return "2x"
    if qstar <= x - 8: return f"base/other q*={qstar}"
    return f"q*={qstar}"

def main():
    log("== HYP-7900 part A: D-graded gate table, D=3..6, N=8..100 ==\n")
    members = {3: [], 4: [], 5: [], 6: []}
    t00 = time.time()
    for D in (3, 4, 5, 6):
        log(f"-- D = {D} (far = {D}(N-1), rung = D/((N+1)D-1), binder 2D-1 = {2*D-1}) --")
        for N in range(8, 101):
            x = D*(N-1)
            fam = [v for v in range(1, N+1) if v != N-1] + [x]
            M, q, a = M_exact_wit(fam)
            rung = F(D, (N+1)*D - 1)
            lo, hi = F(1, N+1), F(2, 2*N+1)
            inwin = lo < M < hi
            isrung = (M == rung)
            if inwin:
                members[D].append(N)
                log(f"   N={N:>3}: M={M} {'(RUNG)' if isrung else '(NOT rung!)'} "
                    f"IN-WINDOW  [q*={q}, a={a}, {branch_label(q, x, D)}]")
            elif N % 10 == 0 or isrung:
                log(f"   N={N:>3}: M={M} loose  [{branch_label(q, x, D)}]"
                    f"{' rung-valued but outside?!' if isrung else ''}")
        log(f"   => D={D} members in [8,100]: {members[D]}\n")
    log(f"== GATE FIT ==  ({time.time()-t00:.0f}s total)")
    g4 = members[4]
    pred = [N for N in range(8, 101) if N % 30 == 1 and N % 7 != 1]
    log(f"  D=4 observed members: {g4}")
    log(f"  conjecture {{N==1 mod 30, N!=1 mod 7}}: {pred}")
    log(f"  MATCH: {g4 == pred}")
    g3 = members[3]
    pred3 = [N for N in range(8, 101) if N % 6 == 1 and (3*N+2) % 5 != 0]
    log(f"  D=3 observed members: {g3}")
    log(f"  HYP-4516 gate {{N==1 mod 6, 5 nmid 3N+2}}: {pred3}")
    log(f"  MATCH: {g3 == pred3}")
    log(f"  D=5 members (binder 9 composite; predicted none): {members[5]}")
    log(f"  D=6 members (predicted none below ~211): {members[6]}")

if __name__ == "__main__":
    main()
