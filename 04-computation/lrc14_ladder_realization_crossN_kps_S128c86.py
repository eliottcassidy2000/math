#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c86 -- HYP-7890: the two-vocabulary unification.

The 2026-07-06 binder-gate vocabulary (HYP-4506/4516/4572: canonical family
F_m(N) = {1..N-2, N, m(N-1)}, mediant/binder competition mod 30) and the
2026-07-19 determinant-slack vocabulary (THM-1230/1235, opus-S396 slack-1
ladder D/(14D-1); boxeph-S123 determinant strata) describe the SAME
near-floor spectrum and never cross-cited.  This script unifies them by
computing the exact realization law of BOTH ladder families across N:

  K-ladder  K_c(N) = {1..N-1} u {cN}       (THM-633 shape, N=12 fully proved)
  F-ladder  F_m(N) = {1..N-2, N, m(N-1)}   (HYP-4572 shape at m=3; opus-S395
                                            ladder m/(12m+5) at N=13)

Sharp questions settled here:
 (Q1) Does THM-633's law transfer to N=13 at c=2, i.e. is M({1..12,26}) = 2/27
      EXACTLY?  If yes, the slack-1 D=2 rung IS realized by a family inside
      opus-S396's own scanned shape {1..12,x}, correcting its "D=2 not found".
      (Hand derivation: at q = 13m+1, a=m, the far element 26 = 2*13 sits at
      distance exactly 2/q since 13m == -1; base {1..12} sits in-band; so
      M >= m-independent cap min(m,2)/q maximized at m=2: 2/27.)
 (Q2) Does F_2(N) TIE the floor 1/(N+1)?  THM-1220's non-rigidity sweep ran
      n=10..18 only; F_2(13) = {1..11,13,24} is the known n=14 tie.  F_2(7) =
      {1..5,7,12} would be an n=8 tie -- a hole in the sweep, and n=8 is a
      PROVED case (Rosenfeld), so a tie there reframes "n=14 is the unique
      non-rigid case" as "the unique non-rigid case IN 10..18".
 (Q3) The general binder law b*(N,m): attained Q minus m(N-1) for the
      F-ladder, extending HYP-4516's m=3 mod-30 gate to all m; and the
      analogous c-law for the K-ladder (b*=1 at N=12 by THM-633 -- where
      else?).

Evaluator: exact all-integer M with pair-sum/diff/2v candidate denominators
(THM-1002 pair-sum lemma; 2v half-turn cusps per MISTAKE-144; diffs for
safety), copied from death-star-S59's gate-verified
lrc_firstgap_crossN_census_deathstar_S59.py (credited).  Gates below must
pass or no conclusion is drawn (MISTAKE-162 positive-control rule).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

# ---------------- exact M (death-star-S59, credited) ----------------
def cand_denoms(S):
    Q = set()
    for v in S: Q.add(2*v)
    for x, y in combinations(S, 2):
        Q.add(x+y); Q.add(abs(x-y))
    Q.discard(0)
    return sorted(Q)

def M_exact_arg(S):
    """Exact M(S) with argmax (a,q) returned (first attaining candidate)."""
    S = sorted(set(S))
    bn, bd = 0, 1
    arg = None
    for q in cand_denoms(S):
        half = q // 2
        for a in range(1, half + 1):
            mn = q
            for v in S:
                r = (v*a) % q
                if r > q - r: r = q - r
                if r < mn:
                    mn = r
                    if mn * bd <= bn * q: break
            if mn * bd > bn * q:
                bn, bd = mn, q
                arg = (a, q)
    return F(bn, bd), arg

# ---------------- gates ----------------
def gates():
    tests = [
        (list(range(1,14)), F(1,14)),
        (list(range(1,13)), F(1,13)),
        (list(range(1,12))+[24], F(2,25)),          # THM-633 c=2 at N=12
        (list(range(1,12))+[36], F(3,37)),          # THM-633 c=3 at N=12
        (list(range(1,12))+[13,24], F(1,14)),       # second tight family THM-1120
        (list(range(1,12))+[13,36], F(3,41)),       # F_3(13) mediant HYP-4506
        ([1,3,4,5,7,13,18], F(3,23)),               # N=7 gap member
        ([1,5,6,11,16,17], F(5,33)),                # N=6 gap member
        (list(range(1,11))+[12,33], F(3,35)),       # F_3(12) trichotomy HYP-4572
    ]
    ok = True
    for S, want in tests:
        got, arg = M_exact_arg(S)
        tag = "OK  " if got == want else "FAIL"
        if got != want: ok = False
        print(f"  GATE {tag} M({S}) = {got} (want {want}) arg={arg}")
    return ok

# ---------------- main ----------------
def main():
    print("== gates ==")
    if not gates():
        print("!! GATES FAILED -- no conclusions below are valid"); return

    print("\n== Part 1: K-ladder K_c(N) = {1..N-1} u {cN} ==")
    print("   rung prediction c/(cN+1); att=Y iff M == rung exactly")
    for N in range(6, 21):
        row = []
        for c in range(1, 7):
            V = list(range(1, N)) + [c*N]
            M, arg = M_exact_arg(V)
            rung = F(c, c*N+1)
            att = "Y" if M == rung else "n"
            row.append(f"c={c}:{M}{'=' if att=='Y' else '!'}r@{arg}")
        print(f"  N={N:2d}: " + "  ".join(row))
        if N == 13:
            V = list(range(1, 13)) + [26]
            M, arg = M_exact_arg(V)
            print(f"  ** Q1 SETTLE: M({{1..12,26}}) = {M} at (a,q)={arg}; "
                  f"2/27 realized: {M == F(2,27)}; slack-1 D=2 "
                  f"{'REALIZED (corrects opus-S396)' if M == F(2,27) else 'not realized here'}")

    print("\n== Part 2: F-ladder F_m(N) = {1..N-2, N, m(N-1)} ==")
    print("   binder b* := attained q at the maximizer minus m(N-1) when the")
    print("   maximizer denominator has that form; floor-tie and gap flags")
    for N in range(6, 21):
        floor = F(1, N+1)
        # first-gap window (1/(N+1), 2/(2N+1))
        top = F(2, 2*N+1)
        row = []
        for m in range(2, 6):
            V = list(range(1, N-1)) + [N, m*(N-1)]
            M, arg = M_exact_arg(V)
            a, q = arg
            bstar = q - m*(N-1)
            flag = ""
            if M == floor: flag = "TIE"
            elif floor < M < top: flag = "GAP"
            row.append(f"m={m}:{M}[b*={bstar}]{flag}@{arg}")
        print(f"  N={N:2d}: " + "  ".join(row))

    print("\n== Part 3: focus answers ==")
    for (lbl, V, cmp) in [
        ("Q2 n=8 tie?  F_2(7)={1..5,7,12}", [1,2,3,4,5,7,12], F(1,8)),
        ("Q2 n=9 hole  F_2(8)={1..6,8,14}", [1,2,3,4,5,6,8,14], F(1,9)),
        ("Q2 known     F_2(13)={1..11,13,24}", list(range(1,12))+[13,24], F(1,14)),
        ("K_2(13)={1..12,26}", list(range(1,13))+[26], F(2,27)),
        ("K_3(13)={1..12,39}", list(range(1,13))+[39], F(3,40)),
        ("K_2(14)={1..13,28}", list(range(1,14))+[28], F(2,29)),
    ]:
        M, arg = M_exact_arg(V)
        print(f"  {lbl}: M = {M} (cmp {cmp}, {'EQUAL' if M == cmp else 'DIFFERS'}) arg={arg}")

if __name__ == "__main__":
    main()
