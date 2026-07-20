#!/usr/bin/env python3
"""lrc_crossN_two_ladders_klein_S322.py -- klein-2026-07-19-S322.

CROSS-N UNIVERSALITY OF THE TWO LADDERS (N = number of speeds, LRC(N+1) case;
first window W_N = (1/(N+1), 2/(2N+1)), mediant 3/(3N+2)).

L-ladder:  L_N(m) = {1..N-2, N} u {(N-1)m}   (L_N(1) = {1..N} = the AP;
           conjecture: L_N(2) = the GW family = the tight locus's second member,
           M = 1/(N+1) plateau; m >= 3: M = m/((N-1)m + 5) when the (5,(N-1)m)
           witness survives the base -- the mediant is EXACTLY rung m = 3)
K-ladder:  K_N(s) = {1..N-1} u {Ns}          (K_N(1) = the AP again;
           Kravitz rungs M = s/(Ns+1))

Questions swept exactly (referee engine, complete-breakpoint exact_M):
 (Q1) Is M(L_N(2)) = 1/(N+1) for all N (the GW/tight-twin identification)?
 (Q2) For which N does L_N(3) attain the mediant 3/(3N+2) -- and does the set
      match THM-1284's populated band {6} u {N == 1 mod 6}?
 (Q3) Do the K-rungs hold s/(Ns+1) universally?
 (Q4) Where the mediant is NOT attained (N = 8..12), what is M(L_N(3)) --
      i.e., HOW does the ladder fail (plateau? off-formula value?).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from lrc14_subgap_referee_klein_S319 import exact_M

def fam_L(N, m):
    base = list(range(1, N-1)) + [N]
    far = (N-1)*m
    return None if far in base else sorted(base + [far])

def fam_K(N, s):
    base = list(range(1, N))
    far = N*s
    return None if far in base else sorted(base + [far])

print(f"{'N':>2} | {'L(2) M':>8} =1/(N+1)? | {'L(3) M':>8} mediant 3/(3N+2)? inW? | K(2), K(3) = s/(Ns+1)?")
for N in range(6, 14):
    floor = Fraction(1, N+1)
    med = Fraction(3, 3*N+2)
    wlo, whi = Fraction(1, N+1), Fraction(2, 2*N+1)
    l2 = exact_M(fam_L(N, 2))
    l3 = exact_M(fam_L(N, 3))
    k2 = exact_M(fam_K(N, 2))
    k3 = exact_M(fam_K(N, 3))
    inw = wlo < l3 < whi
    print(f"{N:>2} | {str(l2):>8} {'YES' if l2==floor else 'NO '}       | {str(l3):>8} "
          f"{'=MED' if l3==med else 'no  '} {'inW' if inw else '---'}      | "
          f"{str(k2):>6}{'=' if k2==Fraction(2,2*N+1) else '!'}, {str(k3):>6}{'=' if k3==Fraction(3,3*N+1) else '!'}")

print()
print("L-ladder full profile per N (m=1..6), value [P]=floor-plateau [F]=formula m/((N-1)m+5) [?]=other:")
for N in range(6, 14):
    floor = Fraction(1, N+1)
    row = []
    for m in range(1, 7):
        V = fam_L(N, m)
        if V is None:
            row.append("dup")
            continue
        Mv = exact_M(V)
        form = Fraction(m, (N-1)*m + 5)
        tag = "P" if Mv == floor else ("F" if Mv == form else "?")
        row.append(f"{Mv}[{tag}]")
    print(f"N={N:>2}: " + "  ".join(row))
print("done.")
