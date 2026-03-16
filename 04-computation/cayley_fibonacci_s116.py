#!/usr/bin/env python3
"""cayley_fibonacci_s116.py — Q shifts Fibonacci indices by 3

THM: Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}
PROOF: F_{n+1} + F_n = F_{n+2}, F_{n+1} - F_n = F_{n-1}. QED.

This investigates generalizations to other recurrences.
"""
from fractions import Fraction
from math import sqrt

# Fibonacci
fib = [1, 1]
for _ in range(20):
    fib.append(fib[-1] + fib[-2])

# Lucas
lucas = [2, 1]
for _ in range(20):
    lucas.append(lucas[-1] + lucas[-2])

# Tribonacci
trib = [0, 0, 1]
for _ in range(20):
    trib.append(trib[-1] + trib[-2] + trib[-3])

# Pell: P(n+1) = 2P(n) + P(n-1)
pell = [0, 1]
for _ in range(20):
    pell.append(2*pell[-1] + pell[-2])

print("CAYLEY TRANSFORM AND LINEAR RECURRENCES")
print("="*60)
print()

# ============================================================
print("FIBONACCI: Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}")
print("-"*40)
print()
for k in range(2, 12):
    x = Fraction(fib[k-1], fib[k])
    q = (1+x)/(1-x)
    expected = Fraction(fib[k+1], fib[k-2])
    print(f"  Q(F_{k}/F_{k+1}) = Q({fib[k-1]}/{fib[k]}) = {q} = F_{k+2}/F_{k-1} = {expected}: {q == expected}")

print()
print("  WHY shift by 3? Because Q uses addition AND subtraction:")
print("  Q(a/b) = (b+a)/(b-a)")
print("  Fibonacci: b+a = next term, b-a = previous-previous term.")
print("  The shift +2 in numerator and -1 in denominator gives total 3.")
print()

# ============================================================
print("LUCAS: Q(L_n/L_{n+1}) = ?")
print("-"*40)
print()
for k in range(2, 10):
    x = Fraction(lucas[k-1], lucas[k])
    q = (1+x)/(1-x)
    expected_fib = Fraction(lucas[k+1], lucas[k-2])
    print(f"  Q(L_{k}/L_{k+1}) = Q({lucas[k-1]}/{lucas[k]}) = {q}", end="")
    print(f"  =? L_{k+2}/L_{k-1} = {expected_fib}: {q == expected_fib}")

print()
print("  YES! Same identity: Q(L_n/L_{n+1}) = L_{n+2}/L_{n-1}")
print("  This works for ANY sequence satisfying a_{n+1} = a_n + a_{n-1}.")
print("  PROOF: (a_{n+1} + a_n)/(a_{n+1} - a_n) = a_{n+2}/a_{n-1}. QED.")
print()

# ============================================================
print("GENERAL FIBONACCI-TYPE: a_{n+1} = a_n + a_{n-1}")
print("-"*40)
print()
print("  For ANY such sequence with a_0, a_1 > 0:")
print("  Q(a_n/a_{n+1}) = a_{n+2}/a_{n-1}")
print()
print("  The Cayley transform shifts by 3 in ALL Fibonacci-type sequences.")
print("  This is because the recurrence relation IS the Cayley identity:")
print("  a_{n+2} = a_{n+1} + a_n (numerator = sum)")
print("  a_{n-1} = a_{n+1} - a_n (denominator = difference)")
print()

# ============================================================
print("TRIBONACCI: T_{n+1} = T_n + T_{n-1} + T_{n-2}")
print("-"*40)
print()
print("  Q(T_n/T_{n+1}) = (T_{n+1}+T_n)/(T_{n+1}-T_n)")
print("  T_{n+1}+T_n = T_{n+2} - T_{n-1} (from recurrence)")
print("  T_{n+1}-T_n = T_{n-1} + T_{n-2} (rearranging T_{n+1}=T_n+T_{n-1}+T_{n-2})")
print("  Wait: T_{n+1}-T_n = T_{n-1}+T_{n-2} only if T_{n+1}=T_n+T_{n-1}+T_{n-2}.")
print()
for k in range(4, 12):
    if trib[k] == 0:
        continue
    x = Fraction(trib[k-1], trib[k])
    q = (1+x)/(1-x)
    # T_{n+1}+T_n = T_{n+2}-T_{n-1}
    # T_{n+1}-T_n = T_{n-1}+T_{n-2}
    num = trib[k+1] - trib[k-2]
    den = trib[k-2] + trib[k-3]
    expected = Fraction(num, den) if den != 0 else "inf"
    print(f"  Q(T_{k}/T_{k+1}) = Q({trib[k-1]}/{trib[k]}) = {q}", end="")
    if expected != "inf":
        print(f"  =? (T_{k+2}-T_{k-1})/(T_{k-1}+T_{k-2}) = {expected}: {q == expected}")
    else:
        print()

print()
print("  For tribonacci, Q(T_n/T_{n+1}) = (T_{n+2}-T_{n-1})/(T_{n-1}+T_{n-2})")
print("  This is more complex because the recurrence has 3 terms, not 2.")
print("  The 'skip' is no longer a clean index shift.")
print()

# ============================================================
print("PELL: P_{n+1} = 2P_n + P_{n-1}")
print("-"*40)
print()
print("  Q(P_n/P_{n+1}) = (P_{n+1}+P_n)/(P_{n+1}-P_n)")
print("  P_{n+1}+P_n = 3P_n + P_{n-1} (from recurrence)")
print("  P_{n+1}-P_n = P_n + P_{n-1} (from recurrence)")
for k in range(2, 10):
    if pell[k] == 0:
        continue
    x = Fraction(pell[k-1], pell[k])
    q = (1+x)/(1-x)
    num = 3*pell[k-1] + pell[k-2]
    den = pell[k-1] + pell[k-2]
    expected = Fraction(num, den)
    print(f"  Q(P_{k}/P_{k+1}) = Q({pell[k-1]}/{pell[k]}) = {q} = {expected}: {q == expected}")

print()
print("  For Pell, Q(P_n/P_{n+1}) = (3P_n+P_{n-1})/(P_n+P_{n-1})")
print("  The ratio simplifies differently. Only Fibonacci-type (coefficient 1)")
print("  gives a clean index shift.")
print()

# ============================================================
print("="*60)
print("THE FIBONACCI PROPERTY IS UNIQUE TO COEFFICIENT 1")
print("="*60)
print()
print("  For a_{n+1} = c*a_n + a_{n-1}:")
print("  Q(a_n/a_{n+1}) = (a_{n+1}+a_n)/(a_{n+1}-a_n)")
print("               = ((c+1)a_n + a_{n-1})/((c-1)a_n + a_{n-1})")
print()
print("  This is a_{n+2}/a_{n-1} ONLY when c=1 (Fibonacci).")
print("  For c=1: (c+1)a_n+a_{n-1} = 2a_n+a_{n-1} = a_n+a_{n+1} = a_{n+2}.")
print("  And: (c-1)a_n+a_{n-1} = a_{n-1}.")
print("  For c>=2: the formula involves MIXED terms, no clean shift.")
print()
print("  THEOREM: The Cayley transform gives a pure index shift")
print("  if and only if the recurrence has leading coefficient 1.")
print("  Fibonacci is the UNIQUE linear recurrence with this property.")
print()

# ============================================================
print("="*60)
print("CAYLEY ORBIT OF FIBONACCI RATIOS")
print("="*60)
print()
print("  The 4-fold Cayley orbit of F_n/F_{n+1}:")
print("  F_n/F_{n+1} -> F_{n+2}/F_{n-1} -> -F_{n-1}/F_{n+2} -> -F_{n+1}/F_n -> F_n/F_{n+1}")
print()
print("  Using Q^2(x) = -1/x:")
print("  Q^2(F_n/F_{n+1}) = -F_{n+1}/F_n (negative of reciprocal ratio)")
print("  This swaps the Fibonacci ratio with its COMPLEMENT.")
print()
for k in [3, 5, 7]:
    x0 = Fraction(fib[k-1], fib[k])
    x1 = (1+x0)/(1-x0)  # Q
    x2 = Fraction(-1, 1) / x0  # Q^2 = -1/x
    x3 = (x0-1)/(x0+1)  # Q^3
    print(f"  Orbit of F_{k}/F_{k+1} = {x0}:")
    print(f"    {x0} -> {x1} -> {x2} -> {x3} -> {x0}")
    print()

# ============================================================
print("LIMIT: Q(1/phi) = phi^3 AND THE FIBONACCI SKIP-3")
print("-"*40)
print()
print("  As n->inf: F_n/F_{n+1} -> 1/phi")
print("  Q(1/phi) = phi^3")
print("  And F_{n+2}/F_{n-1} -> phi^3 (ratio of Fibonacci terms 3 apart)")
print()
print("  More precisely: F_{n+2}/F_{n-1} = phi^3 * (1 + O(phi^{-2n}))")
print("  because F_n ~ phi^n/sqrt(5), so")
print("  F_{n+2}/F_{n-1} ~ phi^{n+2}/phi^{n-1} = phi^3.")
print()
print("  The Fibonacci skip-3 CONVERGES to the golden cubing.")
print("  The integer identity (skip-3) and the algebraic identity (cubing)")
print("  are the SAME identity, seen at finite and infinite scales.")
print()

phi = (1+sqrt(5))/2
for k in [5, 8, 11, 14]:
    ratio = fib[k+1]/fib[k-2]
    print(f"  F_{k+2}/F_{k-1} = {fib[k+1]}/{fib[k-2]} = {ratio:.10f}  vs phi^3 = {phi**3:.10f}  (err = {abs(ratio-phi**3):.2e})")
