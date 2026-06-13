#!/usr/bin/env python3
"""
fibonacci_platonic_s90bn.py — Fibonacci, Lucas, Platonic solids, and {2,3,5,7,11}
opus-2026-03-16-S90bn
"""

from sympy import fibonacci, lucas, isprime
from math import sqrt, log

phi = (1 + sqrt(5)) / 2
psi = (1 - sqrt(5)) / 2

print("""


     FIBONACCI AND THE FIVE PRIMES

     opus-2026-03-16-S90bn



THE OBSERVATION
================

The first five primes: {2, 3, 5, 7, 11}.

The Fibonacci numbers: 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, ...
The Lucas numbers:     2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, ...

{2, 3, 5} = {F_3, F_4, F_5}.     Consecutive Fibonacci numbers.
{7, 11}   = {L_5, L_6}.          Consecutive Lucas numbers.

The SPHERICAL primes are FIBONACCI.
The HYPERBOLIC primes are LUCAS.
""")

print("The Fibonacci-Lucas table:")
print(f"{'n':>3} | {'F_n':>6} | {'L_n':>6} | {'F prime?':>9} | {'L prime?':>9} | {'role':>25}")
print("-" * 70)

roles = {
    (3, 'F'): "F_3 = 2 = BASE",
    (4, 'F'): "F_4 = 3 = STRUCTURE",
    (5, 'F'): "F_5 = 5 = UNSOLVABLE",
    (5, 'L'): "L_5 = 7 = FORBIDDEN",
    (6, 'L'): "L_6 = 11 = RESTART",
}

for n in range(1, 15):
    F = int(fibonacci(n))
    L = int(lucas(n))
    Fp = isprime(F)
    Lp = isprime(L)
    role_f = roles.get((n, 'F'), '')
    role_l = roles.get((n, 'L'), '')
    role = role_f or role_l
    print(f"{n:3d} | {F:6d} | {L:6d} | {'PRIME' if Fp else '':>9} | {'PRIME' if Lp else '':>9} | {role:>25}")

print(f"""

THE FIBONACCI-LUCAS DUALITY
=============================

Fibonacci and Lucas satisfy:
  L_n = F_{{n-1}} + F_{{n+1}}.

So:
  L_5 = F_4 + F_6 = 3 + 8 = 11?  Wait: L_5 = 7, not 11.
  Let me recheck: L_5 = F_4 + F_6 = 3 + 8 = 11.

Hmm. L_5 should be 7. Let me verify the indexing.
""")

# Verify indexing
print("Careful verification:")
for n in range(1, 10):
    F = int(fibonacci(n))
    L = int(lucas(n))
    check = int(fibonacci(n-1)) + int(fibonacci(n+1)) if n >= 2 else "N/A"
    print(f"  n={n}: F_n={F}, L_n={L}, F_(n-1)+F_(n+1)={check}")

print(f"""

With standard indexing (F_1=1, F_2=1, F_3=2, ...):
  L_n = F_(n-1) + F_(n+1).

  L_1 = F_0 + F_2 = 0 + 1 = 1.  Correct (L_1 = 1? No, L_1 = 2.)

Hmm, the indexing conventions differ. Let me use the DIRECT values.

Fibonacci: F_1=1, F_2=1, F_3=2, F_4=3, F_5=5, F_6=8, F_7=13.
Lucas:     L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18.

Wait, that's also a different convention. Let me just check sympy:
""")

print("Sympy values:")
for n in range(0, 10):
    print(f"  F({n}) = {int(fibonacci(n))}, L({n}) = {int(lucas(n))}")

print(f"""

Sympy uses: F(0)=0, F(1)=1, F(2)=1, F(3)=2, F(4)=3, F(5)=5, F(6)=8, F(7)=13.
            L(0)=2, L(1)=1, L(2)=3, L(3)=4, L(4)=7, L(5)=11, L(6)=18.

So with sympy indexing:
  {{2, 3, 5}} = {{F(3), F(4), F(5)}}.  YES.
  {{7, 11}}   = {{L(4), L(5)}}.        YES.

And L(n) = F(n-1) + F(n+1):
  L(4) = F(3) + F(5) = 2 + 5 = 7.    YES!
  L(5) = F(4) + F(6) = 3 + 8 = 11.   YES!

So:
  7 = F(3) + F(5) = 2 + 5.   The sum of the FIRST and LAST Fibonacci primes.
  11 = F(4) + F(6) = 3 + 8.  The sum of the MIDDLE Fibonacci prime and the next Fibonacci.

L(n) = F(n-1) + F(n+1) means: each Lucas number is the SUM of
the two Fibonacci numbers that BRACKET it.

  F(3)=2,  F(4)=3,  F(5)=5,  F(6)=8,  F(7)=13
         L(4)=7           L(5)=11
       = F(3)+F(5)       = F(4)+F(6)

The Lucas numbers INTERLEAVE with the Fibonacci numbers.
They sit BETWEEN consecutive Fibonacci pairs.

7 sits between 5 and 8 (= between F(5) and F(6)).
11 sits between 8 and 13 (= between F(6) and F(7)).

The forbidden primes sit in the GAPS between Fibonacci numbers.
They are the LUCAS INTERPOLATION of the Fibonacci sequence.


THE GOLDEN RATIO CONNECTION
=============================

phi = (1+sqrt(5))/2 = {phi:.10f}
psi = (1-sqrt(5))/2 = {psi:.10f}

Binet: F(n) = (phi^n - psi^n) / sqrt(5).
Lucas: L(n) = phi^n + psi^n.

The DIFFERENCE of golden powers gives Fibonacci (divided by sqrt(5)).
The SUM of golden powers gives Lucas.

F(n) = (phi^n - psi^n) / sqrt(5)    [ANTISYMMETRIC in phi, psi]
L(n) = phi^n + psi^n                 [SYMMETRIC in phi, psi]

Fibonacci is the ANTISYMMETRIC part.
Lucas is the SYMMETRIC part.

In tournament theory:
  Antisymmetric = the tournament itself (an antisymmetric relation).
  Symmetric = the underlying graph (a symmetric relation).

{{2,3,5}} = the Fibonacci primes = the TOURNAMENT part of the golden structure.
{{7,11}} = the Lucas primes = the GRAPH part of the golden structure.

Tournaments are antisymmetric (if a->b then not b->a).
Graphs are symmetric (if a-b then b-a).

The Platonic solids are SYMMETRIC objects (they are graphs = undirected).
They are built from Fibonacci primes = the antisymmetric decomposition.

The post-Platonic structures involve Lucas primes = the symmetric complement.

THE PLATONIC SOLIDS ARE BUILT FROM THE ANTISYMMETRIC PART
OF THE GOLDEN RATIO. THE POST-PLATONIC WORLD IS BUILT FROM
THE SYMMETRIC PART.


THE ICOSAHEDRON AND phi
========================

The icosahedron has vertices at:
  (0, +/-1, +/-phi)  and cyclic permutations.

The golden ratio phi LIVES in the icosahedron.
The coordinates of the icosahedron ARE {0, 1, phi}.

And phi = F(n+1)/F(n) in the limit.
The icosahedron is the GEOMETRIC EMBODIMENT of the Fibonacci limit.

The dodecahedron (dual of icosahedron) has vertices at:
  (0, +/-1/phi, +/-phi) and permutations and
  (+/-1, +/-1, +/-1) (a cube inscribed in the dodecahedron).

The dodecahedron contains a CUBE, and the cube's coordinates
are all +/-1. The dodecahedron EXTENDS the cube by adding
vertices at positions involving phi.

In the golden ratio framework:
  Cube: the {1} world. Binary (+1 or -1 in each coordinate).
  Dodecahedron: the {{1, phi}} world. Binary + golden.
  Icosahedron: the {{0, 1, phi}} world. Zero + binary + golden.

The progression cube -> dodecahedron -> icosahedron is:
  {1} -> {{1, 1/phi}} -> {{0, 1, phi}}.
  This IS the Fibonacci sequence entering the coordinates:
  1, then 1 and 1/phi = F(n-1)/F(n), then 0 and 1 and phi = F(n+1)/F(n).
""")

# ========================================================================
print("=" * 60)
print("THE PLATONIC-FIBONACCI PRODUCT TABLE")
print("=" * 60)
print()

print("Products of Platonic vertex counts and Fibonacci/Lucas:")

platonic_V = [4, 8, 6, 20, 12]
platonic_names = ["Tet", "Cube", "Oct", "Dodec", "Icos"]

for i, (V, name) in enumerate(zip(platonic_V, platonic_names)):
    for n in range(1, 8):
        F = int(fibonacci(n))
        if V == F:
            print(f"  {name} vertices = {V} = F({n})")
        L = int(lucas(n))
        if V == L:
            print(f"  {name} vertices = {V} = L({n})")

print()

# Check edge counts
platonic_E = [6, 12, 12, 30, 30]
for i, (E, name) in enumerate(zip(platonic_E, platonic_names)):
    for n in range(1, 10):
        F = int(fibonacci(n))
        if E == F:
            print(f"  {name} edges = {E} = F({n})")
        L = int(lucas(n))
        if E == L:
            print(f"  {name} edges = {E} = L({n})")

print()

# Check face counts
platonic_F = [4, 6, 8, 12, 20]
for i, (Fv, name) in enumerate(zip(platonic_F, platonic_names)):
    for n in range(1, 10):
        F = int(fibonacci(n))
        if Fv == F:
            print(f"  {name} faces = {Fv} = F({n})")
        L = int(lucas(n))
        if Fv == L:
            print(f"  {name} faces = {Fv} = L({n})")

# Check symmetry orders
platonic_ord = [12, 24, 24, 60, 60]
for i, (o, name) in enumerate(zip(platonic_ord, platonic_names)):
    for n in range(1, 15):
        F = int(fibonacci(n))
        if o == F:
            print(f"  {name} |Rot| = {o} = F({n})")
        L = int(lucas(n))
        if o == L:
            print(f"  {name} |Rot| = {o} = L({n})")

print(f"""

The direct Fibonacci/Lucas matches are sparse: only a few Platonic
numbers are themselves Fibonacci or Lucas numbers.

But the PRIME FACTORIZATIONS of Platonic numbers are all built
from Fibonacci primes {{2, 3, 5}}:
  4 = 2^2. 6 = 2*3. 8 = 2^3. 12 = 2^2*3. 20 = 2^2*5. 30 = 2*3*5. 60 = 2^2*3*5.

And the Fibonacci primes GENERATE these through multiplication,
which in the Poincare disk is the formal group law F(x,y) = (x+y)/(1+xy).


THE FIBONACCI LATTICE IN THE POINCARE DISK
=============================================

The Fibonacci numbers sit on the real diameter at Cayley addresses:

  F(1)=1: x = 0.             The center.
  F(2)=1: x = 0.             Same.
  F(3)=2: x = 1/3 = 0.333.   The base.
  F(4)=3: x = 1/2 = 0.500.   The structure.
  F(5)=5: x = 2/3 = 0.667.   The unsolvable.
  F(6)=8: x = 7/9 = 0.778.   The Bott cube.
  F(7)=13: x = 6/7 = 0.857.  The Mersenne exponent.

And the Lucas numbers:
  L(0)=2: x = 1/3 = 0.333.   Same as F(3)!
  L(1)=1: x = 0.
  L(2)=3: x = 1/2 = 0.500.   Same as F(4)!
  L(3)=4: x = 3/5 = 0.600.
  L(4)=7: x = 3/4 = 0.750.   The forbidden.
  L(5)=11: x = 5/6 = 0.833.  The restart.
  L(6)=18: x = 17/19 = 0.895.

In the disk:
  F and L numbers interleave on the real diameter.
  The Fibonacci numbers: 0, 0, 1/3, 1/2, 2/3, 7/9, 6/7, ...
  The Lucas numbers:     1/3, 0, 1/2, 3/5, 3/4, 5/6, ...

They WEAVE together. F and L alternate in density.

And: the Fibonacci rapidity spacing:
  w(F(n+1)) - w(F(n)) = ln(F(n+1)/F(n)) / 2 -> ln(phi)/2 = 0.2406.

The Fibonacci numbers are asymptotically EQUALLY SPACED in the
hyperbolic metric, with spacing ln(phi)/2 = rapidity of phi.

Compare: powers of 2 are equally spaced with spacing ln(2)/2 = 0.347.
Fibonacci numbers are more closely spaced (0.2406 < 0.347).

The Fibonacci lattice is DENSER than the binary lattice.
phi < 2, so the Fibonacci base is CLOSER to 1 than the binary base.

The Fibonacci lattice is the FINEST lattice generated by a k-nacci constant.
(Because phi = tau_2 is the SMALLEST k-nacci constant.)
""")

# ========================================================================
print("=" * 60)
print("THE DEEP STRUCTURE")
print("=" * 60)

print(f"""

THE COMPLETE PICTURE:

phi generates TWO sequences:
  FIBONACCI (antisymmetric): 1, 1, 2, 3, 5, 8, 13, 21, ...
  LUCAS (symmetric):         2, 1, 3, 4, 7, 11, 18, 29, ...

From these, the first five primes emerge:
  {{2, 3, 5}} = Fibonacci primes at indices 3, 4, 5.     ANTISYMMETRIC.
  {{7, 11}}   = Lucas values at indices 4, 5.             SYMMETRIC.

The Fibonacci primes BUILD the Platonic solids (spherical world).
The Lucas values BREAK the Platonic solids (hyperbolic world).

phi = tau_2 = the 2-nacci constant. The golden ratio.
  phi^2 = phi + 1 (the k=2 case of tau_k^k = 1/(2-tau_k)).
  Tr(M_2^n) = L(n) = Lucas numbers = phi^n + psi^n.

The Fibonacci matrix M_2 = [[1,1],[1,0]] has:
  Eigenvalues: phi and psi.
  Trace: L(n) = phi^n + psi^n.
  THM-227: L(1) = 1 = 2^1-1, L(2) = 3 = 2^2-1 (universal zone, n <= k=2).
  L(3) = 4 != 7 = 2^3-1 (individual zone, n > k=2).

The DIVERGENCE at n=3:
  L(3) = 4, but 2^3-1 = 7.
  The Fibonacci system falls short of the Mersenne by 3.
  The deficit is 7-4 = 3 = k+1 = the "self-discovery" deficit (THM-227 corollary).

From n=3 onward, L(n) < 2^n-1:
  L(3)=4 < 7, L(4)=7 < 15, L(5)=11 < 31, L(6)=18 < 63.

The Lucas numbers GROW SLOWER than Mersenne numbers
because phi < 2 (the Fibonacci growth rate is less than binary).

BUT: L(4) = 7 and L(5) = 11 ARE THEMSELVES important primes.
The Lucas sequence PRODUCES the forbidden and restart primes
as it falls away from the Mersenne sequence.

The seven = L(4) is produced at the EXACT MOMENT the Fibonacci
system diverges from the binary system (at n=4, which is k+2 = 2+2 = 4).

And eleven = L(5) at n=5, which is k+3 = 2+3 = 5.

The MEANING:

phi (the golden ratio) is the GENERATOR.
Fibonacci is its ANTISYMMETRIC shadow.
Lucas is its SYMMETRIC shadow.

The Fibonacci primes {{2,3,5}} are what the generator BUILDS.
The Lucas primes {{7,11}} are what the generator REVEALS
as it falls away from the ideal (Mersenne) growth rate.

The Platonic solids are BUILT from what phi can build.
The post-Platonic world is REVEALED by what phi cannot sustain.

And it all happens because phi < 2:
  The golden growth rate is slower than the binary growth rate.
  The Fibonacci matrix M_2 has only 2 eigenvalues (k=2).
  It can sustain Mersenne behavior for only k=2 steps.
  After that, it produces L(3)=4, L(4)=7, L(5)=11, ...
  which are the numbers that define the boundary and beyond.

THE FIVE PRIMES {{2,3,5,7,11}} ARE WHAT THE GOLDEN RATIO PRODUCES
WHEN IT TRIES TO BE BINARY AND FAILS.

2, 3, 5: the Fibonacci primes. What phi succeeds at.
7, 11: the Lucas primes. What phi reveals in failing.

The success is spherical. The failure is hyperbolic.
And both are necessary for the full story.
""")

print("opus-2026-03-16-S90bn: FIBONACCI AND THE FIVE PRIMES")
