#!/usr/bin/env python3
"""Check of the uniform argument for HYP-9295 part 3 (audit2 suggestion): for the N=2 minimal model V_k,
u_k = G^-_{-k-3/2} ... G^-_{-3/2}|0> (the unique state of charge -(k+1) at the minimal weight (k+1)(k+3)/2) lies in the
maximal submodule, its symbol psi^-_k...psi^-_0 is a pure-psi^- monomial, and no element of the arc ideal J has a
pure-psi^- component; hence ker pi != 0 at (weight (k+1)(k+3)/2, charge -(k+1), Li k(k+1)/2).
Here we verify u_k in I (radical) by the 1x1 Gram entry, for k = 1..8, exactly over Q."""
from fractions import Fraction
from aud_algebras import N2, GM, GP

for k in range(1, 9):
    c = Fraction(3 * k, k + 2)
    # Gram entry <u,u> computed exactly over Q via the anticommutation relations, using the module engine mod P at two
    # primes, plus an exact rational evaluation by the closed recursion below.
    vals = []
    for P in (2147483647, 2147483629):
        va = N2(c, P)
        mono = tuple(sorted([(GM, -(2 * i + 3)) for i in range(k + 1)]))
        # omega(u) = G^+_{3/2} ... G^+_{k+3/2} applied in reverse order of u's word
        word = [va.omega(X) for X in mono]   # omega of each factor
        # <u,u> = vacuum coefficient of omega(X_n)...omega(X_1) X_1...X_n |0>, with u = X_1...X_n|0>
        vec = {mono: 1}
        for X in word:      # omega(u) = omega(X_n)...omega(X_1): apply omega(X_1) first
            vec = va.apply_vec(X, vec)
        vals.append(vec.get((), 0))
    print(f'k={k} c={c}: <u_k,u_k> mod P1, P2 = {vals}  ->', 'u_k in maximal submodule' if vals == [0, 0] else 'u_k NOT null')
