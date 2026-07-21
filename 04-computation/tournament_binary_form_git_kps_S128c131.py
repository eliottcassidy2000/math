#!/usr/bin/env python3
"""tournament_binary_form_git_kps_S128c131.py -- kind-pasteur-2026-07-20-S128c131

BINARY FORMS AND TOURNAMENTS: transitivity is the GIT nullcone.

The characteristic polynomial char_A(x) of the adjacency matrix, homogenized to
char_A(x,y) = y^n char_A(x/y), is a BINARY FORM of degree n -- a point in Sym^n(C^2), on which
SL_2 acts.  The map  T |-> char_A  sends a tournament to a binary form.  Claim:

  * TRANSITIVE  <=>  char_A = x^n  <=>  the GIT NULLCONE / maximally-UNSTABLE binary form
    (a single root of multiplicity n = n > n/2).  [A transitive is strictly triangular =
    nilpotent = spectrum {0}, THM-895's lambda=0 <=> transitive.]
  * The trace moments tr(A^k) are the SL_2-INVARIANTS (power sums = coefficients of the form).
  * GIT STABILITY of char_A stratifies tournaments: unstable (some root mult > n/2),
    semistable, stable (all roots mult < n/2).  This computes the stratification and asks
    which tournaments are unstable, and whether the transitive is the ONLY nullcone member.

By Hilbert-Mumford for binary forms (classical): a degree-n form is UNSTABLE iff it has a root
of multiplicity > n/2; STABLE iff every root has multiplicity < n/2; strictly SEMISTABLE iff
max multiplicity = n/2 exactly.  Root multiplicities of char_A = powers of its irreducible
rational factors.
"""
import sys
import numpy as np
import sympy as sp
from collections import Counter

x = sp.Symbol('x')


def from_bits(bits, n):
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A


def charpoly_int(A, n):
    return sp.Poly(sp.Matrix(A).charpoly(x).as_expr(), x)


def max_root_mult(p):
    """Max multiplicity of any single root = max power of an irreducible rational factor."""
    fl = sp.factor_list(p.as_expr())
    return max((e for _, e in fl[1]), default=1)


def c3(A, n):
    """number of 3-cycles = (1/3) sum over triples of cyclic triangles; via tr contributions."""
    M = np.array(A)
    # cyclic triangles count = trace(A^3)/3
    return int(np.trace(np.linalg.matrix_power(M, 3)) // 3)


print("=" * 84)
print("GIT STABILITY of char_A over tournaments (transitivity = nullcone x^n)")
print("=" * 84)
for n in (3, 4, 5, 6):
    m = n * (n - 1) // 2
    strat = Counter()          # stability class
    multdist = Counter()       # max root multiplicity
    transitive_form = None
    unstable_examples = []
    for bits in range(1 << m):
        A = from_bits(bits, n)
        p = charpoly_int(A, n)
        mm = max_root_mult(p)
        multdist[mm] += 1
        # GIT class
        if mm > n / 2:
            cls = "unstable"
        elif mm == n / 2:
            cls = "semistable(strict)"
        else:
            cls = "stable"
        strat[cls] += 1
        cyc = c3(A, n)
        if cyc == 0:            # transitive (acyclic)
            transitive_form = sp.factor(p.as_expr())
        if cls == "unstable" and cyc > 0 and len(unstable_examples) < 3:
            unstable_examples.append((sp.factor(p.as_expr()), cyc))
    print("  n=%d : transitive char_A = %s   (= x^%d ? %s)"
          % (n, transitive_form, n, transitive_form == x**n))
    print("        GIT strata: %s" % dict(strat))
    print("        max-root-multiplicity distribution: %s" % dict(sorted(multdist.items())))
    if unstable_examples:
        print("        NON-transitive UNSTABLE forms (mult > n/2, but cyclic): %s"
              % [(str(f), "c3=%d" % c) for f, c in unstable_examples])
    else:
        print("        NON-transitive unstable forms: NONE -> transitive is the unique nullcone")
    sys.stdout.flush()

print()
print("=" * 84)
print("THE HALF-DICTIONARY BRIDGE  char_A <-> char_S  via  x -> 2x+1")
print("=" * 84)
print("  S = 2A - J + I, so on regular tournaments char_S(x) = 2^n char_A((x-1)/2) up to the")
print("  all-ones direction; the skew form char_S is EVEN (spectrum +-i*lambda, symmetric),")
print("  the nullcone x^n of char_A maps to (x+1)^n-flavored under the shift. Illustrate n=5:")
for bits in [0, (1 << 10) - 1]:   # transitive-ish and complement
    n = 5
    A = from_bits(bits, n)
    S = (np.array(A) - np.array(A).T)
    pA = sp.factor(sp.Matrix(A).charpoly(x).as_expr())
    pS = sp.factor(sp.Matrix(sp.Matrix(S.tolist())).charpoly(x).as_expr())
    print("  bits=%-4d : char_A = %s" % (bits, pA))
    print("            char_S = %s  (even in x: %s)"
          % (pS, sp.simplify(sp.expand(sp.Matrix(sp.Matrix(S.tolist())).charpoly(x).as_expr()).subs(x, -x)
                             - (-1)**n * sp.expand(sp.Matrix(sp.Matrix(S.tolist())).charpoly(x).as_expr())) == 0))
print()
print("  READING: T -> char_A is a map to degree-n binary forms; transitivity = GIT nullcone")
print("  x^n; trace moments = SL_2 invariants; the fibers are co-spectral classes, which is")
print("  exactly where H (the #P datum) can split (THM-1780, n=6). The characteristic form")
print("  sees the spectral/invariant-theoretic tournament and forgets the permanent.")
