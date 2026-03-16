#!/usr/bin/env python3
"""
solve_via_sa_s90bq.py — Solving problems via the S/A master decomposition
opus-2026-03-16-S90bq

The isomorphism: every problem that can be stated in one S/A costume
can be TRANSLATED to another costume and solved there.

cosh/sinh <-> Lucas/Fibonacci <-> graph/tournament <-> even/odd <-> real/imaginary

Solve in whichever costume makes the problem EASIEST, then translate back.
"""

from math import cosh, sinh, tanh, atanh, log, sqrt, pi, exp, cos, sin
from sympy import fibonacci, lucas, isprime, factorint, symbols, simplify
from sympy import cosh as scosh, sinh as ssinh, tanh as stanh
from sympy import cos as scos, sin as ssin, exp as sexp
from sympy import Rational, sqrt as ssqrt, oo, Sum, Function

x = symbols('x')

print("""


     SOLVING PROBLEMS VIA THE S/A DECOMPOSITION

     opus-2026-03-16-S90bq



PRINCIPLE: Every equation involving cosh/sinh can be solved
as an equation involving Lucas/Fibonacci (or even/odd functions,
or graphs/tournaments, or real/imaginary parts).

The isomorphisms are EXACT. Solving in one costume and
translating back gives the SAME answer.


PROBLEM 1: SOLVE cosh(x) = 3
==============================

In the cosh/sinh costume:
  cosh(x) = 3.
  (e^x + e^{-x})/2 = 3.
  e^x + e^{-x} = 6.
  Let u = e^x: u + 1/u = 6. So u^2 - 6u + 1 = 0.
  u = (6 +/- sqrt(32))/2 = 3 +/- 2*sqrt(2).
  x = ln(3 + 2*sqrt(2)) or x = ln(3 - 2*sqrt(2)) = -ln(3+2*sqrt(2)).
  x = +/- ln(3 + 2*sqrt(2)) = +/- 1.7627...
""")

val = log(3 + 2*sqrt(2))
print(f"  x = +/- {val:.10f}")
print(f"  Verify: cosh({val:.6f}) = {cosh(val):.10f}")

print(f"""

Now solve the SAME problem in the Cayley costume:
  cosh(x) = 3 means the symmetric part of e^x equals 3.
  e^x = cosh(x) + sinh(x). If cosh = 3, then sinh = sqrt(9-1) = sqrt(8) = 2*sqrt(2).
  (Because cosh^2 - sinh^2 = 1.)
  So e^x = 3 + 2*sqrt(2).
  The Cayley transform: n = e^(2w) where w = x/... wait.

  Actually: cosh(x) = 3 means (n + 1/n)/2 = 3 where n = e^x.
  So n + 1/n = 6. This is the SAME equation as before.
  n = 3 + 2*sqrt(2) = (1+sqrt(2))^2.

  The coupling: x_n = (n-1)/(n+1) = (3+2*sqrt(2)-1)/(3+2*sqrt(2)+1)
              = (2+2*sqrt(2))/(4+2*sqrt(2)) = 2(1+sqrt(2))/(2(2+sqrt(2)))
              = (1+sqrt(2))/(2+sqrt(2)).

  Rationalize: multiply by (2-sqrt(2))/(2-sqrt(2)):
  = (1+sqrt(2))(2-sqrt(2)) / ((2+sqrt(2))(2-sqrt(2)))
  = (2-sqrt(2)+2*sqrt(2)-2) / (4-2)
  = sqrt(2) / 2 = 1/sqrt(2).

  So: the coupling is x = 1/sqrt(2) = sqrt(2)/2.
  And tanh(w) = 1/sqrt(2) where 2w = x_original.

  VERIFY: tanh(arctanh(1/sqrt(2))*2) ... hmm, let me just check.
""")

coupling = 1/sqrt(2)
w = atanh(coupling)
print(f"  Coupling x = 1/sqrt(2) = {coupling:.10f}")
print(f"  Rapidity w = arctanh(x) = {w:.10f}")
print(f"  cosh(2w) = {cosh(2*w):.10f}  (should be 3)")
print(f"  n = e^(2w) = {exp(2*w):.10f}  (should be 3+2*sqrt(2) = {3+2*sqrt(2):.10f})")

print(f"""

BEAUTIFUL. cosh(x) = 3 translates to:
  "The coupling is 1/sqrt(2) = the Cayley address of (1+sqrt(2))^2."
  The number (1+sqrt(2))^2 = 3+2*sqrt(2) is the SILVER RATIO squared.

The silver ratio delta = 1+sqrt(2) = 2.414... is the analogue of the
golden ratio for the equation x^2 - 2x - 1 = 0.
  phi: x^2 - x - 1 = 0.   phi = (1+sqrt(5))/2.
  delta: x^2 - 2x - 1 = 0. delta = 1+sqrt(2).

cosh(x) = 3 has the SILVER RATIO as its answer.
Translated to the Cayley line: coupling = 1/sqrt(2).


PROBLEM 2: FIND n SUCH THAT THE S/A RATIO IS EXACTLY 1/2
===========================================================

The S/A ratio at n is: tanh(ln(n)) = (n^2-1)/(n^2+1).
Set this equal to 1/2:
  (n^2-1)/(n^2+1) = 1/2.
  2(n^2-1) = n^2+1.
  2n^2 - 2 = n^2 + 1.
  n^2 = 3.
  n = sqrt(3).

The number with S/A ratio exactly 1/2 is sqrt(3).
""")

n_val = sqrt(3)
ratio = (n_val**2 - 1)/(n_val**2 + 1)
print(f"  n = sqrt(3) = {n_val:.10f}")
print(f"  S/A ratio = (n^2-1)/(n^2+1) = {ratio:.10f}")
print(f"  Coupling of sqrt(3): x = (sqrt(3)-1)/(sqrt(3)+1) = {(n_val-1)/(n_val+1):.10f}")
print(f"  = 2 - sqrt(3) = {2-sqrt(3):.10f}")

print(f"""

The coupling of sqrt(3) is 2-sqrt(3) = 0.2679...
And 2-sqrt(3) = 1/(2+sqrt(3)) = the reciprocal of the silver ratio-ish.

More importantly: sqrt(3) is the MAGNITUDE of the pure-imaginary
coupling at angle 60 degrees: tan(60 deg) = sqrt(3).

So the point where S/A = 1/2 is exactly where the ROTATIONAL content
(imaginary coupling = sqrt(3)) meets the boost content.


PROBLEM 3: TRANSLATE A FIBONACCI IDENTITY TO HYPERBOLIC
=========================================================

Fibonacci identity: F(2n) = F(n) * L(n).

Translate to hyperbolic:
  F(n) = sinh(n*w) / sinh(w)  where w = ln(phi)/2... no.
  Actually, using Binet:
  F(n) = (phi^n - psi^n)/sqrt(5).
  L(n) = phi^n + psi^n.

  F(2n) = (phi^(2n) - psi^(2n))/sqrt(5)
        = (phi^n - psi^n)(phi^n + psi^n)/sqrt(5)   [difference of squares]
        = F(n) * L(n).

  In hyperbolic language, let a = phi^n, b = psi^n. Then:
    F(n) = (a - b)/sqrt(5).  [antisymmetric in a,b]
    L(n) = a + b.            [symmetric in a,b]
    F(2n) = (a^2 - b^2)/sqrt(5) = ((a-b)(a+b))/sqrt(5) = F(n)*L(n).

  The identity F(2n) = F(n)*L(n) says:
    ANTISYMMETRIC(doubled) = ANTISYMMETRIC(original) * SYMMETRIC(original).

  In cosh/sinh language:
    sinh(2x) = 2*sinh(x)*cosh(x).

  This is the DOUBLE-ANGLE formula! And F(2n) = F(n)*L(n) is its
  Fibonacci/Lucas version.

  sinh(2x) = 2*sinh(x)*cosh(x)
  F(2n) = F(n) * L(n)

  These are the SAME identity in different costumes.
  The factor of 2 vs absence of 2 is because F and L are normalized
  differently from sinh and cosh (F uses 1/sqrt(5), L uses 1).
""")

# Verify
print("Verify F(2n) = F(n)*L(n):")
for n in range(1, 10):
    F2n = int(fibonacci(2*n))
    Fn = int(fibonacci(n))
    Ln = int(lucas(n))
    print(f"  n={n}: F({2*n}) = {F2n}, F({n})*L({n}) = {Fn}*{Ln} = {Fn*Ln}, match: {F2n == Fn*Ln}")

print(f"""

PROBLEM 4: TRANSLATE THE ADDITION THEOREM
===========================================

cosh/sinh addition:
  cosh(a+b) = cosh(a)cosh(b) + sinh(a)sinh(b).   [S of sum = S*S + A*A]
  sinh(a+b) = sinh(a)cosh(b) + cosh(a)sinh(b).   [A of sum = A*S + S*A]

Fibonacci/Lucas addition:
  L(m+n) = L(m)L(n) + 5*F(m)F(n)  ... divided by 2.
  Actually: L(m+n) = (L(m)L(n) + 5*F(m)F(n)) / 2.
  And: F(m+n) = (F(m)L(n) + L(m)F(n)) / 2.

  These ARE the cosh/sinh addition formulas!
  S(m+n) = (S(m)S(n) + 5*A(m)A(n)) / 2.    [the 5 comes from the sqrt(5) normalization]
  A(m+n) = (A(m)S(n) + S(m)A(n)) / 2.

Formal group addition:
  F(x_m, x_n) = (x_m + x_n) / (1 + x_m*x_n) = x_(mn).

  Here x = tanh(w) = sinh(w)/cosh(w) = A/S.
  The formal group adds the A/S RATIOS using the tanh addition law.

  tanh(a+b) = (tanh(a) + tanh(b)) / (1 + tanh(a)*tanh(b)).

  This IS F(x_a, x_b) = (x_a + x_b) / (1 + x_a*x_b).

  So: the formal group law is the ADDITION THEOREM for the A/S ratio.

  cosh/sinh: adds the S and A parts separately (4 terms).
  tanh: adds the A/S ratios directly (1 formula).
  Fibonacci/Lucas: adds the F and L values (4 terms, with factor 5).
  Formal group: adds the couplings (1 formula).

  FOUR COSTUMES FOR THE SAME ADDITION THEOREM.
""")

# Verify Fibonacci/Lucas addition
print("Verify F(m+n) = (F(m)L(n) + L(m)F(n)) / 2:")
for m, n in [(2,3), (3,4), (5,7), (3,8)]:
    Fmn = int(fibonacci(m+n))
    Fm = int(fibonacci(m))
    Fn = int(fibonacci(n))
    Lm = int(lucas(m))
    Ln = int(lucas(n))
    pred = (Fm * Ln + Lm * Fn) // 2
    print(f"  F({m}+{n}) = F({m+n}) = {Fmn}, (F({m})L({n})+L({m})F({n}))/2 = ({Fm}*{Ln}+{Lm}*{Fn})/2 = {pred}, match: {Fmn == pred}")

print(f"""

PROBLEM 5: THE EVEN/ODD SIEVE
================================

The sieve of Eratosthenes removes EVEN numbers first, then multiples of 3, etc.

In the S/A framework:
  EVEN numbers = those with x_n having an even numerator (when n is even).
  ODD numbers = those with x_n having an odd numerator (when n is odd).

But more deeply: a number n is EVEN iff it is in the 2-Vitali atom
of an even core. The S/A decomposition of the sieve:

  Step 1: Remove all n with S/A ratio tanh(ln(n)) having the WRONG parity.
    Even n: tanh(ln(n)) = (n^2-1)/(n^2+1). For even n: n^2-1 is odd, n^2+1 is odd.
    Odd n: same formula, but n^2-1 is even, n^2+1 is even.

  Actually, the S/A ratio doesn't cleanly separate even/odd.
  The right decomposition is:

  Every n has: n = 2^a * m where m is odd (the odd core).
    a = the SYMMETRIC exponent (the 2-adic valuation).
    m = the ANTISYMMETRIC core (the odd part).

  Why "symmetric" for the power of 2?
    Because doubling (the Vitali/Bott map) is the SYMMETRIC operation:
    it doesn't change the odd core.
    The odd core is the ANTISYMMETRIC part: it survives all doublings.

  So: n = (symmetric 2-part) * (antisymmetric odd part).
  Primes > 2: the 2-part is 2^0 = 1 (trivial symmetric part).
    ALL the content is in the antisymmetric odd core.
    Primes are PURELY ANTISYMMETRIC (no symmetric doubling).

  Composites: n = 2^a * m. If a > 0: there IS a symmetric part.
    The number has been DOUBLED. It has acquired symmetric content.

  THE SIEVE REMOVES NUMBERS WITH SYMMETRIC CONTENT.
    Removing evens = removing numbers with a > 0 (symmetric 2-part).
    Removing multiples of 3 = removing numbers divisible by 3.
    Each sieve step removes a SYMMETRIC PATTERN from the torus.

  The primes that survive are the numbers with NO symmetric redundancy.
  They are the purely antisymmetric residuals of the sieve.


PROBLEM 6: EVEN AND ODD FUNCTIONS AS S/A IN FUNCTION SPACE
=============================================================

An EVEN function f(x) = f(-x): symmetric under x -> -x.
An ODD function f(x) = -f(-x): antisymmetric under x -> -x.

Every function: f = f_even + f_odd.

Key even functions: cosh, cos, x^2, x^4, ..., ALL CONSTANTS.
Key odd functions: sinh, sin, x, x^3, ..., tanh, arctan, arctanh.

The EVEN functions CANNOT distinguish x from -x.
They see the MAGNITUDE but not the SIGN.
They are the functions that live on the POSITIVE half-line
(because f(x) = f(-x) means they're determined by |x|).

The ODD functions see ONLY the sign (relative to even content).
They vanish at x = 0 (no sign information at the origin).
They are the functions that encode DIRECTION.

tanh is ODD: tanh(-x) = -tanh(x).
  The coupling REVERSES when you reverse the rapidity.
  Reversing the tournament reverses the coupling.

arctanh is ODD: arctanh(-x) = -arctanh(x).
  The rapidity reverses when you reverse the coupling.

cosh is EVEN: cosh(-x) = cosh(x).
  The undirected content is the SAME regardless of orientation.

THE FORMAL GROUP F(x,y) = (x+y)/(1+xy):
  F is ODD in x (fixing y): F(-x,y) = (y-x)/(1-xy) != F(x,y) generally.
  But F(-x,-y) = (-x-y)/(1+xy) = -F(x,y). So F is ODD in BOTH variables simultaneously.

  F is a JOINTLY ODD function: reversing BOTH inputs reverses the output.
  This means: the composition of two reversed couplings gives
  the reverse of the original composition. Consistent!


PROBLEM 7: USING S/A TO FACTOR THE TRANSFER MATRIX
=====================================================

The k-nacci transfer matrix M_k has eigenvalues that come in types:
  Real eigenvalue tau_k: purely SYMMETRIC (or antisymmetric, depending on sign).
  Complex pair lambda, conj(lambda): S part = 2*Re(lambda^n), A part = 2i*Im(lambda^n).

Tr(M_k^n) = tau_k^n + sum 2*Re(lambda_j^n).  [ONLY S parts contribute to trace]

The DETERMINANT det(M_k^n) = det(M_k)^n = (-1)^(n*(k+1)).
  det encodes the PRODUCT of eigenvalues.
  Tr encodes the SUM.

The S/A decomposition of the CHARACTERISTIC POLYNOMIAL:
  p(x) = x^K - x^(k-1) - ... - x - 1.

  Even-degree terms: x^K, -x^(k-2), -x^(k-4), ...  [S part of the polynomial]
  Odd-degree terms: -x^(k-1), -x^(k-3), ...          [A part of the polynomial]

  For k=3: p(x) = x^3 - x^2 - x - 1.
    S part: -x^2 - 1 (even in x).  Wait, x^3 is odd and -x is odd.
    A part: x^3 - x (odd in x).
    S part: -x^2 - 1 (even in x).
    p(x) = (x^3 - x) + (-x^2 - 1) = A(x) + S(x).

    A(x) = x^3 - x = x(x^2-1) = x(x-1)(x+1).
    S(x) = -x^2 - 1 = -(x^2+1).

    The ANTISYMMETRIC part factors as x(x-1)(x+1).
      Roots: 0, 1, -1. The TRIVIAL roots.
    The SYMMETRIC part factors as -(x^2+1).
      Roots: +/- i. The IMAGINARY roots.

    The actual roots of p(x) = A(x) + S(x) are determined by
    the INTERACTION of these two parts.
    The tribonacci roots are NOT at 0, +/-1, or +/-i.
    They are where A(x) = -S(x), i.e., x(x-1)(x+1) = x^2+1.
    x^3 - x = x^2 + 1. This IS the tribonacci equation.

    So: the tribonacci roots are where the ANTISYMMETRIC part
    EQUALS the SYMMETRIC part (up to sign).
    The eigenvalues live at the S/A BALANCE POINTS of the characteristic polynomial.

THIS IS THE DEEPEST APPLICATION:

  The eigenvalues of a matrix are the points where the
  symmetric and antisymmetric parts of its characteristic polynomial
  are in a specific ratio. For the k-nacci matrix, this ratio is -1
  (they must sum to zero: A(x) + S(x) = 0).

  The ROOTS of any polynomial are the S/A balance points.
  Solving an equation = finding where S and A cancel.


THE UNIVERSAL SOLVER
=====================

To solve ANY equation in the S/A framework:

  1. DECOMPOSE the equation into S and A parts.
  2. The solutions are where S + A = 0, i.e., A = -S.
  3. TRANSLATE to whichever costume makes A = -S easiest to solve.

Costumes available:
  cosh/sinh: good for hyperbolic equations, exponentials.
  cos/sin: good for periodic equations, trigonometry.
  Lucas/Fibonacci: good for integer sequences, recursions.
  even/odd polynomials: good for algebraic equations.
  real/imaginary: good for complex analysis.
  graph/tournament: good for combinatorial problems.

EXAMPLE: Solve x^3 - x^2 - x - 1 = 0 (the tribonacci equation).

  A(x) = x^3 - x.  S(x) = -x^2 - 1.
  Need: A(x) = -S(x), i.e., x^3 - x = x^2 + 1.

  In the cosh/sinh costume: let x = 2*cosh(t/3) (a standard substitution for cubics).
  Then: 8*cosh^3(t/3) - 2*cosh(t/3) = 4*cosh^2(t/3) + 1.
  Using cosh(t) = 4*cosh^3(t/3) - 3*cosh(t/3):
    2*cosh(t) + 4*cosh(t/3) = 4*cosh^2(t/3) + 1.
  This doesn't simplify nicely. The tribonacci equation is not
  solvable by a single cosh substitution.

  But in the FORMAL GROUP costume:
  The equation x^3 - x^2 - x - 1 = 0 factors via (x-1):
    (x-1)(x^2) - (x+1) = 0.
    x^2(x-1) = x+1.
    x^2 = (x+1)/(x-1) = Q(... wait, (1+x)/(x-1) is not quite Q.

  Actually: (x+1)/(x-1) = -Q(-x) = -(1+(-x))/(1-(-x)) = -Q(-x).
  Or: (x+1)/(x-1) = -(1-x)^{-1}(1+x) = ... this is getting complicated.

  The SIMPLEST approach: the augmented equation.
  (x-1)(x^3-x^2-x-1) = x^4 - 2x^3 + 1.
  x^4 - 2x^3 + 1 = 0.
  x^3(x - 2) = -1.
  x^3 = -1/(x-2) = 1/(2-x).

  THIS IS tau_k^k = 1/(2-tau_k) from our framework!
  The tribonacci equation, in the augmented form, IS the identity we proved.

  And 1/(2-x) = the Cayley-transform-like object centered at 2.

  Solved by: tau_3 = the unique positive real root.
  tau_3 = 1.8392867552...


THE SENTENCE:

  Every equation is a request to find where symmetric = antisymmetric
  (where they cancel). The S/A decomposition reduces ALL solving
  to finding balance points. The various costumes (cosh/sinh,
  Fibonacci/Lucas, even/odd, real/imaginary) are different LENSES
  for finding these balance points. Use whichever lens makes
  the balance most visible.
""")

print("opus-2026-03-16-S90bq: SOLVING PROBLEMS VIA S/A DECOMPOSITION")
