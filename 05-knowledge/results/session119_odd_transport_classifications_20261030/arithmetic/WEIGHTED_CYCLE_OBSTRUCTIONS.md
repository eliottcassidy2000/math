# One cyclic cokernel behind odd-cycle signs and odd-integer arithmetic

**PROVED elementary matrix and congruence classification.** This gives a precise
shared mechanism, not an identification of all meanings of “odd,” a new Collatz
cycle formula, or progress claiming Collatz convergence. No literature-priority
claim is made for the Smith-normal-form deductions.

## 1. Inheritance and the invariant that must not be forgotten

The repo's [THM-4604, Collatz words and their carry cocycle](../../../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md)
already proves the affine word formula, the cycle point B/(2^A-3^k), and that
traces retain the two diagonal multipliers but lose the ordered carry. It credits
the cycle formula to earlier work. [Lagarias's primary paper on rational cycles](https://matwbn.icm.edu.pl/ksiazki/aa/aa56/aa5614.pdf)
also treats periodic parity sequences and rational cycles, following Böhm–Sontacchi.
We do not rebrand that formula as new. The extension here is the complete cyclic
integer obstruction, its modular fibers, and a single theorem containing the
odd-cycle sign system as a different specialization.

The [negative-cycle shadow packet](../../collatz_negative_cycle_shadow_20261004.md)
warns that an integral endpoint of an open affine word need not supply its exact
odd intermediate cylinder. Below we prove sufficiency only for a *closed cyclic
solution*, whose common denominator is coprime to two. No open-word implication
is silently imported. The exact ordered carry remains indispensable.

## 2. Uniform weighted-cycle theorem

Fix integers k>=1, p>=1, u with gcd(p,u)=1, and a_i>=1 for i modulo k.
Let A=sum_i a_i and D=p^A-u^k. Define the integral cyclic matrix M by

    (Mx)_i = p^(a_i) x_(i+1) - u x_i.

For k=1 both contributions occupy the same entry, so M=[p^(a_0)-u].
For an arbitrary integral right-hand side b=(b_i), define

    N_i(b) = sum_(j=0..k-1) u^(k-1-j) p^(a_i+...+a_(i+j-1)) b_(i+j),

where the empty prefix has exponent zero and indices are cyclic.

**Theorem.**

1. det M = (-1)^(k-1) D, and its Smith normal form is

       diag(1,...,1,|D|),

   with the last entry zero if D=0. Thus the cokernel is cyclic: Z/DZ, with
   the convention Z/0Z=Z.
2. If D!=0, the unique rational solution is x_i=N_i(b)/D. An integral solution
   exists exactly when N_0(b)=0 modulo D. One congruence suffices for all k
   coordinates.
3. D=0 occurs only when p=1 and either u=1, or u=-1 with k even. Then an
   integral solution exists exactly when N_0(b)=0 as an integer, and every
   integer choice of x_0 extends uniquely to an integral cyclic solution.
4. For every modulus m>=1, Mx=b modulo m is solvable exactly when

       N_0(b)=0 modulo gcd(|D|,m).

   When solvable there are exactly gcd(|D|,m) solutions. For D=0 this says
   compatibility modulo m and exactly m solutions.

**Proof of the Smith statement.** For k>=2, the determinant has just the diagonal
and cyclic-shift terms. One (k-1)-minor, deleting the last row and column, is
(-u)^(k-1). Another, taking rows 0,...,k-2 and columns 1,...,k-1, is
p^(a_0+...+a_(k-2)). They are coprime. Therefore the gcd of all (k-1)-minors is
one. By the determinantal-divisor characterization of Smith form, all first
k-1 invariant factors are one and the last is |det M|. For k=1 the statement is
immediate. This also proves rank k-1 when D=0, including the zero 1-by-1 matrix.

Composing the k affine recurrences proves D x_i=N_i(b). Directly, their
numerators obey the useful identity

    p^(a_i) N_(i+1)(b) - u N_i(b) = D b_i.       (1)

When D!=0, gcd(p,D)=gcd(u,D)=1. Hence divisibility of N_0 by D propagates to
all N_i by (1), and the displayed fractions solve every equation. Conversely
any integral solution must satisfy that divisibility.

If D=0, coprimality rules out p>=2. With p=1, u^k=1 gives exactly the stated
cases. The recurrence x_(i+1)=u x_i+b_i closes precisely when N_0(b)=0; all
coefficients are integers and x_0 is free. This proves the singular statement.

The homomorphism b -> N_0(b) modulo D is surjective: its first coefficient
u^(k-1) is a unit modulo D; in the singular case that coefficient is +/-1.
The preceding solvability results identify its kernel with M Z^k. Thus it is
an explicit cokernel coordinate, not just an abstract invariant factor. Reducing
modulo m means quotienting this cyclic cokernel by multiplication by m, giving
exactly the stated compatibility. Smith form gives kernel size gcd(|D|,m), and
every nonempty fiber is a coset of that kernel. This proves all four statements.

**Coprimality matters.** If it is dropped, the scalar determinant need not carry
the full obstruction. For p=u=2 and a=(1,1), M=[[-2,2],[2,-2]] has Smith form
(2,0), not (1,0). The extra factor two is a genuine missing invariant.

## 3. Odd cycles are the unit-weight specialization

Put p=1,u=-1. Then the same equations become

    x_i+x_(i+1)=b_i,     D=1-(-1)^k.

For odd k, D=2. There is a unique rational solution; it is integral exactly
when sum_i b_i is even, since alternating signs disappear modulo two.
The cokernel is one Z/2 obstruction. For even k, D=0: the alternating sum of
b_i must vanish, and one integer parameter remains free.

For b_i=1, this recovers the cycle-bipartition obstruction: odd k forces the
unique rational assignment x_i=1/2, which is not integral or binary; even k
allows alternating x_0,1-x_0 and exactly two binary assignments. For b_i=0,
an alternating-sign function on an odd cycle must vanish over the reals or
integers. Over Z/mZ it can retain precisely the two-torsion permitted by the
modular kernel count. Odd cycle length is controlling a sign monodromy here,
not the parity of every vertex label.

## 4. Odd Collatz states are the coprime-weight specialization

Put p=2,u=3 and every b_i=1. For a prescribed positive exponent word a,

    D=2^A-3^k,    x_i=N_i/D.

D is nonzero and odd. Each N_i is positive and odd: its j=0 summand is odd and
all later terms are even. Its last summand also shows N_i is coprime to three;
D is coprime to three as well.

**Exact cycle criterion.** The word describes a closed cycle of nonzero odd
integers under x -> (3x+1)/2^v_2(3x+1) if and only if D divides N_0. Indeed the
single congruence makes every coordinate integral by the theorem, and division
of odd N_i by odd D makes every x_i odd. Therefore

    3x_i+1 = 2^(a_i) x_(i+1)

has *exactly* valuation a_i, since x_(i+1) is odd. This verifies every
intermediate exponent, not just the final affine endpoint.

The cycle is positive exactly when D>0, because all N_i are positive. If D<0,
all its entries are negative. Thus integrality, positivity and the exact word
are separately accounted for. There is no D=0 exceptional Collatz word.

The least period of an integral cycle equals the least period of its exponent
word. A repeated cycle plainly repeats its exact valuations. Conversely, a
periodic exponent word admits a corresponding rotation of its unique rational
solution, so uniqueness forces the same coordinate period. Words considered up
to cyclic rotation therefore describe cycles without a chosen starting point;
repeated words must still be removed when counting primitive cycles.

The singleton word (2) gives x=1; (1) gives x=-1; (1,2) gives the two-cycle
-5 -> -7 -> -5. The last has even cycle length although both its states are
odd integers. This is a direct warning against conflating those two senses of
“odd.”

## 5. Same determinant and exponent multiset can hide different arithmetic

Take the two words

    w=(1,1,1,2,1,1,4),    v=(1,1,1,1,1,2,4).

They have the same length seven, total exponent eleven, and exponent multiset.
Both have D=-139, so their matrices have the same Smith invariant factors and
same number of solutions over every modulus for a *solvable* right-hand side.
Nevertheless, for b=(1,...,1),

    N_0(w)=2363=-17 D,    N_0(v)=2123, with 139 not dividing 2123.

The first is the actual negative odd cycle

    -17 -> -25 -> -37 -> -55 -> -41 -> -61 -> -91 -> -17.

The second has no integral cyclic solution. Thus the abstract cyclic cokernel,
the determinant, total exponent, and even the exponent multiset do not specify
where the chosen right-hand side lands in that cokernel. The ordered numerator
is the necessary additional coordinate. This sharpens the inherited trace-
blindness warning into an explicit solvability counterexample.

## 6. What is and is not unified

The source object is a weighted cyclic system; the map to Smith form preserves
its abstract integer obstruction and modular kernel sizes. It forgets the image
of b. The explicit N_0 map restores that information. Odd-cycle sign alternation
and odd-integer dynamics fit the same theorem through different weights, but
their obstructions differ: Z/2 in one case, an odd cyclic group in the other.
This explains a connection without treating them as equivalent problems.

A prospective Collatz advance must still exclude the positive integral solutions
with nontrivial primitive exponent words, or prove an independent descent/coverage
theorem. This classification does neither. It instead gives a complete, cheap,
exact obstruction test for each specified word and shows precisely why a scalar
or unordered invariant cannot replace that test.

## 7. Reproduction

Run `python3 check_weighted_cycles.py` and `python3 -O check_weighted_cycles.py`.
Both outputs match `receipt.json`. The standard-library checker uses explicit
exceptions and independent exact determinant elimination, direct modular image
enumeration, and literal accelerated-cycle equations. Its bounded universes are:
2,418 weighted matrices; 7,770 numerator identities; 585 modular systems with
38,250 right-hand sides; and 21,844 Collatz exponent words of lengths at most
seven with each exponent at most four. The latter finite scan is not a theorem
excluding positive cycles beyond that box. The all-parameter results follow from
the proofs above, not extrapolation of the checks.
