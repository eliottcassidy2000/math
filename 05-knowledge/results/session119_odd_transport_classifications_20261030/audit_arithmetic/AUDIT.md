# Independent audit of the weighted-cycle arithmetic bridge

**Verdict: PASS, with the stated coprimality and closed-cycle hypotheses.**
This is an independent mathematical audit and separate finite implementation.
It does not certify a Collatz convergence proof, exclude all positive cycles,
or establish literature priority. No producer implementation is imported.

## 1. Determinant and complete integer obstruction

For k>=2 the weighted cyclic matrix has only two nonzero permutation terms:
the diagonal and the full cyclic shift. Its determinant is therefore
(-1)^(k-1)(p^A-u^k). At k=1 the two contributions share a single matrix entry;
the stated matrix is [p^(a_0)-u], and the formula still holds.

The two displayed (k-1)-minors really are (-u)^(k-1) and a power of p.
Coprimality makes their gcd one. Hence the gcd of *all* (k-1)-minors is one,
forcing every earlier Smith invariant to be one. This proves the cyclic
cokernel, including the rank k-1 case when the determinant vanishes. The
singular cases are precisely p=1,u=1 and p=1,u=-1 with even k. Indeed equality
p^A=u^k and gcd(p,u)=1 force p=1, and then u^k=1. The allowed u=0,p=1 case
is nonsingular and causes no exceptional failure.

The ordered numerator formula follows by composing the affine recurrences;
its direct relation p^(a_i)N_(i+1)-uN_i=D b_i verifies all equations at once.
For nonzero D, p and u are both coprime to D. Thus N_0 divisibility propagates
to every coordinate. The first coefficient u^(k-1) gives a surjective map
onto Z/DZ, with image-kernel equal to M Z^k. At D=0 it is +/-1 and gives
an integer cokernel coordinate. This checks that the specified right-hand
side, not merely the abstract Smith form, is retained.

Modulo m, the obstruction is the image of N_0 in (Z/DZ)/m(Z/DZ), namely
modulo gcd(|D|,m). Every nonempty solution fiber has that same gcd as its size.
This remains true at D=0, giving compatibility modulo m and m solutions;
at m=1 it correctly gives one element in the trivial domain and codomain.

## 2. Closed Collatz words and the ordered hostile

For p=2,u=3,b_i=1, D is odd and nonzero. Each N_i is positive and odd, so an
integral cyclic solution has nonzero odd coordinates and common sign equal to
the sign of D. Each equation then has exactly the prescribed 2-adic valuation,
because its next coordinate is odd. This uses *closed* cyclic divisibility;
it does not turn a merely integral endpoint of an arbitrary open word into
an exact intermediate cylinder.

The least-period assertion also passes: repeated states repeat their exact
valuation sequence. Conversely a rotation preserving the exponent word gives
another solution to the same nonsingular cyclic system, so uniqueness forces
the corresponding coordinate period. Primitive cycles therefore require
primitive words after quotienting by cyclic rotation.

For the two words

    (1,1,1,2,1,1,4), (1,1,1,1,1,2,4),

independent calculation gives D=-139 for both, and N_0=2363 and2123. Only the
first is divisible by D and yields

    -17,-25,-37,-55,-41,-61,-91.

Every accelerated step and its valuation was checked directly. Thus identical
Smith form, length, total exponent and exponent multiset do not determine the
right-hand side's cokernel class. Positivity is separately necessary for the
positive-integer Collatz conjecture; this negative cycle is a hostile to a
proposed information loss, not progress excluding positive cycles.

## 3. Signed-cycle and signed-graph compatibility

The p=1,u=-1 specialization gives x_i+x_(i+1)=b_i, with one factor2 for odd
cycle length and one free coordinate with an alternating RHS condition for
even length. Both are correctly distinguished from parity of integer states.

The separate signed-graph proof also passes. Spanning-tree elimination uses
unit pivots, leaving root coefficients 0 or +/-2 on every chord. A balanced
component has only zero rows and therefore retains all corresponding RHS
compatibility equations. An unbalanced component has one factor2; its other
nonzero rows must have equal RHS values after sign normalization. Hence
n-c unit factors, o factors2, b free kernel directions and e-n+b free cokernel
directions are correct, including isolated vertices and signed loops.

In particular odd target modulus removes the factor-two torsion but does not
remove the free compatibility conditions. The even square with RHS(1,0,0,0)
stays insoluble modulo3. The triangle with RHS(1,0,0) has one solution modulo
odd moduli and none modulo even moduli. These were independently enumerated.
No modular Gram-matrix substitution is used.

## 4. Independent finite universe

`check_independent.py` constructs every matrix literally, including the k=1
collision of diagonal and shift entries. It computes determinants by permutation
expansion and determinantal divisors by gcds of all literal minors, rather than
using either producer's formula as an implementation. It also enumerates actual
modular image fibers and separately evaluates the numerator criterion.

The weighted-matrix universe is k=1,...,4, exponents in{1,2}, p=1,...,4 and
u=-3,...,3 with gcd(p,u)=1. Modular fibers use k<=3,p<=3 and m=1,...,5.
Closed Collatz words use lengths1,...,6 with exponents in{1,2,3}; the two
seven-letter hostile words are added explicitly. Signed RHS controls run
modulo2,...,7. Exact counts are in `normal.json`.

Both normal and optimized Python runs were executed and are byte-identical.
All checks use explicit exceptions. Finite checks corroborate the independent
all-parameter arguments above, not an exclusion of all untested trajectories.
