# A free divisor-sum bit depends on the allowed operations

**PROVED elementary results and exact finite controls.** The general group
facts are standard finite linear algebra. The contribution is their explicit
arithmetic realization and a precise separation of observation symmetries,
named-update symmetries, and source-reconstructing observations. No literature
priority or arithmetic-convergence claim.

Write sigma(n) for the sum of all positive divisors and s(n)=sigma(n)−n for
the aliquot sum. They are different observers. Let [x] denote the positive
rational squareclass of x. For positive integer representatives it can be
recorded by the squarefree kernel.

## 1. Exactly what 3, 5 and19 contribute

For squarefree products of the three named primes, use bits
x=(e_3,e_5,e_19) in F2^3. Multiplicativity gives

    [sigma(3^e_3 5^e_5 19^e_19)] = [6]^e_5 [5]^e_19,

because sigma(3)=4, sigma(5)=6 and sigma(19)=20. The two nonzero classes
[6] and [5] are independent. Thus the sigma-squareclass observer has four
fibers, each containing two integers differing by the prime3 toggle.
The prime3 bit is invisible to this observer, not absent from the integer.

The eight integers are1,3,5,15,19,57,95,285. By contrast all three individual
prime inputs have aliquot sum1. Their equal aliquot values do not imply equal
responses to later attachment.

## 2. Full observational symmetry is much larger than operational symmetry

Let M:F2^m→F2^r be surjective, let K=ker M and k=m−r. Each of the2^r fibers
has2^k elements. Consequently all bijections F preserving the exact observation
MF(x)=Mx form

    product_(y in F2^r) Sym(2^k).                      (1)

This permits a different arbitrary rearrangement in every fiber. It is not
necessarily linear, affine, or compatible with an update rule.

Now retain every **named** bit toggle T_i(x)=x+e_i and demand
F T_i=T_i F for all i. Then

    F(x)=F(0)+x,

by expressing x as a sum of coordinate vectors and commuting those toggles.
Therefore the exact observation-and-update symmetry group is

    {x↦x+c : c in K},  of order2^k.                   (2)

Conversely each such translation preserves M and commutes with all toggles.
This is a complete classification for the declared operation signature.
Toggles are guarded prime attachment/removal on a squarefree integer bank,
not unrestricted multiplication which may create repeated prime factors.

For the three-prime bank above, (1) has order16, while (2) has order2:
identity and the uniform3-toggle. For example

    F(e_3,e_5,e_19)=(e_3+e_5 e_19,e_5,e_19)

preserves every sigma-squareclass but fails to commute with the named5 and19
toggles. It is a legitimate observational symmetry and an invalid symmetry
of the stronger named-update system.

The intermediate affine symmetry count is also exact. Choose a splitting
F2^m=K⊕S. Every invertible affine F preserving M has the form

    F(k,s)=(A k+B s+c,s),
    A in GL(K), B in Hom(S,K), c in K.

Its order is |GL(k,2)|2^(k(r+1)). For m=3,r=2 it has order8, strictly between
16 and2. This classification follows by comparing the two components of MF=M
and imposing invertibility on K. For the zero-augmented Hamming bank, m=8,r=3,
there are32 uniform kernel toggles; the general affine group is much larger.

**Do not conflate port permutation with commuting updates.** A coordinate
permutation sends T_i to T_(pi(i)); it generally normalizes the update family
rather than commuting with each named member. The Fano GL3(2) group also changes
the syndrome basis equivariantly. It is not the same group as (2), which fixes
the observation value and every update name.

## 3. Neutral prime-power states: six objects, several symmetry groups

The arithmetic lane supplies the cited classification that individually
square-sigma prime powers are3,81,343. For fixed u with gcd(u,21)=1, consider

    {3^a 7^b u : a in{0,1,4}, b in{0,3}}.

All six states have the same global sigma-squareclass. If that is the only
observation, all720 permutations preserve it. If we also retain the two named
prime-port partitions (same3-exponent, or same7-exponent), the automorphism
group is exactly S3×C2, order12.

Proof: preserving equality of each coordinate induces a permutation of each
coordinate's equivalence classes, hence a separate permutation of the three
3-states and two7-states. Every such product permutation works. This concerns
allowed exponent-state replacements, not commuting with literal multiplication
by3, nor preserving ordinary magnitude or an aliquot trajectory. It is a
natural port-respecting subgroup, not the whole observational symmetry group.

Global squareclass equality also does not preserve arbitrary character tests
with zeros: sigma(3)=4 and sigma(81)=121 are both squares, but modulo11 the
former has Legendre symbol+1 and the latter0. Quadratic-character comparisons
must retain their unit guard, or add support/valuation data. The squareclass
statement itself remains correct.

## 4. The same3 attachment can reconstruct an integer for a different observer

For coprime positive integers a,n, multiplicativity gives

    s(an)=s(a)n+sigma(a)s(n).                            (3)

This follows by expanding sigma(a)sigma(n)−an with sigma(n)=n+s(n).
If a>=2 then s(a)>0, so

    n=[s(an)−sigma(a)s(n)]/s(a).                        (4)

In particular, on the domain3∤n,

    s(3n)=n+4s(n),        n=s(3n)−4s(n).                (5)

Thus any representation that answers both exact aliquot queries n↦s(n) and
n↦s(3n) must be injective on that domain. It cannot collapse all primes other
than3 to one state merely because each has aliquot sum1. For the requested
primes5 and19, the same attachment gives s(15)=9 and s(57)=23, recovering
5 and19 respectively. This is an exact source-reconstruction theorem, not a
bound on aliquot trajectories.

Meanwhile multiplication by3 preserves sigma-squareclass whenever3∤n,
because sigma(3n)=4sigma(n). The operation is free for that coarse observer
and source-revealing for the paired exact aliquot observer. There is no
contradiction: the observers retain different data.

The coprimality guard is necessary. At n=3, sigma(9)=13 rather than16, and
the right side of (5) would be7 whereas s(9)=4. No update across a shared
prime port can use the coprime formula without its correction factor.

For fixed coprime a, the value pair (n,s(n)) transforms by

    [[a,0],[s(a),sigma(a)]].

Changing coordinates to (n,sigma(n)) diagonalizes this matrix to
 diag(a,sigma(a)). This is the exact addition/multiplication balance, not a
proof that a coarse sum alone closes under future operations. Even exact total
sum probes need not reconstruct n: n=14 and23 both have sigma(n)=24 and
sigma(3n)=96. Their paired aliquot values do distinguish them.

## 5. Rank alone fails a canonical tournament operation

The palette quotient in the tournament lane retains binary rowspaces, not
just their dimensions. Two tiny hostile controls show why. Use the intrinsic
source-first order of a transitive tournament, whose automorphism group is
trivial, and fix any nonzero color v.

* In the linear palette quotient, the two-color lists (v,0) and (v,v) both
  have rank one. Delete their universal source: the remaining colors are
  zero and nonzero, respectively, so the resulting ranks are zero and one.
* In the affine palette quotient, the three-color lists (0,v,v) and (0,0,v)
  both have augmented rowspace dimension two. Delete their universal source:
  the remaining pair is equal in the first case and different in the second,
  giving augmented dimensions one and two.

Thus even a canonical, unmarked source-deletion operation can distinguish
states merged by rank. This is not a lower bound for arbitrary counting
algorithms: it is an exact failure of this proposed operation-preserving
state quotient. The verifier checks the coordinate projections directly.

## 6. Verification and inheritance

The coprime divisor-sum identity and squareclass linearization are elementary
inherited mechanisms, not new facts about sigma. Session107 supplies the code
banks and the warning to distinguish future operations from static decisions.
The present explicit group hierarchy and reconstruction tests declare their
carriers and operations so that a larger symmetry claim can be checked.

`verify_operations.py` enumerates all8! permutations of the three-prime cube,
separates observation-only, affine and named-update-preserving groups, checks
the noncommuting conditional-toggle hostile, enumerates all6! permutations
of the neutral grid, and verifies coprime probe identities over stated finite
ranges. It uses exact arithmetic and explicit exceptions. The general proofs
above, not those finite samples, establish the classifications.

Run `python3 verify_operations.py` and `python3 -O verify_operations.py`.
No third-party dependency is required.
