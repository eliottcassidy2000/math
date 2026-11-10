# Independent audit: exact two-burst arithmetic image

**Verdict: PASS, with one corrected citation request.** Audit date 2026-11-10.
The all-order proofs of the shifted product, injectivity, inverse, dyadic
count, density zero, prime-order criterion, valuation formulas and gcd
identity in `../coverage/PROOF.md` are sound in their explicitly restricted
domains. The fixed-seed divisor criterion is sound for physical positive
unit-start seeds. No assertion of general Chan–Pak coverage follows.

The independent verifier imports no producer code and does not use the
inherited O/Z recurrence. Run from the repository root:

    python3 05-knowledge/results/session241_prime_relations_20261110/coverage_audit/verify.py

It writes `verification.json`, including its source hash. Checks raise explicitly on failure and remain active in optimized Python.
Normal and optimized runs produce identical output.

## Inheritance and scope audit

Read the canonical entries THM-4888,
`01-canon/theorems/THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement.md`,
and THM-4916,
`01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md`.
The O/Z moves, rooted-edge state, general exact source inverse and product
mechanisms are inherited. The present result is a particularly simple
restriction with a binary parser and explicit thin-image count, not the
invention of source inversion or of tournament multiplication. No claim of
external priority was audited or made. Historical directories not available
in the checkout were not silently presumed read.

The concept board used here is physical root conventions, exact arithmetic
inversion, coverage versus local congruences, and multiplicity versus support.
The anchor is exact inverse certification; the niche is a recurrence-free
Hamiltonian count; the wildcard is the exact dyadic boundary term.
Canonical hostile: a modular or shifted-factor relation does not certify
coverage of target integers. Least-used sidecar: actual root-edge labels and
binary runs, rather than residue data alone.

One citation defect was reported to the producer: the unrelated colliding-ID
file `THM-1370-elliptic-jacobian-conjecture-all-dimensions.md` should not be
listed as omission provenance. The relevant distinct slug is
`THM-1370-h-spectrum-omits-7-21-all-n.md`. This does not affect the new proof.

## Independent physical derivation

Let the physical root-edge word be O^A Z O^B Z with A=a+1>=2, B=b>=1.
There are m=A+B+2 spine vertices, ordered so earlier vertices beat later
vertices. The root beats O vertices and loses to Z vertices. Thus the
complete tournament has a+b+4 vertices. The extra initial O represents the
rooted-edge seed; it is not optional.

Every Hamiltonian path is an increasing subsequence of spine vertices,
then the root, then the complementary increasing subsequence. If the root
is first or last, there is exactly one path in each case. For an internal
root, the last left vertex must be Z and the first right vertex must be O.
The two possible last-left indices are k=A+1 and m.

For last-left k, summing over first-right j=1,...,A gives
sum_j 2^(k-j-1)=2^A-1. The additional case j=k+1 gives one, hence 2^A paths.
For last-left m, first-right may be any O position, and the count is

    sum_(j=1,...,A) 2^(m-j-1)
      + sum_(j=k+1,...,m-1) 2^(m-j-1)
    = 2^(A+B+1)-2^B-1.

Adding the two endpoint-root paths proves

    H = 2^(A+B+1)+2^A-2^B+1.

Substitution r=A+1=a+2 and s=B+1=b+1 yields
2H-3=(2^r-1)(2^s+1), independently of any inherited state calculus.
The tournament is strong since the root reaches the first spine vertex,
every spine vertex reaches the final Z, and that Z reaches the root.

**Hostile convention test:** the full physical word OZOZ has H=9, whereas
the construction with a=b=1 has full word OOZOZ and H=19. Any shorthand
O^a Z O^b Z must explicitly mean moves after the rooted-edge seed.

## Exact inverse and counting audit

Write M=(2^r-1)(2^s+1), r>=3, s>=2. Its binary patterns are:

- r<s: 1^r 0^(s-r) 1^r;
- r=s: 1^(2r);
- r>s: 1 0^s 1^(r-s-1) 0 1^s.

These are disjoint types, including the r=s+1 empty middle-one run.
They recover both ordered exponents. A purported parser candidate must still
satisfy the exponent bounds and exact product; those checks reject malformed
strings. This is an all-order uniqueness proof, not evidence extrapolated
from the finite inverse test.

For H<=2^k, every pair with r+s<=k qualifies, no pair with r+s>=k+2
qualifies, and the boundary r+s=k+1 qualifies exactly when r<s.
Consequently the count is

    (k-4)(k-3)/2 + floor(k/2)-2, for k>=4.

The boundary term is necessary: at k=10 the triangle alone gives21, while
the actual answer is24. The reviewer initially omitted this boundary in a
scratch test; that failed immediately and was corrected before reporting.
The leading coefficient 1/2 for (log_2 X)^2 follows by sandwiching arbitrary
X between adjacent dyadic bounds. Density is zero among all integers and
among odd integers. This is about distinct values, justified by injectivity,
not a word-count/multiplicity substitution.

## Prime and fixed-seed sidecars

For a fixed physical seed, the two nonempty O bursts erase E before each
Z. Direct expansion gives

    H=D+2^a U+2^b(2^(a+1)U-d),
    2(H-D)-d=(2^(a+1)U-d)(2^(b+1)+1).

Since physical unit-start seeds have U>=2d>0, the shifted target on this
family is positive. Enumerating positive factors B of the prescribed
2^(b+1)+1 shape and requiring ((2(H-D)-d)/B+d)/U=2^(a+1) with a>=1
is necessary and sufficient. It need not be injective for arbitrary seeds.

The prime-order criterion is valid because -1 is the unique element of
order2 in a prime field. The valuation identities are ordinary odd-prime
LTE, retaining the base valuation; they do not assume non-Wieferich primes.
For the gcd proof, work modulo the entire common divisor, not merely its
prime radical. Put x=2^g, R=r/g, S=s/g. If R is odd, x^RS is both1 and-1,
so the odd common divisor is1. If R is even, S is odd and a Bezout relation
uR+vS=1 has v odd; x=(-1)^v=-1 modulo the entire common divisor. This proves
it divides x+1, including multiplicities. Reverse divisibility is immediate.
These statements preserve the factorization of 2H-3, not primality of H:
H=19 and H=33 give shifted products35 and63, respectively.

## Finite-exact verification universe

- Literal full adjacency subset DP for every 1<=a,b<=4 (16 tournaments,
  orders6 through12), with converse invariance as a separate control.
- All vertex permutations independently counted whenever order<=8.
- Transitive and directed-triangle controls have H=1 and H=3.
- Missing-seed hostile has H=9 instead of19.
- Binary patterns checked on 14,400 ordered exponent pairs.
- Gcd identity checked on 10,000 ordered exponent pairs.
- Independently enumerate inverse candidates by exact divisibility by
  2^r-1, then test whether the quotient minus1 is a power of2. This is not
  the producer's run parser. Every integer target from1 through100000 is
  checked, including even targets; its image has91 elements.
- Exact counts verified at k=10,20,40,80,160,320, giving respectively
  24,144,684,2964,12324,50244.

All tests passed. They are finite regression evidence, separate from the
all-order arguments above. They do not audit every possible seed graph,
every O/Z word, or every tournament.

## Actual consequence

There is an explicit infinite strong-tournament family carrying a uniquely
recoverable ordered pair of exponents and a Mersenne/generalized-Fermat
product in its shifted Hamiltonian count. The exact image is sparse, with
only quadratic growth in bit length. This is a rigorous arithmetic bridge
and a quantified limitation of this particular language. It neither proves
nor disproves Chan–Pak, excludes no new value for all tournaments, and does
not turn shifted-product compositeness into compositeness of H. The universal
omission of7 and21, if used, remains an inherited theorem, not a conclusion
of these image tests.
