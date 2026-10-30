# Odd character symbols on every odd modulus: complete support classification

**PROVED elementary classification with exact finite controls.** No literature
priority, new Ramsey bound or major-conjecture resolution is asserted. The
prime-power valuation-frame classification is inherited from
[Odd-order arithmetic, signed fibres, and odd-cover descent, §3](../../odd_prime_function_cycle_20261018d.md).
[LEM-004](../../../../01-canon/theorems/LEM-004-tournaments-are-odd-functions.md)
provides the odd-sign/circulant-tournament dictionary. The extension here retains
all composite-modulus strata and classifies their allowed support, rather than
silently applying a prime-field character to every nonzero residue.

## 1. Definitions and classification theorem

Let N>1 be odd and U(N)=(Z/NZ)^×. A sign character is a homomorphism
χ:U(N)→{±1}. A partial odd symbol is a function σ:Z/NZ→{0,±1} such that

    σ(0)=0,  σ(-d)=-σ(d),  σ(ud)=χ(u)σ(d) for every u in U(N).

Values and negation are ordinary integers, not elements of characteristic two.
For a divisor m>1 of N, let S_m consist of the residues of additive order m:

    S_m={ (N/m)t mod N : t in U(m) }.

Write N=∏p^{a_p}. Every sign character has a unique expression

    χ(u)=∏_{p in E} (u/p),

where E is a subset of the distinct prime divisors, and (u/p) is the Legendre
sign. Put f=∏_{p in E}p; the empty product is one. This f is the least divisor
of N through which χ factors, its conductor for the current finite unit action.
The character is odd exactly when the number of primes p≡3 mod4 in E is odd.

**Theorem 1.** If χ(-1)=+1, the only partial odd symbol is zero. If χ(-1)=-1,
then the following is a complete and unique parametrization:

- On every stratum with f not dividing m, σ is forced to vanish.
- On every stratum with f dividing m, choose independently ε_m in {0,±1}, and set

    σ((N/m)t)=ε_m ∏_{p in E}(t/p),  t in U(m).

There are exactly 3^D partial odd symbols, where

    D=#{m|N:f|m}=∏_{p in E}a_p ∏_{p not in E}(a_p+1).

Here f>1 because χ is odd, so no m=1 term has been counted. Requiring nonzero
values on exactly a prescribed set of strata A gives 2^|A| symbols when every
m in A is divisible by f; otherwise it gives none.

## 2. Proof, including stabilizers and conductor

Each sign character kills squares. On U(p^a), the kernel of reduction to U(p)
has odd order p^{a-1}, so every homomorphism to {±1} kills it. The square map on
U(p) has kernel {±1}, hence its image has index two; the only nontrivial sign
character is the quotient by squares, the Legendre sign. CRT therefore gives
exactly the displayed products, with independently chosen prime factors.
Reduction U(N)→U(m) is surjective for m|N: lift a unit at each prime dividing m,
and choose a unit freely at omitted prime factors. A local Legendre factor at p
survives this reduction exactly when p divides m. This proves the conductor claim
and χ(-1)=∏_{p in E}(-1)^((p-1)/2).

The unit action on S_m is transitive. The stabilizer of the representative N/m
is precisely

    H_m={u in U(N):u≡1 mod m}.

A nonzero value at that representative extends equivariantly if and only if χ is
trivial on H_m. Otherwise a stabilizer element with character -1 forces that
value to equal its negative, hence to be zero. Triviality on H_m is equivalent
to factoring through U(m), hence to f|m. In the compatible case every point is
u(N/m), and the rule χ(u)ε_m is well-defined, with one independent frame value.
This proves every assertion, including necessity; it is not merely a construction.
If χ is even, covariance at u=-1 and oddness force every value to be zero.

## 3. Exact support and tournament consequences

Let z count nonzero *residue differences*. For fixed odd χ the complete support
polynomial is

    Σ_σ z^{|supp σ|}=∏_{m|N, f|m}(1+2 z^{φ(m)}).

Every supported stratum supplies one of its two signs and has φ(m) residues;
choosing zero contributes one. Thus the formula counts the whole family, not
just its total size.

Use σ to orient the pair x,y when σ(y-x)≠0, with x→y exactly when σ(y-x)=+1.
Oddness makes this a well-defined partial orientation of the complete graph.
Its degree is uniform: outdegree and indegree both equal |supp σ|/2. The number
of oriented unordered pairs is N|supp σ|/2, so the corresponding edge-count
polynomial replaces each exponent φ(m) by Nφ(m)/2. These are integers since
m>1 is odd and φ(m) is even.

**Corollary 2: full support classification.** A nonzero sign on every nonzero
residue with a *single global unit sign character* exists exactly when

    N=p^a,  p≡3 mod4,  χ=(·/p),

and then there are 2^a labelled symbols. To prove necessity, full support includes
S_p for every prime p|N, so f divides every such p. Oddness gives f>1. Therefore
N has only one distinct prime divisor, f=p, and p≡3 mod4. Sufficiency and the count
are the inherited valuation-frame formula. This does not say composite odd N
admits no circulant tournament: arbitrary odd symbols always do. It excludes only
the additional single-character covariance constraint.

For any nonzero partial symbol, its stabilizer within the unit multipliers is
exactly ker χ. Indeed covariance says a multiplier either preserves σ or changes
it to -σ, and those functions differ when σ≠0. Translations preserve the oriented
difference relation as well. Thus (Z/NZ)⋊ker χ is an explicit orientation-preserving
affine subgroup of order Nφ(N)/2. The other affine multipliers reverse all arcs.
No claim about the full automorphism group, which may be larger, is made.

Distinct characters cannot share a nonzero symbol: pick a nonzero σ(d); then
χ(u)=σ(ud)/σ(d) recovers every character value. Consequently, if one wishes to
count partial odd symbols covariant under *some* sign character, their exact total is

    1+Σ_{χ odd}(3^{D_χ}-1).

The zero symbol is counted only once. The theorem itself concerns a fixed χ.

## 4. Minimal composite hostile

Take N=15 and χ(u)=(u/3), an odd character of conductor three. The allowed strata
have additive orders 3 and 15; the additive-order-five stratum must be zero.
For a literal obstruction, multiplier 11 fixes residue 3 modulo15, but χ(11)=-1.
Thus covariance would impose σ(3)=-σ(3), impossible for a sign. A naive
nonzero valuation-frame extension from prime powers fails at this fixed point.
There are 3^2=9 partial symbols for this χ, but no full tournament symbol with it.

This is an actual map from arithmetic character support to a circulant relation.
It does not imply a Ramsey extremal value, a Hamiltonian-path maximization theorem,
or an unrestricted Boolean circuit lower bound.

## 5. Reproduction

`check_characters.py` enumerates every sign character on each odd N=3,...,151 by
assigning signs to the literal square cosets and testing their multiplication
table. It independently compares signed constraint components on all nonzero
residues against the divisor/conductor prediction, checks support polynomials,
and verifies exact unit stabilizers on one nonzero symbol per eligible stratum.
For N<=15 it additionally enumerates all literal odd {0,±1} functions and checks
all covariance equations. The all-N proof is the stabilizer/CRT argument above.
Normal and optimized runs agree with `receipt.json`; explicit exceptions survive -O.
