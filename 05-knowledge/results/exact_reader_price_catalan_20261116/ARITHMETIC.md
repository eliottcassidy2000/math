# Additive versus multiplicative target selection

## A wheel equivalence, not a primality equivalence

For an odd modulus W and gcd(M,W)=1, the map
x=−2p/M modulo W is an affine permutation and x+2=2(M−p)/M.
Thus the local unit conditions for Goldbach (p,M−p) are exactly the
local unit conditions for a twin pair (x,x+2), with parity handled
separately. If M is divisible by an odd prime ℓ, the two forbidden
residues coalesce, so the admissible count modulo ℓ changes from ℓ−2
to ℓ−1. By CRT these local counts multiply for a squarefree wheel.
The affine map is target dependent and does not send actual primes
to primes.

For M=8, the wheel 3·5·7 has 1·3·5=15 formal channels, though the
small Goldbach pair (3,5) uses exceptional wheel primes. For M=9,
parity forces (2,7) if both summands are prime. These decompositions
partition the first four primes into (3,5) and (2,7), two fundamentally
different roles rather than two symmetry-related witness pairs.

Multiplicatively, 8+1=9 yields 2⁶−1=(8−1)(8+1)=7·9, exactly the
Bang–Zsigmondy exception at exponent six; the Catalan link is already
in `math/05-knowledge/results/collatz_mod6_20260917_synthesis.md`,
lines 85 and 397. The two absent tournament totals are the first two
values Φ₃(2^k) for k≥1: 7 and 21. Their next value 73 is realized by
an apex over a seven-vertex ordered spine: exactly vertices 2,6,7 beat
the apex; it beats all other spine vertices. This can be checked by the subset DP in verify.py (zero-based mask 98). The product (2^k−1)Φ₃(2^k)=2^{3k}−1 gives
3·21=63 at k=2, but a tournament with the product count is not a
witness for either factor. By ord₇(2)=3 and LTE,
v₇(Φ₃(2^k))=0 if 3|k, and 1+v₇(k) otherwise.
If Φ₃(2^k) is prime, k must be a power of three:
write k=3^a d with gcd(3,d)=1 and use
Φ_{3^{a+1}}(2) as a proper factor if d>1.

## Additive/multiplicative division of labor

For primes p,q≥5, in a proposed tournament join H=p+q−1,
the modulo-six types of the ports are forced: H≡1 requires (+,+),
H≡3 requires (−,−), H≡5 requires (+,−) in either order.
Here + and − mean residues +1 and −1 mod six.
Thus a local opposite-type seam cannot act as a single universal −1
join even under Goldbach.

First reduce multiplicatively via SCC ordinal sums. One needs to
construct primes other than 7, products 7p with p an odd prime other
than 3 and 7, and 49,63,343; the factorization argument is immediate
by separating the 7-adic valuation and using 49,343 to generate each
exponent ≥2. Past a finite cutoff, these targets are ±1 modulo six.
Suppose one can build physical offset seams H=p+q±1, each for opposite
prime types with known ready port witnesses, and sufficiently large
multiples M of six have a Goldbach pair among the required ready
ports. Then for target N≡δ mod six, δ=±1, choose M=N−δ and that
offset. This covers only the required generators, leaving the finite
cutoff and finite atoms to verify. This is conditional both on exact
physical seam construction and ready/compatible Goldbach witnesses.
Quantitatively, if ordered Goldbach representations have a uniform
lower bound cM/log²M and each role has o(M/log²M) unavailable primes,
the union bound preserves one pair. Relative density one alone
does not supply that little-o bound.

## A fixed prime reservoir cannot finish

Let the left bank be any finite set of (prime r, offset δ), δ=±1,
r≥5, with intended right prime q=N−r−δ. Assign a distinct prime
ℓ>r+1 to each bank element. By CRT one may choose N≡r+δ modulo
its ℓ for all elements simultaneously while prescribing N≡1
or 5 modulo six. The combined progression is primitive, so Dirichlet
gives infinitely many prime N in it. For sufficiently large such N,
every corresponding q is a proper multiple of its assigned ℓ,
hence no bank element serves. The same argument for 7 times a prime
uses N=7p and p≡(r+δ)/7 modulo each assigned ℓ≠7.
This rules out a static finite left bank only when the right side
is required to be prime. It neither obstructs target-dependent banks
nor rules out composite right profiles.

Archive threads to consult: `arithmetic_braids_20260917_collatz.md`,
§4, already refutes an early Collatz-core→Goldbach shortcut;
`paley_h31_exact_rate_and_apery_wheel_20261014.md`, lines 33–38,
points to the inherited twin/Goldbach parity synthesis
`d19_parity_shift_synthesis_20261014.md` (THM-4650/4651).
