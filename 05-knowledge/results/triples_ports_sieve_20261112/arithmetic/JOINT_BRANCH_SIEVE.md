# A joint sieve for three formal ternary counts

Status: **PROVED** arithmetic normal form, density, minimal period and infinite
composite progression; **FINITE-EXACT** local receipts and the prime seed.
The scalar triple below is formal. This note does not establish a simultaneous
tournament source, a preserved target, or return closure for its three entries.

## Inheritance and the new question

The inherited ready family has count `A_k=99*2^k+23`, for `k>=1`.
[WHEEL.md](WHEEL.md) classifies its first-25-prime sieve. The ternary affine
branches suggest testing the stronger arithmetic predicate on

\[
(A_k,B_k,C_k)=(99\,2^k+23,\ 297\,2^k+71,\ 297\,2^k+73)
             =(A_k,3A_k+2,3A_k+4).
\]

The source is an exponent; the target is this scalar triple; the map is the
displayed affine formula. The preserved predicate here is avoidance of a
specified finite set of prime divisors. All graph, root, endpoint and reset
information is absent. A graph interpretation needs those sidecars separately.

The cheap hostile question is whether all three branches can even survive
the same small-prime sieve. They can, but the common exponent phase becomes
much more rigid than for `A_k` alone.

## J1. Exact normal form for all positive exponents

Let `J(k)` mean that none of the first 25 primes, from 2 through 97, divides
any of `A_k,B_k,C_k`. Then `J(k)` holds if and only if

\[
k=6+12t,\qquad t\ge0,
\]

and `t` avoids every class in this table:

| responsible prime | modulus for t | forbidden residues |
|:---|---:|:---|
|23|11|5, 6|
|31|5|1|
|47|23|22|
|53|13|1|
|59|29|11|
|71|35|18|
|79|13|3|
|83|41|13, 16|

Proof. Prime 5 excludes `k=3,1,0 mod 4` through the respective branches
`A,B,C`, so `k=2 mod 4`. Prime 13 excludes `k=2 mod 12` through `C`.
Prime 7 excludes `k=1 mod 3` through `B`, removing `k=10 mod 12`.
Exactly `k=6 mod 12` remains. Substitute `k=6+12t`.

For each odd prime `p`, powers of 2 have period `d=ord_p(2)`. After the
substitution their period divides `d/gcd(d,12)`. Complete enumeration of
that small period gives precisely the table; the receipt retains both the
original branch-specific forbidden classes and the substituted classes.
This is an exhaustive local proof for every exponent, not an extrapolation
from the first 10,000 exponents. All three counts are odd for `k>=1`.

## J2. Density, CRT decomposition and minimal period

The two conditions modulo 5 and 35 combine into 27 allowed residues modulo
35: remove the seven residues `1 mod 5` and the additional residue 18.
Combine the two modulus-13 conditions as well. The independent CRT clocks
are then

\[
\begin{array}{c|rrrrrr}
\text{modulus}&35&11&23&13&29&41\\
\text{allowed residues}&27&9&22&11&28&39.
\end{array}
\]

Consequently the relative density in `t` is

\[
\frac{27}{35}\frac9{11}\frac{22}{23}\frac{11}{13}
\frac{28}{29}\frac{39}{41}=\frac{64152}{136735}.
\]

The exact density among positive exponents `k` is therefore

\[
\boxed{\frac{5346}{136735}}.
\]

A period is

\[
P=12(35\cdot11\cdot23\cdot13\cdot29\cdot41)
 =1642460820,
\]

with exactly `64216152` accepted residues. This period is minimal.
For every prime `q` dividing `P`, the following witness has
`J(k) != J(k+P/q)`:

| q | k | P/q | J(k), J(k+P/q) |
|---:|---:|---:|:---|
|2|6|821230410|1, 0|
|3|6|547486940|1, 0|
|5|54|328492164|1, 0|
|7|222|234637260|0, 1|
|11|30|149314620|1, 0|
|13|30|126343140|1, 0|
|23|270|71411340|0, 1|
|29|834|56636580|0, 1|
|41|162|40060020|0, 1|

Indeed, periods of a nonconstant periodic Boolean function form a subgroup
of the integers. Its least positive period divides `P`; any proper divisor
divides some `P/q`, which these witnesses rule out. Periodicity on positive
integers extends uniquely to this two-sided periodic function. The least
period in `t` is correspondingly `136871735`.

The sum of branchwise densities cannot produce this result: all three
counts share the same exponent. The necessary phase `6 mod 12` is the
first exact example of that dependence. Nor can the three branches choose
their own CRT phases independently.

## J3. A prime seed and an infinite all-composite hostile

At `k=6`, the triple is

\[
(6359,19079,19081).
\]

All three entries are prime, by complete trial division through their
square roots. The last two differ by 2. This is one finite prime seed.
No infinitude of prime triples follows from the positive sieve density.

The same seed produces an explicit infinite hostile progression. Their
respective orders of 2 are `3179,9539,3180`. Set

\[
M=\operatorname{lcm}(P,3179,9539,3180)=239978082932247660.
\]

For every integer `u>=1`, at exponent `k=6+M*u` all three counts pass the
first-25-prime sieve, while `6359|A_k`, `19079|B_k`, and `19081|C_k`.
Each count is larger than its fixed divisor, so **all three are composite**.
The proof uses periodicity and strict growth; the receipt also checks 100
members by modular arithmetic without constructing their huge values.

This repairs the tempting implication from compatible finite prime filters
to prime return counts. What survives is a completely classified finite
sieve and a genuine prime test point. Actual graph return capability and
prime infinitude remain separate questions.

## J4. An inherited sharp aliquot interface

Let `M(n)=A135244(n)`, the largest integer whose proper divisors sum to `n`,
or 0 if none exists. The inherited
[sharp equality theorem, equation (4)](../../uniform_return_metagraph_20261111/arithmetic/PROOF.md)
gives, for every `n>=2`,

\[
M(n)\le(n-1)^2,\qquad
M(n)=(n-1)^2\ \Longleftrightarrow\ n-1\text{ is prime}.
\]

For completeness, any source `m` for such an `n` is composite. If `p` is
its least prime divisor, its largest proper divisor is `m/p>=sqrt(m)`.
The distinct proper divisors 1 and `m/p` give `n>=1+sqrt(m)`. Equality is
possible exactly for `m=p^2`. Conversely a prime square has proper-divisor
sum `p+1`, attaining the bound. This is the previously proved mechanism,
now applied to the joint scalar counts.

The prime seed therefore yields exact identities

\[
M(6360)=40436881,\quad
M(19080)=364008241,\quad
M(19082)=364084561.
\]

All three targets are already within the inherited million-target prefix;
these are not new enumeration terms or a range extension. For every
`k=6+239978082932247660*u`, `u>=1`, J3 instead proves the strict inequalities

\[
M(A_k+1)<A_k^2,\quad M(B_k+1)<B_k^2,\quad M(C_k+1)<C_k^2.
\]

Thus the map from primality to extremal aliquot equality is genuinely
predicate-preserving, while merely passing the finite wheel does not
preserve that equality. The receipt checks the seed proper-divisor sums
directly and the source-bound equality boundary on sources 2 through 5000.

## Reproduction and controls

Run `python -B joint_branch_sieve.py` and
`python -O -B joint_branch_sieve.py` from this directory. All checks remain
active under optimization. Independent count paths are the displayed
coprime-clock product and generalized CRT inclusion-exclusion. Every local
phase word is exhausted; an additional 10,000-exponent literal modular
comparison checks the compressed predicate. Minimality witnesses are
checked against both the compressed predicate and the literal 25-prime
tests. Seed primality uses complete trial division.

The finite universe is the first 25 prime divisors and three displayed
formal affine counts. The periodicity theorem quantifies over all positive
exponents. It does not classify other branch words, all primes, or the
tournament metagraph.
