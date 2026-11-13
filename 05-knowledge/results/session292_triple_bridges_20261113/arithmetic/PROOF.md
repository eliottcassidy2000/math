# Conditional CRT fibers and inexhaustible hidden factor births

**Status: PROVED in the scopes below.** This refines, rather than rediscovers,
[E7 of the endpoint compiler](../../endpoint_coverage_compiler_20261111/CONNECTIONS.md).
E7 already constructs actual tournaments with H=99·2^k+23 for all k≥1,
finds all eleven active small primes, computes the survivor density, and
constructs infinitely many composite survivors using 419. Here we minimize
its wheel, give a 360-cell factored decoder, and strengthen its composite
hostile to arbitrarily many hidden factors inside *every* surviving cylinder.
Neither arithmetic survival nor rejection changes those physical realizers.

## 1. The least wheel

Put F(k)=99·2^k+23, k≥1. For the first 25 primes (exactly those <100),
p divides F(k) iff k belongs to the following active class:

|p|r|d=ord_p(2)|
|---:|---:|---:|
|5|3|4|
|13|1|12|
|17|1|8|
|19|9|18|
|29|27|28|
|37|23|36|
|53|25|52|
|59|35|58|
|61|0|60|
|67|56|66|
|83|15|82|

The derivation is inherited: solve 2^k=−23/99 in the cyclic subgroup
<2> of F_p^*, with direct exceptions when p divides 2·99·23. No other
small prime divides a term. The 29 and 37 classes are contained in the
5 class: both periods are multiples of 4 and both residues are 3 mod4.
Delete those redundant clauses. The lcm of the remaining nine periods is

    M=61,209,720=2^3·3^2·5·11·13·29·41.

This is the **least** period of the accepted/rejected exponent sequence.
For every prime q|M, the following k and k+M/q have opposite acceptance:

|q|k|the unique rejecting prime on the rejected side|
|---:|---:|---:|
|2|5|17|
|3|21|19|
|5|0|61|
|11|56|67|
|13|29|53|
|29|93|59|
|41|69|83|

These are finite modular certificates, reproduced without constructing F(k).
The k=0 row is harmless: replace both k values by k+M to remain in k≥1.
A least period L of an M-periodic one-sided sequence divides M: periodic
extension to Z gives gcd(L,M) as a period. If L is a proper divisor of M,
it divides M/q for some q|M; that shift would preserve every value, contrary
to the table. This proves minimality.

The earlier period 428,468,040=7M remains a valid period. Its recorded
226,472,946 survivors are seven copies of the new minimal wheel; no inherited
count or density is retracted.

## 2. A conditional product, not independent local probabilities

Keep the shared core C=360 for primes 5,13,17,19,61. The four other clauses
have periods d=66,52,58,82 with new factors q=11,13,29,41 respectively.
Writing k=b+360t, b∈[0,359], the p-clause is possible iff

    b=r_p mod g_p,       g_p=gcd(360,d_p).

If possible it excludes exactly one t class mod q_p=d_p/g_p; if not,
it excludes none. The four q_p are distinct primes and coprime to 360.
Consequently the exact number of lifts of b to Z/MZ is

    N_b = 0, if a core clause rejects b;
    N_b = product_p [q_p − 1_(b=r_p mod g_p)], otherwise.

**Proof.** Solving 360t=r_p−b mod d_p first imposes the g_p condition;
dividing by g_p leaves an invertible coefficient mod q_p and one forbidden
residue. CRT makes these four *conditional* coordinates independent. Their
sum and distribution modulo 8 are

    sum_b N_b = 32,353,278;
    density   = 5,392,213 / 10,201,620;
    mod8: (6909279,0,7419360,0,6909279,3696000,7419360,0).

The script proves its finite arithmetic inputs directly, verifies all
61,209,720 residues by a separate Boolean sieve, and also checks the count
by generalized-CRT inclusion/exclusion. These are certificates for a symbolic
periodic theorem, not a new primality census.

**General fiber principle.** For arbitrary finite congruence constraints,
retain a core C containing the overlaps. If d_i=g_i q_i, g_i|C,
gcd(C,q_i)=1 and q_i pairwise coprime, the same conditional product applies
(with multiple forbidden residues replaced by their distinct lifted residue
sets). Without these hypotheses, use compatible fiber products rather than
multiply counts. Pairwise nonempty intersections do not suffice for a shared
choice: {0,1}, {1,2}, {0,2} in a three-element state space are the minimal
three-guard hostile.

## 3. Every cylinder hides arbitrarily many new primes

**Theorem (hidden-factor refinement).** Let a,c be positive integers, b≥2,
gcd(b,c)=1, and f(k)=a b^k+c. Fix any infinite progression k=r mod M,
k≥k_min≥0, and any finite forbidden set of primes P that divides no f(k)
on this progression. For every h there is a subprogression on which every
f(k) is divisible by at least h distinct primes outside P. The primes and
subprogression can be found by finite search. In particular every surviving
25-prime cylinder for F contains such subcylinders with h primes >97.

**Proof.** No prime dividing b divides any f(k). Suppose S is the finite
list of primes already forced on the current progression. Choose one
exponent k0 there. For p∈S let e_p=v_p(f(k0)). Refine the progression by
requiring k=k0 modulo ord_(p^(e_p+1))(b), for every p∈S, as well as its
old modulus. All these constraints have the common solution k0. On the
refined progression f(k)=f(k0) mod p^(e_p+1), so every S-valuation is fixed
at e_p, not merely bounded below. Since f(k) grows without bound, eventually
f(k)>product_(p∈S) p^e_p. Therefore some prime q∉S divides this larger f(k1).
It cannot belong to P and cannot divide b. Refine again by
k=k1 mod ord_q(b), retaining all previous conditions. Now q is forced as
well as S. Repeat h times. Each growth/search step terminates, because the
progression is infinite and f is unbounded. □

This is an elementary exponential-sequence application of valuation pinning,
not a claim of a new general theorem in prime-distribution theory. No
unproved prime infinitude in this sequence is needed. It strengthens the
inherited single-419 hostile into all-cylinder, unbounded-complexity failure
of any fixed finite sieve to certify primality. The explicit inherited
419 cylinder can already be shortened to k=2+1,162,984,680t, t≥1.
An independent concrete four-factor control is

    F(34)=149·347·2111·15583=1,700,807,049,239.

Their base-2 orders are respectively148,346,1055,1113. Consequently every
exponent k=34+30,670,763,342,539,320t, t≥0, passes the small-prime wheel
and has these four distinct factors. This is a certificate of the refinement
mechanism, not the proof for arbitrary h (which is the induction above).

## 4. What this gives the tournament task

The exact 99 ray is physical and source-valid. Its successor arithmetic
F(k+1)=2F(k)−23 already comes from adding the prescribed tail vertex.
The sieve tells which exponents satisfy a finite divisibility guard; it
neither chooses an arbitrary target nor supplies new endpoint interfaces.

The transferable tool is **condition before composing**. A candidate
constructor with several modular guards can store the shared clock b and
multiply only the independent residual fibers. The necessary extension
is to replace b by the full retained endpoint/occupancy state and prove
that each permitted fiber corresponds to an actual graph construction.
Session291's twelve formal port choices have a true divisor G(m), but their
formal sum is not yet a physical addition. There is no map here that replaces
that missing operation by CRT.

The triple hostile above and the hidden-factor theorem identify different
failures: shared local guards can already be incompatible; even a complete
finite guard can omit arbitrarily complex hidden structure. The positive
counterpart is the marked finite-carrier identity theorem in the manuscript
lane: when the **entire response family** has a proved dimension bound,
finite observations really can determine all continuation. Seek that bound
for source-valid endpoint transitions, not merely for residues of H.
