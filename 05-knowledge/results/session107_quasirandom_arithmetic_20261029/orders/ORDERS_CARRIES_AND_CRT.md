# Wieferich delay, two-clock slopes, and exact Collatz CRT charts

Session107, 2026-10-29. **PROVED elementary assembly; FINITE-EXACT checks;
OPEN all-source Collatz convergence.** No literature-wide novelty claim.

## Findings first

The primes 1093 and 3511 delay the binary clock at the first prime-square
lift, but do **not** delay the ternary clock. Their base-3 Fermat quotients
are respectively 312 and 7, both nonzero. Consequently the first lifting
character of a Collatz slope `3^J/2^D` forgets D but still detects J.
This gives an exact local obstruction to a naive binary-clock mixing model,
not an obstruction to uniform source residues: inside every fixed exact
valuation cylinder, the source and every endpoint separately are perfectly
uniform modulo every modulus coprime to 6. They are jointly correlated,
not independent. Carry and source mark remain indispensable.

The phrase “3n=1 collatz” is ambiguous. We treat the established 3n+1 map
as primary and state the same finite-word algebra for 3n−1 separately.
The literal equation 3n=1 has no positive integer solution and is not a
Collatz iteration. No convergence claim for either signed map is made.

## 1. Recovered provenance and external facts

The local repository packet
`math/05-knowledge/results/fermat_wieferich_order_20261014e.md`
already proves the full split/grow tower and its source-bound phase compiler.
Its inverse edges are for **pn±1** at the marked prime p; substituting
p=1093 does not turn that into ordinary 3n+1. Its numerical rows, reproduced
here, are inherited. `THM-4644` also converts the orbit splitting into counts
of invariant circulant tournaments. In particular, at 1093, −1 belongs to
<2>, so a doubling-invariant tournament is impossible; that theorem uses
the odd-order part of doubling. This is not a Collatz result.

Session106 `recurrence/SOURCE_INTERVAL_TRANSFER.md` supplies the exact
integer gauge n=r+2^(D+1)t, y=h+2·3^J t and interval survival guards.
The present assembly adds the two-clock quotient readout and explains
precisely which residue correlations these inherited charts do and do not
control.

For external verification, [Dorais–Klyve, JIS 14 (2011)](https://cs.uwaterloo.ca/journals/JIS/VOL14/Klyve/klyve3.html)
define base-2 Wieferich primes by 2^(p−1)=1 modulo p² and report their
search through 6.7·10^15. We use this as a dated primary computational
reference, not a claim that 1093 and 3511 exhaust all such primes or a
current search record. [Keith Conrad's expository note](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/wieferich-primes.pdf)
records the distinction between Fermat's little theorem and the historical
necessary Wieferich condition for a hypothetical first-case counterexample
to Fermat's last theorem. Neither theorem proves Collatz convergence.
All exact facts used below are proved or independently computed here.

## 2. Inherited order lifting, with exact correlation consequences

Let p be an odd prime, a≥2 an integer with p∤a, d=ord_p(a),
and s=v_p(a^d−1).
Then

    L_k=ord_(p^k)(a)=d p^max(0,k−s).

Proof: a returning exponent must be dt. Writing a^d=1+p^s u, p∤u,
the binomial expansion gives v_p((1+p^s u)^t−1)=s+v_p(t).
For p∤t the leading term has valuation s and every other term has larger
valuation. Raising to p increases the valuation by precisely one: the
first term has valuation s+1, the intermediate terms at least 2s+1,
and the last term ps>s+1 for odd p except s=1,p=2, which is excluded.
Iterate and combine with the prime-to-p case. The least t then gives L_k.
Since d divides p−1, this also proves v_p(a^(p−1)−1)=s.
Thus Fermat's little theorem guarantees s≥1, while Wieferich means s≥2.

Multiplication by a on the units modulo p^k has (p−1)p^(k−1)/L_k
cycles. At a lift k→k+1 the pL_k points over any cycle split into p cycles
if k<s, and form one growing cycle otherwise. This is already in the
recovered packet; its following probability interpretation is useful.

For U,V independent uniform exponents modulo L_k,

    P[a^U=a^V mod p^k]=1/L_k.

For every fixed integer lag h the periodic equality correlation is exactly

    average_t 1[a^(t+h)=a^t mod p^k] = 1[L_k divides h].

These follow respectively from injectivity on one period and cancellation
of a^t. Hence a delay keeps both collision rate and old full-period
correlation unchanged, where an ordinary lift decreases the collision rate
by p and destroys that old-period correlation. These are statements about
an exponent ensemble, not about the sequence of valuations of an orbit.

For coprime prime powers the joint clock period is the **lcm**, not the
product, of the local periods, because a single exponent drives every
coordinate. The CRT product alone does not make those clocks independent.

|p|ord_p(2)|v_p(2^ord−1)|ord_(p²)(2)|q_p(3)|ord_p(3)|ord_(p²)(3)|
|---:|---:|---:|---:|---:|---:|---:|
|1093|364|2|364|312|7|7651|
|3511|1755|2|1755|7|270|947970|

Here q_p(a)=(a^(p−1)−1)/p modulo p. In particular 3^7−1=2·1093,
so the small ternary order at 1093 has valuation exactly one. A small
order and a delayed lift are different phenomena.

## 3. Two-clock Fermat-quotient character

For a unit x modulo p² define q_p(x) by the same quotient. It is
well-defined modulo p and

    q_p(xy)=q_p(x)+q_p(y),   q_p(x^h)=h q_p(x)

for all integer h. Multiplying x^(p−1)=1+p q_p(x) modulo p² proves
the first identity; inversion and induction prove the second.
For p≠2,3 and the exact affine slope α=3^J 2^(−D) modulo p²,

    q_p(α)=J q_p(3)−D q_p(2).                         (1)

The additive character exp(2πi q_p(x)/p) is therefore identically one
on <2> at a base-2 Wieferich prime: it is an exact unmixed observable.
For 1093 and 3511, equation (1) reduces to J q_p(3), not zero in
general. Over a formal uniform clock square (J,D) modulo p, this readout
has p points in each fiber, since q_p(3)≠0. The Wieferich defect aligns
its kernel with the D direction; it does not lower the rank of the full
two-clock map. No uniform-clock law is asserted for Collatz trajectories.

A slope congruent to 1 modulo p satisfies α=1+pu modulo p² and
q_p(α)=−u modulo p. Thus it lifts to 1 modulo p² iff q_p(α)=0.
For the two primes in the table, a slope already equal to 1 modulo p
lifts iff J=0 modulo p. This condition is necessary for a full affine
identity lift, but is not sufficient without the carry condition below.

### Joint subgroup repair is exact, but does not randomize a trajectory

Let H_k=<2,3> inside the units modulo p^k. For each of the two marked
primes, reduction H_2→H_1 is surjective and its kernel has order exactly
p: the ambient reduction kernel has order p, and 3^ord_p(3) is a
nonidentity element of that kernel, by the base-3 valuation-one check.
Thus |H_2|=p|H_1|, despite |<2> mod p²|=|<2> mod p|. Moreover

|p|size of H_1|size of H_2|index of H_2 in units mod p²|
|---:|---:|---:|---:|
|1093|364|397852|3|
|3511|3510|12323610|1|

Indeed in the cyclic group of nonzero residues modulo p the subgroup
generated by two elements has size lcm of their orders, giving
lcm(364,7)=364 and lcm(1755,270)=3510. This classical cyclicity can also
be obtained directly: the exponent of a finite subgroup of a field's
multiplicative group is attained by combining commuting prime-power
components; all its elements solve X^e−1=0, so its size is at most e,
and the element of order e makes equality hold. The kernel count proves
the displayed prime-square sizes without requiring a primitive-root
calculation at p². At every higher level the same base-3 valuation-one
argument supplies the entire reduction kernel, so
|H_k|=|H_1|p^(k−1) for all k≥1.

This is a precise repair of the *available slope subgroup*, not a proof
that actual clock pairs visit it uniformly or even exhaust it.

## 4. Carries and the affine identity-lift test

Fix ε in {+1,−1}, an exact positive valuation word a_1,...,a_J,
D_i=sum_(j≤i) a_j, D=D_J, and

    C=sum_(i=1)^J 3^(J−i) 2^(D_(i−1)).

Direct substitution gives the signed odd-map endpoint

    y=(3^J n+εC)/2^D.                                (2)

Modulo an odd modulus this is the affine map F(n)=αn+β,
α=3^J 2^(−D), β=εC 2^(−D). For p≠2,3 it reduces to the identity
modulo p iff 3^J=2^D modulo p and C=0 modulo p. Under those two
conditions write α=1+pu and β=pv modulo p². Then

    F^t(n)=n+t p(un+v) mod p²                         (3)

for every nonnegative t. Indeed composition adds u,v, because every
cross-term has p² as a factor. A point is fixed iff un+v=0 modulo p;
otherwise its period in this formal affine permutation is p. All points
are fixed iff u=v=0. Hence the identity lifts precisely when BOTH

    J q_p(3)−D q_p(2)=0 mod p,    C=0 mod p².          (4)

The sign changes v, not these divisibility criteria. At a Wieferich prime
the D term vanishes but the carry test remains. Formula (3) is iteration
of the fixed affine word map on a finite ring; its iterates are **not**
asserted to be repeated realizations of that valuation word. Repeated
native execution needs every successive exact 2-adic guard.

## 5. Exact CRT uniformity is not orbit mixing

The word (2) has the exact source cylinder

    n=r mod 2^(D+1),
    r=(2^D−εC) (3^J)^(-1) mod 2^(D+1),
    n=r+2^(D+1)t,   y=h+2·3^J t,
    h=(3^J r+εC)/2^D.                               (5)

Proof of exactness: the congruence makes the final quotient odd. Let C_i
be the prefix carry. The identity C−3^(J−i)C_i=0 modulo 2^(D_i)
and invertibility of 3 modulo powers of 2 imply that every prefix
quotient n_i=(3^i n+εC_i)/2^(D_i) is integral. For i<J, its next
recurrence 3n_i+ε=2^(a_(i+1))n_(i+1) makes n_i odd. Together with
the imposed oddness of n_J this proves every valuation is exactly a_i.
Conversely any exact realization gives the final congruence. Since 3 is
odd it specifies one odd residue class modulo 2^(D+1). Equivalently one
can append letters by Session106's unique congruence split. For ε=−1, positive odd inputs
stay positive (3n−1≥2), so the same source range t≥0 is valid.

For any m coprime to 6, take t through any m consecutive integers.
Both n and y separately run once through every residue modulo m,
because their respective slopes 2^(D+1) and 2·3^J are units.
The same is true for every prefix endpoint: its slope in the final-word
parameter is 3^i 2^(D+1−D_i), again a unit. This covers p^k at both
Wieferich primes at every k, with no defect whatsoever.

But the pair (n,y) has just m possible residues, lying on the graph of
(2), each with frequency 1/m, instead of m² pairs with frequency 1/m².
Thus all one-coordinate marginals can be uniform while joint randomness
fails maximally. Exact endpoint/source information must not be discarded.
A finite interval of N consecutive parameters has each residue count
floor(N/m) or ceil(N/m), giving normalized discrepancy at most 1/N.
This survives Session106's interval guards inside a fixed word; changing
word weights or selecting one infinite orbit is a different operation.

## 6. Scope and next genuine obligation

The useful transfer is to track three separate coordinates: the dyadic
exact-word address, the odd-modulus clock phase, and the affine carry.
Wieferich delays constrain the second coordinate, not the first, and do
not erase the third. Neither Fermat's theorem nor CRT supplies temporal
independence along one orbit, uniform survival weights across words, or
source repayment. A usable quasirandomness claim for Collatz must name its
ensemble, observable, conditioning and error, then prove it without
replacing these exact carriers by independent random clocks.

Run `python3 check_orders.py` and `python3 -O check_orders.py`.
The checker uses only the standard library and explicit exceptions.
Its finite scope is stated in the JSON receipt; these checks support but
do not replace the all-parameter proofs above.
