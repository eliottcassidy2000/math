# Phase-clock symmetry, an exact drift sum, and an eleven-clock construction

**PROVED IN NOTE:** the exact reversal criterion for integral coefficients;
dihedral invariance and a drift sum rule; complete uniformity boundaries
within the Fermat-shaped and Mersenne-shaped families; and an actual
infinite H25/H3 paired-rebase family using clock11. **FINITE-EXACT:**
the catalog of all odd primes through43 and the declared controls.
No classification of all successful prime clocks is claimed.

This follows the [four-phase construction](../actuators/PROOF.md) and
[sevenfold realization](PERIODIC.md). It separates two facts that the
small clocks3,5,7 happen to share: integral coefficients and uniform
additive drift. Neither property is equivalent to primality.

## 1. The exact phase object

Let p>=3 be odd, k=ord_p(2), Q=2^k−1, h=Q/p. The formulas in this
section also apply when p is composite. For1<=i<p put

    r_i=h i,             β_i=i mod2,
    a_i=p(Q−rev_k(r_i))/Q,
    d_i=1 if 2i<p, otherwise0,
    P_i(n)=p g(n)−a_i n+p d_i(n mod2).             (1)

The inherited exact block formula for g gives

    P_i(2^k n+r_i)=P_i(n)+κ_i,
    κ_i=p(g(r_i)−1+d_i β_i)−a_i r_i.              (2)

Although a_i may be rational, κ_i is always an integer: the product
`a_i r_i=i(Q−rev_k(r_i))` is an integer. Thus failure of integrality
in the nonconstant coefficient is not failure of the drift identity.

**Exact integrality criterion.** All a_i are integers iff reversal of
k-bit words preserves the set

    {h,2h,...,(p−1)h}.                             (3)

Indeed `a_i=p−rev_k(hi)/h`. Since reversal is an involution, divisibility
of all the reversed values by h makes it a permutation of (3). It is
enough to check one representative of every doubling orbit, because
reversal converts cyclic rotation into the opposite cyclic rotation.
This is a concrete arithmetic predicate, not a claim that every prime
has an integral phase bank.

## 2. Which symmetry forces common drift?

Use indices modulo p in1,...,p−1. Left rotation of the k-bit word gives

    a_(2i)=(a_i+p d_i)/2,
    a_(p−i)=p−a_i,       d_(p−i)=1−d_i.            (4)

Consequently, when the coefficients are integral,

    (2i)a_(2i)=i a_i mod p,
    (p−i)a_(p−i)=i a_i mod p.                      (5)

Thus i a_i is constant modulo p on every orbit under doubling and
negation. If p is prime and `<2,−1>=F_p*`, there is one such orbit and

    a_i = K i^(−1) mod p,    K=a_1.                (6)

The residues a_i lie strictly between0 and p, so (6) specifies their
integer representatives. Integrality remains a separate hypothesis.

The drift has the same orbit invariance even without that hypothesis.
Define `c_i=β_i(p d_i−a_i)`. Pulling back the last bit of the phase word
gives

    P_i(2n+β_i)=P_(i/2)(n)+c_i.                    (7)

After k such steps, κ_i is the sum of c over that cyclic itinerary.
It is therefore unchanged by doubling. For negation, set
`X_i=p d_i−a_i`. Then

    c_i−c_(p−i)=X_i.

Summing (4) around a doubling orbit gives `Σ p d_i=Σ a_i`, so the sum
of X_i vanishes. Hence

    κ_(2i)=κ_i,       κ_(p−i)=κ_i.                 (8)

For a composite p an itinerary may repeat a shorter orbit; the same
argument with k terms still holds. For prime p, every nonzero phase
has period k.

**Uniformity criterion.** If p is prime and doubling together with
negation is transitive on F_p*, then all phases have the same drift.
In general, common drift is equivalent to equality of the values on
representatives of those orbits. Transitivity is sufficient, not claimed
necessary. This explains why p17 has room for two distinct drifts while
p11 is forced to have just one.

## 3. A telescoping law determines the total drift

Put `χ4(p)=1` when p=1 mod4 and `−1` when p=3 mod4. Then

    Σ_(i=1)^(p−1) κ_i = k(p−χ4(p))²/8.             (9)

This holds for every odd p, with k as above. In particular a uniform
prime clock must have

    κ = k(p−χ4(p))² / [8(p−1)].                   (10)

**Proof.** Complementation in (4) first gives
`Σ a_i=p(p−1)/2`. Let j=2i modulo p, represented in1,...,p−1, so
`j=2i−p(1−d_i)`. Multiplying the first identity in (4) by j gives

    j a_j = i a_i + p i d_i −(p/2)(1−d_i)a_i.

Sum over all i; the left and first right sums cancel by permutation.
Writing m=(p−1)/2 gives

    Σ_(i>p/2) a_i = 2Σ_(i<=p/2) i = m(m+1).

Doubling maps the upper half bijectively to the odd residues, and on
that half `a_(2i)=a_i/2`. Hence

    Σ_(i odd) a_i=(p²−1)/8.

It follows that

    Σ_i c_i
      =p·#{odd i<=p/2}−(p²−1)/8
      =(p−χ4(p))²/8.

Each c_i is counted k times when all k-step phase itineraries are
summed, proving (9). Uniformity then gives (10). QED.

For p11, k10, equation (10) forces drift18. For p13 it also forces18,
even though its a_i have denominator5. For p17, (9) forces only the
average16, consistent with its two values2 and30. The sum rule therefore
retains both the successful common clock and the failure boundary.

## 4. Two complete structural families

### Fermat-shaped moduli

For every m>=1 put p=2^m+1. Primality is unnecessary. Its order k is2m:
2^m=-1 modulo p, and a smaller possible order would be at most m, where
2^k−1<p. Now h=2^m−1, and the two m-bit halves of r_i are
`i−1` and its complement. Reversing them proves

    a_i=1+rev_m(i−1).                              (11)

Thus every coefficient is integral at every size in this family.
However, κ_1=2, whereas for m>=3,

    κ_3=6·2^(m−2)+6.                               (12)

To obtain (12), use `r_3=3(2^m−1)`,
`g(r_3)=3·2^(m−2)+3`, `a_3=1+2^(m−2)`, and d_3=β_3=1 in (2).
For m1 and m2 the clocks3 and5 have one doubling orbit and uniform
drift2. Therefore the **only uniform members of this entire shaped
family** are3 and5. In particular a later Fermat prime cannot restore
uniformity merely by being prime; p17 already exhibits the proved split.

### Mersenne-shaped moduli

For every k>=2 put p=2^k−1. The order is k, h=1, and

    a_i=p−rev_k(i)

is always integral. The first phase has

    κ_1=2^(k−1),

while for k>=3 the third phase has

    κ_3=2^(k−2)+2.                                 (13)

This follows from r_1=1, r_3=3, g(1)=g(3)=1, and their explicit reversals.
The two displayed values differ for every k>=4. At k2 and k3, doubling
and negation are transitive, giving uniform clocks3 and7. Hence these
are the **only uniform members of the Mersenne-shaped family**.

These two classifications explain the small3/5/7 pattern without
asserting that all integral uniform clocks belong to those families.
The next section supplies an explicit counterexample to that assertion.

## 5. Clock11 and its actual H25/H3 consumer

For p11, Q1023 and h93. Direct reversal gives
`rev_10(93)=744=8·93`, so a_1=3. Since2 has order10 modulo11, reversal
closure of that one orbit proves every coefficient integral. By (6),

    a_i=3 i^(−1) mod11,
    P_i(1024n+93i)=P_i(n)+18.                     (14)

Use source phase i1 and multiplier9. Its target phase is i9 and integer
carry is0, with a_1=3,d_1=1,a_9=4,d_9=0. Therefore

    F(n)=[25P_1(n)−3P_9(9n)]/11
        =25g(n)−3g(9n)+3n+25(n mod2)

has the exact loop

    F(1024n+93)=F(n)+36.                           (15)

The integer coefficients come from25=3 mod11; the drift is
`(25−3)·18/11=36`. This targets the deletion counts25 and3 already
present in the actual paired tournament constructions.

The physical realization retains the constant. Let the H25 pair at odd
address B be

    M(B)=25g(B)+122B+16(B mod2)+48,      M(B)+2,

with the two six-vertex deletion cores and cuts from
[the paired-reset construction](../../h13_paired_descent_20261111/modular_cycles/PROOF.md).
Let N(C)=3g(C)+12C+10 be the C3 unit pair used in [PERIODIC.md](PERIODIC.md).
Then

    E11(n)=M(64n+31)−N(576n+46)
          =F(n)+3206,
    E11(1024n+93)=E11(n)+36.                       (16)

For example, the source block31 gives
`g(64n+31)=g(n)+n+(n mod2)`, and the target block46 gives
`g(64x+46)=g(x)+34x+26`. Substitution at x=9n proves (16).
Both the odd H25 address and the actual source-dependent constant3206
are part of this construction.

### An explicit infinite family cancels the constant

For every k>=7 put

    q_k=50·2^k−7,
    m_k=125·2^(k−7)−88,
    n_k=1024^m_k q_k+93(1024^m_k−1)/1023.          (17)

Then m_k>0 and

    E11(q_k)=3168−2250·2^(k−6)=−36m_k,
    E11(n_k)=0.                                   (18)

Indeed g(49)=23 and g(449)=219, and the relevant trailing-bit identities
are

    g(q_k)=g(49)+3·2^(k−3)·50−2,
    g(9q_k)=g(449)+279·2^(k−6)·50−30.

Inserting them in (16) proves the first equality in (18); iterating
(15) proves the second. Thus for

    B_k=64n_k+31,       C_k=576n_k+46,

both H25 counts equal the literal C3 pair at right addresses
`2C_k+1` and `4C_k+3`. The first seed is q7=6393 with defect−1332;
37 copies of the ten-bit word0001011101 cancel it. This is an infinite
actual paired family, not merely an abstract gap of the desired size.
No arbitrary-prefix preparation is claimed for this eleven-clock family.

## 6. Scoped catalog and reproduction

The following is the **complete finite catalog of odd primes3..43**.
“Integral” means every a_i, and distinct drifts are listed without
multiplicity. This is not an exhaustive classification of primes.

|p|ord_p(2)|Integral?|Distinct drifts|
|---:|---:|:---:|---|
|3|2|yes|2|
|5|4|yes|2|
|7|3|yes|4|
|11|10|yes|18|
|13|12|no|18|
|17|8|yes|2,30|
|19|18|no|50|
|23|11|no|36|
|29|28|no|98|
|31|5|yes|10,16,38|
|37|36|no|162|
|41|20|no|30,170|
|43|14|yes|54,66,122|

```bash
python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/prime_layers/phase_clocks.py
```

The exact script checks the symmetry identities, drift sum and integrality
at every listed phase; Fermat-shaped moduli for m1..8 and Mersenne-shaped
moduli for k2..8; the eleven-clock actuator at n0..127; three members of
the infinite family; and the four actual H25/C3 charts at n0,1 by subset
path counting. The source compiler is pinned by its Git-blob hash.
The large family members are certified by the all-address identities,
not by enumerating their exponentially many Hamilton paths.

The structural boundary is now precise. Reversal closure decides integer
coefficients. Doubling and negation determine which drifts must agree.
The telescoping law fixes their total. A physical dyadic address solution
retains the otherwise missing constant, and (17) pays that constant.
None of these steps asserts that every original target enters the family.
