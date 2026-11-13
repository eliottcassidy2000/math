# Finite leading and trailing digits with the physical ports retained

**PROVED IN NOTES.** Every finite base11 leading word and every unit
suffix can be imposed simultaneously in the actual module-repair family,
even while requiring exact seven-adic valuation one. These are infinitely
many **finite tournament counts with specified finite digit windows**.
The middle digits and total length are unrestricted. This is not selection
of an arbitrary complete integer, or a claim of prime cofactors.

Fix n≥1 and put

    A=15·100^n, B=11^n, Y_m=A·2^m−B, m≥0.

[selection/PROOF §§2–5](selection/PROOF.md) constructs both strong repair
witnesses, of orders5n+m+5 and5n+m+6, with the common primitive rooted
state(Y_m,5B,5B,Y_m+6B). Thus every exponent selected here retains those
actual graphs, the one-arc neutral repair and its fresh-sealed two-port
interface. No infinite graph or merely formal count is used.

## 1. The finite-digit theorem

Let k≥1, let c be a unit modulo11^k, and let w be the integer represented
by a prescribed length-d base11 word whose first digit is nonzero:

    11^(d−1) ≤ w < 11^d.

Then infinitely many m≥0 satisfy all three conditions:

1. Y_m≡c modulo11^k;
2. v_7(Y_m)=1;
3. the first d base11 digits of Y_m are exactly the word w.

Condition2 may be omitted. Neither an upper bound on total digit length
nor any condition on the intervening digits is asserted. In particular,
conflicting overlapping digit windows are avoided by allowing arbitrarily
large counts, rather than declaring a specified short integer attainable.

## 2. An exponent progression fixes both arithmetic conditions

Here is the elementary arithmetic behind
[RESIDUE_SELECTOR §§1–2](RESIDUE_SELECTOR.md). Two has order10 modulo11:
2^5≡−1 and2^2≠1. Also

    2^10=1+93·11,  93≢0 modulo11.

Binomial expansion inductively gives

    2^(10·11^(j−1))≡1+93·11^j modulo11^(j+1), j≥1.

Consequently the order of2 modulo11^k is exactly

    L=10·11^(k−1).

Its powers list all units. Multiplication by A preserves that set, as
does subtraction of B because n≥1. Thus a unique m0 in0,...,L−1 solves
the chosen suffix, and every m=m0+Lt retains it.

Modulo7 the actual count is2^(n+m)−2^(2n), so it is divisible by7
exactly when m≡n modulo3. Write m=n+3b in that case. In arithmetic
modulo49, division by11^n is legitimate and

    Y_m/11^n = 15·(200/11)^n·8^b−1
             ≡ 7(2+5n+b) modulo49,

using200/11≡1+35,15≡1+14 and8≡1+7. This remains true for negative b:
these are units, and their inverses have the same first-order binomial
expansion. Hence v_7(Y_m)=1 exactly when2+5n+b≠0 modulo7.

Now gcd(L,21)=1. Among the21 exponent classes m0+Lt modulo21L,
exactly seven satisfy m≡n modulo3. Moving between consecutive such
choices adds3L to m and hence the unit L to b modulo7. Exactly one of
the seven fails the valuation-one test. Choose any of the six survivors,
represented by m*≥0. Then **every** exponent

    m=m*+21L·ell, ell≥0,                            (1)

has both the prescribed11-adic suffix and exact factor7. L is even;
only coprimality with21 is needed here. No coprimality with2 is claimed.

## 3. Irrational rotation gives leading digits on that same progression

We use only the elementary fact that log_11(2) is irrational. Otherwise
log_11(2)=a/b with positive integers a,b would imply2^b=11^a,
contradicting unique prime factorization. No transcendence theorem is
required. Multiplying this irrational number by the positive integer21L
preserves irrationality.

For completeness, every forward orbit of an irrational rotation is dense,
and so is every tail of that orbit. To see this, partition the unit circle
into N equal intervals and compare0,alpha,...,N alpha. Two lie in one
interval, so some positive q has0<||q alpha||<1/N. Write that distance as
epsilon. The multiples j q alpha for0≤j≤floor(1/epsilon) form an
epsilon-net of the circle: they step by epsilon in one of the two
orientations, with a final gap at most epsilon. Starting from any given
orbit point simply translates this net. Taking N arbitrarily large proves
that every open interval is visited after every prescribed starting time,
and therefore is visited infinitely often.

The leading word w corresponds to the nonempty logarithmic interval

    I=[log_11(w)−(d−1), log_11(w+1)−(d−1)) ⊆ [0,1).

For a positive integer with at least d digits, membership of the
fractional part of its base11 logarithm in I is exactly the desired
leading-word condition. Select a closed subinterval J strictly inside I,
with positive length and positive distance from both endpoints. This
also avoids wraparound when I touches0 or1.

Along(1), the fractional parts of

    log_11(A)+m log_11(2)

are dense, with infinitely many visits to the interior of J. The actual
count has logarithm

    log_11(Y_m)=log_11(A)+m log_11(2)
                  +log_11(1−B/(A·2^m)).             (2)

The last term is negative and tends to zero. Once m is sufficiently large,
its absolute value is smaller than the distance from J to the endpoints
of I. Every subsequent visit to J therefore puts the actual logarithm
inside I. Also Y_m eventually has at least d digits. Thus infinitely many
actual counts satisfy the leading word while retaining both arithmetic
conditions. This proves the theorem.

There is an exact terminating search: compile m*, then enumerate(1) and
test the leading digits using integer division by powers of11. The proof
guarantees eventual success. It supplies no efficient bound for the first
successful exponent; numerical logarithms are unnecessary in the test.

The same argument can prescribe leading and unit trailing digits of the
cofactor Y_m/7: impose count suffix7c, and replace the constant phase
log_11(A) by log_11(A/7). Exact valuation one remains. This does not impose
primality of that cofactor.

## 4. Why this is still sparse rather than an exact target compiler

For fixed n, the number of actual values Y_m≤X is exactly

    max(0, 1+floor(log_2((X+B)/A))),

hence O(log X). Any progression selected above is a subset of these values.
The gaps Y_(m+1)−Y_m=A·2^m grow without bound. Every Y_m is also4 modulo5.
Thus finite-window prefix/suffix feasibility, even with an exact7-adic
condition and a primitive physical interface, cannot imply attainment of
every odd integer. It permits an arbitrarily long unspecified middle;
a full target fixes that middle and the total length as well.

No new numerical search or external literature claim is needed for this
note. Its physical realization and finite modulus certificates are the
linked selection and residue packets; the additional real-phase argument
is the self-contained pigeonhole proof above.
