# Multiply-perfect numbers: a real inherited tournament interface

Session65. **RECOVERED:** HYP-3546, part(A), is the Euclid–Euler theorem
restated through the dimension of the tournament arc cube. Its other
topological and LRC speculation is not imported. **PROVED elementary
extension:** the complete multiperfect test below on that specified family.
**PROVED classical necessary condition:** the odd-factor valuation budget.
No odd multiperfect nonexistence theorem or new priority claim is made.

## 1. Exact arc-family transport

At tournament order v=2^r, r>=2, the number of independently oriented pairs is

    D_r=binom(2^r,2)=2^(r-1)*(2^r-1).

Put M_r=2^r-1. Multiplicativity of the divisor sum, and coprimality of M_r
with 2, give the exact identity

    sigma(D_r)=sigma(2^(r-1))*sigma(M_r)=M_r*sigma(M_r),
    sigma(D_r)/D_r=sigma(M_r)/2^(r-1).

Therefore D_r is multiply perfect if and only if

    v_2(sigma(M_r)) >= r-1,

and its integer multiplier is sigma(M_r)/2^(r-1). The multiplier is at
least2 because sigma(M_r)>=M_r+1=2^r. Equality2 is equivalent to M_r being
prime: the two divisors1,M_r already sum to2^r, and any additional divisor
strictly raises the sum. Thus the old even-perfect equivalence is recovered
without assuming the composite-apex cases are never multiperfect.

The decisive example is r=4: order16 has120 arcs; M=15 and sigma(15)=24,
so sigma(120)/120=3. It is triply perfect, not perfect. Another tested
composite-apex example is r=10: sigma(1023)=1536, giving the triply perfect
arc dimension523776. These examples are not an exhaustive classification.
Orders4,8,32 give
6,28,496 with multiplier2. The exact test is for this family, not all
triangular numbers and certainly not all integers.

Connection contract: map the tournament order to its complete arc-set size,
then to the Mersenne apex. The scalar arithmetic identity preserves the
multiperfect predicate on D_r. It destroys every tournament orientation,
path, module, and native move. The missing divisor factorization/valuation
must be retained to evaluate the predicate; no tournament structure proves
the primality of the apex. For r>=2, D_r is even, so this family cannot
represent a hypothetical odd multiperfect witness. This last observation
blocks a purported reduction of the odd nonexistence problem through this
restricted family.

## 2. The exact odd valuation ledger, not a sufficient existence test

Let n be odd and n=product p_i^a_i. Then

    v_2(sigma(n))
      = sum_(a_i odd) [v_2(p_i+1)+v_2(a_i+1)-1].

For even a_i the geometric sum has an odd number of odd terms, so contributes
zero. For odd a_i, factor p_i^(a_i+1)-1 and use the elementary two-adic
lifting identity for an even exponent; subtracting v_2(p_i-1) yields the
displayed term. Multiplicativity adds these valuations.

If sigma(n)=k*n, this ledger must equal v_2(k). Every odd exponent costs at
least one unit. In particular:

* odd k forces all a_i even, hence n is a square;
* k=2 forces exactly one odd exponent a, with p=a=1 modulo4. This is the
  classical Euler form, with the remaining prime powers a square;
* in general the number of odd exponents is at most v_2(k).

These are necessary, not sufficient. The odd integer21 has sigma(21)=32,
but21 does not divide32. Its large binary valuation does not supply the
missing odd divisibility. The cancellation in section1 was special to D_r.

For primary literature on this ledger and its extensions, see Chen–Luo,
*Odd multiperfect numbers*, Bulletin of the Australian Mathematical Society
88 (2013), 56–63, DOI10.1017/S0004972712000858; the accessible version is
https://arxiv.org/html/1102.4396v1, section2. The proof used here is the
elementary valuation calculation above, not an import of every claim in
that paper. The inherited HYP-3546 snapshot is kept with its mixed status.

## 3. Where a Turing reduction actually applies

For each specified n, sigma(n) is computable by finite divisor enumeration.
A machine can scan odd n>1 and halt when n divides sigma(n). Such halting is
equivalent to the existence of an odd multiperfect number, not to a finite
tournament-graph isomorphism. The witness has finite positive verification;
absence of all such witnesses is a universal statement. This describes the
semantics of a search, not verification of a particular small machine table.

An odd perfect witness is also an odd multiperfect witness (k=2). Hence
nonexistence of all odd multiperfect numbers implies nonexistence of odd
perfect numbers. No converse or mutual completeness reduction is established.

## 4. Exact finite controls

    python -B check.py > normal.json
    python -O -B check.py > optimized.json
    cmp normal.json optimized.json

The standard-library verifier checks all arc-family exponents2..24,
independently matches factored sigma against literal divisor enumeration
for n=1..10,000, and checks the odd valuation/parity laws on all5,000 odd
inputs. It includes120 and21 as the scope hostiles. These finite tests do
not establish an exhaustive list of multiply-perfect arc dimensions or any
nonexistence result for odd numbers.
