# Independent review of the valuation and fixed-support compiler

**Verdict: PROVED within the stated scopes.** I reviewed A1--A4 in
[`../arithmetic/PROOF.md`](../arithmetic/PROOF.md), inspected the producer's
`verify.py`, and wrote a separate checker importing none of that code.
The exact row criterion, exponent bound, and hand exclusions pass. No
odd-perfect-number nonexistence result or historical novelty claim follows.

## Exact valuation rows

For distinct odd primes p,r the formula is correct. If `d=ord_r(p)>1`,
then `r` divides neither d nor `p-1`; on `d|n` odd-prime lifting gives
`v_r(sigma(p^(n-1)))=v_r(p^d-1)+v_r(n)`. Off that divisibility condition the
valuation is zero. If `d=1`, cancellation of `p-1` leaves exactly `v_r(n)`.
The order-one exception is essential, rather than a cosmetic replacement.

Every requested positive valuation imposes both an activation divisor and
an exact r-adic exponent in n. Taking their LCM L with the role requirement
does three things at once: it makes all positives active, detects valuation
overflow caused by other targets' orders, and detects forced activation of
a requested zero coordinate. These tests are necessary because `L|n`.
They are sufficient because n=L realizes them all whenever L>1.

The role check is exact: ordinary lengths are odd, and Euler-role lengths
have 2-adic valuation exactly one. Orders with too much 2-part cannot
activate; those remaining cannot force the wrong role when combined by LCM.
For an all-zero ordinary row L=1, a sufficiently large odd prime length
avoids all finitely many activation divisors while making the exponent
positive. This covers the subtle zero-exponent boundary. An empty target
set causes no additional failure. Local Euler role is an exponent condition;
the prime itself must additionally be 1 modulo4 when used as an actual
odd-perfect Euler prime.

The construction controls only named targets, not the entire factorization
of the geometric sum. Distinct source rows have independent exponents;
an imposed shared exponent needs another compatibility step. The proof
retains both restrictions explicitly.

## The joint order obstruction

The data `ord_11(3)=5`, baseline2, and `ord_23(3)=11`, baseline1, give the
claimed strict improvement over a support-only constraint. Requiring both
edges forces `55|a+1`, hence at least three factors of11. The exact vector
(2,1) is impossible although each requirement alone is realizable; length55
realizes (3,1). Since11 is3 modulo4, it cannot be the Euler prime. Its
exponent in a hypothetical odd perfect number is even, so the conditional
conclusion is `11^4|N`. It requires the two actual transfers, not merely
three primes occurring somewhere in the support.

## Finite exponent bound

For an **exact specified complete support** S and fixed Euler role q,
every possible supplier p is included in the budget equation for r.
Each allowed supplier contributes at most its order baseline plus
`v_r(a_p+1)`; each disallowed supplier contributes zero. With
`A=max a_p`, this proves

    A <= C+(|S|-1)log_3(A+1).

The producer's choice of B is a genuine certificate. For h=|S|-1 and
integer k>=2h, the ratio of successive values of
`3^(k-C)/(k+1)^h` is

    3*(1-1/(k+2))^h >= 3*(1-h/(k+2)) > 3/2.

Thus a strict crossing at B persists for every larger integer. The same
ratio proves such a B exists. There is no numerical-logarithm rounding
assumption, and no inference from a finite search cutoff to the infinite
tail. All exponents satisfying the budget equations lie strictly below B.

The dependence on the actual support primes through C is indispensable.
This proves termination for each specified support, not a uniform finite
list of possible supports. Missing suppliers from exterior primes would
invalidate the budget upper bound, so it cannot be applied to an induced
subgraph or an arbitrary partial support. Conversely, even exact internal
budgets are used only as necessary tests; the final full sigma equality
retains possible exterior output factors. The added empty-support,
singleton-support, and no-possible-Euler-prime paragraph correctly covers
the trivial decision cases outside the stated s>=2 theorem.

## Hand exclusion and independent tests

On S={3,5,7,11,13}, the maximal role-sensitive graphs are strongly
connected for both candidate Euler primes, and the abundancy ceiling exceeds
2. The stronger calculation is valid:

- q=5 gives C=3 and B=13. Only source3 can supply target13; within the
  bound its supply is zero or one, less than target13's positive even
  exponent.
- q=13 gives C=4 and B=14. Only source11 can supply target5; its length
  is at most13, so its order-one supply of5 has valuation at most one,
  less than target5's positive even exponent.

Both contradictions cover all exponent branches because the bound has
already proved those branches finite. They do not rely on extending a
square-component factorization to larger exponents.

`arithmetic_check.py` independently evaluates literal integer geometric
sums, rather than importing the producer's LTE evaluator. It checks all
9,072 exponent vectors in the two boxes, verifies their final integer sigma
equalities are false, verifies their internal budgets fail, and verifies
the stated hand contradiction on every vector. The vector counts are
3,888 and5,184. Each candidate has exactly one factor2 in its divisor sum,
providing a nonvacuous Euler-parity control.

For the row criterion it also enumerates complete modular periods for
sources3,7,13, targets5,11, and exact target valuations0 through2 in both
roles. The period is665,500 for each source. Incremental geometric sums
modulo5^3 and11^3 test the full accepted and rejected valuation sets.
The full state (partial sums and powers) returns to its initial value;
the period is divisible by4, so it preserves both roles. A period
representative of length1 still has positive-exponent representatives
after adding a whole period. These checks cover 1,996,500 lengths.

The checker separately tests the sharp length55, the two insufficient
single antecedents, the order-one denominator, and exponent nonmonotonicity
`sigma(3^2)=13`, `sigma(3^4)=121`. All require calls remain active under
`-O`. Normal and optimized receipts are byte-identical; the source is
hashed at startup and checked again before output.

From the repository root:

```text
python3 05-knowledge/results/proof_principles_20261103/rank_audit/arithmetic_check.py --output /temporary/path/arithmetic.json
python3 -O 05-knowledge/results/proof_principles_20261103/rank_audit/arithmetic_check.py --output /temporary/path/arithmetic.optimized.json
```

Compare with `arithmetic_receipt.json` and
`arithmetic_receipt.optimized.json`. These finite computations corroborate
the implementation and the explicit consequences; the unbounded statements
are established by the separate arguments above.
