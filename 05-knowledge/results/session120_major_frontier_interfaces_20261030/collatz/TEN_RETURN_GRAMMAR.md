# One non-one core: the remaining E12 candidates cannot accumulate 10-returns

**PROVED conditional on the same inherited certificate as
[EIGHTH_DEFECT_TRANSFER.md](EIGHTH_DEFECT_TRANSFER.md).** This is a bounded
structural consequence, not another E12 family closure. The source already
uses return decomposition and the fee below; the explicit small-defect
one-core grammar is the additional deduction here.

**Successor status:** the E12 obstruction classes below are now excluded by
[precision_lifts_ramsey_deficits_20261030/COLLATZ.md](../../precision_lifts_ramsey_deficits_20261030/COLLATZ.md).
The native grammar remains a valid conditional structural statement.

Write gamma=5878847/19531250. Section5, equation(C2), of
`terminal_fees_phase_covers_20261030/COLLATZ.md` proves that every non-one
native10 excursion has price at slope lambda11 at least

    1+gamma=25410097/19531250.

The all-ones native10 excursion u has clock T10=3906250 and price1/5.
For a non-one primitive11 excursion w, cut at every actual10-return. There
is no unfinished suffix, because the endpoint is an11-return. Let a be the
number of all-ones10 blocks and b the number of non-one10 blocks. Necessarily
b>=1. Exactly one11-return occurs, at the endpoint; hence

    S11(w)=sum_blocks S10(block)+1,
    delta11(w)=sum_blocks delta10(block)-1+Delta A,
    Delta=lambda12-lambda11>0.

Consequently

    delta11(w)>=b(1+gamma)+a/5-1+Delta A.           (R)

No intermediate state is erased at higher precision: the decomposition is
at actual10-returns and the single11-return is explicitly charged.

## Low-price rigidity

If delta11(w)<gamma+1/5, then (R) forces b=1 and a=0. Therefore w itself
is a primitive10 excursion whose endpoint also returns modulo5^11.

More generally, put j=5A-8S11(w). In every remaining E12 obstruction from
this packet, 3<=j<=9. Since

    delta11(w)=j/8-A/(5D)<j/8<=9/8,
    D=78125000,

and two non-one10 blocks would give a price strictly above
1+2gamma>9/8, **b=1 for every surviving candidate**. Also

    a < 5(j/8-gamma).

Moreover the inherited core fee implies j10(c)/8>1+gamma, hence
j10(c)>=11. The score identity proved below then strengthens the bound to
a<=floor((j-3)/2). Combining both restrictions gives:

| j of native11 candidate | maximum all-ones10 blocks a |
|---|---:|
|3|0|
|4|0|
|5|1|
|6|1|
|7|2|
|8|2|
|9|3|

In particular every j=3 OR j=4 candidate is already primitive10. This covers
all words in the narrow E12 families `(a12,b11)=(3,1)` and `(1,2)`, since
they all have j=3. Directly, their required prices are respectively less
than2/5 and less than4/5-eta, both less than gamma+1/5.

## Exact one-core normal form and retained next tangent

Every surviving candidate has the form

    w = u^l c u^(a-l),    0<=l<=a<=3,

where c is a non-one primitive10 excursion. Let its tangent at precision11
be Y in F_5^2. The endpoint condition is Y=-a(4,4). Conversely, any such c
and a<=4 yield a primitive11 excursion in every displayed placement:
before c the cumulative tangent is kv for0<k<=l<5, hence nonzero; after c
it is (k-a)v and becomes zero only at the last block. Within any native10
block there is no10-return, and consequently no11-return. This proves both
necessity and sufficiency of the return-level grammar.

The placement also preserves A, S11, and the native11 tangent at precision12.
To see the last assertion without dropping a digit, write q=5^10 and retain

    F_c(X)=(1+q alpha_c)X+q beta_c mod5^12,
    F_u(X)=(1+q alpha_u)X+q beta_u mod5^12.

Let C=(alpha_c,beta_c), U=(alpha_u,beta_u) in (Z/25Z)^2. Since q^2 is
zero modulo5^12, composition ADDS these two-digit pairs, independently of
placement. Therefore

    C+aU = 0 mod5,
    X11 = (C+aU)/5 mod5.

The quotient is well defined by the first condition. The checker computes
U=(4,4) modulo25 for this particular T10 and independently verifies all
3125 pair/a compositions. No arbitrary entry translation is substituted:
this is an exact block rearrangement preserving the stated invariants,
not a claim that the reordered words have identical full trajectories.

Finally, if j10(c)=5A_c-8S10(c), then the exact score and clock identities give

    j11(w)=j10(c)+2a-8.

Indeed each u has eighth-defect1/4 at level10, and the final11-return adds
one to the score. Thus j10(c)=j+8-2a. Exact possible values appear in
`ten_return_grammar.json`. The remaining analytic problem is now to exclude
appropriate low-defect *primitive10 cores* with their retained two-digit
boundary, not to control accumulation of arbitrarily many10-return blocks
inside the candidate native11 excursions.

## Quantized core fees further restrict whole E12 families

Apply the same quantized transfer lemma to the native10 core inequality
S10(c)<=lambda11 A_c-(1+gamma). Here the first allowed eighth-defect is
11/8, so at slope lambda12 the core costs at least1+eta. The all-ones10
block costs exactly

    lambda12 T10-S10(u)=1/5+Delta T10=6/25.

Subtract the unique11-return. The exact grammar therefore gives the tiered
uniform bound

    delta11(w)>=eta+(6/25)a10.                    (T)

For an unpaid E12 family with a12 all-ones11 blocks and b11 non-one11
blocks, let K be the TOTAL number of internal all-ones10 blocks in those
non-one11 blocks. Then

    a12/5+b11 eta+(6/25)K < 1.

This yields the following additional global restrictions:

| (a12,b11) | maximum K |
|---|---:|
|(0,1)|2|
|(1,1)|1|
|(2,1)|0|
|(3,1)|0|
|(0,2)|1|
|(1,2)|0|

In particular every unpaid candidate in the `(2,1)` family, including its
j=5 portion, must use a primitive10 core with no internal10-return. These are conditional
obstruction restrictions: they do not assert that any surviving core exists.
The earlier 21 defect partitions remain necessary but now have these extra
internal-clock constraints; none is falsely claimed eliminated wholesale.

## Reproduction and scope

    python3 check_ten_return_grammar.py > /tmp/ten.json
    python3 -O check_ten_return_grammar.py > /tmp/ten-O.json
    cmp /tmp/ten.json /tmp/ten-O.json

The checker uses exact fractions and explicit exceptions. It checks the
finite grammar arithmetic and modular composition, not the existence of
these native cores. E12 and Collatz convergence remain OPEN.
