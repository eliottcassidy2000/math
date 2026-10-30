# Paid first-return kernels: a sharp tenth layer and nine next obligations

**PROVED relative to the independently audited ninth-layer theorem:** for
every finite positive valuation word w of integer clock A,

    S_10(w) <= (1220703/1953125) A.                         (C1)

The coefficient is sharp. If w contains a non-one letter, the following
weaker-than-two-clock, but positive, extra margin also holds:

    S_10(w) <= (1220703/1953125) A - 488281/1953125.         (C2)

**PROVED reduction:** the eleventh-layer sharp bound is equivalent to nine
families of short native excursion-price inequalities described below.
**FINITE-EXACT:** complete non-one first-return prices at levels1–4, for all
word lengths and all positive letters, via independently checked potentials
and literal upper-bound witnesses. **SUPERSEDED frontier:** the nine
unbounded-word families at level10 are now closed by the successor below.
The all-layer content bound and universal Collatz remain **OPEN**.

**Successor,2026-10-30:** [terminal fees with two-layer epoch clocks](../terminal_fees_phase_covers_20261030/COLLATZ.md) closes all nine E11 families by an independently checked F2 certificate. It proves the sharp E11 bound and positive non-one slack, then reduces E12 to seven open families. The OPEN E11 statements below describe this packet's former frontier and are superseded by that proof; its E10 theorem and exact-kernel criterion remain in force.

No new certificate for Session113's *overcounting envelope* is claimed.
This proof controls actual content through a different, source-preserving
first-return decomposition.

## 1. Inheritance and reformulation

For a finite positive word x, use

    F_x(X)=(P_x X+B_x)/Q_x,
    P_x=3^len(x), Q_x=2^A_x, A_x=sum(x),
    c(x)=v5(gcd(|P_x-Q_x|,B_x)),
    S_e(w)=sum_(nonempty prefixes x of w) min(c(x),e),
    lambda_e=(5/8)(1-5^(-e)).

Concatenation is chronological, F_(uv)=F_v F_u. The current state is an
invertible integral five-adic affine map; c(x)>=e exactly when F_x is
identity modulo5^e. No source is resampled or changed in the decomposition.

The inherited [Session113 E9 certificate](../session113_symmetry_consumers_20261030/uniform_content/REPAIRED_CERTIFICATE.md)
proves

    S_9(w)<=lambda_9 A,
    S_9(w)<=lambda_9(A-2) if w contains a letter>1.        (1)

That second inequality is the load-bearing input to C1. The prior
[suffix-entry theorem](../entry_cosets_resonant_echoes_20261030/COLLATZ.md)
optimizes fixed-controller means and gives an error depending on the entry
clock. Here controller changes themselves are allowed: cuts occur only at
actual identity returns in a declared quotient, where the lower-level state
really can be treated as identity without forgetting its next digit.

## 2. Primitive returns and their two-digit boundary

An **e-excursion** is a nonempty positive word b such that F_b is identity
modulo5^e, with no proper nonempty prefix having that property. Cut any
actual word at every such return. This gives complete e-excursions followed
by a possibly empty unfinished suffix.

On each complete excursion and the suffix, the original lower-layer charge
is exactly S_e of that block from identity: the actual preceding map is
identity modulo5^e, so right multiplication by it does not change any
truncated identity test. At precision e+1 the preceding map is NOT erased.
Write an e-return map as

    F(X)=(1+5^e alpha)X+5^e beta modulo5^(e+1).

Its retained boundary is (alpha,beta) in F_5^2. For e>=1, chronological
composition adds these pairs, because5^(2e) is divisible by5^(e+1).
Thus the actual next-layer return condition is that the accumulated pair
is(0,0). This is a25-state boundary, not a pair of independent averages.

Define the price of a complete excursion by

    delta_e(b)=lambda_(e+1) A_b-S_e(b).                   (2)

Every new next-layer return incurs exactly one additional unit of charge.
The endpoint fee in(2) includes the level-e charge at the end of b.

### The exceptional all-ones excursion

The map F_1=(3X+1)/2 has order T_e=2*5^(e-1) modulo5^e. Its unique
all-ones e-excursion is1^T_e. Direct geometric-sum/LTE arithmetic gives

    S_e(1^T_e)=(5^e-1)/4,
    delta_e(1^T_e)=1/5,
    (alpha,beta)=(4,4) at its endpoint.                  (3)

For the last identity, F_(1^T)(X)=mX+(m-1), and
((3/2)^T_e-1)/5^e=4 modulo5. The same vector has additive order5 at every
level. This supplies the normalization for the paid return clock.

### Exact repetition checksum: factorial valuation

If b is an e-excursion and r=c(b)>=e, every proper prefix of b has height
less than e. Repeating b cannot change those heights: the preceding block
power is identity modulo5^e. At the Nth endpoint, the elementary geometric-
sum/LTE formula gives c(b^N)=r+v5(N). Consequently, for every N>=1,

    S_infinity(b^N)
       =N(S_e(b)-e+r)+v5(N!)
       =N(S_e(b)-e+r+1/4)-s_5(N)/4,                  (3a)

where s_5 is the base5 digit sum. This is the inherited periodic-mean
mechanism specialized to native first-return blocks; the additional
observation is that every proper phase is already below the return depth.
It gives an exact finite-word checksum, including its clock N A_b, and
connects the new prices to the earlier fixed-controller framework:

    delta_e(b)=lambda_(e+1)A_b-R(empty,b)+(r-e)+1/4.

Here R(empty,b)=S_e(b)-e+r+1/4 is the limiting charge per complete
repetition, not the charge divided by its clock. No uniform bound on R is
assumed in this identity. The digit-sum remainder
is favorable for exact repetitions; arbitrary switching is handled by the
paid kernel instead of pretending to share this factorial formula.

### Concurrent completion grammar: a usable bound with its missing term retained

Incoming Session116, commit `b435f833bc`, proves the complete maximizing-entry
language for each primitive controller and reduces the nonprimitive case to
its primitive root. See its [phase-completion theorem](../session116_dimension_owner_towers_20261030/collatz/PHASE_COMPLETION_THEOREM.md),
which credits the earlier completion/rotation result in our predecessor.
For any native e-excursion b, its finite row-score algorithm gives an exact
R_max(b)>=R(empty,b). Therefore(3a) immediately supplies the valid lower bound

    delta_e(b)>=lambda_(e+1) A_b-R_max(b)+(r-e)+1/4.    (3b)

This is a concrete bridge between the fixed-controller classifier and the
switching kernel. It preserves the original word, its clock, and its actual
endpoint height r; entry optimization is used only as an upper bound on the
subtracted quantity. It never replaces the source or inserts a maximizing
entry into the actual path.

The bridge does not yet settle the nine families below. The controller b
still ranges over unbounded native words, and even a hypothetical universal
periodic estimate R_max(b)<=5A_b/8 would yield only

    delta_e(b)>=(r-e)+1/4-A_b/(8*5^e),                 (3c)

not a uniform one-unit fee. The missing clock-dependent term identifies
exactly why a periodic-mean bound alone is insufficient here. The fixed-
controller optimizer grammar can help certify selected excursion families;
the paid finite-word kernel is still needed for unrestricted switching.

## 3. A five-position potential that permits arbitrary switching — PROVED

On the25 boundary states set

    Phi(alpha,beta)=((-alpha) mod5)/5,                   (4)

where the residue representative is0,...,4. Phi>=0 and Phi(0,0)=0.
Let r(y)=1 when y=(0,0) and0 otherwise.

For an all-ones excursion y=x+(4,4),

    r(y)+Phi(y)-Phi(x) <= 1/5.                           (5)

If the multiplier digit does not wrap, Phi increases by1/5 and there is
no identity return. At a wrap it decreases by4/5; an actual identity return
can add at most1. Translation remains retained for the actual return test.
For an arbitrary other transition x->y,

    r(y)+Phi(y)-Phi(x) <= 1.                             (6)

If y=0 this is1-Phi(x)<=1. Otherwise it is at most4/5.
These are exact inequalities for all states and transitions.

**One-fee lifting theorem.** Suppose S_e(w)<=lambda_e A for all words, and
suppose every e-excursion containing a non-one letter has delta_e>=1.
Then S_(e+1)(w)<=lambda_(e+1)A for every word.

**Proof.** Apply(5) and(6) to every complete excursion in its actual order.
Their fees pay every next-layer return plus the final nonnegative Phi;
the initial state is0. The unfinished suffix has no level-e return and
therefore no next-layer return. Its remaining cost is nonnegative by the
assumed lower-layer bound. Summing the unchanged block clocks proves the
claim. Neither the number of controller changes nor the entry length enters
a separately charged overhead. QED.

A uniform proof that every non-one e-excursion has price at least1, for
all e>=1, would therefore prove every sharp finite-layer inequality by
induction, and then S_infinity<=5A/8 for every finite word. This **OPEN
sufficient hypothesis** is stronger than necessary; the exact kernel test
in section5 remains available if a cheaper non-one excursion is found.

The elementary base S_1<=A/2 follows because one letter cannot be identity
modulo5 (its translation is a unit), so successive identity returns use at
least two positive clock units, including the first one.

### A quantitative slack version

Suppose the stronger input

    S_e(w)<=lambda_e A-gamma whenever w has a non-one,
    gamma>=1

is known. Put Delta=lambda_(e+1)-lambda_e=1/(2*5^e).
A non-one complete excursion has price at least gamma+Delta A_b.
Subtracting its at-most-one right side in(6) leaves at least
`gamma-1+2 Delta`, since a word containing a non-one has clock at least2.
If the only non-one occurs in the unfinished suffix, its cost is at least
gamma+2 Delta instead. Consequently

    S_(e+1)(w)<=lambda_(e+1)A-(gamma-1+2 Delta)           (7)

for every word containing a non-one. No sharpness of this new margin is
asserted. A complete non-one excursion in fact has clock at least3, but
that improvement is not needed for the uniform statement(7).

## 4. Tenth-layer closure and its exact boundary

Apply the preceding theorem at e=9 with gamma=2 lambda_9>1 from(1).
This proves C1. Formula(7) gives gamma_10=2 lambda_10-1=488281/1953125,
which proves C2. In exact values,

    lambda_9=488281/781250,
    lambda_10=1220703/1953125.

Sharpness follows from3,906,250 consecutive ones: their tenth-layer charge
is2,441,406, exactly lambda_10 times their original clock. This is a finite
word and is realizable by a positive source; no infinite all-ones positive
orbit is inferred.

The spare margin has fallen below1. It is therefore invalid simply to
repeat this one-fee induction at level10. The next section retains more of
the boundary arithmetic to reduce that remaining problem exactly.

## 5. Exact return-kernel criterion, and nine eleventh-layer families

For each tangent X in F_5^2 define g_e(X) as the minimum of delta_e(b) over
all non-one e-excursions ending at X, with infinity for an empty class.
Under the lower-layer bound, every excursion costs at least Delta A_b>0;
there are only finitely many words below a fixed cost. Therefore any finite
minimum is attained. Include the additional all-ones edge v=(4,4) of price1/5.

**Exact lifting criterion.** Given the level-e all-word bound, its sharp
level-e+1 counterpart holds iff every nonempty zero-sum sequence of these
native excursion increments has total price at least1.

**Proof.** Cut at every next-layer return. Each resulting complete group
is a nonempty zero-sum sequence of e-excursions and needs exactly one unit
for its terminal return; if it contains an earlier next-layer return, cut
there too. The final incomplete group and unfinished suffix have
nonnegative cost. Conversely, a zero-sum sequence costing less than1 can
be realized by its minimizing positive words. Their chronological
concatenation has at least one next-layer return but price less than1,
violating the all-word inequality. No artificial kernel increment is allowed
without an actual excursion witness. QED.

### Why the next test has only nine forms

At e=10, C2 gives every non-one primitive excursion of clock A_b>=3 the
price bound

    delta_10(b)>=gamma_10+(lambda_11-lambda_10)A_b
                  >1/4.                              (8)

The final strict inequality is exact: its excess at A_b=3 is1/(4*5^10).
Every all-ones excursion costs1/5. Hence five or more excursions already
pay at least1, and four or more non-one excursions pay strictly more than1.
A nonempty zero-sum sequence of at most four all-ones excursions is impossible.

It remains precisely to check

    sum_(i=1)^b g_10(X_i) >= 1-a/5
    whenever sum_i X_i = -a v,                         (9)

for

    b=1, a=0,1,2,3;
    b=2, a=0,1,2;
    b=3, a=0,1.

These nine families are **necessary and sufficient** for the sharp E11
bound, using the already proved E10 bound and(8). Their boundary has only25
states, but each g_10 still minimizes over arbitrarily long native words.
This is the remaining arithmetic obligation, not a claim that those nine
families have been checked at level10.

Equivalently compute the first three min-plus convolutions of g_10 and
inspect the indicated diagonal targets. If a proposed simpler fee bound
fails, this exact criterion remains available; one need not discard the
return decomposition or revert to independent controller errors.

## 6. Exact small-level price tables and hostile controls

The producer and independent auditor determine all25 non-one prices at each
e=1,2,3,4, over **all lengths and all positive valuation letters**. Their
minimum values are

| e | min_X g_e(X) | A of a minimizing word | tangent |
|---|---:|---:|---|
|1|13/5|6|(1,3)|
|2|74/25|8|(3,0)|
|3|653/125|26|(2,0)|
|4|2953/625|126|(4,1)|

The fourth value is smaller than the third: monotone growth of the fee must
not be extrapolated. All four exceed1, but that is not an all-e proof.
These are new complete first-return cost classifications, not a replacement
for the stronger inherited E9 all-word verification or the analytic E10 lift.

The full modular graph supplies only the *first-return price*. The audited
25-state kernel then supplies a paid potential for arbitrary sequences of
those excursions. Both endpoints and the non-one flag are retained. In
particular, reducing a letter modulo ord_(5^(e+1))(2) must not turn a
non-one letter into the exceptional letter1. The all-letter tail inequality
in the auditor handles that distinction without assuming a bounded alphabet.

See AUDIT.md for exact certificate universes, lower-bound edge inequalities,
100 literal integer witnesses, all-letter tails, corruption control and the
independent650 universal potential inequalities. General Collatz convergence
requires additional links beyond a finite-layer content bound and is not
claimed by C1 or by a possible future all-layer proof.
