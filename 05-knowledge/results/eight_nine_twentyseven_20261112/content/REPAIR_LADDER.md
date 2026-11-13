# A count-preserving content ladder from neutral block replacement

**PROVED.** This combines Session295's weighted evaluation identity with
the actual native reset in this packet. Every graph below is literal,
has the stated order and Hamilton count, and has the stated exact raw
and uniform-context contents. It does not give prescribed-target coverage.

## 1. The two inherited physical blocks

Let C and N be the two strong eight-vertex tournaments in
[constructors/PROOF.md, Section3](../constructors/PROOF.md):

    C=(56,117,129,22,100,76,141,59),
    N=(212,137,186,97,170,195,22,72).

The native operation complements every arc of C except the arcs of the
named Hamilton spine `(1,0,4,5,6,2,7,3)`. It produces N with the same
Hamilton count513. Their slot polynomials are

    F_C=x^8+34x^6+232x^4+246x^2,
    F_N=x^8+36x^6+238x^4+238x^2.

Their raw contents are9 and3, respectively. The relevant exact evaluations
are small enough to expose both the common minimum and the early tariff:

| j | F_C(j) | v_3(F_C(j)) | F_N(j) | v_3(F_N(j)) |
|---|---:|---:|---:|---:|
| 1 | 513 | 3 | 513 | 3 |
| 2 | 7128 | 4 | 7320 | 1 |
| 3 | 52353 | 3 | 54225 | 2 |
| 4 | 268128 | 2 | 277728 | 1 |

Thus j=4 attains both raw minima simultaneously. This common minimizing
argument matters: raw fixed divisors do not generally multiply for
different factors.

## 2. Exact content at every repair stage

Choose any word W of length m>=1 over the physical blocks C,N, and form
its ordinal join in that order. Let r be the number of C blocks and s
the number of N blocks, so r+s=m. Then

    |W|=8m,                   H(W)=513^m,
    h(W)=3^(2r+s),
    g(W)=3^(2r+s+1) if r>=1,
    g(W)=3^s       if r=0.                              (1)

Here h and g mean raw ordered-cover content and the exact uniform-context
content, respectively.

*Proof.* The slot polynomial is `F_W=F_C^r F_N^s`. For every integer j,
its3-adic valuation is at least2r+s, with equality at j=4. Thus the raw
3-adic valuation is exactly2r+s. Every other possible raw prime divides
`F_W(1)=513^m`, so only19 needs exclusion. Both F_C(2) and F_N(2) are
units modulo19, excluding19. This proves the raw formula.

[Session295, ARITHMETIC.md](../../session295_construction_atlas_20261113/ARITHMETIC.md)
proves the exact weighted change of basis

    g(W)=gcd_(j>=1) (j−1)! F_W(j).                     (2)

For completeness, evaluating and taking finite differences are mutually
inverse integral triangular transforms. Multiplying coordinate j by
`(j−1)!` keeps both transforms integral, because their lower-triangular
entries only acquire factorial ratios with larger numerator index. Thus
the generated weighted ideals agree. Evaluations beyond the degree add
no smaller ideal, by the same Newton expansion.

Suppose r>=1. At j=1,2,3 the valuations in (2) are respectively

    3r+3s,  4r+s,  3r+2s.

Each is at least2r+s+1. At every j>=4, the factorial contributes at least
one additional factor3 to the raw lower bound. At j=4 this is exactly
one, so the context valuation is2r+s+1. The context and raw contents have
the same prime support by [PROOF.md, TheoremA](PROOF.md), completing the
case r>=1. If r=0, j=2 already attains the raw value s and its factorial
is a unit. Hence g=h=3^s. QED.

Consequently start with m copies of C and replace one block at a time
by its actual native reset N. At every stage the order remains8m and
the Hamilton count remains513^m. The exact context exponents follow

    2m+1, 2m, ..., m+2, m.

Each replacement except the last removes one factor3; the last removes
two. The raw exponent decreases by one at every replacement. For m=1
the sequence is simply3,1. The jump is caused by the first minimizing
evaluator moving from j=4, where the excursion tariff is1, to j=2,
where it is0. The exponent m+1 is absent from this specified family;
no assertion is made that other same-count sources cannot attain it.

## 3. What the quotient forgets

C and N are strong and nonisomorphic: their complete path-cover profiles
differ. The ordered strong-component decomposition therefore recovers the
entire word W. The2^m words give2^m pairwise nonisomorphic tournaments of
order8m with Hamilton count513^m. For fixed r there are `binom(m,r)` such
classes, all with the same full slot polynomial and context response.

This is a concrete loss of information: the polynomial observer remembers
block multiplicities but forgets their order. The source graph retains
that order, while the literal block replacement acts on a chosen named
component. This is an exact family contribution, not a formula for the
total multiplicity at Hamilton count513^m.

The incoming Session294 exceptional-block/triangle families and Session295
lower envelopes supply the arithmetic language. The present addition is
the count-preserving physical repair path and its exact final tariff
change, using the newly supplied native reset.

## 4. A bounded hostile probe for complete normalization

The stronger statement

    min {g(T): H(T)=N}=gcd(N,3)

is **OPEN** here, with the minimum over all tournament orders and N an
attainable positive odd count. The right side is a universal lower bound:
the parity argument in the constructor proof forces3|g whenever3|H.
The neutral N block attains this lower bound at H513; its powers need
not, as (1) shows.

The complete inherited bank through order7 gives77 distinct H fibres.
Exactly73 attain the proposed lower bound within that order range.
The four **bounded-order** misses are:

| H | minimum g through order7 | proposed unrestricted minimum |
|---:|---:|---:|
| 55 | 5 | 1 |
| 171 | 9 | 3 |
| 175 | 7 | 1 |
| 189 | 63 | 3 |

All four minimizing witnesses have order7. These are not global
counterexamples: the finite bank does not decide larger realizations
at the same H. Every fibre's minimum and a literal witness are saved
in `receipt.json`, under `bounded_H_fibres` and
`bounded_normalization_misses`.

There is useful retained structure in the last three misses. Their
literal minimum witnesses are respectively

    T171=(42,28,25,112,97,70,7),
    T175=(14,28,56,112,97,67,7),
    T189=(22,44,88,49,98,69,11).

All three are regular and have identical ordered-cover coordinates
`(C2,...,C7)=(1638,6678,15120,20160,15120,5040)`. Their polynomials are

    x^7+28x^5+112x^3+c*x,  c=30,34,48.

The common weighted tail gcd is126, so their exact content is
`gcd(H,126)`. In every fixed uniform exterior E their actual responses
therefore obey

    H(E[T175])−H(E[T171])=4H(E[TT1]),
    H(E[T189])−H(E[T171])=18H(E[TT1]).

This produces literal linear response directions from the hostile
probe. It supplies no count-preserving normalization by itself.

## 5. Reproduction and scope

The literal source profiles and neutral reset are independently computed
by `../constructors/escape.py` and recorded in `../constructors/escape_receipt.json`.
The displayed four evaluations are direct substitutions into their complete
polynomials. The all-m result is the proof above and requires no additional
census or claim of exhaustive source replacement. Incoming attribution is
pinned to `origin/main` commit `9138a603e2`.

The bounded normalization table is generated by the same321-representative
pass in `check.py`, with532 ordinary classes after converse multiplicities
and transitive orders1,2. No extra tournament enumeration is used. Run
`python3 -B 05-knowledge/results/eight_nine_twentyseven_20261112/content/check.py`
from the repository root to refresh its exact receipt.
