# Complete near-square spectrum at branching defect at most six

**PROVED using the complete source classification and exact finite
certificates below.** The source family has infinitely many orders. The
result is a complete intersection with that family, not a classification
of every tournament or a proof of Landau or Chan–Pak.

## 1. Exact answer and quantifiers

Let T be a strong tournament with a marked root r, d>=2 in-neighbours,
and positive defect b=F_r(2)-H(T)<=6. Then the values of H(T) that are
one more than an integer square are EXACTLY

    {4^m+1 : m>=1}
       union {37,101,145,485,3845,4357,10001,15877,35345}.       (1)

Each value has an explicit rooted source in the inherited grammar.
The corresponding A002496 values are exactly

    {F_j : j>=1 and F_j is prime} union {37,101,4357,15877},     (2)
    F_j=2^(2^j)+1.

The conditional phrase in the set is essential: we have not decided
which later Fermat candidates are prime. The four additional primalities
are certified here. Formula(2) is a complete structural classification of
the intersection, with the remaining prime predicate plainly retained.
The exceptional A002496 value2 has no tournament source, since H is odd.
Roots with d=1 are outside this theorem.

For b=6 alone, the complete near-square spectrum is

    {37,101,145,10001,35345}.

There are six rooted ordinary isomorphism classes in this intersection:
two at H37 and one at each other value. Thus its A002496 intersection
is exactly{37,101}.

## 2. The source inventory and exact prefix transfer

The input is the complete classification in
[the fixed-defect proof](../defect/PROOF.md), with the two-inport
transitive-residual languages in `defect/two_inport.json`. For b2,b4,b6
these contain respectively5,16,31 neutral tails

    3^k W, k>=1,

and there are two isolated words at b4, two at b5, and four at b6.
There is no b3 source. The contact letter3 means the backbone vertex
beats both ordered inports. The word is intrinsic because the inports
have a unique orientation and the residual transitive order is unique;
there is no unaccounted word permutation or converse quotient.

The other b6 sectors have no near-square targets:
[MINIMUM_SQUARES](MINIMUM_SQUARES.md) excludes both the d3 collars and
the d2 residual-C3 families. Higher d cannot occur at this budget because
b>=2^d-2. This makes the inventory exhaustive for the present theorem.

Here is a general exact transfer that avoids fitting a dyadic ray. For
the fixed tail kernel V consisting of W, its two inports and root, let
H=H(V), q=H(V-r), s be its number of Hamilton paths starting at r, and

    F=sum_(S containing r) e_V(S,r) H(V\S).

Insert any nonempty tournament P whose vertices lose to r and beat every
other vertex of V. Put Q_P(2)=sum_X H(P[X])H(P\X). A Hamilton path can
visit P only at its beginning and immediately after r, so at most twice.
If both P segments are nonempty, the fixed-kernel contribution is F-q:
the first V segment ends at r and cannot consist only of r. If only the
second P segment is nonempty the contribution is F, while if only the
first is nonempty it is H-s. Thus the exact total is

    (F-q)(Q_P(2)-2H(P))+(F+H-s)H(P)
      = E Q_P(2)+B H(P),
    E=F-q,    B=H+q-s-E.                                  (3)

For P=TT_k this becomes E*2^k+B for every k>=1. The nonempty condition
matters: substituting k0 need not recover the original kernel when q!=s.
The program computes H,F,q,s by literal first-vertex and reverse path
tables for every tail, then independently checks104 actual prefix graphs.
All52 coefficient pairs and the source receipt hash are saved.

## 3. A finite modular certificate that applies to all exponents

For every tail, n^2=E*2^k+B-1 is tested modulo

    16, 3, 5, 7, 13, 17, 31.

For k>=4 the power of2 is0 modulo16. For every other displayed modulus,
2^120=1; the script checks these identities explicitly. Therefore the
necessary square-residue conditions depend only on k modulo120. For each
of all120 classes the receipt records a first obstructing modulus, or0
if none obstructs. A tail with no surviving residue is ruled out for ALL
k>=4. Values k1,2,3 are checked directly.

Exactly the following nine tails survive this modular certificate:

| Defect | Tail W | E | B | Exact square-plus-one exponents |
|---:|---|---:|---:|---|
|4|122|25|5|none|
|4|13|16|1|every even k>=2|
|4|132|31|5|k9|
|4|201|25|5|none|
|4|212|30|5|k1,4,7|
|4|221|34|5|k7|
|4|231|34|5|k7|
|6|1122|69|17|k9|
|6|2211|78|17|k7|

The final column is proved next; survival of congruence checks alone is
not treated as a solution. The nonsurviving tails' small-exponent hits are
b2: H5; b4: H17; b6: H101 and145. The isolated near-square words are
b4: `1`, H5; b6: `221` and `231`, each H37. The other isolated values
are15,13,13,11,15, all excluded directly.

## 4. Exact finite factor equations for the survivors

The B1 survivor is n^2=16*2^k, giving exactly even k>=2. Its values are
4^m+1 for every m>=3; the finite H5,H17 supply m1,m2.

Every other survivor has

    n^2-c^2=C*2^k, c=2 or4.                               (4)

Write c=2^t, C=2^v C_0 with C_0 odd, and K=k+v. The finitely many K
through2t+2 are checked directly. For larger K, n=2^t u with u odd.
Since u^2-1 is divisible by8, K>=2t+3. The two factors n-c,n+c have
gcd2c, and exactly one has2-adic valuation t+1. Therefore their odd
parts a,b satisfy a b=C_0 and

    |2^(K-2t-2) b-a|=1.                                 (5)

Conversely every ordered odd factorization C_0=a b and sign for which
(a+/-1)/b is a power of2 at least2 gives exactly one candidate K through
(5), and the resulting n is checked in(4). There are finitely many such
factor pairs; this proves completeness without any bound on k.

The exact lists are:

| c | C | All(k,n) |
|---:|---:|---|
|2|25|none|
|2|31|(9,126)|
|2|30|(1,8),(4,22),(7,62)|
|2|34|(7,66)|
|4|69|(9,188)|
|4|78|(7,100)|

The receipt saves every successful odd factor pair, sign and power of2,
as well as the small-exponent cases. This gives the final column above
and proves(1). It is a finite certificate of an infinite classification.

## 5. Prime values and the connection contract

The additional values145,485,3845,35345 are divisible by5, and
10001=73*137. The other four exceptional values37,101,4357,15877 receive
full-factorization Lucas certificates using n=6,10,66,126. Finally,
4^m+1 can be prime only when m is a power of2, by the elementary
factorization of x^u+1 for odd u>1. This proves(2).

The source-to-tail map preserves a literal rooted source, the exact defect,
and its future endpoint interface. The map to(E,B) preserves every count
along that prefix family but loses the tail word and source identity.
The further map to a square or prime target loses the available return
operations. Thus this classification improves target-guided enumeration
and proves whole-source-family exclusions; it does not supply a uniform
return constructor. Its arithmetic predicate can be composed with the
exact A135244 extremal interface from [PROOF](PROOF.md#5-an-exact-a002496a135244-extremal-interface),
but the prime guard is still required there.

Run `python -B tail_squares.py`. The code uses explicit checks, literal
kernel path counts, exact modular arithmetic, exhaustive factor pairs and
primality certificates. `tail_squares.json` is the complete certificate.
The source-language proof and the elementary arithmetic proof play distinct
roles; a finite list of square hits by itself would not prove this result.

## 6. Weak sources and the complete certified prime table

The independently reviewed weak SCC factorization in
[RESPONSE_FAMILIES](../metagraph/RESPONSE_FAMILIES.md#4-weak-graphs-have-a-canonical-composition-law)
extends the SAME prime intersection(2) to all rooted tournaments with
at least two inports and defect at most6, whether strong or weak.
Indeed an upstream component before the root's SCC forces H<=b; the
complete small contributions have values1,3,5. A root in the first SCC
has a strong branching source followed by a suffix of Hamilton count t
dividing b. The only new nontrivial suffix case through6 is t3 after a
strong defect2 source, which makes H divisible by3. Transitive suffixes
leave H unchanged. Thus no new square-plus-one value with positive square
root occurs. If the square root0 is allowed, weak H1 is an additional
value; it is not prime. This is a scoped corollary of the audited weak
classification, not a general extension of every strong-source theorem.

Applying(2) to the certified A002496 table through square root1,000,000
admits exactly eight prime values:

    5,17,37,101,257,4357,15877,65537.

The table contains54,110 primes including2. Every one of the other54,101
odd prime values is excluded from ALL source orders at branching defect
at most6. Equivalently, every rooted tournament representing one of those
values and having at least two inports has defect at least7. This is an
all-order exclusion applied to a finite prime table, not an assertion
that7 is attainable or optimal.

`prime_source_filter.py` saves the exact eight admitted source graphs and
the counts in `low_defect_prime_filter.json`. The graphs are selected from
the already certified source families; no new global minimum-order claim
is made.
