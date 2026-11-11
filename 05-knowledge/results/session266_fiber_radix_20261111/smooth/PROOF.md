# Complete leading-growth classification for every 3,5,7-smooth count

**PROVED using inherited all-order exclusions of7,21 and complete finite strong
inventories at3,5,35,49,63,343.** This is an infinite-family classification of
polynomial degree and leading coefficient, not all lower coefficients. No
external priority claim. General Chan–Pak remains OPEN.

## 1. Six scalar atoms suffice

Write A_n(h) for ordinary tournament isomorphism classes with h directed
Hamiltonian paths. The inherited SCC law says h multiplies across the ordered
strong components, and the maximal number r(h) of nontrivial strong components
is exactly the eventual polynomial degree of A_n(h).

Inside the multiplicative monoid of 3,5,7-smooth positive odd integers excluding
7 and21, the irreducible nonunit elements are exactly

    3, 5, 35, 49, 63, 343.                            (S1)

Here irreducible refers to this **realizable-count monoid**, not integer
primality or strong connectivity. A strong tournament can have a reducible
scalar count; it then does not contribute a maximal-length core word.

*Proof of (S1) and realizability.* Let h=7^a3^b5^c.
If a=0 use3 and5. If a is positive even use49^(a/2), then3 and5. If a is odd
and c>0 use35 once,49^((a-1)/2), then the remaining3 and5. If a>=3 is odd and
c=0 use343 once,49^((a-3)/2), then3. The remaining case a=1,c=0 is excluded
for b=0,1 and uses63 then3 for b>=2. Each of the six displayed values is
irreducible because every proper integer factorization contains7 or21, or
because it is prime. Every one has an inherited strong realization. Hence all
these arithmetic factorizations are actual tournament SCC constructions, and
there are no other excluded 3,5,7-smooth counts. QED.

This is why the proof does not assume the full odd-spectrum conjecture. It
closes an entire smooth multiplicative submonoid using six known inputs.

The numbers of strong isomorphism classes at these six counts are respectively

    1, 1, 2, 8, 2, 103.                              (S2)

The last value is recovered from the inherited complete polynomial

    P_343(z)=12z^8+67z^9+16z^10+8z^11,

not a fresh census. Source: opus S50's
`04-computation/experiments/h_spectrum_multiplicities_20261021/strong_ears_375.out`,
with the completeness argument in
`05-knowledge/results/h_spectrum_multiplicities_20261021.md`.
The other five entries are already in the frozen Session24 inventory. All
sources are pinned in `SOURCES.json`.

## 2. Exact degree for the entire smooth family

Except for h=7,21, put h=7^a3^b5^c, a,b,c>=0. Then

    r(h)=b+c+floor(a/2),                              (S3)

unless a=1,c=0,b>=2, in which case

    r(7*3^b)=b-1.                                    (S4)

The h=1 convention is r=0, A_n(1)=1. The exceptional two missing values have
A_n=0 and no growth degree.

*Proof and equality cases.* A maximal core word can only use the six scalar
atoms: any other core count can be replaced by at least two actual strong
atoms, increasing word length while preserving h. Let u,v,w be the numbers
of35,63,343 atoms, and let d be the number of49 atoms. Then

    a=u+v+3w+2d,
    r=b+c+d-v+w
      = b+c+a/2-(u+3v+w)/2.                         (S5)

For even a, the minimum penalty is0, attained uniquely by u=v=w=0,d=a/2.
For odd a, the parity constraint makes the penalty at least1. Equality occurs
in exactly two patterns:

    (u,v,w)=(1,0,0), available when c>=1;
    (u,v,w)=(0,0,1), available when a>=3.

If neither is available, a=1,c=0. Its only legal pattern is v=1, requiring
b>=2 and penalty3, giving(S4). These list all maximal multisets, not just
one witness. QED.

For example the numerical square63² has degree5, although juxtaposing two
H63 cores gives only two nontrivial components. The dominant construction
uses49 and four triangles. Squaring an irreducible count need not double its
maximal core length.

## 3. Closed leading coefficients

Write A_n(h)=L(h)n^r+O_h(n^(r-1)) when h>1 is attained. Each individual
ordered core word of length r contributes n^r/r!. If a multiset uses k_q
components of scalar type q, there are r!/product k_q! orders and
product P_q(1)^k_q choices of internal strong class. Therefore its contribution
to L is product P_q(1)^k_q/k_q!.

Let a=2m be even. There is only the49/3/5 maximal pattern, so

    L(7^(2m)3^b5^c)=8^m/(m! b! c!).                  (S6)

Let a=2m+1 be odd and suppose c>=1 or m>=1. The two available maximal
patterns in(S5) give

    L(7^(2m+1)3^b5^c)
      = [2c*8^m/m! + 1_(m>=1)*103*8^(m-1)/(m-1)!] / (b! c!).  (S7)

The exceptional ray has

    L(7*3^b)=2/(b-2)!  (b>=2).                       (S8)

These formulas completely specify the leading coefficient at every attained
3,5,7-smooth count. Component vertex sizes affect lower polynomial terms but
not this leading coefficient. Their exact contribution remains available by
using the full P_q(z) rather than P_q(1).

Examples:

| h | degree r | leading L |
|---|---:|---:|
|63|1|2|
|189|2|2|
|49|1|8|
|343|1|103|
|16807=7^5|2|824|
|1225=35²|3|4|
|3969=63²|5|1/3|
|1715=7³*5|2|119|

Cumulative ordinary multiplicity through N consequently has leading term
L(h)N^(r+1)/(r+1). Converse-merged multiplicities have half the ordinary
leading coefficient for h>1, because the inherited self-converse sector has
strictly smaller degree. These statements do not identify labeled frequencies
with ordinary frequencies; automorphism-weighted strong masses would replace
(S2) in that channel.

## 4. One exponential encodes every leading coefficient

Use formal variables X,Y,Z for valuations at3,5,7. The coefficient of
X^b Y^c Z^a t^r in

    exp(t*(X+Y+2YZ+8Z²+2X²Z+103Z³))                  (S9)

is the sum of leading contributions from r-component multisets. Selecting the
largest nonzero t degree recovers(S3)-(S8). This exponential is not counting
unordered SCCs: the r! orders have already canceled the r! in the asymptotic
binomial coefficient. Keeping that cancellation explicit prevents a false
commutativity claim about the actual tournament construction.

## 5. A toric loss-of-state bridge to radix carries

The map from six core symbols to the monomials X,Y,YZ,Z²,X²Z,Z³ forgets which
factorization produced the same scalar count. It has concrete positive
relations

    35²=49*5²,
    63²=49*3⁴,
    343²=49³,
    35*3²=63*5.

Every side has actual SCC realizations. The sides have different numbers of
nontrivial components, hence different singleton-padding dimensions and
polynomial growth orders. Equality of scalar H does NOT identify their
recursive growth states or provide a bijection between their padded families.

This is a precise parallel with the money-plant map p(x)↦p(B): carry relations
identify distinct coefficient states after scalar evaluation, but their daily
trajectories differ. Here the evaluation is a monomial map in prime valuations,
and the lost state is factorization length, vertex size, order and internal
strong identity. The missing sidecar is the multigraded core word. Both maps
are quotient maps with an explicit kernel mechanism, not a shared numeral.

On the signed Collatz atlas, ordinary product of h codes is the edge-code
product derived in `../collatz/PROOF.md`. Thus these same relations are exact
collisions of composed signed-edge codes, with different multiplicity-growth
orders attached. None gives a Collatz convergence theorem or a physical
operation sending H to its next Collatz code.

## 6. Exact factorization elasticity: the maximum length ratio is5/2

Restrict to factorizations into the six irreducible scalar atoms(S1), rather
than arbitrary strong components with reducible counts. Assign h=7^a3^b5^c
the additive weight ell(h)=a/2+b+c. On the six atoms its values are

    1,1,3/2,1,5/2,3/2.

Thus every atom-factorization of length k satisfies
k<=ell(h)<=(5/2)k. The longest and shortest atom-factorizations of the same h
therefore have ratio at most5/2. Equality is attained by

    63² = 63*63 = 49*3*3*3*3,

and by every power of this relation. The exact elasticity of this scalar
factorization monoid is5/2. This measures how much singleton-padding dimension
can change under scalar equality while every factor remains irreducible.
It is not a bound for arbitrary SCC words: a strong tournament can have a
reducible H and hence create much shorter non-atom words.

Both sides are actual tournament constructors. Two H63 cores need16 vertices;
a minimal H49 core and four triangles need19. Padding the former makes equal
orders available as well. These are equal-count families, not isomorphic
objects or a size-preserving move between individual tournaments.

## 7. Beyond smooth numbers: a conditional full degree law

If the complete attainable monoid were all odd positives except7,21, its
irreducibles would be: every odd prime p!=7; every7p for an odd prime
p!=3,7; and49,63,343. For h=7^a3^b m with gcd(m,21)=1 and c=Omega(m), the
same proof gives r=b+c+floor(a/2), except the ray a=1,m=1 where r=b-1.
This last paragraph is CONDITIONAL on the full support statement, or locally
on the needed prime and7p witnesses. It is not used in the unconditional
3,5,7-smooth theorem. The generally unknown atom masses would still be needed
for leading coefficients and lower polynomial terms.
