# Proofs for the arithmetic and incidence bridges

**PROVED elementary statements and reformulations.** Classical structures
and inherited results are identified at their point of use. There is no
external priority claim.

## A Square divisor sums define a binary dependence system

Fix distinct primes `P={p_1,...,p_k}`. Let L contain all prime divisors of
the numbers `p_i+1`. Form the binary matrix

```text
A[l,i]=v_l(p_i+1) mod 2.
```

For a subset represented by `x in F2^k`, let `n_x=product_(x_i=1) p_i`.
Then

```text
sigma(n_x)=product_(x_i=1)(p_i+1),
Q(sigma(n_x))=A x.
```

Consequently square divisor sum is equivalent to `x in ker A`. Every such
subset yields an actual aliquot total `n_x+s(n_x)=sigma(n_x)`, with the
usual base convention s(1)=0. This is a statement about squarefree sources;
repeating a source prime changes its local factor to `1+p+...+p^a`.

The columns define a binary vector matroid. Nonempty inclusion-minimal
kernel supports are its circuits. This is the classical vector-matroid
construction, applied to the inherited square-class transport identity
from [Session83](../session83_aliquot_divisor_grammar_20261028/CARRIER_AND_144_BRIDGE.md).

For P={2,3,5,7,11,13}, on row coordinates (2,3,7), the columns are

```text
2:010, 3:000, 5:110, 7:100, 11:010, 13:101.
```

Here strings display coordinates in the stated order, rather than binary
integer significance. The only minimal zero supports are

```text
{3}, {2,11}, {2,5,7}, {5,7,11}.
```

The 7 coordinate occurs only at source prime 13, so no kernel support uses
13. The other five columns have rank two, hence eight kernel supports;
their source products are `1,3,22,66,70,210,385,1155`. All are checked by
both parity and independent divisor addition.

The prime 3 is the only possible zero column among ALL prime inputs.
Indeed p+1=t^2 gives p=(t-1)(t+1); primality forces t-1=1, hence p=3.
No uniqueness is claimed for parallel pairs or three-element circuits.

### A1 An explicit arithmetic Fano representation

For the primes `(7,2,19,5,89,59,29)`, the squarefree parts of their
successors are `(2,3,5,6,10,15,30)`. On row primes (2,3,5) these are
precisely the seven nonzero vectors of F2^3. This checks all prime and
square-factor conditions explicitly; it does not rely on a prime-producing
conjecture.

The matrix has rank three, so its kernel has dimension four and 16 elements.
There are no weight-one or weight-two words because the columns are
nonzero and distinct. Each unordered pair of columns determines a third
by addition; the 21 pairs group into seven weight-three words. The sum of
all seven columns is zero, so complementing gives seven weight-four words
and the all-one weight-seven word. Including zero accounts for all 16:

```text
W(z)=1+7z^3+7z^4+z^7.
```

This is the usual parity-check construction of the [7,4,3] Hamming code.
Its weight-three supports are the Fano lines; its weight-four supports are
their complements and are also minimal dependencies. The arithmetic map
is from a subset of the seven named primes to the product source and the
square-class vector of its divisor sum. It does not turn an arbitrary
Fano orientation into a tournament, nor assert that a Hamming code already
solves a Hamiltonian-path problem.

This matroid is not graphic. A representing graph would have neither loops
nor parallel edges, since the matroid has neither one- nor two-element
circuits. If the component ranks are r_i, then sum r_i=3 and a simple
component has at most binom(r_i+1,2) edges. Hence the whole graph has at
most binom(4,2)=6 edges, whereas this matroid has seven. Contradiction.
The positive multigraph model of the six-prime bank therefore has a real
boundary, rather than being a universal hidden-graph explanation.

## B Persistent equal totals with different exceptional layers

The starting equality is

```text
sigma(33)=sigma(35)=sigma(47)=48.
```

All three integers are odd. By multiplicativity, for every a>=0,

```text
sigma(2^a*m)=48(2^(a+1)-1),     m in {33,35,47}.
```

This is the strongest immediately transported predicate: **equal sigma**
survives the entire common doubling operation. It says nothing by itself
about equality of sources, aliquot sums, divisor shapes, or prime supports.

### B1 Exactly one square-total layer

Put k=a+1. If `48(2^k-1)` is square, then `3(2^k-1)` is square, hence
`2^k-1=3u^2` for an integer u. For k>=3 the left side is 7 modulo 8,
whereas `3u^2` is in {0,3,4} modulo 8. This is impossible. At k=1 the
total is 48 and not square. At k=2 it is 144. Therefore square total
occurs iff a=1, at sources 66,70,94. QED.

### B2 Exactly one divisor-balance layer on the two triple towers

For m=33 or 35 and a>=2, the prime exponents are (a,1,1), giving

```text
F=4(a+1)-2=4a+2,   U=3,   S=7.
```

Thus F=U+S iff a=2. At a=1 the source is squarefree, with (F,U,S)=(6,3,6);
at a=0 it has (2,2,2). Neither satisfies the equation. The unique balance
layer is therefore 132 and 140, one step after the unique square layer.

For the third tower m=47, a>=2 gives (F,U,S)=(2a,2,3), so equality would
require 2a=5. The cases a=0 and a=1 are respectively the prime base state
47, which does satisfy F=U+S trivially, and the nonbalanced semiprime 94.
Thus it has no composite F=U+S layer. This distinguishes the prime base
exception from the p^2qr branch.

The complete fibres are **FINITE-EXACT with a proved bound**:

```text
sigma^-1(48)  = {33,35,47},
sigma^-1(144) = {66,70,94,115,119},
sigma^-1(336) = {132,140,182,188,195,249,287,299}.
```

Since sigma(n)>=n+1 for n>1, target t requires 2<=n<=t-1; n=1 is checked
separately. Both factorization and direct divisor addition exhaust these
ranges. The only F=U+S members in the last fibre are 132 and 140.

The general F=U+S classification is inherited, not reproved as new:
prime bases, prime cubes, and p^2qr. The equality of cardinalities does
not identify U and S as disjoint classes. The 923 sorted exponent profiles
in the check are an independent finite control of that inherited formula.

## C Divisor heights are right adjoints to taking powers

For k>=1 put `R_k(N)=product_p p^floor(v_p(N)/k)`. For positive b,N,

```text
b^k divides N  iff  b divides R_k(N).
```

This follows coordinate by coordinate from `k v_p(b)<=v_p(N)`. In the
divisibility poset it says that the power map is left adjoint to R_k.
For fixed b>=2, `h_N(b)>=k` is the same condition. Therefore the number
of nonunit bases of height at least k is `tau(R_k(N))-1`. This recovers
Session43's all-divisor formula and its exact 143-element cycle basin.

The same exponent box underlies F,U,S, but the operations differ:
squarefree truncation uses min(e,1), while R_k uses floor(e/k).
Neither should be substituted for the other merely because both compress
an exponent. Under multiplication the omitted remainder gives an exact
carry:

```text
R_k(NM)=R_k(N)R_k(M) product_p p^floor((r_p+s_p)/k),
r_p=v_p(N) mod k,  s_p=v_p(M) mod k.
```

For example R_2(2)=1 but R_2(4)=2. Retaining the remainders repairs the
product law; coprime multiplication has no carry.

At `N=2^5 3^7 5^3 7^4 11 13^2` the repeated-divisibility implementation
and exponent-minimum implementation agree on all 5759 nonunit divisors.
The cycle is exactly `4,2,5,3,7`; 13 maps into it at 2. All outputs are
at most 7, and 6 maps to 5, so height>=2 means entry into that cycle.
The remaining 5616 bases have height one and end at the terminal 1.
There is no finite maximum defining h_N(1); its terminal status is explicit.

## D Integral edge-sum reconstruction

For a connected finite simple graph G=(V,E), set
`B:Z^V -> Z^E`, `(Bx)_(uv)=x_u+x_v`. Choose a spanning tree and root.
Propagating along its edges gives

```text
x_v=epsilon_v t+c_v,   epsilon_v in {+1,-1}, c_v integral,
```

for integral edge data q and root parameter t. Each remaining edge uv
imposes `(epsilon_u+epsilon_v)t=q_uv-c_u-c_v`. Opposite signs give a
compatibility condition; equal signs fix 2t. Thus:

* In the bipartite case, all compatible q have a free root parameter;
  choosing integral t gives integral labels.
* In the nonbipartite case, a single odd cycle fixes t. Compatibility
  may still hold over Q while t is half-integral.

This is the inherited rational kernel law in the E atlas. Its integral
refinement follows directly from minors. Deleting one root column from
the spanning-tree rows gives determinant +/-1, so the first n-1 Smith
factors are 1. For bipartite G the rational rank is n-1, completing the
list. Otherwise the rank is n. Every n-row minor is even, because mod 2
the all-one vertex vector lies in the kernel. A spanning odd-unicyclic
subgraph has determinant +/-2: peel tree leaves, then expand the odd
cycle determinant. Such a subgraph exists by taking a spanning tree and
an edge whose endpoints have the same tree parity. Therefore the final
Smith factor is exactly 2. QED.

The resulting cokernels are `Z^(m-n+1)` for bipartite G and
`Z^(m-n) plus Z/2` for nonbipartite G. These are classical signed-incidence
facts, with the proof included to fix the integral scope and transpose
convention. See [Zaslavsky, Theorem 4.1 and section 4.6](https://people.math.binghamton.edu/zaslav/Oldcourses/510.S18/mts.pdf).

### D1 Reconstruction is weaker than realization

On a triangle, assigning the square 1 to all three edges gives the unique
rational labels (1/2,1/2,1/2). All rational compatibility conditions hold,
but integer and distinctness conditions fail. By contrast (6,19,30) gives
distinct positive labels and square edge sums 25,36,49.

For bipartite G, the sign change `x_v=epsilon_v f_v` turns B into an
oriented difference map. A tree has no cycle compatibility constraints,
but the path labels (0,1,2) still repeat the absolute difference 1. Thus
unimodularity does not imply gracefulness. The familiar additional
obstruction is the squared-difference discriminant of THM-2761.

The examples K3,3, K5 and Petersen have the Smith data in README. K3,3
is graceful: put labels (0,1,2) on one shore and (3,6,9) on the other;
its nine absolute differences are precisely 1,...,9. It also has the
inherited positive square-sum realization with shores (1,193,673) and
(3,168,483). These are different labellings for different predicates.

Petersen's standard presentation has outer edges i--(i+1), inner edges
(i+5)--((i+2 mod 5)+5), and spokes i--(i+5), for i mod 5. Contract each
spoke. The two cycle step sizes 1 and 2 give every pair of five branch
sets, so the result is K5. Petersen itself has degree three and no
4-cycle, excluding literal K5 and K3,3 subgraphs. Subgraph exclusion
therefore cannot be used as minor exclusion or planarity.

Finally, the same triangle has Laplacian Smith factors (1,3,0), in contrast
to the signless-incidence factors (1,1,2). The graph alone does not select
the relation module: the operator and coefficient ring are essential.

## E Descent and continuation retain different information

For coprime opposite-parity Euclid pairs m>n>0, the Berggren branches are

```text
U(m,n)=(2m-n,m), A(m,n)=(2m+n,m), D(m,n)=(m+2n,n).
```

Their slope images are respectively (1,2), (2,3), and (3,infinity).
The only valid pair at slope 2 is (2,1), while slope 3 has the wrong
parity. Therefore each nonroot pair has a unique inverse branch:

```text
m<2n:        (n,2n-m),
2n<m<3n:     (n,m-2n),
m>3n:        (m-2n,n).
```

Each is again positive, ordered, coprime and of opposite parity; m+n
strictly decreases. Descent terminates at (2,1), and reversing it proves
unique generation. This is an elementary proof of the classical pair
version of the Berggren tree, not a new classification claim.

The determinants of U,A,D are +1,-1,+1. Hence chi is a word character,
analogous to lambda on the commutative monoid of positive integers.
The analogy preserves sign transport only. U followed by A and A followed
by U generally have different endpoints despite the same letter counts.
The parent rule reconstructs both the word and its sign in the check.

The inherited every-modulus witness is the pair
`(16M+2,4M+1)`. It is primitive of opposite parity and its Pythagorean
triple is congruent to (3,4,5) modulo M, while its branch character is -1.
The root has character +1. THM-4808 proves this at every M; the check
independently follows inverse branches for M=1,...,200.

The sign alone still has a two-state streaming update. Residue non-descent
does not prove that the whole branch word is necessary to retain this sign;
the full endpoint distinguishes branch words, a different observation.

For comparison, a strict decreasing height on every edge of the divisor
five-cycle would imply `h(4)>h(2)>h(5)>h(3)>h(7)>h(4)`, impossible.
One can classify that finite dynamics by cycle and basin, but cannot
reuse Berggren's well-founded descent without changing the operation.

### E1 The abstract continuation quotient

Let X carry specified partial operations and an observation o. Add an
invalid state bottom and make every operation total by sending illegal
inputs, and bottom itself, to bottom. Give bottom its own observation.
Define `x~y` iff `o(w x)=o(w y)` for every finite word w of operations.

This relation preserves current observation (empty word) and is stable
under every operation (prepend that operation to each continuation).
Conversely every observation-preserving, operation-stable equivalence is
contained in it, by induction on word length. Thus it is the coarsest
lawful deterministic quotient for the declared observation and operations.
This is the standard continuation/Myhill–Nerode mechanism, including
legality as an observable; it is not automatically a finite quotient.

For multiplication, the repo already classifies the resulting mu and
lambda quotients. For tournament composition, the retained child colours
and ports depend on which substitutions and moves are allowed. For square
sums, neither an unlabelled graph nor the square/nonsquare bit guarantees
the same future attachments. The abstract quotient gives a precise target
for research, but computing it can be as hard as the original problem.

### E2 The inherited even-semiprime ray, with its primality guard

This is the even case of the incoming
[Session84 bijection](../session84_guarded_unification_20261028/BERGGREN_SEMIPRIME_BRIDGE.md),
included with its short proof to fix the map and failure condition.
For an odd prime q, sigma(2q)=3(q+1) is square iff
`q=12k^2-1` for a positive integer k. Indeed a square root z is divisible
by 3; writing z=3u gives q+1=3u^2. Parity of q forces u even, say u=2k.
The converse is immediate provided q is prime.

Then

```text
(q-2,2z,q+4)=(12k^2-3,12k,12k^2+3)=3 Phi(2k,1).
```

The pair (2k,1) is primitive and of opposite parity, and D increments k.
The scale 3 is exact because Phi(2k,1) is primitive. At k=1 this gives
q=11, source 22 and square sigma 36. Since q=-1 modulo 3, multiplication
of 2q by 3 is coprime and multiplies sigma by 4, yielding the source 66
and total 144 at the first point. These candidate-ray steps are not aliquot
iterations. At k=5, q=299=13*23, and actual sigma(598)=3*14*24=1008;
the formal product 3(q+1)=900 no longer computes sigma.
