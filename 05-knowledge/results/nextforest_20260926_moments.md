# Finite jets, exact rational observers, and the carry forest

**Status: PROVED elementary observer and obstruction results; FINITE-EXACT controls.**
No Collatz convergence theorem or literature-priority claim. This note continues
[the ordered forest](forest_20260926_excursions.md) and
[the Bernoulli carry graph](forest_20260926_bernoulli.md).

## 1. Inheritance and recovered connections

Anchor: how much of an actual integer survives an observer of its ordered
carry polynomial. Niche: the repo's nonlinear GMC channel recovery and complete
Hermite banks. Wildcard: the Prouhet product as a source-preserving hostile.
The board is exact source, observer point, finite jet, sparsity, block separation,
and bounded-correction rank.

Closest proved mechanism: the forest identity
`z^L D_y=3^a D_x+B_w-(2-z)A_I`, with exact nonnegative carry polynomial A_I.
Canonical hostile: two distinct arithmetic carriers can have the same observer.
Corrected near miss: neither a finite number of scalar features nor a finite
jet order alone implies an unbounded fixed fibre. Least-used sidecar: the
observer's exact denominator, or an adaptive sparsity certificate.

Three older repo mechanisms are relevant, with their scopes intact:

* [THM-2631 / homogeneous Wick-channel decoder](../../01-canon/theorems/THM-2631-homogeneous-wick-channel-linear-decoder-and-private-row-no-go.md)
  forbids a particular **linear** channel decoder; its proof does not forbid
  nonlinear recovery.
* [THM-2639 / equal-mass two-rung certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md)
  uses a nonlinear ideal identity to recover nonvanishing on a collided face.
  Its free-semigroup and coefficient-torus hypotheses do not transfer to Collatz.
* [THM-4443 / arbitrary-jet precision](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md)
  treats a complete Hermite observer on a **degree-bounded polynomial space**.
  Its exact inverse needs the complete bank and unit data. Unbounded digit
  degree falls outside that theorem's fixed bank.

The legitimate transfer is to specify the observer, its admissible source set,
its nonlinear reconstruction rule, and its precision cost. No GMC or Hermite
noncancellation theorem is imported into Collatz.

Write `D_n(z)=sum_s b_s(n)z^s` and
`C_n(z)=sum_(s>=1)c_s(n)z^s`, where
`c_0=1`, `c_(s+1)=floor((c_s+3b_s)/2)`.
A jet of order r here means all Hasse derivatives of orders 0 through r-1
at z=1. Equality of these jets is equivalent to equality of the corresponding
first r ordinary power moments of the bit/carry positions.

## 2. Prouhet twins along actual finite Collatz prefixes

**Theorem P.** Fix an odd positive source a, a finite shortcut horizon H,
a jet order r>=1, finitely many primes p, and finitely many integer polynomials
P_i that are nonzero at every `a_t=T^t(a)`, 0<=t<=H. There are pairs of positive
odd integers `(n_E(L),n_O(L))`, for arbitrarily large integer L, such that:

1. both sources have the same first H parity steps as a;
2. for every t<=H, the r-jets of both D and C at their actual time-t states agree;
3. every selected `v_p(P_i(T^t(n)))` equals `v_p(P_i(a_t))` for both sources;
4. `max(n_E,n_O)/min(n_E,n_O)` tends to infinity with L.

**Construction and proof.** Choose Q>a, divisible by a sufficiently high power
of 2 to preserve H parity steps, and by sufficiently high powers of every
selected prime. More precisely, take

`v_2(Q)>=H+1+max_(i,t) v_2(P_i(a_t))`

when 2 is selected, and at least H+1 otherwise. At odd selected p take
`v_p(Q)>max_(i,t) v_p(P_i(a_t))`.
Let e_t be the number of odd steps before t and set

`Q_t=3^(e_t) Q / 2^t`.

This is an integer and retains strictly more p-adic divisibility than the
selected polynomial values at a_t. For any integer Z, the source `a+QZ`
has that parity prefix and its t-th state is `a_t+Q_t Z`. Polynomial
congruence then preserves each specified valuation exactly.

Partition j=0,...,2^r-1 by parity of its binary digit count. Write the two
sets as E and O, and use

`n_E(L)=a+Q sum_(j in E) 2^(L(j+1))`,
`n_O(L)=a+Q sum_(j in O) 2^(L(j+1))`.

Take L larger than the binary lengths of every `3a_t+1` and every `3Q_t`.
All the blocks are disjoint, including their multiplication-by-three carries.
Let `H_q(z)` be the carry polynomial for multiplication 3q with initial
carry zero, as opposed to the initial carry one defining C_q. Then

`D_(T^t n_E)-D_(T^t n_O) = D_(Q_t)(z) z^L P_r(z^L)`,
`C_(T^t n_E)-C_(T^t n_O) = H_(Q_t)(z) z^L P_r(z^L)`,

where the exact product identity is

`P_r(X)=product_(i=0)^(r-1)(1-X^(2^i))`
`      =sum_(j=0)^(2^r-1) (-1)^(s_2(j)) X^j`.

The product is divisible by `(X-1)^r`, proving both jet equalities. This is
an identity of actual binary blocks and carry streams, not a formal pair
of polynomials without integer realizations.

The largest block index `2^r-1` and the next index `2^r-2` have opposite
parity. Thus one source's leading block is L bits above the other's.
For all sufficiently large L their ratio exceeds `2^(L-1)` and is
asymptotic to `2^L`. This proves (4).

**Exactly what fails.** No function of this fixed jet/valuation bank can
approximate log(n) with one uniform bounded additive error on all positive
sources, or n with one uniform multiplicative factor. It would assign the
two members the same estimate, while their logarithmic gap is unbounded.
The theorem also applies to a fixed finite transcript of these observers,
not just to one time slice.

**What survives.** When r>=2, the first digit-position moment
`M(n)=D'_n(1)` bounds the highest set-bit position:
`floor(log_2 n)<=M(n)` for n>=2. Therefore each exact fixed M-fibre is finite.
The shared moment values in Theorem P vary with L. The theorem does **not**
produce unbounded integers in one fixed moment fibre, and does not refute
arbitrary nonlinear ranks of the moments or an adaptive jet order.
Nor are the twin sources two points on one orbit, and finite H never
becomes an infinite common orbit by this construction.

Connection ledger: source is an actual finite orbit pair; target is the
joint jet transcript; map is evaluation of Hasse derivatives at 1;
preserved predicates are the declared valuations, prefix and moments;
destroyed information is the separated block allocation; restoration
requires more jets, exact digit data, or another injective observer. The
cheapest hostile is r=1 with one high block on each side; higher r uses the
same product rather than a new probabilistic assumption.

## 3. A single exact rational observer is already lossless

**Theorem R.** For every positive rational z=p/q in lowest terms with z!=1,
the map `n -> D_n(z)` is injective on nonnegative integers. If q>1 and n>0,
the reduced denominator is exactly

`denominator(D_n(p/q)) = q^floor(log_2 n)`.

For the denominator statement put h=floor(log_2 n). Multiplying by q^h
gives an integer congruent to p^h modulo every prime divisor of q; hence
it is coprime to q. For injectivity, a collision gives a nonzero polynomial
with coefficients in {-1,0,1} and leading coefficient +/-1 having root p/q.
The rational-root argument forces q=1. If p>=2, the leading monomial
strictly dominates the sum of all smaller ones, so it cannot vanish.
(The excluded p=q=1 is digit count and is highly noninjective.)

There is also an explicit decoder when p>=2. For an actual digit value x,
its denominator is coprime to p and

`b_0 = numerator(x) * denominator(x)^(-1) mod p`

is 0 or 1. Replace `x` by `(q/p)(x-b_0)` and repeat. This successively
recovers every bit and terminates. For p=1, q>=2, the denominator gives
h and the numerator is a base-q digit string with the binary bits reversed.

Thus z=3/2 lies inside the forest's positive-carry regime z<2 while
retaining the exact source. Its denominator already encodes binary height.
An inequality that treats this value only as a real magnitude may discard
that arithmetic precision, even though the exact value is lossless.
There is no free compression: the denominator has h+1 binary digits at
z=3/2, and a decoder still reconstructs the original integer.

Observer point matters. At z=2 the value is n itself. At z=1 finite jets
have Theorem P collisions. At the golden ratio, `D_4(z)=z^2=1+z=D_3(z)`.
A dimension argument on the full polynomial vector space cannot decide
injectivity on the constrained {0,1}-coefficient source set.

Connection ledger: target is one exact rational, map is evaluation,
preserved predicate is the complete integer source, lost data is none,
sidecar is exact numerator/denominator arithmetic; rounding destroys the
conclusion. This is a positive control against overextending any finite-
feature obstruction.

## 4. A nonlinear adaptive moment bank also reconstructs the source

**Theorem N.** If `S=s_2(n)`, then S and the S power sums
`p_j=sum_(b_s=1) s^j`, 1<=j<=S, determine n uniquely.

Let the unknown set-bit positions be s_1,...,s_S. Newton's recursion

`e_0=1`,
`k e_k=sum_(j=1)^k (-1)^(j-1) e_(k-j) p_j`

recovers the coefficients of
`Q_n(X)=X^S-e_1 X^(S-1)+...+(-1)^S e_S`.
Its roots, with multiplicity, are precisely the positions, so n is the
sum of their powers of two. Equivalently an order S+1 Hasse jet at 1
suffices. The same statement holds for C_n using total carry K and K
moments, treating a coefficient 2 as a repeated position.

This is the valid nonlinear analogue of the older GMC lesson: collided
linear observations need not prevent nonlinear recovery after the right
structural hypothesis is retained. Here that hypothesis is the number
of atoms, and the operation is Newton reconstruction; it is not the
GMC two-ray Bezout identity. The observer size grows with sparsity.

An **OPEN** constructive target is to propagate the support polynomial
Q_n, or a smaller certified nonlinear invariant of it, across forest
resets with a quantitative contraction. Reconstructing n and then
checking its orbit is only a recoding. A useful result must bound a
coefficient, root geometry, or carry interaction without presupposing
termination. Prouhet twins are required hostile controls for fixed-order
truncations; rational evaluations are required positive controls for
claims about finite observers.

## 4b. Exact nonlinear support-polynomial transport

The adaptive reconstruction has a concrete transport law. Define

`Q_n(X)=product_(b_s(n)=1)(X-s)`,
`R_n(X)=product_(s>=1)(X-s)^(c_s(n))`.

For odd n, put `m=U(n)` and `a=v_2(3n+1)`. Then

`Q_n(X)^3 * X * R_n(X) = Q_m(X-a) * R_n(X+1)^2`.       (N1)

**Proof.** The exact digit/carry relation is

`3D_n(z)+1+C_n(z)=z^a D_m(z)+2C_n(z)/z`.

Every coefficient on both sides is a nonnegative integer. Regard a
coefficient k at z^s as k copies of the position s, and send that multiset
to the monic polynomial having those roots. Addition becomes multiplication;
shifting positions by a becomes replacing X by X-a. This gives (N1)
coefficient by coefficient, preserving every root and its multiplicity.
For an even shortcut step the corresponding identity is simply
`Q_(Tn)(X)=Q_n(X+1)`.

Degrees in (N1) give the inherited digit budget
`S(m)=3S(n)+1-K(n)`. Taking the logarithmic derivative as a formal rational
function gives a different exact observer:

`3 Q'_n/Q_n + 1/X + R'_n/R_n`
`    = Q'_m(X-a)/Q_m(X-a) + 2 R'_n(X+1)/R_n(X+1)`.

Expansion at infinity recovers every positional power-moment identity.
The poles and their multiplicities retain the positions that finite jets
can forget. There is no analytic claim at a root; these are identities
of polynomials and rational functions.

For example n=m=1 and a=2 give
`Q_1=X`, `R_1=(X-1)^2(X-2)`; both sides of (N1) are
`X^4 (X-1)^2 (X-2)`. This terminating control prevents mistaking a
nontrivial support factorization for strict progress.

The map connects the nonlinear Newton reconstruction to the exact forest
carrier: it trades spatial coefficient transport for integer-root
multiplicity transport. It preserves the full source and carry streams,
not a new arithmetic invariant proven to decrease. A proposed polynomial
height, resultant, root-gap, or logarithmic-derivative inequality must
respect the translated target and repeated carry roots. Forgetting those
shifts would lose the very endpoint information being sought.
A fixed observer `log|Q_n(xi)|`, at a nonintegral real xi, is still
`sum_s b_s log|xi-s|` and is excluded by the
[additive-digit circuit theorem](nextforest_20260926_circuit.md); the
new rank would need interactions, a source-dependent observer, additional
carry data, or selected returns.

## 5. A bounded-correction obstruction to polynomial position weights

The previous forest hostile handled nonnegative combinations of D_n(z)
for z>1. Signed finite combinations do not rescue a lower-bounded rank.
After combining repeated evaluation points, let the largest z with
nonzero coefficient have coefficient a. If a<0, the odd spikes
`n=2^h+1` make the potential tend to minus infinity, even after adding
`c log n`, a bounded correction, and an arbitrary function of digit count
(the latter is constant on these spikes). If a>0, the forest's padded
51->77 family has a positive increment tending to infinity dominated
by that largest z, while digit count is the same at the endpoints.
Consequently no nonzero finite signed combination of evaluations z>1,
plus those terms, is both lower bounded and globally nonincreasing.
This conclusion concerns that specified form of rank, not general
nonlinear functions of the injective rational evaluations.

There is a direct polynomial-position counterpart with a different hostile.
For R>=1 define

`n_R=51+sum_(j=1)^R 150*2^(10j)`,
`m_R=77+sum_(j=1)^R 225*2^(10j)`.

Since `3*150/2=225` and all blocks are disjoint,
`U(n_R)=T(n_R)=m_R>n_R`. Their digit counts both equal `4+4R`, while

`M(m_R)-M(n_R)=1+4R`.

For any real polynomial weight w(s) of degree d>=1 with positive leading
coefficient a_d, set `W(n)=sum_(b_s(n)=1) w(s)`. The high block increment
at position t is

`delta_w(t)=w(t)+w(t+5)+w(t+6)-w(t+1)-w(t+2)-w(t+4)`.

It has degree d-1 and leading coefficient `4d a_d`. Summing t=10j gives

`W(m_R)-W(n_R) = 4 a_d 10^(d-1) R^d + O(R^(d-1))`.

The low block contributes only a constant. Therefore
`W(n)+c log n+f(s_2(n))+g(n)`, with arbitrary real c, arbitrary f, and
bounded g, cannot be nonincreasing on every sufficiently large odd edge.
The logarithmic increment is bounded and f cancels; W's increment diverges.
This includes an arbitrarily signed lower-degree polynomial beneath a
positive leading term. It does not exclude nonlinear moment ranks,
selected-return inequalities, or unbounded corrections retaining more data.

## 6. Reproduction and next use

Run:

```text
python -X utf8 -B 04-computation/experiments/nextforest_20260926_moments.py
python -X utf8 -B -O 04-computation/experiments/nextforest_20260926_moments.py
```

The [script](../../04-computation/experiments/nextforest_20260926_moments.py)
uses explicit checks, no assertions and no imported producer code. Its universe
is every n=1,...,1023 for rational decoding and Newton reconstruction;
60 Prouhet constructions with horizons 0,3,8,16, orders 1,...,5 and three
guard spacings each; actual sources up to 2252 bits; and 80 growing
block-family edges. It checks 465 actual paired time slices, full literal
carry-block recomposition, all 2/3/5 valuations of x,x+1,x^2+1, and all
claimed jets. Dense/sparse source rows, the golden-ratio collision and
signed lower-order polynomial weights are explicit controls.
It also expands all coefficients of (N1) for every odd n<1024 and the
even translation identity for every positive even n<1024. Both normal
and optimized runs pass 179980 gates with identical output.
The [frozen output](nextforest_20260926_moments.out) records the finite counts.
The unbounded statements above rest on their proofs, not those counts.
