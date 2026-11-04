# Arbitrary overlapping intervals: exact scalar classification

Session 188. **Proved:** simultaneous XOR reversals of intervals in a transitive
block quotient have an exact scalar-preservation criterion, and their quotient
action is universal. As a consequence, a named transitive module cannot be
replaced by any different tournament on the same names invisibly to the full
two-port response, in any fixed ambient quotient.

The distinction between simultaneous original-module XOR and a sequence of
reversals of currently valid modules is essential. The former can create a
directed triangle from a transitive triple. No assertion below silently
reclassifies such an operation as a legal sequential module reversal.

## 1. Setup

Let A_1,...,A_k be nonempty tournaments on disjoint fixed names, and let

    R=A_1 -> A_2 -> ... -> A_k.

Substitute R at a fixed root z of any tournament Q. Every union of consecutive
blocks I_[a,b]=A_a union ... union A_b is a module of Q[R]. Choose an arbitrary
family of intervals a<b, each used once modulo two, and simultaneously toggle
every edge internal to each chosen original interval. Denote the result by T'.

There is a quotient tournament K on the k block names and bits epsilon_i such
that the new root module is

    R'=K[A_1^(epsilon_1),...,A_k^(epsilon_k)],              (1)

where exponent one means whole converse of that block. Indeed each internal
edge in A_i is toggled the same number of times, and each cross-block edge
between i and j is toggled according to the intervals containing both.

Singleton block intervals, if also allowed, change only epsilon_i. They do
not affect the quotient classification below, but they are excluded from the
count of interval patterns because they introduce extra internal choices.

For any tournament B, let q_h(B) count ordered spanning covers by h nonempty
directed paths. In particular q_1(B)=H(B). Equivalently

    q_h(B)=sum_(ordered nonempty partition S_1,...,S_h)
                       product_j H(B[S_j]).                (2)

## 2. Scalar classification, uniformly in the ambient quotient

**Theorem 1.** For every h,

    q_h(R') >= q_h(R).

Equality holds simultaneously for all h if and only if K is transitive. If K
is nontransitive, then

    H(Q[R'])-H(Q[R]) >= 2 H(Q) product_i H(A_i) > 0.        (3)

Hence scalar equality of the ambient Hamilton counts is equivalent to
transitivity of the new block quotient K.

**Proof.** Fix any named subset S of the root module. Originally R[S] is an
ordered join of its nonempty pieces A_i[S], so

    H(R[S])=product_(nonempty pieces) H(A_i[S]).

In R'[S], each piece is unchanged or conversed, which preserves its Hamilton
count. Choose any Hamilton path in the induced quotient on the occupied block
names, and concatenate one Hamilton path in each piece in that order. This
gives an injection from the product counted above into Hamilton paths of
R'[S]. Thus H(R'[S])>=H(R[S]) for every S. Applying this separately to every
part of every ordered partition in (2) proves coefficientwise domination.

If K is transitive, R' is an ordered join of the same blocks, permuted and
possibly conversed. Every induced subset then has exactly the same product
Hamilton count as before; hence all q_h agree.

If K is nontransitive, it has at least three Hamilton paths. One elementary
proof starts with a directed triangle's three paths and inserts each remaining
vertex into each path by any fixed deterministic valid insertion rule.
Deleting the inserted vertex recovers the original path, so this preserves
distinctness. Concatenating full block paths in each of those quotient orders
gives H(R')>=3 product_i H(A_i), while H(R)=product_i H(A_i).

Finally the occupied-run expansion is

    H(Q[B])=sum_h c_h(Q,z) q_h(B),

with nonnegative integer coefficients and c_1(Q,z)=H(Q)>0. Compressing the
h nonempty B-runs in a Hamilton path gives this bijection. Combine the
coefficientwise domination with the strict h=1 bound to obtain (3). QED.

The scalar query is an actual full-response entry: when a=c=0 both ports
precede every old vertex in the order x->y, so F_T(0,0)=H(T).

## 3. Every quotient tournament occurs, with a unique interval pattern

Write eta_ab for the bit selecting interval [a,b], a<b. For i<j let delta_ij
be the bit saying that the original arc i->j is reversed. Then

    delta_ij = sum_(a<=i, b>=j) eta_ab modulo 2.            (4)

Extend delta by zero when its first index is 0 or its second index is k+1
(the displayed indices run from 1 to k). Two finite differences invert (4):

    eta_ab = delta_ab + delta_(a-1,b)
                        + delta_(a,b+1) + delta_(a-1,b+1) mod 2.  (5)

All four entries in (5) have valid increasing indices except these explicitly
zero boundary entries. Summing the finite differences telescopes back to (4).

Thus the 2^(binomial(k,2)) interval parity patterns are in bijection with all
named quotient tournaments K. Exactly k! patterns have transitive K, and
exactly those preserve the scalar Hamilton count in every fixed ambient
quotient and for every choice of the nonempty blocks.

This universal quotient action is why one cannot apply the scalar invariance
of sequential module converses to arbitrary overlapping original intervals.
The nontransitive quotient cases are nevertheless separated immediately by
Theorem 1. The scalar-preserving cases can still have different full named
responses; scalar equality is not asserted to be a collision.

## 4. Complete no-collision result for transitive source modules

**Theorem 2.** Let M be a fixed physical set, pi a linear order on M, and U any
tournament on M different from TT_pi. For every fixed rooted quotient Q,

    F_(Q[TT_pi]) != F_(Q[U])

as full named two-port functions.

**Proof.** If U is nontransitive, then H(U)>=3. Every induced nonempty subset
of U has at least one Hamilton path, whereas every induced subset of TT_pi
has exactly one. Equation (2) gives q_h(U)>=q_h(TT_pi) for all h, with strict
difference at least two at h=1. The occupied-run expansion therefore gives

    F_(Q[U])(0,0)-F_(Q[TT_pi])(0,0) >= 2H(Q)>0.

If U is transitive, it is a different named linear order. S187's universal
named transitive-order theorem supplies a physical two-port query that
separates the two orders in every Q. Its proof uses coefficientwise path-cover
monotonicity of the rank queries and c_1=H(Q)>0; it does not assume a common
homogeneous pair or an ordered ambient context. These two cases exhaust U.
QED.

Consequently every nonempty simultaneous XOR family of nontrivial original
intervals of a transitive module changes its full named response in every
ambient quotient. Formula (5) ensures that such a family changes at least one
core edge. This covers overlapping and nested intervals of every size, not
only two crossing intervals or disjoint reversals.

It also covers arbitrary replacement of that transitive module, not only
replacement obtained by reversals. It does **not** cover arbitrary nested
reversals inside a nontransitive source module. In Theorem 1 with nontransitive
blocks A_i, the scalar-preserving transitive quotient permutations and child
converses remain a separate full-response problem.

## 5. Full coefficient controls and provenance

No proof above drops cubic or quartic coordinates. Strictness at one physical
query is enough to exclude equality of the full function. For concrete
controls the verifier retains the complete Boolean difference polynomial for
both a crossing-interval example and a nested-interval example inside a strong
ambient quotient. Each whole-converse control is also retained separately.
Unlisted coefficients are exactly zero in the full Boolean basis, not omitted
unknown coordinates. Arbitrary contrasts are never assumed cubic.

The inherited transitive-order result was checked directly at
[S187 collision/PROOF.md](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session187_named_response_20261104/collision/PROOF.md),
blob f8136d492a089b63415e098dcfee9812cd4b1302. The occupied-run identity is the
same positive path-cover substitution used there. The present additions are
the arbitrary-interval scalar classification and the extension from named
transitive reorderings to arbitrary same-name replacements. Historical
novelty is not asserted.

Run `python verify.py` and `python -O verify.py`. The checker uses explicit
exception gates and standard-library arithmetic. Its finite tests corroborate
the implementation; the path-cover injections establish the unbounded claims.
