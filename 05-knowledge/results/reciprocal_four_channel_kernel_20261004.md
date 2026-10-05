# Four reciprocal channels explain the square and its correction

2026-10-04. **PROVED:** the reciprocal pair identity, the determinant-scaled
trace identity and its three-plus-one representation split, the carry loss
of these characters, and a lossless typed trace/anchor decoder for positive
valuation words. **FINITE-EXACT:** the declared controls. Four channels do
not by themselves define a tournament. No Collatz coverage or convergence
claim follows.

Artifacts: [program](../../04-computation/experiments/reciprocal_four_channel_kernel_20261004.py)
and [saved output](reciprocal_four_channel_kernel_20261004.out).

## 1. Recovery and the carrier being squared

[The quadratic-rank note, section6](quadratic_escape_rank_atlas_20261004.md#6-the-first-and-third-quadratics-act-on-repeated-word-presentations)
already proves that doubling an ordered valuation word sends its multiplier
`lambda` to `lambda^2`, while `J=lambda+1/lambda` goes to `J^2-2`.
[The affine-axis note, section3](geometry_collatz_drift_carriers_20261004.md#3-same-axis-means-one-primitive-pattern-different-axes-retain-order)
proves the ordered carry decoder and the axis/multiplier composition law.
The present addition identifies the four-channel source of the correction,
extends it to unequal factors, and states exactly what this compression loses.

Closest proved mechanism: those two notes' word matrices and decoder.
Canonical hostile: words12 and21 have identical diagonal data but different
carries. Corrected near miss: two neutral eigenweights are not two globally
invariant trivial channels in the full tensor representation. The useful
sidecar is the affine anchor, which records the coordinate of the eigenspace.
The live board is **reciprocal weights / ordered pairs / determinant /
representation character / ordered carry / exact source guard**.

The source object is a two-dimensional matrix or a two-element weighted
alphabet. The target is its tensor character. Squaring this character is
an operation on an observable. Squaring a matrix retains dimension two;
taking a tensor square produces dimension four. Neither operation has been
identified with a tournament-size operation. No historical undefined f/g
law is needed here.

## 2. The exact four-channel kernel

For a positive scalar a define `J(a)=a+a^(-1)`. The Cartesian product of
the two weighted alphabets `{a,a^(-1)}` and `{b,b^(-1)}` has four ordered
positions:

| Position | Weight |
|---|---:|
|(+,+)|ab|
|(+,-)|a/b|
|(-,+)|b/a|
|(-,-)|1/(ab)|

Summing the two diagonal and two mixed positions separately proves

    J(a)J(b)=J(ab)+J(a/b).                                  (1)

This identity is valid for any nonzero scalars; positivity gives the exact
equality boundary `J(a/b)>=2`, with equality iff a=b. Thus for a repeated
factor lambda,

    J(lambda)^2=J(lambda^2)+2.                              (2)

The +2 is the combined weight of the two mixed positions, each of weight
one. With unequal positive factors those positions contribute the retained
relative mode `J(a/b)`, not a constant. For example a=3/2 and b=9/8 give
the mixed contribution25/12.

Taking `a=b=sqrt(lambda)` gives a second, compatible identity:

    X=sqrt(lambda)+1/sqrt(lambda),     X^2=J(lambda)+2.       (3)

The two square maps therefore act at different coordinates. If
`Y=J+2=X^2`, word doubling acts by

    lambda -> lambda^2,
    J -> J^2-2,
    X -> X^2-2,
    Y -> (Y-2)^2.                                          (4)

These are linked coordinate rules; (3) does not assert that `J -> J+2`
is another word-doubling operation. Nor does the +2 count two newly added
tournament vertices. The four positions have no specified pairwise
orientation. In the repeated positive case the two mixed weights tie,
so ordering weights alone would not give a tournament anyway.

## 3. The matrix identity and the actual three-plus-one split

Let M be a two-by-two matrix over the rationals or reals, with trace t and
determinant delta. Direct expansion gives

    tr(M tensor M)=t^2,
    tr(M^2)=t^2-2delta,
    tr(M tensor M)=tr(M^2)+2delta.                           (5)

No nonsingularity is needed for (5). If delta>0, normalize
`A=M/sqrt(delta)`, so `det(A)=1`. Its trace X then satisfies

    X^2=tr(A^2)+2.

For a positive triangular carrier with diagonal entries P,Q, put
`lambda=P/Q`. The eigenvalues of A are `sqrt(lambda),1/sqrt(lambda)`;
thus `tr(A^2)=J(lambda)`. The determinant normalization is precisely why
the correction in (3) is two rather than `2PQ`.

**PROVED two-coordinate iteration.** Without normalization, actual matrix
squaring updates trace and determinant together:

    (t,delta) -> (t^2-2delta, delta^2).                     (5a)

Thus the first trace step is `t^2`, `t^2-1`, or `t^2-2` when the initial
determinant is0,1/2, or1, respectively. This provides one matrix framework
for the earlier quadratic trio, with a crucial distinction: only the
determinant fibers0 and1 remain constant under iteration. Indeed
`delta^2=delta` has just those two real solutions. The middle fiber1/2
becomes1/4, then1/16; freezing its correction at one would change the map.
For the exact hostile

    M=(0,-1/2;1,0),

the first three traces of M, M^2, M^4 are `0,-1,1/2`, whereas iteration
of `t^2-1` would give `0,-1,0`. Retaining the determinant is the required
extra coordinate. This does not make three unrelated integer dynamics
conjugate, nor does it supply an actual Collatz-word operation for `t^2-1`.

There is a canonical representation decomposition in characteristic zero:

    V tensor V = Sym^2(V) direct-sum exterior^2(V),          (6)

of dimensions three and one. To prove it, the factor-swap involution
commutes with every `M tensor M`. Its +1 eigenspace has basis
`e00,e01+e10,e11`, and its -1 eigenspace has basis `e01-e10`.
The latter is multiplied by det(M). Consequently

    tr(Sym^2 M)=t^2-delta=tr(M^2)+delta.

For normalized A the two characters are `J+1` and1, totaling `J+2`.
The script checks the full change-of-basis identity, not only its trace.

In a diagonal eigenbasis the four tensor weights are `lambda,1,1,lambda^-1`.
The two weight-one positions are individually neutral for that diagonal
matrix. They are not two canonical trivial summands for all matrices at
once. One lies inside the three-dimensional symmetric representation and
mixes with its other coordinates when the matrix is not diagonal. The
global decomposition is (6), not a representation `A^2 direct-sum1
direct-sum1`: the assignment `A -> A^2` is not multiplicative for
noncommuting matrices.

The elementary word witness is

    M_(1212)=(81,85;0,64),
    M_(1122)=(81,73;0,64).                                 (7)

The first is the square of the carrier for12. The second separately
squares the carriers for1 and2 before composing. Their diagonal characters
agree, but their ordered carries do not.

## 4. Characters forget exactly the coordinate we still need

For a nonempty positive valuation word w, let

    P=3^r, Q=2^A, M_w=(P,B_w;0,Q),
    F_w(n)=(P*n+B_w)/Q,
    alpha_w=B_w/(Q-P).

Here P and Q are unequal by prime factorization. With the shear
`T_alpha=(1,alpha;0,1)`, one has

    M_w=T_alpha diag(P,Q) T_alpha^(-1).                    (8)

Thus conjugacy characters of an individual rational carrier cannot detect
its anchor. The tensor, symmetric and exterior characters in this note all
ignore B. More explicitly, the kth symmetric matrix is triangular with
diagonal entries `P^(k-i)Q^i`, i=0..k, so its trace is their sum. Tensor
products and exterior powers likewise retain only products of diagonal
entries. Even a trace of an ordered product of these upper-triangular
carriers depends only on the products of their diagonals. Changing their
order can change the carry without changing any such character.

This statement concerns the natural matrix representation and its algebraic
tensor constructions over characteristic zero. It is not a statement about
all representations of every finite affine group.

**PROVED positive repair: retain the full marked operator.** The tensor
matrix itself is not carry-blind. In its ordered-pair basis,

    (M tensor M)_(00,00)=P^2,
    (M tensor M)_(00,01)=PB,
    (M tensor M)_(11,11)=Q^2.

The positive square roots recover P,Q, and division by P recovers B. The
remaining entries can then be checked against the reconstructed tensor.
Furthermore

    (M_v tensor M_v)(M_u tensor M_u)
        =(M_v M_u) tensor (M_v M_u).

This is a genuine recursive four-channel carrier that retains order and
carry. Its numerical entries and marked basis contain the information;
it is not a four-valued label. Over unrestricted matrices, M and -M have
the same tensor square; the positive-diagonal type removes that ambiguity.
The loss occurs on passing to characters, diagonal weights, or an
unmarked conjugacy class, not on retaining the full marked tensor.

For example the two chronological compositions of letters1 and2 are
`(9n+5)/8` and `(9n+7)/8`. Composing the first with the inverse of the
second gives the nontrivial formal translation `n -> n-1/4`. Every kth
symmetric character of its two-by-two matrix is k+1, the same as the
identity. This formal group operation is not asserted to be a legal
integer Collatz path.

### A smallest total-cost obstruction to one four-state label

At fixed r and A, there are `binomial(A-1,r-1)` positive valuation words.
The carry decoder below shows that their carries are distinct. Thus a
lossless additional finite label, after r,A are known, needs at least this
many values, or at least `ceil(log2 binomial(A-1,r-1))` bits.

The first total cost for which a fibre exceeds four is A=5, r=3. All six
words have `J=1753/864`:

| Word | Carry B | Anchor alpha | Exact source residue modulo64 |
|---|---:|---:|---:|
|113|19|19/5|55|
|122|23|23/5|43|
|131|31|31/5|19|
|212|29|29/5|57|
|221|37|37/5|33|
|311|49|49/5|61|

For A<=4 every fixed-r fibre has at most three words. The six-word example
therefore gives the stated minimal total-cost test. A single four-state
label cannot retain this fibre alongside r,A. A recursive tree using a
four-symbol alphabet can still encode it and arbitrarily larger fibres:
the retained address then carries unbounded information. This is not an
obstruction to recursive finite-alphabet compression.

## 5. Typed trace plus one rational anchor repairs the formal carrier

**PROVED.** On the type of nonempty positive valuation words, `(J,alpha)`
is a lossless encoding of the marked word. This strengthens the earlier
axis/multiplier statement only by replacing the multiplier with its typed
trace; the carry decoder itself is inherited.

Indeed

    J=(P^2+Q^2)/(PQ)

is reduced because P and Q are coprime. Its denominator is exactly
`3^r*2^A`, recovering r,A and the ordered type `P=3^r,Q=2^A`.
Although an ambient J identifies lambda and its reciprocal, that reciprocal
cannot have this positive-word prime type. Next

    B=(Q-P)alpha.

For r>=2 split off the first letter:

    B=3^(r-1)+2^a_1 B_tail,

where B_tail is positive and odd. Hence
`a_1=v2(B-3^(r-1))`. Divide to recover the tail carry and repeat. The
last carry must be1 and the last exponent is recovered from total cost A.
This both proves injectivity and gives a validator: arbitrary rational
inputs that fail these conditions are rejected by the script.

The repaired encoding composes exactly. If u is followed by v, then

    lambda_uv=lambda_v lambda_u,
    alpha_uv=[lambda_v(1-lambda_u)alpha_u
              +(1-lambda_v)alpha_v]/(1-lambda_v lambda_u).   (9)

The denominator cannot vanish on positive valuation words. Equation (9)
comes from composing `F_w(n)=lambda_w n+(1-lambda_w)alpha_w`; it retains
order. Repeating the same word fixes alpha, as the inherited note proves.
The coefficients in (9) are not generally positive, so it is not a convex
averaging formula.

The word's exact arithmetic source cylinder is also recoverable:

    n=(Q-B)P^(-1) mod2Q.                                   (10)

To see sufficiency, `P*n+B=Q mod2Q` makes the final formal endpoint odd.
The carry recursion successively forces the listed divisions and odd
intermediate states; necessity follows by composing the actual steps.
An equivalent direct induction splits off the first letter and applies
the same congruence to its tail. Positive sources remain positive.

This repairs the formal map and its source guard, not its first-hit or
home status. At27, word12 is actual but1212 is not. At5, word42 is actual
but its first step has already reached1; exporting a first-hit certificate
must remove that root padding. A supplied source and its relevant proof
obligation remain separate from the trace/anchor code.

The anchor is an unbounded rational sidecar, so this is not a constant-bit
storage claim. What is fixed is the number of rational coordinates.

## 6. Exact controls

Run either command; normal and optimized output agree:

```text
python 04-computation/experiments/reciprocal_four_channel_kernel_20261004.py
python -O 04-computation/experiments/reciprocal_four_channel_kernel_20261004.py
```

The script is independent of previous implementations and uses explicit
exceptions instead of assertions that disappear under optimization.

* All23 positive rational values a/b with1<=a,b<=6 give529 ordered kernel
  checks, including the exact equality boundary of the relative mode.
* All81 matrices with entries in{-1,0,1}, including singular matrices,
  satisfy (5) and the full three-plus-one change of basis in (6).
* All625 matrices with entries in{-1,-1/2,0,1/2,1} are squared four times,
  giving2500 exact trace/determinant updates. The determinant1/2 hostile
  checks the first failure of freezing the middle quadratic's correction.
* Every positive valuation word of total cost at most9 is included:
  511 trace/anchor round trips,1022 actual source replays, and3577
  symmetric-character checks at degrees0..6. All511 full marked tensor
  operators independently recover P,Q,B in both list and tuple encodings.
* Every ordered pair of words of total cost at most5 is composed:
  961 matrix, full-tensor, trace and decoded-word controls for (9).
* Hostiles retain the six-word collision, the unequal carries in (7),
  the character-invisible nonzero translation, the illegal repeat at27,
  the overall-sign ambiguity before positive typing, and the first-hit
  failure at5. Six invalid exact trace/anchor inputs and four invalid
  marked tensors are rejected.

The useful bridge is therefore an exact character identity with an explicit
repair for its lost coordinate. It does not identify four weighted positions
with four tournament isomorphism classes or transport an unproved coverage
claim through the quotient.
