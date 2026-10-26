# A common method: certify what an observation can recover

**PROVED elementary observation criteria; SYNTHESIS of the linked audited
examples.** These criteria are not new general theorems. The contribution
is to identify the exact observer, repair and remaining obligation in each
lane, including cases where an apparent loss is inferable after all.

## 1. Recovery is a factorization question

Let X be a specified source class, O:X→Y an observation, and f:X→Z the
target. A function g on O(X) with f=g∘O exists if and only if

    O(x)=O(x') implies f(x)=f(x').

Necessity follows by applying g; for sufficiency define g(y) to be f(x)
for any x in the fibre. This criterion does not imply an effective,
efficient, continuous or algebraic decoder without further hypotheses.
The source class X and the image O(X) are part of the statement.

Consequently a collision is useful only when it changes the actual target.
Equal conic endpoints with different branch signs obstruct recovering the
sign; equal radical with different signs obstructs a different decoder.
Equal path counts without equal tournaments obstruct classification by H.
None of these automatically obstructs another target or another observer.

There is a basic information bound. For a fibre over y, let
`m(y)=|f(O^-1(y))|`. A sidecar taking at most K values can recover f only
if every m(y)<=K. Conversely, if every fibre has at most K target values,
one can index those values within each fibre to obtain such a sidecar
(as a set-theoretic assertion, with no effectiveness promise).

This is not an algorithmic lower bound. A one-bit target can always be
supplied as an oracle bit. The prime-sieve obstruction concerns a specified
family of finite congruence observations; it does not show that primality
needs unbounded memory in every algorithm. In the conic example, however,
the one bit is supplied by an explicit composable character of the branch
matrices, and the two states over each pole prove that bit necessary for
recovering the oriented tournament from its pole.

## 2. Continuation requires the operations and their guards

For an observation O and named partial operations, introduce an absorbing
failure state with a distinct observable value. Define the future profile

    P_O(x)=(O(w x)) over all finite named operation words w,

including the empty word. Equality of these profiles is the coarsest
operation-stable equivalence refining equality of O. Applying one operation
and then a word is another word, proving stability. Any stable refinement
must preserve every word by induction, proving minimality. This is the
inherited operation-kernel criterion of [THM-840 / hamming-five-continuation](../../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md)
and [strong-atom C5](../strong_atom_continuation_vitali_20261026/README.md).

On a finite universe, repeated partition refinement terminates. On an
unbounded universe, neither termination at one fixed depth nor a bounded
number of retained bits follows. Changing named operations to freshly
sampled or independently relabelled operations changes P_O itself; it is
not merely an optimization of the same theorem.

The path reconstruction now makes one particular all-order statement
possible: the one-step legal permutation guards already distinguish all
fixed-path tournament presentations. Hence the closed state count on X_n
is exactly `2^binom(n−1,2)`. This is a state-distinction theorem, not a
lower bound on the time required by every tournament algorithm.

## 3. Three different repair mechanisms

**Structural inference.** The observer U(T) initially appears to forget
orientations. The separator theorem and deletion-reachability lemma recover
them from U(T) itself. Its fibres on the tournament source class are
singletons. Adding an orientation sidecar would store inferable data.
This depends on the source promise: an arbitrary oriented graph need not be
the complete Hamilton-path support of a tournament. The procedure and
proof do not silently assume that recognition problem solved.

**Finite state completion.** A conic point has two possible transported
orientations. The determinant character gives their explicit two-sheet
cover, with all sixteen states reachable. It is stable under all three
Berggren branches and sufficient for this orientation target. It does not
recover the integer triangle, its primality, its abc height, or the joint
field degree retained in the separate strong-atom problem.

**Height-dependent exhaustion.** The corrected prime sieve provides

    for every H, a B(H), for every 3<=x<=H:
        accepted_at_B(H)(x) iff x,x+2 are prime.

One may take B(H)=floor(sqrt(H+2)). The prime-equality exceptions are
essential to this statement. A fixed stage has composite false positives,
and its periodic survivor count alone does not supply an unbounded prime
family. The quantifiers for twin infinitude would instead require

    for every H, some x>H that passes every corrected prime test.

These assertions are different. The coherent class−1 gives an explicit
limit point with no positive bounded lift. The exact finite recurrence
for pi_2 does not by itself prove that its values are unbounded.

## 4. Weights, order and image membership remain separate

The F1 language makes a second distinction explicit. Passing from branch
words to a commutative monoid preserves the determinant character but loses
order; augmentation then loses the character too. For prime-divisor
polynomials, augmentation counts the proper-divisor monomials, whereas
evaluation at the primes computes their weighted sum. Neither scalar generally
recovers the shape or prime labels. Their kernels differ, as the explicit
aliquot collision shows. Calling all of them an underlying geometry would
erase the actual maps that explain the successes and failures.

Likewise, a bijection is not an existence proof at unbounded height. The
twin-radius triangle map retains primality on its marked image; counting
all primitive triangles instead changes that image. The faithful map
T→U(T) retains H and Aut, so it transports the existing7/21 exclusions;
it does not make an arbitrary partial digraph admissible or realize a new
H value. The field cover supplies actual Paley tournaments, all with the
same H189, rather than a variable family realizing the full spectrum.

## 5. A research procedure supported by these examples

1. Specify the source class, target and operations before choosing a model.
2. Recover the closest exact map and its known hostile example.
3. Test whether apparent missing data are inferable from structure before
   storing them or declaring a no-go theorem.
4. If recovery fails, exhibit a target-changing collision and determine
   whether a composable finite state, a retained weight/label, or a growing
   horizon repairs this particular loss.
5. Recheck the repaired map under every permitted operation and at its
   first boundary; the next prime square and the five-vertex unused arc
   are examples of boundaries that expose the mechanism.
6. State the remaining image, uniformity, positivity or termination
   obligation separately from the encoding theorem.

This is a cohesive account of the audited examples, not an empirical claim
that every historical repository advance has been inspected or is correct.
Shared logical forms or repeated small integers motivate tests; an asserted
equivalence between open problems still needs actual reductions preserving
their predicates and quantifiers. The connections proved here are valuable
precisely because those maps and their failure boundaries are explicit.
