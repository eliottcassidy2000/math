# Independent review: all-greedy endpoint descent

19 October 2026.

## Verdict and publication disposition

**PASS for both stated descent theorems**, conditional on the previously
audited matched-block height theorem and its cited explicit logarithmic-form
input. The reverse-injective code lemma, finite repeated-factor argument,
all-m V/F result, and arbitrary fixed-template U result for m>=2 are valid.
No mathematical theorem repair is required.

The publication includes the two scope/wording clarifications supplied in the reviewed
overlay. These preserve every formula and proved theorem:

1. The all-greedy theorem closes the endpoint-descent goal behind the earlier
   typed-left-special question. Its guaranteed pair is source-anchored OR
   has different incoming types. It does not establish the original stricter
   occurrence hypothesis requiring different incoming types in every large
   prefix.
2. In the general-template m=2 proof, the second matching U begins at y_1,
   immediately after iteration 1 and before iteration 2. It does not begin
   at y_2 after both iterations.

The frozen source is unchanged. The complete reviewed text is
ALL_GREEDY_ENDPOINT_DESCENT.scope-reviewed.md, accompanied by the exact
SCOPE_ONLY_OVERLAY.patch and SCOPE_OVERLAY.json.

    Frozen source SHA-256:
    d4fed9827b5de52240486490389cbdd299d90fb830eb64a3ad707418a1b20558

    Scope-reviewed text SHA-256:
    a65677f28afb8a736f99985375d7e8e44740c7949386906dc5b82ebaba4f1b28

## 1. Frozen material and executable checks

All hashes and sizes in the source manifest verify, as does SHA256SUMS. The
assigned submitted verifier has hash
054389e778111ceb97577d0b3cb5c0d4902603ddaa6d8a4440e38b6bf241208d.
Both submitted ordinary and optimized runs reproduce the supplied receipt
byte for byte: **40,345 explicit checks**.

The independent verifier imports no submitted implementation. Its ordinary
and optimized outputs agree byte for byte: **156,750 explicit checks**.
Its independent controls include:

- All 3,138 nonempty binary words through length 10 and ternary words through
  length 6, at every factor length.
- All 67,284 repeated-factor pairs in those controls, including 9,552
  overlapping pairs, with 38,136 constructive left shifts checked.
- 107 exact lawful whole-word controls over 34 fixed expansion templates,
  including all K/V templates of macro-length at most four and larger
  templates K^32, V^20, K^128 and (KV)^7.
- Multiple-reset cases with floor(theta_U)=35, hence potential code choices
  of 35 and 36 conditional resets, as well as zero/one-reset cases.
- Thirteen complete literal Terras orbit replays, macrobranch guards for
  every whole-word control, exact checkpoint types, source anchors,
  pairwise distinctness and direct endpoint descent.
- Independent logarithmic enclosures using 192 terms for log 2 and
  log(3/2), to cross-check all mechanical floor counts.
- Synthetic affine countercontrols showing why separated output images
  and the noncrossing source anchor cannot be omitted from the abstract
  reverse-decoding argument.

The universal results rest on the proofs below. A finite verification count
does not supply an all-m or all-U quantifier. Both submitted and independent
checks use explicit exceptions, so optimization does not remove them.

## 2. The arithmetic image and injectivity hypotheses

The prior frozen image-guard packet has separately passed independent review.
Its exact relevant identities are also reconstructed independently here.
For a positive lawful input, V outputs have finite type v_3(y-1)=4 and
residue 82 modulo 243; VF outputs have finite type 5 and residue 1 modulo
243. These are disjoint output images, independent of the entering integer.

Each fixed chronological code word has an affine branch
(3^h x+b)/2^n with strictly positive slope. Thus it is injective on its
lawful domain, even though the whole Collatz map is not globally injective.
The argument only cancels a branch after identifying its code symbol from
its actual output-image class.

For a general nonempty expansion template U over K,V, its terminal letter
fixes its output type t: 3 if it ends in K and 4 if it ends in V. This does
not depend on the earlier letters or input. Since

    F^j(z)-1=(3/4)^j(z-1),

the image of UF^j has type t+j. At the coefficient level the independent
implementation checks the exact identity

    b_(UF^j)-M_(UF^j)=3^j(b_U-M_U).

Since h_U>t, the term containing the entering integer cannot cancel this
valuation. In particular the two branches UF^a and UF^(a+1) have disjoint
finite types t+a and t+a+1 for arbitrary a>=0. No bounded-type assumption
is present. The exceptional fixed point F(1)=1 cannot intervene because
U has a finite output type and hence never outputs 1.

Literal prefix-freeness is not required: V is itself a literal prefix of
VF. The proof keeps the actual specified code parse and uses image-class
separation, rather than attempting to infer code boundaries from the raw
parity string.

## 3. Reverse-injective checkpoint lemma

Write y_0=x and y_i for the actual state after i complete code images.
Every such completed prefix is before the extra final coefficient crossing,
and has cumulative multiplier S_i>=1 and positive carry e_i. Therefore

    y_i=S_i x+e_i>x for every i>=1.

Suppose y_i=y_j for 0<=i<j<=m. If i=0, this contradicts the displayed strict
inequality. If i>0, disjoint output images identify the same preceding code
symbol on both sides. Its affine injectivity gives
y_(i-1)=y_(j-1). Repeating exactly i times gives y_0=y_(j-i), again impossible.

All cancellations use identical actual lawful branches. All preceding
states exist within the given finite trajectory. The proof neither
reorders iterations nor removes an unmatched reset, and no inverse for the
global Collatz map is assumed. The conclusion is that all completed
iteration checkpoints are pairwise distinct, independently of whether the
final endpoint descends.

The synthetic controls illustrate both essential hypotheses:

- E(x)=2x+1 and C(x)=(x+3)/2 are injective affine maps with positive carry.
  Starting at 1, E gives 3 and C fixes 3; the cumulative slopes 2 and 1
  both satisfy the noncrossing bound. Their output images overlap, so an
  internal repeat is possible without image separation.
- D(x)=(x+2)/2 on inputs 2 modulo 4 has even output image, disjoint from
  E's odd image, but D fixes 2. Its slope is below one, so it violates the
  source-prefix multiplier anchor. Separated images alone are insufficient.

These are abstract affine countermodels, not claimed Collatz bank words.
The actual greedy families satisfy both required hypotheses.

## 4. Minimal-left-occurrence alternative

For a repeated length-k symbol factor with zero-based starts i<j, if i>0
and the immediately preceding symbols agree, shifting both starts one place
left preserves a length-k match. The first symbols agree by the predecessor
assumption; the next k-1 agree by the old match. The shifted factors remain
inside the same finite word, including when the occurrences overlap or k=1.

Thus a pair with globally minimal earlier start must either have i=0 or
different predecessors. Equivalently, any repeated pair can be shifted left
until one of those alternatives holds. The independent controls test this
constructive version for every repeated pair, in addition to testing the
globally minimal pair.

For the actual coded trajectory, an anchored pair cannot be equal by the
positive-carry source argument. A pair with different predecessors cannot
be equal because the two image classes are disjoint. This provides a direct
way to eliminate the matched-block return alternative. It is consistent
with, and weaker than, the all-checkpoints-distinct reverse-decoding lemma.

## 5. Mechanical complexity and variable-length encoding

For V, let theta=log(s_V)/log(4/3). Exact rational multiplier comparisons
give 0<theta<1 and the completed-iteration invariant

    1<=S_n<4/3,
    f_n=floor(n theta)-floor((n-1)theta).

A length-k increment factor is a function of the starting phase in [0,1).
Its possible discontinuities occur at the k points {-j theta}, 1<=j<=k.
Together with the fixed endpoint zero, they partition [0,1) into at most
k+1 half-open intervals, on each of which the full k-symbol factor is
constant. The floor convention assigns boundary values to the correct
half-open intervals. Consequently the finite prefix has at most k+1
different length-k symbol factors.

This uses no recurrence theorem and makes no assertion that every infinite
rotation factor occurs in the current prefix. Irrationality follows from
prime factorization as stated; the factor upper bound itself does not need
any partial-quotient estimate.

For m>=3, k=ceil(m/100) satisfies m>2k. Therefore the m-k+1 actual occurrence
positions exceed k+1, forcing a repeated factor wholly in f_1...f_m.

The corresponding code images have lengths 122 or 124. Two equal symbol
factors encode exactly the same literal parity word, of length r>=122k.
They may have different absolute bit offsets in the full trajectory, but
that does not affect literal equality. Both complete encoded futures fit
before the final extra F because their symbol endpoints are at most m.
Thus 100r>=122m=L is established at actual macroblock checkpoints.

There is no incorrect inference that arbitrary literal parity factors have
complexity at most k+1. The proof counts symbol factors and then encodes one
matching pair at its known code boundaries. This is sufficient for the
audited matched-block theorem and avoids the extra offset counts that a
global literal-factor bound would require.

Since the selected checkpoint values are unequal, that theorem proves
strict descent at the full actual first-crossing endpoint. The final F is
retained throughout; it is not removed by a symbolic simplification.

## 6. All-m V family and the small cases

The argument above covers every m>=3. For m=1 and m=2, the actual V family
has lengths 124 and 246 respectively, so the manuscript's conservative
bound N<=250 is correct. Both are in the audited expansion-first
short-range theorem's N<=2^52 range. Thus the statement for every m>=1 is
valid.

All finite parity words retain their unique lawful source residue modulo
2^N and infinitely many positive lifts. The result therefore concerns real
positive cylinders and excludes nonnegative endpoint displacement on them;
it is not a claim that these cylinders are empty.

The reset counts are unbounded with m, while the previously audited bound
on consecutive V runs gives eventual avoidance of the one-percent
pure-expansion repetition criterion, including the stated reparsing scope.
The new proof closes the whole greedy family's endpoint realizability and
does not rely on increasing a numerical episode cutoff.

## 7. Arbitrary fixed U and the m=2 boundary

For any fixed nonempty U over K,V, every letter expands. Starting from a
completed-iteration slope in [1,4/3), appending U followed by the maximal
number of noncrossing F resets returns the slope to [1,4/3). The additional
final F then gives the first explicit-reset crossing. Intermediate U
macroblocks only increase the slope, and every conditional reset was
chosen to retain slope at least one. Hence the complete word remains in
the audited bank scope.

Let theta_U=log(s_U)/log(4/3)=a+beta, with a=floor(theta_U) and 0<beta<1.
The reset counts satisfy

    j_n=floor(n theta_U)-floor((n-1)theta_U)
       =a+floor(n beta)-floor((n-1)beta).

Unique prime factorization excludes rational theta_U. The two code images
are UF^a and UF^(a+1), and the preceding image calculation gives disjoint
types t+a and t+a+1. All checkpoints are therefore distinct by reverse
decoding. The binary choices j_n-a have complexity at most k+1 by the same
phase partition argument.

For m>=3, equal k-symbol futures have at least k|U| bits, while total
expansion length in the actual parse is m|U|. Thus the same one-percent
matching condition holds, independently of how long U is or how many
conditional resets each iteration needs.

For m=2, take the literal U at the original source y_0 and the literal U
at y_1, immediately after the first complete iteration. Both are wholly
present. Their common length is |U|=L/2. The prefix giving y_1 has multiplier
at least one and positive carry, so y_1>y_0. The matched-block theorem
proves strict endpoint descent without a factor-complexity argument. The
overlay merely makes this intended start explicit; the submitted checker
already uses the correct states[1] checkpoint.

The corollary properly stops at m>=2. A general U can be arbitrarily long,
so its one-copy word cannot automatically be assigned to the short range.
The independent controls invoke an m=1 theorem only for U=V. No arbitrary-U
m=1 statement is inferred from the finite controls or the V special case.

The template is fixed within a word. Allowing unrelated U_n would require
fresh proofs of output separation and factor complexity; neither is
silently supplied by this result.

## 8. Exact scope of the earlier typed-factor question

The frozen earlier question sought, eventually in every prefix, a matching
pair with different incoming types. The new proof provides the disjunction
of different incoming types OR a pair anchored at the original source.
The anchor's strict positive-carry inequality supplies the missing state
inequality even without a type assigned to the original source.

This is enough to prove the entire endpoint-descent goal. It does not prove
the stricter different-types-only occurrence statement, and the reviewed
Section 7 says so explicitly. This is the only substantive scope wording
clarification; it changes no step of either descent proof.

The reusable terminal-image criterion is valid with its explicit source
anchor, separated output classes, actual factor occurrences and matched-block
scope hypotheses. F does not have a constant output type on arbitrary inputs,
so the entire bank alphabet cannot be assigned artificial fixed labels to
apply it. The manuscript explicitly retains that obstruction.

## 9. Dependencies and retained open problems

The endpoint step uses the independently audited matched-block height and
short/long-overlap theorem. Its cited explicit logarithmic-form input remains
a dependency; it is not independently re-proved by this finite certificate.
The required first-crossing, lawful-source, positive-integer and literal
factor hypotheses are verified above.

Neither theorem supplies entry into this bank language for an arbitrary
integer, covers all first-crossing bank families, excludes all possible
Collatz cycles, or determines the entire realized S_d/2-adic boundary.
The original different-types-only finite occurrence question also remains
distinct from the proved endpoint goal. These limits are preserved in the
reviewed publication text.

Reproduction:

    python3 review/independent_all_greedy_controls.py
    python3 -O review/independent_all_greedy_controls.py

Compare with the independent receipt and optimized output. The two submitted
reruns are preserved separately. The overlay metadata records both source
and reviewed hashes for publication provenance.
