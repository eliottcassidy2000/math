# Scoped independent review: local rotation-defect ideal

Date: 20 October 2026.

Verdict: **PASS for the all-word theorem, its stated arithmetic corollaries,
and the critical-slope correction.** The numerator/vertex terminology
clarification requested during review is now applied. This is a proof review
of the supplied mechanism, not a
global Collatz audit, a novelty/priority certification, or a convergence proof.
The repository checkout was read only. No Git operation was performed.

## Proof findings

1. **Connected graph content.** For a nonconstant binary word, at least two
   distinct rotations exist. Carry injectivity makes every edge difference
   nonzero. If the full-edge gcd is h, graph connectivity gives one common
   residue for all rotation carries. A zero transition gives 2r=r, hence r=0
   modulo h without requiring 2 to be a unit at this stage. A one transition
   then gives h dividing D. Conversely, gcd(D,B) is odd and propagates through
   the exact rotation recurrence. Thus the full-edge content is exactly g.

2. **Context stripping.** Each full edge is 2^s a^t times its local defect.
   Since D=2^q-a^p is coprime to 2a for odd a, the established h dividing D
   permits cancellation of each context. The local gcd both divides h and is
   divisible by h, so it equals h. Maximal stripping is unnecessary. The proof
   applies to every connected graph and to powers of words, not only to
   primitive necklaces or coprime shapes.

3. **Injectivity and width bound.** The first bit is the carry modulo 2; after
   subtracting a^(p-1) when that bit is 1 and dividing by 2, the suffix is
   recovered. Known length and weight supply the claimed injectivity. Moving
   10 to 01 strictly increases the carry by its positive context factor. This
   proves the unique extrema and the bound
   (2^(ell-h)-1)(a^h-2^h)/(a-2). Distinct equal-shape cores necessarily have
   0<h<ell, so their bound is positive. Adjacent-swap cores have defect 1.

4. **Cycle gates and signs.** The reduced denominator is |D|/g. Scaling by
   an odd nonzero additive parameter b gives actual cycle vertices bB_i/D,
   with the prescribed parity word on odd-denominator rationals. All vertices
   are integral exactly when |D|/g divides b. Since B_i>0, positivity is
   equivalent to b/D>0. Oddness of b is needed for the asserted parity map;
   the numerical scaling formulas themselves do not need it. A primitive
   word and its powers produce the same rational orbit, with the power
   factor H>0 multiplying both D and B, and hence g.

5. **Scaled content and lattice interpretation.** The integer affine-carry
   list bB_i has full-edge content |b|g by direct scaling. If b shares a
   factor with 2a, one should not try to cancel the context using |b|g as a
   unit; apply the unscaled theorem first. A spanning tree generates all
   full pair differences by paths, and a nonzero one-row integer list has
   its gcd as its single nonzero Smith invariant. The localization-and-
   recovery explanation is valid because the unscaled generator is coprime
   to 2a. No primitive-shape assumption is hidden here.

Terminology clarification applied: the latest note calls b B_a the
numerator-carry list and explicitly gives the actual ax+b cycle vertices as
b B_a/D. Equivalently, b B_a are free-cycle vertices for the additive parameter
bD. This removes the earlier ambiguity between numerators and actual cycle
coordinates. The pinned main-note SHA-256 begins `4e9ac5a9`; the earlier
`1a6d808d` revision had the terminology ambiguity.

The finite-bank norm paragraph is also valid under its explicit fixed-norm,
uniform-comparability and every-block monotonicity assumptions: c*y<=C*x is
incompatible with y/x>lambda^L for arbitrarily large L. An existing exact
expanding witness is the block K in
`math/05-knowledge/results/collatz_unbounded_return_centers_20261019.md`,
with multiplier 3^19/2^30>1, positive carry, and arbitrarily large positive
sources in each finite repetition cylinder. This does not obstruct
nonuniform gauges, unbounded states or selected-return descent. A precise
source link/name for the note's J would improve provenance; it does not affect
the new ideal theorem.

## Scoped critical-slope correction

The actual THM-4717 source was read, including its setting, statements,
proof, status and related-entry gloss. Its YAML title calls log_3 2 “the
slope of every integer cycle”; Statement 2 explicitly says “the critical
slope is exactly the slope of every integer 3x+1 cycle.” Its surrounding
discussion of limits does not qualify that finite-cycle equality. The
related-entry gloss also omits the conditional convergent hypothesis.

THM-4633's setting uses k odd members and K total shortcut steps, so the
translation to the present notation is p=k and q=K. Its product identity
gives

    q log 2 - p log 3 = sum log(1+1/(3n))
                       in (0, p log(1+1/(3m))].

Consequently the correction's strict logarithmic upper error and the bounds

    log 2 / log(3+1/m) <= p/q < log_3 2

are correct. The lower equality is attained by the trivial cycle. The
upper equality is impossible for every finite word by unique factorization.
The minimum-to-infinity limit follows by squeezing, conditional on a family
of such cycles existing. THM-4633 Statement 3 explicitly needs
k^2 <= (3/2)m log 2 for its convergent conclusion. The proposed correction
is therefore approved on the stated scope. No endorsement or audit of
THM-4717's topology, exceptional-point limits, or Cantor/conjugacy claims is
included.

## Independent finite controls

`independent_controls.py` imports no producer or repository mathematics.
It was run in normal and optimized Python, with byte-identical receipts.
The all-length result rests on the proof above; these controls test exact
arithmetic and adversarial cases.

- All 1,004 nonconstant words of lengths 2 through 9 at each of bases 3, 5,
  and 9, plus the negative-17 word at base 3: 3,013 word/base cases.
- A binary-parent spanning tree, with independently chosen amounts of
  context stripping: 71,707 allowed context/core identities and bounds.
- Exact Fraction replay for b=-5,-3,-1,1,3,5: 18,078 rational-cycle closures,
  with parity, integrality, sign, least-period and scaled-content checks.
- Injectivity and unique extrema for all 108 length/weight/base shapes.
- 36 word-power cases including the positive nonunit, sparse unbalanced,
  trivial, and signed nonunit examples.

The producer's two supplied receipts were compared and are byte-identical.
The producer's verifier was inspected, including its explicit RuntimeError
checks and finite universes, but its 5.6-million-check run was not repeated.

Reproduce the independent controls from the workspace root:

    python3 session11-collatz-review/independent_controls.py
    python3 -O session11-collatz-review/independent_controls.py

Both outputs should match `independent-receipt.json`. Audited source and
artifact SHA-256 values are in `AUDIT_PINS.json`.
