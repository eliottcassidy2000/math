# Runner clearance curvature: exact recurrence, signed knots, and an unbounded-order hostile

Status: **PROVED scoped**, with exact reproducible controls. No external literature novelty claim. LRC(14) remains **OPEN**. This note completes the delta-ladder idea proposed historically in `01-canon/theorems/THM-1820-lrc-is-a-moment-nullcone-problem-relation-lattice-pairing.md`, without inheriting its experimental interpretations.

## Result first

For every finite nonempty family of positive integer speeds, the complete distribution of its clearance function has an exact finite signed height-curvature representation. Its rational moment generating function has a noncanceling pole at the reciprocal of the actual loneliness maximum, including when the target is attained only on a measure-zero set. The representation is invariant under common speed dilation. Its recurrence order, however, is unbounded even for primitive **two-speed** rows: `(1,N)` with even N has exactly **N-1** nonzero positive-height poles. Thus fixed runner count supplies no universal finite recurrence order for this exact representation.

The positivity hypothesis fails already at `(1,4)`: the positive-height atoms are `5/2 at 1/5`, `-3/2 at 1/3`, and `5/2 at 2/5`. Negative knots are interior changes in the inverse speed of a moving excursion boundary, not component births or deaths.

## 1. Objects and the exact map

Put `f_v(t)=min_i ||v_i t||` on the unit circle and let `M=max f_v`. Partition the circle into maximal affine segments of f, retaining actual endpoint times and active signed runner slopes. Every slope s is a nonzero integer. For a segment with endpoint heights l<h, add `1/|s|` to c_h and subtract `1/|s|` from c_l. Aggregate equal heights, discard zero coefficients, and write

`ν_v = Σ_h c_h δ_h`.

All heights and coefficients are rational. For δ≥0, let `μ_v(δ)=Leb{t:f_v(t)≥δ}`. A segment contributes its full length below l, `(h-δ)/|s|` between l and h, and zero above h. Consequently

`μ_v(δ)=Σ_h c_h (h-δ)_+`.                                      (1)

This identity is valid at every δ≥0, including the maximum; f has no constant segment, so finite equality level sets have zero Lebesgue measure. Equivalently ν on positive heights is the distributional second derivative of μ. Its atom at zero is a normalization sidecar, with `Σ_h c_h=0`, `Σ_h c_h h=1`.

**Preserved predicate:** `M≥a` iff the largest surviving positive-height atom is at least a. **Lost coordinates:** phase location, runner ownership, ties, ordering of excursions, and chronology. **Required sidecar for an actual witness:** an affine endpoint time attaining the maximal height, with its tied active slopes. The measure alone is a maximum certificate, not an arrival/owner transport certificate.

**Dilation:** the covering map t↦dt preserves uniform circle measure, so μ_dv=μ_v; uniqueness of distributional derivatives gives ν_dv=ν_v. The segment description independently explains this: d times as many copies with slopes d times as large exactly cancel. Primitive normalization is therefore respected.

## 2. Why the final atom cannot cancel

M>0: finitely many runner zeros do not cover the circle. No affine segment has zero slope, so every attainment of M is an endpoint of positive-length segments strictly below it in their interiors. Every incident segment contributes a positive `1/|s|` at M. No segment has lower endpoint M, hence no negative contribution occurs there. Thus c_M>0 and all atoms have height at most M.

This is the indispensable endpoint theorem. A fixed-threshold mass test misses an isolated equality witness; the *whole height distribution* still detects it because its support closes at M. No assertion that one fixed mean proves equality is made.

## 3. Exact rational moments and minimal recurrence

For k≥1 define `A_k=Σ_h c_h h^k`. Direct integration on each affine segment gives

`A_k = k ∫_0^1 f_v(t)^(k-1) dt`.

Let K be the set of nonzero positive-height atoms and r=|K|. Then

`R_v(z)=Σ_(j≥0) A_(j+1) z^j = Σ_(h∈K) c_h h/(1-hz)`.          (2)

Every displayed residue is nonzero and the poles are distinct, so the reduced denominator is exactly `Π_(h∈K)(1-hz)`. In particular the nearest positive pole is at `1/M`, and the exponential sequence A_k has minimal constant-coefficient recurrence order exactly r. Proof: any annihilating polynomial P satisfies `Σ_h c_h h^k P(h)=0` for all k≥1; the Vandermonde matrix at the r distinct nonzero nodes forces P(h)=0 at every h. Conversely their product annihilates. This is a statement about this particular recurrence representation, not an impossibility theorem for all algorithms or certificates.

At most 2r exact moments, with r known as an upper bound, suffice for standard rational reconstruction, but this is an instance-dependent method. Signed residues prohibit treating ν as a positive atomic measure. The *original* pushforward distribution f_*Leb is positive and has the ordinary moment inequalities; these two measures must not be conflated.

## 4. An all-height hostile at fixed runner count

Let N=2m≥2 and v=(1,N). Then the complete atom formula is

`c_0 = -4 + 2/N`,
`c_(k/(N+1)) = 2 + 2/N` for `1≤k≤m`,
`c_(k/(N-1)) = -2 + 2/N` for `1≤k≤m-1`.                       (3)

There are no other atoms.

**Proof.** Reflect across t=1/2 and work on [0,1/2], where the slow runner distance is t. In each fast tooth, the rising branch is Nt-k and the falling branch is k+1-Nt. The lower envelope with t changes rising slope N→1 at `t=k/(N-1)`, for k=1,...,m-1, creating an interior curvature coefficient `1/N-1`. It changes rising slope 1→-N at `t=k/(N+1)`, for k=1,...,m, creating a maximum coefficient `1+1/N`. Reflection doubles both. At zero, normalization gives the stated c0; equivalently the zero branches contribute it directly.

Since gcd(N-1,N+1)=1 for even N, none of the displayed positive and negative heights coincide: equality requires k to be a multiple of N+1, impossible in its range. Every coefficient is nonzero. There are m+(m-1)=N-1 positive-height poles. Therefore r=N-1 is unbounded on this primitive two-speed family. Its exact maximum is `m/(N+1)`, with terminal atom `2+2/N`.

The same formula explains the first hostile N=4 and the positive control N=2. It also shows why component-death data alone cannot reconstruct the mass profile: the negative intermediate knots are genuine and invisible to a barcode listing just component maxima. The barcode *does* preserve M; no false stronger obstruction is asserted.

## 5. Two hypothesis/test/refinement cycles

1. **Bold analogy:** exact rational generating functions in THM-3181 might turn clearance levels into a positive finite spectral measure. **Test:** exhaustive rational lower envelopes for all 162 nonempty subsets of {1,...,8} of size at most four. **Failure:** (1,4), with its negative 1/3 atom. **Repair:** signed height curvature with endpoint positivity only; the preserved target is the maximum, not positivity of all coefficients.
2. **Bold analogy:** fixed runner count might bound a recurrence order, making exact finite moment reconstruction uniformly small. **Test:** continue primitive (1,N), not just a bounded random census. **Failure mechanism:** new inverse-boundary knots appear at every fast tooth. **Repair:** prove formula (3), establishing exact unbounded minimal order N-1, while retaining an instance-dependent complete instrument.

The least-used sidecar recovered from the archival sample is the actual endpoint/owner slope, not a tournament orientation. No intrinsic pairwise orientation was forced.

## 6. Reproducible computation

`probe.py` constructs the actual minimum of affine runner branches by exact interval clipping. It tests all 162 nonempty subsets of {1,...,8} of size≤4, five moment identities, normalization, terminal positivity, and common doubling for every row of size≤3. Controls include the tight 13-speed AP and the historical deep well `{1,...,12,182}`. The latter yields M=14/183 but does **not** assert that this is its only lonely time.

`family_check.py` compares (3) against exact envelopes for every even N=2,...,80. A separate fixed-threshold wall sweep at five rational thresholds for each N checks the positive-part mass formula directly against the original inequalities, without using the computed affine envelope. Both scripts use explicit exceptions, not assertions; normal and optimized outputs byte-match.

Reproduce from this directory:

    python3 probe.py > probe-output.json
    python3 family_check.py > family-output.json
    python3 -O probe.py > probe-optimized.json
    python3 -O family_check.py > family-optimized.json
    cmp probe-output.json probe-optimized.json
    cmp family-output.json family-optimized.json

No finite computation proves a uniform LRC statement. No lower bound for the maximum beyond existing LRC facts is obtained here. The new usable output is the exact transform, terminal noncancellation mechanism, and unbounded recurrence-order obstruction.
