# Independent audit of the moving-centre braid certificate

**Status: independent analytic/source audit PASS; normal, optimized, and
frozen-output replay PASS. Audit frozen.** The general sufficient tube
criterion is proved analytically. The named coarsening and synthetic
controls are FINITE-EXACT. No old certifier or witness was edited.

The primary proof is
[planar_jc48_sep08_moving_tubes.md](planar_jc48_sep08_moving_tubes.md).
I read it and the full standalone source, including the coefficient
expansion, exact separation calculation, adaptive selection, and
crossing extraction. The final source's additional literal-table and
shape guards are included in this acceptance.

## 1. Complete affine substitution and uniform Rouché bound

Write the Taylor expansion at a starting base point and centre as

    F(u0+U,z0+Z)=sum a_jk U^j Z^k.

The source computes this by the ordinary binomial expansion of every
monomial of the literal polynomial. Substituting
`U=theta du` and `Z=theta dz+w` gives the exact coefficient rule

    C_(j+k-l,l) += a_jk du^j binom(k,l) dz^(k-l),
    0<=l<=k.

This is precisely its `moving` loop. Contributions to the same power
are added before absolute values are taken. In particular all
cancellations due to common root motion survive. There is no Taylor
truncation: all available base and root powers are retained.

For a proposed radius `R`, the left side
`B(C_01)R`, where `B=max(|Re|,|Im|)`, is a lower bound for
`|C_01 w|` on `|w|=R`. On `0<=theta<=1`, every other term is bounded
above by `A(C_jk)R^k`, where `A=|Re|+|Im|`. Thus a positive
stored margin is a strict Rouché inequality uniformly on the entire
segment. It includes the constant residual, every higher power of
theta, every nonlinear power of w, and the varying linear terms
with indices `(j,1), j>0`. Comparing to `C_01 w` gives exactly one
counted zero in the disk; hence that zero is simple.

As a separate algebra check, I reconstructed the actual quartic from
the resultant of

    U(T)=T^4-2T^2,
    V(T)=T^6-3T^4/2+T^3/3-T.

Its full coefficient dictionary equals the prototype's literal table.
I then used independent symbolic affine substitution in two named
exact complex controls and compared all **38** resulting coefficients
with the source expansion. Every coefficient agreed. This check
does not use the old continuation implementation.

## 2. Exact separation, endpoint gluing, and the actual homotopy

For a pair of centres put `d(theta)=d0+theta delta`. The function

    B(d(theta))=max(Re d,-Re d,Im d,-Im d)

is the maximum of four real affine functions. All changes of its
active affine function occur when a coordinate is zero or when
`Re d=Im d` or `Re d=-Im d`. Between such points it is affine,
so its minimum occurs at a listed breakpoint or an endpoint; a
flat minimum is covered as well. The source enumerates exactly
these rational candidates, retaining only those in `[0,1]`.
Consequently the returned quantity is the exact minimum L-infinity
distance, not a numerical sample or a guessed lower bound.

Euclidean separation is at least this quantity. Each radius is
chosen as one sixteenth of the smallest relevant pair minimum,
and every pair's strict radius-sum inequality is also tested.
Positive separation therefore gives pairwise disjoint Euclidean
disks at every segment parameter. Four one-root disks exhaust the
four roots of the verified monic quartic; there is no fifth root
or unaccounted multiplicity.

At a joint between consecutive segments, the centre for each label
is exactly the same rational number on both sides. Its two endpoint
disks are concentric. The smaller disk contains one root and is
contained in the larger disk, which also contains exactly one root.
Thus these are the same root. This replaces the old tiny-endpoint
isolation test by a valid stronger geometric circumstance; it does
not assume that arbitrary nearby endpoint proposals have matching
labels.

Each actual root and its affine centre lie in the same moving
convex disk. Their straight contraction stays in that disk, while
different disks remain disjoint. These contractions give a labelled
configuration homotopy, agree at every joint, and agree at the
beginning and end after matching the identical unordered centre
sets. The actual base configuration is identified with the rational
centre configuration by this fixed contraction path.

For several new loops, identical ordered initial centres give the
same contraction and hence one common based identification. The
method does not license choosing unrelated conjugations for
different loops. This is the relevant condition when transferring
the accepted tube method to a future multiple-loop proof.

## 3. The complete inherited input and independent word extraction

The input is exactly the old smooth-loop witness, with 9,554 rows
and 9,553 segments, based at `4+3i`. Its six edge endpoints are

    4+3i, 1/32, i/32, -1/32, -i/32, 1/32, 4+3i.

The source verifies that all input points advance strictly along
these edges. It selects new endpoints within one edge at a time,
so each coarsened base segment follows the original path exactly.
There is no shortcut across a corner or across the interior of the
loop. I independently checked all 9,554 input rows have exactly
two base coordinates and four two-coordinate root centres; the
final source checks these shapes as well.

The old source is used only through its raw hash and safely parsed
literal coefficient table. The new table is checked equal to it,
and the full monic degree-four condition is checked separately.
No old continuation or braid-extraction routine is imported. The
old gzip and raw JSON are hash-pinned. Their old root-path claims
are not a substitute for the new moving-disk inequalities.

The greedy doubling and bisection need not find a maximal feasible
segment: feasibility has not been proved monotone, and the primary
does not claim otherwise. Every chosen segment is directly accepted
by its own exact margin and separation tests. Its endpoint index
strictly advances. If necessary the adjacent original endpoint is
tested and a failure aborts; the old subdivision is not automatically
accepted just because it was previously certified in another model.

For the new centre polygon, multiplication by `1+i/4` is an
orientation-preserving complex projection. The code verifies
nonzero real pair differences at all projected endpoints and
computes each crossing time by an exact rational formula. It
rejects simultaneous crossings, nonadjacent crossing strands, and
zero imaginary separation. For each segment an affine real pair
difference has at most one crossing, so the complete event list
is exhausted by the tested sign changes. The current order is
updated at every event.

The initially left strand passing below gives the positive
counterclockwise half twist in the declared convention. Only
adjacent inverse letters are cancelled. The resulting word is
computed independently and then compared with the stored named
word; that comparison alone is not treated as the proof of the
braid. The tube homotopy above is what relates the centre polygon
to the actual roots.

## 4. Positive and hostile controls

The common-translation example has

    F(u,z)=product_(a=0)^3 (z-u-a).

When `dz=du`, substitution leaves `F` independent of theta, so
every positive-theta coefficient vanishes for **any** common
translation. The large literal translation `10^6(1+i)` checks
this cancellation in the code. Constant separations and the four
positive Rouché margins then certify the whole segment. The
unbounded statement is this algebraic cancellation, not a census
of translation sizes.

Swapped endpoint labels give an actual affine-centre collision
halfway along the segment, despite separated endpoint sets. The
separation calculation correctly returns zero and rejects the
segment. The additional exact examples cover an interior
`Re=Im` minimum and a flat minimum interval.

The curvature hostile is

    F=(z-u^2(1-u))(z-2)(z-3)(z-4).

Its small root starts and ends at zero and has vanishing first
path residual, but at `u=2/3` it equals `4/27>1/16`. Therefore
the proposed constant-centre disk is actually inadequate; this is
not merely a failed sufficient inequality. The full higher-term
Rouché bound rejects it. This pays the distinction between a true
moving tube and a first-order cancellation test.

## 5. Independent replay and accepted scope

Independent runs of

    python3 -B 04-computation/planar_jc48_sep08_moving_tubes.py
    python3 -B -O 04-computation/planar_jc48_sep08_moving_tubes.py

both pass **10,324 always-active gates**. Their complete outputs
are byte-identical to the frozen primary output: **374 bytes**.
They reconstruct **92** accepted moving segments in **654** exact
attempts from the complete 9,553-segment input path. The full reduced
word is exactly `[3,1,-3]`. The semantic digest records all selected
endpoint indices, radii, and positive margins. Wall time is not
part of the deterministic output and is not a universal claim.

| Artifact | SHA256 |
|---|---|
| New source, 11,199 bytes | `735e928d0d39bc43d6f1ec3e3d8a4965b8181813cf6094d6e8b50bb8fcf59990` |
| Frozen output and both independent replays | `14606ab9af4e2b6d25f26b4b6c7c7638c9577c51c4eb8f652bf87cf717d12ea5` |
| Chosen-segment semantic record | `77c95276801a4096438b0c2c018124b2bf35fb3323b7f572689f99a6752f73e9` |
| Inherited polynomial producer | `d5f0d1ce16b439c13d936f35ce38f558b1c5210f892a11218ab7f4f4695198e1` |
| Input gzip witness | `3ae1dfc0ada5d2fe148cb4da1b4643fcfd10280e66ac8a73ad783b6a007b22b9` |
| Input raw JSON | `06d050671046d9145ce23c82fa7fbbbeaa2d8a73ba62fd664b74cdc8046820be` |

No mathematical correction remains. This accepts the general
sufficient moving-tube method and the one named exact coarsening.
It supports using the method in **new** production. Any new curve
still requires its own literal polynomial, complete path data,
uniform certificates, common based gauge, and downstream topology
or group proof. No unproduced braid conclusion or modification of
an older pinned certificate is included in this acceptance.
