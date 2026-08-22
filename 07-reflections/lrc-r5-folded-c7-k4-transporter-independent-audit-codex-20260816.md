# The folded heptagon obstruction is sound, but its missing coordinate must be source-side

**Status: independent audit supporting PROVED THM-3524 + FINITE-EXACT.**  The
submitted script, output, and reflection had
LF SHA-256 values

```text
5fe3f696f122869462ff73b1b4ebdb957fa8ca7ee3692c25ef94f0f7efae81cf,
b97db9c8b0582d868c8386be5f5d8dff72968062a9b5ae79956678cd86c0b98a,
e8de54937b19871cb391e41b7efbfd5d52570bfad35edead68ba5e92befc00e1.
```

They were inspected and replayed without modification.  Normal and optimized
executions are byte-identical to the submitted 25-line output.  A disjoint
verifier reconstructs the source and target by a different route and accepts
the mathematical claims in the scope stated below.

## Verdict

The package is sound as a **marked split-field transporter classification and
fixed-linear-class obstruction**.  It proves all of the following.

1. The negation-orbit quotient of the Paley tournament on seven vertices has
   four orbit blocks and off-diagonal `K4` support.
2. No one of the `48`, `192`, or `2401` specified labelled source-mode
   allocations admits a common scalar drift convolution, even after
   independent amplitude and cyclic-origin calibration.
3. More generally, no frequency-independent channel mixing followed by one
   common circulant can send the rank-three THM-2594 source bank to the
   rank-four U_full Walsh bank.
4. Four channel-dependent invertible circulants do give a formal split-field
   transport.  This is the sharp survivor, not a common-base transplant.

It does **not** construct a U_full ancestry relation, Boolean coupling,
physical current, absolute graph `H^1` class, grouped coefficient,
bispectrum, scalar-row exclusion, or LRC(14) theorem.

## 1. Folded Paley and Walsh algebra

For the Paley orientation `x -> y` when `y-x` belongs to `{1,2,4}`, folding
under negation gives

```text
O0={0}, O1={1,6}, O2={2,5}, O3={3,4}
```

and the independently reconstructed directed arc-count matrix

```text
0 1 1 1
1 1 2 2
1 2 1 2
1 2 2 1.
```

Because `-1` is a nonresidue modulo seven, negation reverses each arc; hence
the off-diagonal matrix is symmetric.  The exact quotient is therefore a
**weighted looped symmetric four-block digraph**.  Its off-diagonal support
is a bidirected or undirected `K4`; the three diagonal ones retain the
internal arcs of the nontrivial two-point blocks.  Promotion should keep this
qualification whenever the slogan “the fold is K4” is used.

Multiplication by two cycles `O1 -> O2 -> O3 -> O1`.  Independently, the
three nontrivial characters of `F_2^2` cut its four vertices into

```text
01|23, 02|13, 03|12.
```

This identifies the common three-matching grammar, but it supplies no
canonical map between the two carriers.

The descent statement also checks exactly.  The inversion-fixed rational
group algebra is

```text
Q[C7]^inv = Q x Q(zeta7+zeta7^-1),
```

where the cubic factor has minimal polynomial
`x^3+x^2-2x-1`, discriminant `49`, and cyclic Galois group.  It is not
`Q[F_2^2]=Q^4`: the former has two rational simple factors and the latter
four.  In the certified split field, the three cubic roots are independently
recovered as

```text
424661381881782662835,
169426766391843876494,
550417624219391222408.
```

Thus the marked identification becomes available only after splitting and
forgets rational descent.  This is not the `S3` cubic carrier of the fixed
sporadic Keller coordinate equations.

## 2. Completeness of the three allocation searches

The independent verifier first generates all `7^4=2401` ordered labelled
mode choices and then filters by predicates, rather than reusing the
candidate's permutation construction.

| census | exhaustive universe | count | exact common | gauge/amplitude common | kernel classes |
|---|---|---:|---:|---:|---:|
| trivial-preserving | first target receives mode `0`; fold classes are exactly `{0,1,2,3}` | 48 | 0 | 0 | always 4 |
| arbitrary folded bijection | four fold classes are exactly `{0,1,2,3}` in any target order | 192 | 0 | 0 | always 4 |
| unrestricted selection | every function from four target channels to seven labelled modes, repetitions allowed | 2401 | 0 | 0 | always 4 |

Amplitude is removed by normalizing at frequency zero.  Cyclic drift-origin
gauge is removed by taking a canonical representative over all thirteen
character phases.  Every allocation retains four distinct canonical ratio
rows.

These universes are complete for what their names say.  The `48` and `192`
objects are marked split-field lifts, not rational folded-algebra maps.  The
`2401` objects are coordinate selections, not arbitrary linear combinations.
The separate rank theorem below handles arbitrary fixed linear channel
mixing.

## 3. Rank obstruction by an independent linear system

The raw THM-2594 aggregate has only the three columns `theta=0,1,2`.  An
integer `3 x 3` minor on rows `(1,3,4)` is

```text
48132400318030844942585467137428417511789621695304247254245705613312000,
```

so its characteristic-zero rank is exactly three.  Its two-dimensional
Fourier transform is nonzero at all `91` coordinates.  Independent reductions
at `547`, `1093`, and `2003` each reproduce support `91/91` and rank three.

The target was rebuilt from the THM-3514 independent atom engine by direct
`tau` slices, without importing the submitted transporter or its U_full
parent.  Its four Walsh rows have all `52` Fourier coordinates nonzero and
rank four.  The first four drift columns already have determinant

```text
154930576847978146912 mod 572252886246508880869.
```

For each of the eight antipodal sign lifts, the candidate eliminates the
frequency multipliers and obtains projective equation rank twelve.  The
independent audit instead keeps the multipliers.  With sixteen entries of a
fixed `4 x 4` channel matrix `M` and thirteen scalars `lambda_k`, it solves

```text
M x_k = lambda_k y_k,       k in F_13.
```

The resulting `52 x 29` system has rank `25`, nullity `4`, in all eight
cases.  Every nullspace basis vector has zero projection to all thirteen
`lambda_k`; its `M` component annihilates the entire source span.  Eliminating
the multipliers independently gives rank twelve.  This proves, rather than
dimension-counting heuristically, that the projective solution space is
annihilator-only.

The finite systems are redundant with the general obstruction.  For all
seven source profiles at once, any fixed channel mixer and any common drift
operator have

```text
Y = M X C,
rank(Y) <= rank(X) = 3 < 4.
```

The conclusion remains true even if `C` is an arbitrary `13 x 13` matrix;
circulant and invertible are not needed for the negative rank inequality.

The positive formal survivor also reproduces exactly: the four full kernels
and four augmentation kernels have support `13/13`, rank four, and digests

```text
9f03166d1b8d730673aeccf7bff9196d2de42e0ff9fa7c08680d38a772dc041f,
351976b6d1ac51915d686b32e1cc0f3f3dfa53d1f33f7c8831f7a74a495e3025.
```

They are channel-dependent.  Their existence is automatic for pairs of
cyclic vectors and does not pay the ancestry debt.

## 4. THM-3518 does not supply the fourth source channel

The phase and cut-cycle certificates of
`THM-3518-ufull-all-role-phase-normalized-dual-contraction-and-cut-cycle-certificate.md`
are valuable target-side typing data, but neither is a fourth THM-2594 source
response.

- The factor `zeta^(q_t a_L)` is defined on the U_full endpoint guard-sheet
  contraction.  No map identifies `a_L` with a THM-2594 common-base root or
  deep coordinate.  If one common phase/right action were supplied on the
  source, it would be an invertible column operation and preserve rank three.
  If phases are chosen separately by target channel, one has left the
  fixed-mixing/common-circulant class and merely reintroduced the four
  channel-dependent transporters.
- THM-3518's edge arrays are gradients.  The incidence rank is seven, the
  cycle rank is six, and every one of the six cycle coordinates is zero.
  Appending a zero absolute-cycle packet adds no independent response.  The
  nonzero matrix-tree polynomial depends on a chosen cut representative; it
  is not a hidden `H^1` channel.
- Most importantly, THM-3518 starts after U_full endpoint marginalization and
  explicitly retains no common ancestry base.  A lawful fourth source channel
  must be defined on the THM-2594 Boolean fibre before that loss.

Therefore the answer to the requested sidecar question is **no**.  A fourth
independent common-base source observable is necessary for this linear class,
but not sufficient for a physical U_full transplant.

A concurrent, still-provisional exact probe makes the distinction sharper.
THM-2594's genuine fixed-absolute-root reindex is source-side and raises the
union of the slaved and absolute profile planes from rank three to rank four.
Nevertheless its full fourteen-channel projective system, all twelve affine
torsor dilations, all 56 one-named-sidecar choices, and all 127 binary
absolute-sidecar choices still have zero excess beyond source annihilators.
Thus the newly found fourth coordinate does not by itself solve the common-
operator transplant.  This incoming result is not part of the present audited
package and should be independently promoted, if at all, on its own evidence.

## Promotion outcome

THM-3524 promotes the package as a finite-exact theorem with all four scope
clauses visible:

```text
PROVED FINITE-EXACT over the certified split field;
48/192 are marked folded bijections and 2401 are labelled selections;
the general no-go is only fixed channel mixing plus one common right operator;
no U_full common ancestry/current/H1/bispectrum/LRC(14) consequence.
```

Use “weighted looped quotient with off-diagonal K4 support” for the exact
Paley fold, and phrase the missing fourth channel as necessary, not sufficient.
No concrete mathematical repair to the submitted files is required.

## Independent artifact

Run

```text
python -B 04-computation/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_folded_c7_k4_transporter_independent_audit_20260816.py
```

The independent script/output/semantic SHA-256 values are

```text
0bf56db1acfd83d6f161dd95b5d428df90f8cb5e6ec24e729081f24a86f287f2,
c20c59264590cd42b108363c91736df8c9e921ac94b09c1f5fbb5d977f300859,
a6d6ad7a9891af8bddcba6fe7cf9f2554214731043844098f8556d1003885222.
```
