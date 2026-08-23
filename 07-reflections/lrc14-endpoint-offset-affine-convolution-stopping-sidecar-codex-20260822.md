# Endpoint fibres do not linearly recover the deep-offset profile

**Status:** FINITE-EXACT stopping sidecar on the thirty typed THM-3672
positive non-cover charts.  This has no theorem ID, makes no genuine-cover
claim, and proves no LRC(14) conclusion.

## The proposed connection and its exact failure

THM-3710 and THM-3713 both expose a cyclic thirteen-slot object on the same
thirty charts, but the common `C_13` syntax does not by itself define a map.
The exact comparison uses

```text
source: (F^+_(k,l),F^-_(k,l)) in Z^13 x Z^13,
        the two signed full-grid successor-endpoint fibres of THM-3710;

target: Delta_(k,l) in Z^13,
        the deep-offset three-site defect of THM-3713.
```

For every multiplicative reindexing `alpha in F_13^*`, test the most general
two-input cyclic convolution

```text
T_(alpha,K)(F)_u
 =sum_j K^+_j F^+_(alpha u-j)
  +sum_j K^-_j F^-_(alpha u-j),                    (1)
```

and its affine version `c+T_(alpha,K)`.  The free kernel index absorbs every
additive affine shift `beta`.  Thus (1) contains direct affine reindexing,
arbitrary mixtures of the two signs, and every fixed cyclic finite-difference
refinement.  If it existed, it would preserve the exact chartwise offset
coefficients and the cyclic translation law, up to the declared automorphism
`u -> alpha u`.

There are `30*13=390` equations.  The linear system has 26 kernel unknowns;
the affine system has 27.  For every one of the twelve values of `alpha`, the
coefficient and augmented ranks are

```text
                         mod 1000003       mod 1000033
linear system               (26,27)           (26,27)
with constant term          (27,28)           (27,28).            (2)
```

Either prime already proves nonexistence over `Q`: a full-rank minor modulo
the prime is a nonzero integer minor, so the rational coefficient ranks are
26 and 27 while the corresponding augmented ranks are at least 27 and 28.
The second prime is an independent hostile control.  Hence no map of the
form (1), affine or linear, recovers the deep-offset profile on even this
single finite control.

The characteristic-thirteen shadows disagree even earlier.  All thirty raw
offset numerator profiles are the zero vector modulo thirteen, whereas the
signed endpoint fibres form thirteen distinct profiles.  Dividing each
offset profile by its own integer content produces nonzero primitive shapes,
but that chart-adaptive nonlinear saturation cannot repair a rational linear
intertwiner.

## What survives and what was lost

The exact chart-variation ranks over `Q` are

```text
rank(endpoint)=6, rank(offset)=5, rank(joint)=6.                  (3)
```

Consequently an affine linear fit from the 26 endpoint coordinates to the
13 offset coordinates exists on these thirty points.  Equations (2) show
that every such fit breaks all affine `C_13` covariance.  This is a
finite-dataset dependence, not a canonical connection, and it should not be
transported to a cover row.

The loss is geometric.  The THM-3710 source keeps signed boundary
multiplicity only on the two 182-grid cosets.  It discards the complementary
endpoint phases, the pairing and ancestry of interval boundaries, and the
interior mass seen by the deep window.  The THM-3713 target integrates the
whole marked interval set against thirteen translated deep windows, but
forgets which endpoint and owner created that mass.  Neither 26-vector
contains the data needed to reconstruct the other equivariantly.

The smallest plausible repair is therefore not another reindexing or Hasse
jet.  It is the full signed endpoint ledger together with interval
pairing/owner ancestry and the deep-window integration kernel; equivalently,
retain the five-colour offset profile itself as the sidecar.  Scalar Hasse or
moment coincidences remain separately testable, but (2) rules out using any
fixed convolution of the two endpoint fibres to recover the complete offset
profile.

## Connection ledger

| field | exact content |
|---|---|
| source | THM-3710 two-sign full 182-grid endpoint fibres on thirty typed charts |
| target | THM-3713 thirteen-colour deep-offset defect on the same charts |
| tested map | all twelve automorphism-twisted two-input convolutions (1), with and without an intercept |
| preserved if successful | exact coefficients, chart identity, sign input, and affine cyclic covariance |
| destroyed upstream | non-182 endpoints, boundary pairing/ancestry, interval interiors, deep-window mass, and owner semantics |
| needed sidecar | full endpoint ledger plus pairing/owner and window kernel, or the offset profile itself |
| cheapest decisive test | one 390-row modular rank test; `(26,27)` or `(27,28)` is already a rational obstruction |

## Reproduction

```bash
python3 -B 04-computation/lrc_endpoint_offset_affine_convolution_stopping_audit.py
python3 -B -O 04-computation/lrc_endpoint_offset_affine_convolution_stopping_audit.py
```

The normal and optimized transcripts must equal the stored output byte for
byte.  The script pins the THM-3672 engine, its stored transcript, and the
THM-2334 interval referee; it also pins the complete reconstructed endpoint
and offset ledger by digest.  Raw-LF hashes are

```text
script  50fdda055620684a7306539d7a79be43645375493daf244136e2effcd34c8b39
output  5ceeb3afe5d6e1c6b68de92cf776307122fcfc3cac205b6ca8024fccf8572eba
```
