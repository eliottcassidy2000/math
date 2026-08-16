# The fourth, sixth, eighth, and ninth source directions all miss the same endpoint connection

**Status: INDEPENDENT AUDIT ACCEPTS THE SCOPED FINITE-EXACT PACKAGE;
promotion recommended only as a split-field representation theorem.**  The
independent companion is
`04-computation/lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py`.
It imports neither submitted candidate script and proves no U_full ancestry,
physical current, absolute graph `H^1`, bispectrum, scalar-row exclusion, or
LRC(14) conclusion.

## Verdict

The two incoming source-side results are sound in their stated scope.  The
fixed absolute-root marginal supplies exactly one new direction beyond the
theta-slaved marginal, while the correctly typed root-difference and owner
refinements raise the source ranks further:

```text
bank                                 rank
-----------------------------------------
theta-slaved                           3
fixed absolute root                    3
their union                            4
root difference s=u-q                  6
theta-resolved root difference         8
owner k=0 or k=1                       6
owner k=0,1 union                      9.
```

Every one of these banks nevertheless has the same sharp obstruction.  No
frequency-independent channel map followed by one common `C_13` operator
maps its spectral curve projectively to the four U_full Walsh channels.
Every formal solution annihilates the source curve.

## 1. Clean-room reconstruction

The audit reruns the hash-pinned THM-2594 primary common-base construction,
then forms all source marginals directly from the joint Boolean table

```text
N(u,q,ell,theta).
```

It separately rebuilds the endpoint target from the independent THM-3514
atom engine, not from either candidate probe.  The resulting source table
digests are

```text
slaved:          666200f5e79e48c02e7beaeec90ebfedccb708bd1e383d1a01997d0abf5073f7
absolute root:   9353c6e645b570c659a5d42491a1dd7c79825d4bd338a476aabbd0afe10e69f2
root difference: ad63d008da3dfa8525db7ef4e724cf6e3a971b3f6efa3d6da9791fbd8003ed91.
```

Their raw supports are respectively `12`, `18`, and `72` of the `91`
septimal-by-drift cells, but every one has all `91` Fourier coordinates
nonzero.  This independently confirms that spectral support, row rank, and
connection compatibility are genuinely distinct invariants.

## 2. The hostile system retains the multipliers

For source vectors `X_k`, target vectors `Y_k`, a fixed channel map `M`, and
one common circulant with Fourier multipliers `lambda_k`, the required
relation is

```text
M X_k = lambda_k Y_k,       k in F_13.                 (1)
```

Rather than trusting the candidates' eliminated wedge equations, the audit
keeps all thirteen `lambda_k` as unknowns in one augmented linear system.
For every tested bank, the projection of the full nullspace to the thirteen
multiplier coordinates has rank zero.  Its `M` components also annihilate
every source vector.  A separately coded wedge system gives the identical
annihilator excess zero.

The complete certificate ledger is

```text
bank             source  augmented  augmented  lambda projection
                 rank    rank       nullity    rank
---------------------------------------------------------------
slaved 7            3       25         16             0
absolute 7          3       25         16             0
union 14            4       29         40             0
difference 7        6       37          4             0
resolved 21         8       45         52             0
owner k=0           6       37          4             0
owner k=1           6       37          4             0
owner difference    6       37          4             0
owner union 14      9       49         20             0.
```

Thus the rank-four source boundary was necessary but not sufficient.  The
rank-six, rank-eight, and rank-nine enrichments fail for the same structural
reason: none carries the endpoint's common connection.

## 3. Exhausted gauges and marked sidecars

The audit independently exhausts:

- all twelve nonzero dilations of the `F_13` torsor for the rank `4`, `6`,
  `8`, and `9` banks;
- all eight sections of the folded `C_7` antipodal marking;
- all `56` named absolute-root sidecars;
- all `127` nonempty binary sums of absolute-root rows; and
- all `7^4=2401` root-difference allocations.

Every dilation has the same annihilator-only certificate.  Every folded
four-row root-difference section has augmented rank `29` and nullity zero.
All `56` named sidecars have rank four.  Of the binary sidecars, exactly the
all-rows mask `127` remains rank three; the other `126` have rank four, and
all `127` still fail.  None of the `2401` allocations has an exact or
amplitude/shift-gauge common kernel.

## 4. Why the tournament does not repair the map

The exact Paley fold is a weighted looped symmetric quotient whose
off-diagonal support is bidirected `K4`.  Choosing one representative from
each antipodal pair produces eight tournament sections, but those sections
are gauge choices, not extra response channels.

The independent affine-code calculation gives

```text
C_1(K4;F_2) = Z_1(K4;F_2) direct-sum G_section,         (2)
```

where `G_section` is the translated three-bit section code.  Hence the
tournament marking is a complement to the cycle space, not a hidden
absolute `H^1` class.  Its skew tournament matrices may be invertible as
left mixers, but invertible left mixing cannot manufacture the missing
source-to-target connection in (1).

This resolves the apparent “tournament of size four” shortcut precisely:
the carrier has four vertices and six pair observables, while the exact
quotient has ties in both directions, endpoint tree responses lie in
`B^1`, and the eight orientations merely choose sections of the fold.

## 5. The strongest survivor

The failed implication is now sharply located.  Each marginal discards the
support incidence that could identify source cells with endpoint chambers.
The smallest surviving source object is

```text
(u,s,ell,theta),       q=u-s,                           (3)
```

and the smallest endpoint object is

```text
(a,d,C,D),             b=a+d.                          (4)
```

The only typed bridge currently justified is the affine skeleton

```text
u=a+c,       s=+/-d,                                   (5)
```

together with the translated owner phase.  No theorem yet identifies the
source word/cell/deep support with the endpoint chamber predicates `C,D`.
The next decisive experiment must therefore live on the joint address
relation

```text
(u,s,ell,theta; a,d,C,D)                               (6)
```

before summing over `u`, `a`, `theta`, or chamber labels.  A nonlinear
Boolean support map, a cospan with a retained owner switch, or a typed
connection with independently derived holonomy remains possible.  Another
post-marginal fixed linear fit is excluded by the present audit.

## Scope boundary

This is a finite exact theorem about representations over the certified
split field.  It validates the rank and no-go ledgers of the two candidate
probes and identifies their common stopping mechanism.  It does not show
that the joint support relation (6) exists, does not transport THM-2471
ancestry, and does not create a physical word current or nonzero absolute
graph cohomology class.  LRC(14) remains open.

## Reproduction

Run

```text
python -B 04-computation/lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_common_base_connection_obstruction_independent_audit_20260816.py
```

Normal and optimized outputs are byte-identical.  The pinned semantic digest
is

```text
544b07a84c6806ea63c48f5227b78d74844f466dd0b446b6a320ea8560238895.
```
