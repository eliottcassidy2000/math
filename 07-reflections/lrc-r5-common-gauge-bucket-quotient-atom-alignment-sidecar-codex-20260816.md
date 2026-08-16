# The endpoint bucket tensor forgets all thirteen common-gauge pullback coordinates

**Status: VERIFIED FINITE-EXACT FOLLOW-UP SIDECAR; the common-gauge source
package is independently accepted in a disjoint clean-room audit.**  The exact companion is
`04-computation/lrc_r5_common_gauge_bucket_quotient_hostile_probe_20260816.py`.
This is a quotient-loss theorem over the certified split field.  Its hostile
pair functions are not asserted to be Boolean, positive, or physically
realizable endpoint factors.  No current, grouped coefficient, row exclusion,
or LRC(14) conclusion follows.

## Verdict

The new common-gauge positive depends essentially on the alignment of the
absolute guard sheets inside each chamber/drift bucket.  The complete 117-row
endpoint bucket table does not determine even one of its thirteen owner
pullbacks.

More strongly, perturbations confined to the 48 owner-active punctured
buckets can change the entire thirteen-vector of pullbacks arbitrarily while
leaving every endpoint bucket sum fixed.  Hence all K4 ranks, Walsh spectra,
weighted-tree factors, and full/restricted bridge values can remain unchanged
while the common-gauge pullback profile is prescribed freely.

The actual endpoint pair function remains important: it lies at one special
atom-level point and gives all thirteen pullbacks nonzero.  The present result
proves that this is new atom-alignment information, not a consequence of the
previous endpoint spectral closure.

## 1. The quotient map

Let

```text
A=F_13 x {L,M,R},             |A|=39,
V=F_p^(A x A),                dim V=1521.              (1)
```

An endpoint pair function is an element `E in V`.  The bucket map remembers
only

```text
B_(C,D,d)(E)
 =sum_(a in F_13) E((a,C),(a+d,D)),                    (2)
```

for the `3*3*13=117` chamber/drift types.  Every bucket contains thirteen
absolute left-sheet positions.  Therefore

```text
dim ker B=117*(13-1)=1404.                             (3)
```

The owner-active punctured part has four chamber corners and twelve nonzero
drifts, so its within-bucket kernel already has dimension

```text
4*12*(13-1)=576.                                      (4)
```

Every endpoint `K4 x F_13` bucket/Walsh/tree statistic factors through `B`.
It is blind to (3) and, in particular, to (4).

## 2. The thirteen ancestry functionals

The common-gauge source measure retains the common offset

```text
c=a-u=b-q.                                            (5)
```

For each owner frequency `k`, its atom-pair functional is

```text
m_k(omega,nu)
 =sum_c M(omega,nu,c) zeta_13^(-kc).                  (6)
```

The endpoint-weight pullback is the pairing

```text
P_k(E)=<m_k,E>.                                       (7)
```

The thirteen vectors `m_k in V*` have rank thirteen.  None is constant on
the fibres of (2): every `m_k`, including `m_0`, varies inside exactly the 72
chamber/drift types realized by the common-gauge ancestry service.

Thus no `P_k` descends to the bucket quotient.  The failure at `k=0` is
important: the missing datum is not merely the twelve primitive Fourier
phases.  The whole regular `C_13` sheet fibre, including its trivial
coordinate, is atom-alignment data.

## 3. Exact surjectivity inside the 48 live buckets

For each owner-active bucket, use the twelve elementary zero-sum pair
functions

```text
e_(a,C;a+d,D)-e_(0,C;d,D),       a=1,...,12.           (8)
```

These are a basis of the 576-dimensional space in (4).  Pairing (8) with the
thirteen `m_k` gives a `13x576` matrix.  Its rank is exactly thirteen.

One explicit nonzero minor uses only the `LL` buckets at drifts `1`, `2`, and
`6`.  Its columns are

```text
(d,a)=(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),
      (2,2),(2,3),(2,4),(2,5),(2,6),(2,7),
      (6,2),                                            (9)
```

and its determinant is

```text
89297302780816677760 mod 572252886246508880869.       (10)
```

Consequently the restriction map

```text
ker(B)|_(48 live buckets) -> F_p^13,
delta E |-> (P_0(delta E),...,P_12(delta E))           (11)
```

is surjective.  For any desired change vector—and in particular the negative
of the actual thirteen-vector—there is a split-field pair-function
perturbation preserving all 117 bucket sums and producing that change.

This is a linear non-implication theorem.  It does not say that the resulting
perturbation is a valid endpoint geometry.

## 4. A two-pair minimal witness

Even a single basis vector is visible to every owner character.  In the
`LL,d=1` bucket, compare absolute left sheets `a=2` and `a=0`:

```text
delta E=e_((2,L),(3,L))-e_((0,L),(1,L)).              (12)
```

Every bucket sum in (2) is unchanged.  Yet the thirteen pullback changes are
all nonzero; their exact tuple is frozen in the companion output.  This is
the minimal support witness: two atom pairs inside one already live bucket.

It pinpoints the first information lost by the endpoint factorization:

```text
atom pair (a,C;b,D)
  -> chamber/drift bucket (C,D,b-a)
  -> K4/Walsh/tree spectra.                            (13)
```

The first arrow forgets precisely the coordinate needed by the common-gauge
owner pullback.

## 5. Tournament and cohomology boundary

The four chambers form the K4 carrier, and choosing orientations on its six
edges gives tournament-section gauge data.  Neither operation sees the
absolute left sheet inside (2).  The missing `F_13` fibre is therefore not a
hidden K4 cycle and not an absolute graph `H^1` class.  It is a separate
regular-torsor sidecar over each chamber/drift bucket.

Equivalently, the useful address object is not just

```text
(C,D,d),
```

but

```text
(a,C; a+d,D; c=a-u=b-q).                              (14)
```

The tournament and tree layers act on the chamber coordinates; the owner
characters act on the retained absolute-sheet/common-offset coordinates.
Collapsing either layer before the pairing destroys different information.

## 6. What the next theorem must preserve

The new actual endpoint pair function passes the atom-level test: its
thirteen values `P_k(E)` are all nonzero.  Equation (11) shows why a proof
cannot replace that function by its bucket table, even if every bucket,
Fourier mode, chart factor, and spanning-tree determinant is retained.

A temporal/Fubini transplant must preserve the exact atom-pair alignment in
(14) and show that the endpoint factors defining `E` live on, or push forward
from, the same chronological THM-2471 relation current.  A synthetic finite
fibre product over the common atom labels always represents the scalar
pairing (7); that formal construction alone cannot supply the missing
timeline or physical endpoint identity.

## Scope boundary

The source common-gauge package is now independently accepted as scoped
finite-exact; this script pins the independent audit hash and semantic digest.
The quotient rank, minor, and two-pair witness here are therefore verified
finite-exact follow-up consequences.  The perturbations lie in a split field
and need not preserve positivity, Boolean structure, pair-level support, or
endpoint-factor provenance.  They prove only that bucket-level data do not
force pullback nonvanishing.  The actual all-owner positive is not refuted.
LRC(14) remains open.

## Reproduction

```text
python -B 04-computation/lrc_r5_common_gauge_bucket_quotient_hostile_probe_20260816.py
python -B -O 04-computation/lrc_r5_common_gauge_bucket_quotient_hostile_probe_20260816.py
```

Normal and optimized outputs are byte-identical.  The semantic digest is

```text
cc535a60a93f5652b1980268a274f7a7280cc722e3c707bdae78467f6d367ddf.
```
