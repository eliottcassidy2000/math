# HYP-2231: LRC n=14 Uses a Sin/Cos Even/Odd Wall Carrier

**Status:** OPEN LRC proof-use hypothesis with exact S655 finite scout.

## Claim

The even/odd connection helps LRC `n=14` only after the scalar diagonal is
lifted back to its retained wall carriers.  The user's sine/cosine prompt gives
the type discipline:

```text
sin  -> odd boundary carrier: zeros and active floor constraints
cos  -> even derivative/slack carrier: adjacent response at those walls
cot  -> log-derivative pole ledger: scalar ratio that forgets side channels
```

The visible Goldbach/Lemoine duplicate branch is

```text
7 + 7 = 14
7 + 2*7 = 21.
```

But the LRC floor wall is not carried by the duplicate center `(7,7)`.  It is
carried by the off-diagonal odd complement pairs around `7`:

```text
(1,13), (3,11), (5,9).
```

Each pair sums to the even target `14`.  The odd target `21` is therefore a
scalar shadow of the same central prime, not the active wall certificate.  The
usable LRC object is:

```text
odd complement wall pairs + even derivative slack + C=27 gcd/carry labels.
```

This sharpens HYP-2116's parity tower: at `n=14`, the odd layer is the active
boundary and the even layer is the doubled prime-7 slack/derivative layer.
It is downstream of HYP-2230/S654: HYP-2230 identifies the carry coordinate
`k` in `v=r+27k` as the shared parity and `mod 14` obstruction, while this
hypothesis identifies which odd wall pairs and derivative slack labels must be
retained before that carry coordinate is useful.

## S655 Evidence

`04-computation/lrc_sin_cos_parity_carrier_s655.py` computes AP and `Vstar` at
`n=14`.  Both rows have exact floor

```text
M = 1/14
```

and the same best times:

```text
1/14, 3/14, 5/14, 9/14, 11/14, 13/14.
```

At those times, the active odd complement pairs are exactly

```text
(1,13), (5,9), (3,11), (3,11), (5,9), (1,13).
```

At AP's `t=1/14`, the depth minima are

```text
v2=0: 1/14 = 1*delta
v2=1: 1/7  = 2*delta
v2=2: 1/7  = 2*delta
v2=3: 3/7  = 6*delta.
```

The active slopes at `t=1/14` are `[1,-13]`, and the other floor times carry
the corresponding odd-pair slopes `[5,-9]`, `[3,-11]`, `[-3,11]`, `[-5,9]`,
and `[-1,13]`.

The AP and `Vstar` rows also keep the same `C=27=3^3` gcd-shell counts:

```text
{1: 9, 3: 3, 9: 1}.
```

The AP single-swap census has

```text
rows = 195
below = 0
tight = 1
loose = 194.
```

The only tight non-AP single swap is the known `Vstar` move:

```text
12 -> 24.
```

It keeps the same gcd-shell counts and the same odd complement wall pairs.  No
single-swap row falls below the floor.

S655 also checks the cotangent half-product

```text
prod_{k=1}^{(C-1)/2} cot(pi*k/C)^2 = 1/C
```

for odd `C=15,21,27,29`.  This is the warning label: the scalar is clean, but
for composite `C=27` it erases the gcd-shell side channel that AP and Vstar
need.

## Repo Interpretation

This refines the S642/S643/S644 even/odd work.  Those sessions showed that
same-pair Goldbach/Lemoine projections have an invertible pair carrier and that
the duplicate branch at `p=7` gives the visible `(14,21)` shadow.

S655 says that this shadow is too coarse for LRC.  The active floor constraints
are not the diagonal pair `(7,7)`, but the odd boundary pairs symmetric around
`7`.  This matches the S653 Basel lesson: scalar identities become useful only
after the packet/product side channel is retained.

For LRC `n=14`, the side channels are:

- odd boundary sine carrier: active wall pairs;
- even derivative/cosine slack: local response at the boundary;
- pair-sum pinch oracle: the THM-369/THM-401 floor witness clock;
- `C=27` gcd shell carrier: the HYP-2222 fixed pocket;
- owner/carry conservativity: the HYP-2167/S611/S612 lift obligation.

The cotangent scalar belongs in this list only as a diagnostic.  It says
"there is a product/log-derivative ledger," but it does not decide which gcd
stratum owns a wall or which lift/carry label is conservative.

## Proposed LRC Move

The next proof target is a no-leak lemma:

```text
If an n=14 row preserves the odd complement wall pairs and the C=27 gcd-shell
carrier, then pair-pinch and owner/carry labels force AP, Vstar, or strict
looseness above 1/14.
```

This is narrower than classifying all rows.  It says the even/odd bridge should
be used as a wall-pair filter, while the known `C=27` machinery handles the
remaining composite-shell lift.

The immediate computational extension is to move beyond single swaps:

```text
enumerate fixed odd-wall-pair rows -> refine by gcd-shell mass -> apply
S611/S612 carry visibility -> search for any residual floor candidate.
```

If the residual set collapses to AP and `Vstar`, HYP-2231 becomes a concrete
route to the LRC `n=14` conservativity seam.

## Assumption Challenge

The tournament vertices in S655 are not runners.  The useful vertices are proof
carriers:

```text
odd boundary, even derivative slack, cotangent pole ledger, pair-sum pinch,
C=27 gcd shell, owner/carry conservativity, raw parity numerology.
```

The quotient preserves floor/tight evidence, odd wall-pair labels, and shell
carrier data.  It destroys raw runner order and exact phase data, which must be
reattached only when the carry/lift proof needs it.

The challenged assumption is that the even/odd bridge should identify `14` and
`21` directly.  The better bridge is:

```text
14 = odd complement wall sum
21 = diagonal scalar shadow
27 = composite carry/gcd clock.
```

Tournament fingerprints for the proof-carrier majority graph:

```text
score_hist={0:1,1:1,2:1,3:1,5:3}
directed_3cycles=1
scc_sizes=[3,1,1,1,1]
hamiltonian_paths=3
```

The top tied carriers are `pair_sum_pinch_oracle`,
`odd_boundary_sin_carrier`, and `C27_gcd_shell_carrier`.

**See also:** HYP-2230, HYP-2229, HYP-2228, HYP-2222, HYP-2219, HYP-2218,
HYP-2167, HYP-2164, THM-401, THM-369,
`04-computation/lrc_sin_cos_parity_carrier_s655.py`,
`05-knowledge/results/lrc_sin_cos_parity_carrier_s655.out`,
`07-reflections/lrc14-sin-cos-even-odd-wall-carrier-s655.md`.
