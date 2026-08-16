# The U_full owner gate exposes a full-spectrum K4 x F13 boundary carrier

**Status: FINITE-EXACT ENDPOINT FACTORIZATION; NOT INDEPENDENTLY AUDITED.**
This realizes the actual U_full H-guard address in the frozen endpoint bank
and proves full Walsh-by-drift support for the resulting `K4 x F_13` table.
The endpoint pairing is still Cartesian.  No common-ancestry support
relation, physical current, row exclusion, or LRC(14) consequence is claimed.

## 1. Inheritance pass

The closest proved mechanisms are THM-3479's refined endpoint bank and
THM-3505's rule that a translated sheet coordinate must be retained before
Fourier summation.  The canonical hostile is the checkerboard kernel of the
minimum joint-address gate.  The corrected near miss is that equality of
circle points or interval components is not THM-2471 ancestry.  The
least-used actual sidecar is the interaction between the H guard of speed
one and the owner `in` factor of speed thirteen.

For

```text
ell=(tau,-alpha,-beta,0,0,0,0,alpha,beta),              (1)
```

only the H guard depends on `tau`.  Its unsafe window on the normalized
`91`-circle translates by `7 tau`.

## 2. The 39 guard atoms collapse to 26 owner-boundary atoms

Write

```text
s=91t=7a+r,       a in F_13,       0<=r<7.              (2)
```

The common guard partition has chambers

```text
C_L=[0,1),       C_M=[1,6),       C_R=[6,7).            (3)
```

The previous guard-sheet probe proves that on `(a,C)` the tau-dependent
factor is `g_C(a+tau)`, and that every single and pair guard kernel is
nonzero at every thirteenth-root frequency.

The new observation is set-theoretic.  In `PATTERN_E`, the owner coordinate
has speed `13`, zero character shift, and mode `in`.  Hence every point of
the E-set satisfies

```text
distance(13t,Z)<1/14.
```

Using (2), `13t=a+r/7` modulo one.  Therefore, away from the irrelevant
choice of representatives at two interval endpoints,

```text
r in [0,1/2) union (13/2,7).                            (4)
```

The endpoint engine uses the corresponding half-open representative.  In
particular, the middle chamber is empty before any Fourier transform or
numerical cancellation.  The actual E support uses only

```text
F_13 x {C_L,C_R},                                       (5)
```

so there are `26` active guard atoms.

## 3. Exact atomization of the endpoint bank

For fixed `(alpha,beta)`, remove only the H guard and leave every other
Boolean factor unchanged.  Split the resulting E intervals by (2)--(3).
Let

```text
A_(alpha,beta)(a,C) = the AX endpoint sum on that atom,
B_(alpha,beta)(a,C) = the BY endpoint sum on that atom.  (6)
```

Linearity gives exact reconstruction

```text
AX_(alpha,beta,tau)
  =sum_(a,C) g_C(a+tau) A_(alpha,beta)(a,C),

BY_(alpha,beta,tau)
  =sum_(a,C) g_C(a+tau) B_(alpha,beta)(a,C).             (7)
```

The companion processes `7,107,008` unguarded intervals and `7,108,460`
atom-split intervals.  It checks (7) directly against the parent engine at
all thirteen tau values for four hostile `(alpha,beta)` pairs.  The
reconstructed `13^3` bank has the inherited exact digest

```text
1fabc5cfdbaa1455e10cd6bf9264488133616a7b0ff381623d729b4b4bfa9682,
```

and recovers

```text
q_H  =320618948602619577408,
q_q5 =503604956476841920373,
q_H-q_q5=389266878372286537904                           (8)
```

in the certified split field.

## 4. Pair drift leaves exactly 52 supported types

For a pair of active atoms put

```text
d=a_R-a_L in F_13.                                      (9)
```

The tau sum depends on the common left sheet through its character phase
and otherwise only on `(C_L,C_R,d)`.  Before summing the common sheet there
are `26^2=676` pair addresses.  After the covariance factorization, there
are

```text
2 x 2 x 13=52                                           (10)
```

supported drift types.  Every one is nonzero separately for `q_H`, `q_q5`,
and their difference.  All `65=5*13` types involving the middle chamber
vanish because of (4), although all of their abstract guard kernels are
nonzero.  Thus the zero pattern is support geometry, not Fourier
cancellation.

## 5. The precise K4 carrier and its full spectrum

Encode `C_L=0`, `C_R=1`.  At every drift, the four chamber-pair states are

```text
(0,0), (0,1), (1,0), (1,1) in F_2^2.                  (11)
```

They are the four vertices of the canonical `K4` carrier.  If `B_ij(d)` is
the bridge bucket, its Walsh channels are

```text
W_0 = B_00+B_01+B_10+B_11,
W_L = B_00+B_01-B_10-B_11,
W_R = B_00-B_01+B_10-B_11,
W_X = B_00-B_01-B_10+B_11.                             (12)
```

The last row is the mixed/checkerboard coordinate isolated abstractly by
the minimum joint-address gate.  Exact reduction gives:

1. the `4 x 13` corner table has row rank `4`;
2. every `W_0,W_L,W_R,W_X` is nonzero at all `13` drifts; and
3. the `F_13` Fourier transform of every Walsh row is nonzero at all `13`
   frequencies.

Thus this endpoint table has complete support in the product dual

```text
widehat(F_2^2) x widehat(F_13).                         (13)
```

This is an exact tetrahedral spectral-closure sidecar.  It is not the LRC
`7 x 13` bispectrum: there is no sevenfold carrier here, and the endpoint
pairs have not been restricted to lawful ancestry.

The four states in (11) also do not automatically form a tournament.  They
give the vertices and Walsh/matching representation of `K4`; an intrinsic
pairwise orientation would still have to come from chronology, flux, or
another proved observable.  This is the same representation type as the
depth-two binary-word carrier in THM-3510, not a cross-domain identification.

## 6. Connection and loss contract

| field | exact content |
|---|---|
| source | pre-merge U_full AX and BY endpoint factors with the H guard removed only |
| target | `F_13 x F_2^2`, represented by drift and two boundary-chamber bits |
| map | `s=7a+r`, owner support (4), and `d=a_R-a_L` |
| preserved | the full frozen character bank, both role values, bridge, all guard characters, and endpoint phases |
| destroyed | the common sheet after its Fourier phase is summed; all root/owner/word/source ancestry fields are absent |
| hostile | all abstract middle kernels are nonzero, but actual middle endpoint support is empty |
| needed sidecar | one character-independent THM-2471 common-base/root/owner/word/source/horizon support predicate |
| cheapest next test | restrict the `676` pair addresses by such a predicate before multiplication and retest (12)--(13) |

The decisive advance is that the missing coordinate is no longer an
unspecified parity bit.  The actual endpoint geometry supplies a fully
supported `K4 x F_13` address table.  What remains missing is semantic and
sharply typed: a lawful ancestry support relation on its pair addresses.

## 7. Reproduction and hashes

```text
python -B 04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py
python -B -O 04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py
```

The semantic ledger is

```text
c7b1f92cf7d8d1c09f387ff609b258480f121e39a464b9287332c32c8491d796.
```

The script and output LF SHA-256 hashes are

```text
7593c216294fbf39d14654627620f4ce22ac7c706f4ab9848d43abfcf372e61b
6830c18bd6413fea7c347e65c168664880621679897a630a8982446517f391aa
```

No grouped current, all-unit projector, physical chronology, scalar-row
exclusion, LRC(14), D5 flux map, or Jacobian-conjecture consequence follows.
