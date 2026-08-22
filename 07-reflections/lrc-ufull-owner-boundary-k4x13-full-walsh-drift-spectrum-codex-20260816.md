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

Using (2), `13t=a+r/7` modulo one.  Therefore, with the endpoint engine's
exact half-open representative,

```text
r in [0,1/2) union [13/2,7).                            (4)
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

Three tempting equality restrictions fail in the inherited bucket
normalization.  Their `(q_H,q_q5,bridge)` values are

```text
same sheet:
  (118992333014385604245,
   366746771947440878083,
   324498447313453607031),

same chamber:
  (237441608630540813143,
   468104475565795325224,
   341590019311254368788),

same guard atom:
  (319186094366758452423,
   268974005044605905796,
    50212089322152546627).
```

None has bridge `389266878372286537904`.  This refutes those three named
restrictions as reconstructions of the frozen Cartesian bridge.  It does not
rule out a differently weighted or genuinely ancestry-typed coupling.

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

Under the independently audited relation-role typing, these weighted role
differences remain in the exact coboundary image `B^1`, not in a nonzero
`H^1` flux class.  Full bucket support, Fourier support, or matrix rank does
not promote that `B^1` carrier to a physical current or an ancestry relation.

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
| map | `s=7a+r`, owner support (4), `d=a_R-a_L`, and—after one gauge choice—the equivariant sheet-label map `u=a+c` |
| preserved | the full frozen character bank, both role values, bridge, all guard characters, and endpoint phases |
| destroyed | the common sheet after its Fourier phase is summed; physical root nodes and all owner/word/source/horizon ancestry fields are absent |
| hostile | all abstract middle kernels are nonzero, but actual middle endpoint support is empty |
| needed sidecar | one character-independent THM-2471 common-base/owner/word/source/horizon support predicate using one common root gauge on both legs |
| cheapest next test | restrict the `676` pair addresses by such a predicate before multiplication and retest (12)--(13) |

The decisive advance is that the missing coordinate is no longer an
unspecified parity bit.  The actual endpoint geometry supplies a fully
supported `K4 x F_13` address table.  What remains missing is semantic and
sharply typed: a lawful ancestry support relation on its pair addresses.

There is also an exact reduction of the root-label debt.  The guard pullback
is `g_C(a+tau)`.  In THM-2471's chart,

```text
phi(y,u)+tau/13=phi(y,u+tau),                         (14)
```

and THM-2538 uses the same physical translation `x -> x+tau/13`.  Thus both
sheet and root labels carry the covariant `+tau` action.  Every equivariant
bijection between the two thirteen-element torsors has the form

```text
u=a+c,                 c in F_13.                     (15)
```

There are thirteen choices and no canonical `c`.  In the Fourier convention
of the guard probe, changing `c` multiplies a primitive frequency-`k`
coefficient by `zeta^(kc)`, so it cannot change nonvanishing.  If one common
`c` is used on both endpoint legs, then

```text
u_R-u_L=a_R-a_L=d.                                   (16)
```

Independent gauges would instead add `c_R-c_L`; the desired common ancestry
base must therefore carry one common gauge.  Equations (15)--(16) align only
the thirteen sheet labels.  They are not a bijection from the `39` guard
atoms (or `26` active atoms) to roots—the chamber bit remains a sidecar—and
they do not identify any guard sheet with a physical THM-2471 root.  Absolute
root-label origin is therefore no longer a nonvanishing obstruction;
common-base support, common gauge, horizons, and source-versus-arrival typing
remain the first missing map.

## 7. Independent branch-transplant audit

The verdict is negative and exact: **no currently proved object supplies a
lawful character-independent common-ancestry support predicate on the actual
U_full guard atoms or endpoint pairs.**  The candidate ledger is:

| proved object | what it actually supplies | why it does not supply the predicate |
|---|---|---|
| `THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary.md` | a genuine finite Boolean stalk with linked nodes `w_u`, `X_(u,a)`, `Y_(q,e')` and base/root/sheet data | no map from an actual U_full endpoint cell to those stalk coordinates is constructed |
| `THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy.md` | a signed ANOVA defect and its cut transform | a signed contraction is not a Boolean support relation, and the theorem explicitly leaves common owner/arrival/deep ancestry open |
| `THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary.md` | exact recovery *if* predecessor and arrival fields are already on one ancestry base | common ancestry is a hypothesis of the recovery formula, not an output of it |
| `THM-2594-realized-theta-slaved-contraction-at-the-r5-window.md` | one genuine linked-node common-base contraction on the canonical r=5 row | it is a different interval table and explicitly is not the THM-2334 U_full endpoint bank or a generic transplant |
| `THM-3479-literal-half-twist-relation-current-two-transplant-certificate.md` | the actual U_full refined endpoint bank | its worker returns separately summed `AX` and `BY` scalars and no shared ancestry key |
| `THM-3505-r5-projective-slope-pencil-and-anchored-haar-basis.md` | a formal binary-label-to-Haar-rectangle map | that map explicitly loses circle-cell measure, factor lineage, common-base identity, sheets, horizons, and the address |

The failed equality hostiles are also type-revealing.  The prior exact atom
probe gives bridge `540653486701996040250` for either equal cut segment or
equal maximal cyclic component and `167726070588785644466` for equal circle
point, versus the frozen bridge in (8).  Moreover, an interval-component key
is indexed by the character-dependent final E decomposition.  Point equality
retains one circle coordinate, whereas THM-2471 ancestry links generally
different nodes over one base; MISTAKE-293 forbids identifying those notions.
The new same-sheet/same-atom hostiles above show that adding the honest guard
address does not repair the type mismatch.

The corresponding transplant contract is therefore:

| field | audit verdict |
|---|---|
| source | actual pre-merge U_full left/right endpoint cells, refined by `(a,C)` |
| target | a THM-2471-style finite Boolean ancestry fibre carrying both endpoint factors |
| available map | guard projection `s -> (a,C)`, pair drift `d=a_R-a_L`, and the noncanonical equivariant label torsor `u=a+c` |
| preserved predicate | H-guard safety and its exact tau covariance; the full Cartesian endpoint character bank |
| information destroyed | physical root realization, owner sheet, word sheet, source sheet, horizons, chronology, and the lawful left/right pairing |
| missing sidecar | a character-independent relation `R subset Omega_L x Omega_R` whose records carry those fields and one common root gauge before either endpoint sum |
| cheapest decisive test | on the bridge classes only, restrict the `26^2=676` active atom pairs by `R`, insert both factors before marginalization, and require inverse DFT to recover (8); only then test K4 factors or chronology |

This audit follows MISTAKE-281/283/293/300/313: nonzero Fourier factors,
matching row/word labels, one outer base, shared selector vocabulary, and
matching partial constructors are each weaker than an exhibited fibre
product with every factor present.

## 8. Reproduction and hashes

```text
python -B 04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py
python -B -O 04-computation/lrc_ufull_guard_sheet_drift_bucket_bridge_probe_20260816.py
```

The semantic ledger is

```text
0b31a992ba23cecd05f28ae353133531f41cc6d84a4a935c34a12d77fd3db590.
```

The script and output LF SHA-256 hashes are

```text
a1d4b667812949001fc863ba881ff7409bbae3c568a6bf7bc24c9dc88b2766b1
ed1ef68e733684d7a015af314bc196fdc71e95b2a5a8d1d98d933553163a2e95
```

No grouped current, all-unit projector, physical chronology, scalar-row
exclusion, LRC(14), D5 flux map, or Jacobian-conjecture consequence follows.
