# Independent hostile audit: the current-leg inverse-owner digit is nonredundant

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED CURRENT-LEG
`r_owner` SIDECAR.  LRC(14) remains OPEN.**  The submitted implementation is
not imported.  The audit rebuilds the thirteen digit profiles from the
independently audited Boolean-square source and common-base endpoint object.

## Radix type and exact reconstruction

At THM-2471 depth `d=13^5`, write

```text
X_(u,a)=(w_u+a)/d,       a=r_owner+13c,
r_owner=a mod 13.
```

Multiplication by `13^4` gives

```text
13^4 X_(u,a)=(w_u+r_owner)/13 mod 1.
```

Thus `r_owner` is the final current-leg inverse digit immediately before the
collision-root coordinate `w_u`.  The exhaustive radix census has `13^4`
values of `a` in each digit.  Under `a -> 13^5-1-a`, the digit and upper word
map to `r -> 12-r` and `c -> 13^4-1-c`.

The independent source construction first builds the packet
`1_Q P_169 1_E`, folds it by `13^4`, pulls the thirteen digit windows, and
then pulls the thirteen collision-root windows.  It does not import the
candidate or relabel a post-integrated table.  On every one of the `33` joint
open intervals and for each root `u`,

```text
sum_r U_(u,r)(y)=U_u(y)
```

holds pointwise.  The decomposition has `72` upper-fold pieces, digit-window
piece counts

```text
(12,5,5,6,5,5,7,5,5,6,5,5,12),
```

and `34` common boundaries.  Its exact digest
`9e2a61ea8ec00830aa93ec4daa42574a9179f5aa80a09e3388dd9af0964f4cb1`
agrees with the submitted candidate.

## Reflection and square typing

The sheet identity

```text
1-X_(u,a)(y)=X_(12-u,13^5-1-a)(1-y)
```

induces

```text
(u,r_owner,state) -> (12-u,12-r_owner,state XOR 2).
```

The source-state rule is the image of each active root under `u -> 12-u`:
`{0}` goes to `{12}` and `{0,6}` goes to `{6,12}`.  It is not three-root cut
complementation, which is the distinct square involution `state XOR 3`.
Reflection preserves singleton/doubleton multiplicity and reverses only the
left/right bit.  The audit checks this on every joint profile interval and
recovers the candidate's palindromic branch-mass vector.

This correction matters: the stronger provisional statement that reflection
complements the spine cut fails already at the left endpoint.  The exact
root-image rule is both weaker and sufficient for the claimed XOR-2 profile
symmetry.

## Common-base, phase, and marginal gates

On each endpoint segment, the audit retains `r_owner` in the selected left
root-service value before multiplying by the selected right value and the
exact `Q(13y)` harmonic jump.  Same-root products vanish pointwise.  The
literal guard is restored directly at `(alpha,beta,tau)=(0,0,0),(1,0,6),`
and `(6,6,12)`.

The response phase is `zeta_13^beta`.  Inversion uses

```text
13^(-3) zeta_13^(-(alpha+tau*t)).
```

Summing the digit coordinate recovers the independently audited Boolean
square twice: every one of the `2197` phased character rows before inversion,
and every `4 x 13` table entry after inversion.  These two marginal checks are
the arrival/source-time, Fubini, sign, and normalization gates for this scoped
one-base object.

## Nonredundancy and flat hostile

Let `T(state,r,t)` be the inverse tensor.  The flattening ranks in axis order
`(state,r_owner,relation)` are

```text
actual:               (4,4,6)
digit-flat hostile:   (4,1,4).
```

The hostile is the unique uniform digit lift of the same audited square
marginal.  After sequentially centering state, digit, and relation, the ranks
are

```text
actual pure interaction:   (3,4,6)
flat pure interaction:     (0,0,0).
```

The actual pure tensor has all `3*12*12=432` allowable triple Fourier modes;
the flat hostile has none.  The full actual tensor has all `676` coefficients,
whereas the flat lift has only `52` lower-order coefficients.

At fixed relation `(1,0,6)`, the actual `4 x 13` matrix has rank `4`, all
`52/52` entries are nonzero, and its pure three-way slice has rank `3`.  The
flat matrix has rank `1`, also has full entry support, and its centered slice
is zero.  Entry support therefore does not establish nonredundancy; the rank
and centered hostile do.

All source, gamma, tensor, flat, spectrum, ANOVA, and fixed-slice digests agree
with the submitted candidate.  Normal and optimized transcripts are
byte-identical after line-ending normalization.  The independent semantic
digest is
`7063720d0e0e4847ce752102de83274ea47d7740fc435a64bae425dbd7100121`.

## Coordinate ledger and next gate

The audited thirteen-valued coordinates remain distinct:

| coordinate | stage | meaning |
|---|---|---|
| `b_source=floor(n/13)` | before source `P_169` transfer | high source-time inverse digit |
| `s=u-q` | source root pair | collision-root difference |
| `r_owner=a mod 13` | current THM-2471 `P_(13^5)` leg | final inverse-owner digit |

The present audit proves no intertwiner among them.  The matching rank-six and
`432/432` signatures for `b_source` and `r_owner` are evidence for a possible
shared response ceiling, not evidence that the coordinates coincide.

The cheapest lawful next experiment is a genuine pre-integration joint table
retaining `b_source` and `r_owner` over the Boolean square.  It must recover
both independently audited one-coordinate tensors as exact marginals.  A
formal tensor product after either marginal would not preserve ancestry or
source/current timing.

## Connection contract

| field | exact answer |
|---|---|
| source | THM-2471 current-sheet inverse index `a` over one audited owner base |
| target | `V_4(state) x F_13(r_owner) x F_13(relation)` |
| map | factor `13^5=13*13^4`, retain `a mod 13`, pull collision root, then couple endpoint factors and invert |
| preserved | current digit, collision root, common coordinate, literal guard order, endpoint phase, Boolean state |
| destroyed | four upper current digits, source-time digit, root difference, pointed tail, exact address, chronology |
| mandatory marginal | sum over `r_owner` recovers the audited square before and after inversion |
| decisive hostile | digit-flat lift with the same marginal but branch rank one and zero pure interaction |
| survivor | actual ranks `(4,4,6)` / centered `(3,4,6)` and `432/432` triple modes |

This is not the source-time `P_169` digit, source root difference, a deep
label, full inverse ancestry, THM-2334 exact address, `U_clock` chronology,
physical current, row theorem, row exclusion, or LRC(14).

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_independent_audit_20260816.py
```

LF SHA-256 for script/output:
`4b0bd05ffa6195ff484433329e334d771bc27e7cd380136b50b45e7248bb98ba` /
`d1b7aef51b58c28afc40ec0d08319bdaf0fae7a3f5e681dddab3d0f1d2f1a543`.
