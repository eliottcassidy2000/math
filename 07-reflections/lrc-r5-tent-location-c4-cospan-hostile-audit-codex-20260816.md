# The tent odd part specifies a marked C4 interface, not a current

**Status: FINITE-EXACT TYPING SIDECAR; USEFUL FORMAL COSPAN CANDIDATE;
BLOCKED BEFORE PHYSICAL CLOSURE AND COEFFICIENT TRANSPORT.  NO RESPONSE
AMPLITUDE, CLOCK, CHRONOLOGY, PHYSICAL CURRENT, D5 CLASS, ROW EXCLUSION, OR
LRC(14) CONCLUSION.**

The accepted exceptional-location section has the exact arc-reversal split

```text
h       = (12,12,9,3,0,0),
h_even  = (12,12,6,6,0,0),
h_odd   = ( 0, 0,3,-3,0,0) = 3*j_mid,                (1)
j_mid   = e_(B->C)-e_(C->B).
```

Equation (1) is useful, but its type is decisive: `h_i` is the `r0` location
of a source-proportionality exception on pointed arc `i`.  It is not a
diagonal response coefficient, an edge flow, or an `F13` word-current.  The
same six-coordinate tuple can be written in all those spaces only after a
map has been supplied.

The exact verdict is:

> (1) gives a canonical middle-coordinate **candidate interface**.  If one
> formally adds the missing Boolean edge and marks the normalized odd basis,
> the coefficient `3` has a unique divergence-free `C4` completion and a
> mod-13 seam `12=-1`.  The present data construct neither the closure edge
> nor the map from exception location to response/current coefficient.

Thus the connection is worth retaining as a pre-registered cospan, but it is
blocked before coefficient transport and cannot be promoted to a current.

## Inheritance and the corrected near miss

The closest exact source is the independently audited pointed-six carrier:

```text
A <-> B <-> C <-> D,
```

with `A=0,B=1,C=3,D=2`, and the alternating incidence tree `P7`.  Its six
directed-arc response lines are injective, while the static path has
`H1=0`.  The closest closure template is the existing formal addition of
`D--A`, yielding the consistently oriented cycle

```text
A->B->C->D->A.                                        (2)
```

The canonical hostile is
[THM-3496](../01-canon/theorems/THM-3496-marked-graph-kummer-degree-square-and-finite-coefficient-frobenius-flux-extinction.md):
its graph-to-Kummer map exists only after orientation, generator,
uniformizer, meridian, and exponent-one markings, and it proves no additive
map to the characteristic-zero response module.  MISTAKE-417 is the nearby
warning that formal Fourier support or a shared coordinate tuple does not
create a coupled observable.

The corrected near miss is sharper than “location is not amplitude.”  There
are two natural normalizations of `j_mid`, and they give different exact
seams.

## 1. The normalization invoice

Let the directed arcs themselves be literal oriented chains.  Relative to
the chosen orientation `B->C`,

```text
e_(B->C) |->  [BC],
e_(C->B) |-> -[BC].                                  (3)
```

Therefore the literal chain realization sends

```text
j_mid |-> 2[BC],             3*j_mid |-> 6[BC].      (4)
```

This agrees with the earlier exact boundary

```text
partial(j_mid)=2(C-B).
```

By contrast, calling the coefficient of the odd-basis vector `j_mid` the
single-edge coefficient uses the normalized marking

```text
nu(j_mid)=[BC].                                       (5)
```

Over `F13`, (5) differs from (3) by the unit `2`; it is legitimate after it
is declared.  It is not available in characteristic two.  Consequently:

| convention | middle coefficient | unique constant `C4` cycle | mod-13 cochain seam |
|---|---:|---|---:|
| marked odd-basis normalization (5) | `3` | `(3,3,3,3)` | `12=-1` |
| literal directed-arc chain (3) | `6` | `(6,6,6,6)` | `24=11=-2` |

Both seams are nonzero, but their scalar differs.  Since THM-3496's
exponent-one normalization is load-bearing, the claimed seam `-1` must be
labelled as the first convention, not as an unmarked consequence of (1).

## 2. What the unique completion does and does not extend

For the oriented cycle (2), evaluation on the middle edge is an isomorphism

```text
ev_BC: Z_1(C4;R) -> R,                                (6)
```

because the four divergence equations force all edge coefficients equal.
Thus, in the normalized convention, middle coefficient `3` has the unique
formal completion

```text
J_formal=3(1,1,1,1).                                  (7)
```

However, restricting (7) back to the three path edges gives

```text
res_P4(J_formal)=(3,3,3),
```

whereas (1) supplies only

```text
odd-location coefficients=(0,3,0).                   (8)
```

Hence (7) is not an extension of the whole odd location section.  It is the
unique divergence-free completion of one selected coordinate after the two
outer coordinates have been discarded.  The even location baseline in (1)
is discarded as well.

There is also no typed edge `D--A` in the current tensor.  The four owner
states form a Boolean square abstractly, but states `D` and `A` share no
retained root, and no same-temporal-copy `U_clock` return has been built.
Thus (2), (6), and (7) remain formal before coefficient transport is even
asked for.

## 3. Chain current and H1 cochain are different types

The tuple (7) is divergence-free as a chain, so it represents

```text
3*[C4] in H_1(C4;Z).                                  (9)
```

After separately marking the edge-coordinate pairing, the same tuple may be
read as a one-cochain.  Its integral period is

```text
s_Z=3+3+3+3=12.                                      (10)
```

The lawful coefficient comparison is the integral base-change cospan

```text
H^1(C4;Z) -> H^1(C4;F13),
H^1(C4;Z) -> H^1(C4;F2).                             (11)
```

Applied to (10), it gives

```text
s_F13=12=-1 != 0,
s_F2 = 0.                                             (12)
```

The Boolean support cochain `(1,1,1,1)` is indeed exact:

```text
(1,1,1,1)=delta(0,1,0,1) in C^1(C4;F2).              (13)
```

But as a **chain**, `(1,1,1,1)` is the nonzero generator of
`H_1(C4;F2)`; a graph has no two-cells whose boundary could kill it.  Thus
“the Boolean all-one mask is exact” is correct only in the cochain/`H^1`
sense.

Nor is (12) a reduction map from the mod-13 class to the mod-2 class.  There
is no unital ring map `F13 -> F2`: the relation `13*1=0` would map to
`13*1=1`.  The two classes are separate specializations of the integral
class (10).  Arc-reversal `+/-` eigenspaces also merge in characteristic two,
so the odd decomposition itself has no mod-2 descent.

The actual finite response table lives in a different split field of
characteristic

```text
755373809845391722745761,
```

and the characteristic-zero response is additively torsion-free.  There is
no nonzero additive coefficient map from `F13` to either response module.
In particular, the notation `12=-1` belongs to the formal `F13` cochain; in
the split response field, the integer images of `12` and `-1` are different.

## 4. The exact response hostile blocks scalar transport

The support coincidence is real but does not transport the coefficient.  I
independently rebuilt the accepted diagonal profile and formed the actual
middle-arc odd response

```text
d_r0(r1)=k_(B->C)(r0,r1)-k_(C->B)(r0,r1)             (14)
```

on the two reflected exceptional rows `r0=3,9` and the ten non-root digits
`r1 in {1,2,3,4,5,7,8,9,10,11}`.  All twenty entries are nonzero, but each
row contains eight distinct response amplitudes, and the two row vectors
have rank two.  They obey the exact reflection

```text
d_9(r1)=-d_3(12-r1).                                 (15)
```

The reconstructed profile and this twenty-entry ledger have SHA-256

```text
d1c7e561538ac6abb7c631c642d547aea7ca8dcfa33296e849a22a5061d2f595
7ee80b176cc42b9c900605645de1da8aff803d12e2a7b4c91df2965e3ad87ce3,
```

and the first values are
`359051715977130886352290` and `728944175358539724970630`.

Thus `h_odd` correctly singles out the support of the middle exception, but
the location coefficient `3` cannot map through any `r1`-blind scalar
transport to the actual response.  This is a direct hostile, independent of
the characteristic obstruction above.  An `r1`-dependent section is not
refuted, but retaining and typing that section is precisely additional data
not present in `h`; no such map is currently constructed.

## 5. Incoming exact sidecars sharpen, but do not repair, the map

The concurrent [K4 XOR closure atlas](k4-xor-reversal-mask-and-characteristic13-closure-boundary-codex-20260816.md)
gives an independent Boolean explanation of (13).  After adding a
**chosen** transitive tournament gauge, the hypothetical `DA` edge changes
the path reversal mask into the cut

```text
AB+BC+CD+DA=delta({A,C}) in C^1(K4;F2).
```

Multiplying its four entries by three has no effect over `F2`, while over
`F13` it multiplies the oriented seam `4` to `12=-1`.  This is exact
compatibility between two base changes of the same marked integral cochain.
It does not manufacture `DA`: the tournament orientation is extra gauge
data and the XOR map forgets amplitudes, roots, address, and time.

The concurrent [common ten-space extraction](lrc-r5-third-digit-common-ten-quotient-not-transfer-codex-20260816.md)
tests the most tempting repair of the `r1` dependence.  On the next-digit
pointed bundle, every child surjects onto the same quotient `Q` of dimension
ten, but there is no common ten-dimensional section and no child-to-child
transfer.  More pointedly, arc reversal does not act on `Q`: adjoining the
reversed parent space costs exactly two dimensions, both localized in the
middle arc pair.  Passing to the repaired eight-dimensional quotient erases
that middle orientation defect.

Therefore the new quotient cannot be used as a canonical target for
`h_odd`.  It either retains noncanonical child-dependent lift data or kills
the two coordinates where the odd middle contrast lives.  This does not
exclude every future address-dependent transport, but it closes the obvious
quotient shortcut.

## 6. The useful broken cospan

Keeping the semantic coefficient sorts explicit gives

```text
Loc_exception^-
      | coefficient of j_mid
      v
Z_address  - - - tau missing - - ->  Z_current  <-ev_BC-  Z_1(C4;Z_current).
                                                               (16)
```

If one forgets units and identifies the two copies of `Z`, (16) becomes the
formal shadow cospan selecting (7).  That shadow preserves:

- middle support under arc reversal;
- the marked middle coefficient;
- divergence zero after formal completion; and
- the integral seam and its two separate base changes.

It destroys:

- the meaning of `h` as an exception-location section;
- the even baseline and outer location coordinates;
- source cells and endpoint response values; and
- clock, temporal copy, arrival, and physical-current data.

Four realization gates remain, in order:

1. construct the same-copy `D--A` clock/address edge;
2. construct `tau`, taking address-location contrast to a response/current
   coefficient with a declared normalization;
3. construct a chain-to-cochain constitutive marking; and
4. map this owner-state `C4` semantically to THM-3496's seven-chart cycle.

Only after all four gates would the normalized candidate enter the existing
marked D5 cospan

```text
H^1(C4;F13) --seam--> F13 <--seam-- H^1(C7;F13)
                                      -> Kummer H^1.    (17)
```

Under the normalization (5), its pre-registered image in (17) would be

```text
-kappa_lambda.                                        (18)
```

Equation (18) is a conditional prediction, not a constructed class.  Even a
future proof of all four gates would inherit THM-3496's theorem that this
Kummer exponent is not an additive Hamiltonian response flux.

## Connection contract

| field | exact answer |
|---|---|
| source | six exception locations `h`, not response amplitudes |
| first map | arc-reversal half-difference in the location module |
| exact image | `3*j_mid` |
| formal completion | middle evaluation inverse on a hypothetical oriented `C4` |
| normalized cycle | `3(1,1,1,1)` |
| literal-arc alternative | `6(1,1,1,1)` |
| normalized seam | integral `12`, mod-13 `-1`, mod-2 zero |
| Boolean hostile | cochain exact, chain nonzero; no `F13 -> F2` map |
| first blocker | no lawful same-copy closure edge |
| second blocker | no location-to-response coefficient transport; exact response hostile rules out every `r1`-blind scalar version |
| quotient hostile | the common ten-space has no common section or arc-reversal action; stabilization kills its two middle orientation defects |
| strongest survivor | exact support alignment, explicit marked interface, and conditional `-kappa_lambda` prediction; an address-dependent section remains untyped |
| scope | no response amplitude, physical current, D5 realization, row exclusion, or LRC(14) |

## Reproduction

```text
python -B 04-computation/lrc_r5_tent_location_c4_cospan_hostile_audit_20260816.py
python -B -O 04-computation/lrc_r5_tent_location_c4_cospan_hostile_audit_20260816.py
```

Normal and optimized output byte-match the stored transcript.  Script,
output, and semantic SHA-256 are respectively

```text
856eda8328665ab1a060c37dc889506a097c1c91c75036bd1b22c2763b8adc54
bf83ef8043625a87e8a56c8478bf5c0ab0f502bf2747e72a1da074b23d774dea
d662e9ffb12e35f60b34fcf86d5cc60788aded7de28fdcc132cfe1fc693290c3.
```
