# Independent hostile audit: the source branch repairs the atom clock, not the current

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS THE SCOPED SOURCE-ALIGNED
PACKAGE.  LRC(14) remains open.**  This audit was reconstructed from the
hash-pinned THM-2594 profile engine and the independent THM-3514 endpoint-atom
engine.  Neither submitted source-aligned script nor any of its bucket arrays
was imported.  Candidate artifacts were read only after the clean-room run had
fixed its own digests and exact profiles.

## Verdict

The source refinement

```text
f_omega^src = 1_Q P^K(e P_omega),    e_nu = e P_nu
```

does retain the endpoint/source guard address before transfer.  Its retained
inverse branch recovers the erased sheet as `floor(n/13)`, and the resulting
39-atom source labels can be paired on one THM-2471 first-collision base.  The
finite endpoint-weighted sum is therefore type-correctly rewritable as one
simple-kernel integral by finite linearity.

That statement has a strict boundary.  THM-3514's `AX/BY` entries have already
been integrated over a separate endpoint variable before they become scalar
coefficients of the outer-base kernel.  The calculation does not identify a
physical relation current at a common endpoint node.  It also does not produce
a grouped coefficient, row exclusion, THM-2512 bridge, or LRC(14).

## Independent reconstruction

The raw `T^2` transition was rebuilt directly as all 1521 intersections
between the 39 guard atoms and their transfer preimages.  Every entry is
nonzero, yet its exact rank is three (checked independently modulo three
primes), and all thirteen source-sheet rows agree within each chamber.  Thus
full support masks total loss of every nontrivial source-sheet character.  On
each of the 169 inverse branches, however, the high base-13 digit recovers the
source sheet exactly:

```text
sheet = floor(n/13),   with 13 branches per sheet.
```

Splitting `e` before transfer, applying `P^2` atomwise, restricting by `Q`, and
then applying `P_(13^5)` gives twenty nonempty source atoms on each leg.
Atomwise recombination restores the unsplit `P^2 e`, `U`, and `V` profiles
exactly.  Before gauge restriction all `20 x 20 = 400` root pairs occur; the
common gauge `a-u=b-q`, equivalently `s=-d`, leaves 362 pairs and all 48 live
owner-active chamber/drift types.

Keeping the absolute common offset through Fourier analysis gives rank four
at every owner character.  At owner character zero the support is exactly

```text
(13,12,12,13),
```

while each of the twelve nonzero owner characters has full
`(13,13,13,13)` support.  This zero-mode exception was recovered by the audit
before candidate comparison and is a useful guard against overphrasing the
spectral statement.

## Endpoint bank and the fixed relation residue

The complete THM-3514 `(1,0,t)`, `t in F_13`, endpoint pair bank was rebuilt
from the independent atom tables.  Its digest is

```text
c28119c8b54f47e5b7a46f1508fbba604b0e3997eaadb05b03ad28edd9aed468.
```

All 169 residue-by-owner pullbacks are nonzero.  The same is true after
double-centering over all 39 endpoint addresses.  As an additional hostile,
the audit also double-centred over the twenty live source addresses; all 169
entries remain nonzero there as well.

For the fixed relation's refined class `(1,0,6)`, the endpoint total is

```text
225010624370142818572 mod 572252886246508880869,
```

and both the full and all-address doubly-centred source-aligned profiles have
support `13/13`.  These exact values agree entry-for-entry with the candidate.
The interpretation is only the THM-2334 residue pushforward attached to the
class `(1,0,6)`.  It is not the exact-address orbit sum `C(a;X,m)`: the residue
map has already merged all exact addresses in that class.

## Candidate comparison

After the independent semantic surface was fixed, comparison with commit
`2a069ab1b` gave exact agreement on:

- raw-transition support, rank, row identity, and digest;
- inverse-branch sheet recovery;
- 400 all-root pairs, 362 common-gauge pairs, and 48 active types;
- every owner-character rank and Walsh--drift support count;
- the full endpoint-bank digest and all thirteen endpoint totals;
- all thirteen fixed-class full and doubly-centred profile entries; and
- the bridge full and doubly-centred profiles.

The audit semantic digest is
`822957b53ef3bcaa8509fa68bd2c0d3080dec7632db3eca063e2d947162bd4b9`.
Normal and `python -O` runs are byte-identical; the pinned script LF hash is
`827644be815bf3152507af9ff82c48b76e843126b285edeb76121bf2ed347396`.

## Surviving obstruction

The positive result lives on the `U_full` side of the existing two-transplant
boundary: it has the complete refined endpoint bank.  The `U_clock` side has
the delayed same-clock realization.  Nothing here identifies those two
objects or supplies a single typed row having both properties.  The remaining
work is therefore structural, not spectral: isolate an exact address inside
the nonzero residue class and realize the preintegrated endpoint factors as a
lawful physical current without losing the source branch or common clock.
