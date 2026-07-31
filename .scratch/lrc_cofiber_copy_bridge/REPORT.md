# THM-2806 / THM-2771 cofiber-copy bridge audit

**Status:** `FINITE-EXACT`, scratch only.  Both companions replay
byte-identically under ordinary and optimized Python.  No theorem, row
exclusion, root-deck correction, or LRC(14) conclusion is asserted here.

## Verdict

There is a real geometric explanation of the numerical factor `2` in the
rail-eight cell `(s,t,e)=(0,4,1)`, and in fact of every primitive multiplier
`2,121,265,254` in THM-2771:

> after all inherited coefficient filters, each primitive multiplier is the
> literal number of positive, equal-weight, equal-content interval copies.

There is no signed cancellation.  In the exceptional target column `t=12`,
the larger raw wing is cut by the delayed terminal selector into alternating
live and zero copies; the live counts are exactly `121,265,254`.

This positive result does **not** give the desired physical cospan from the
THM-2806 common sheet to THM-2771's right cofiber.  At `(0,4,1)` both
cofiber copies already fail the native `E3` factor, one also fails native
`c2`, neither has source carrier-twist support, and their endpoint masks do
not transport from the common sheet.  Equal interval content survives only
after forgetting precisely the sidecars needed for an endpoint/root/Cech
map.

The intrinsic coefficient Bockstein is additive on the literal copy split.
This is a lawful coefficient-level decomposition, not a realization of
THM-2771's target convolution or target-role-to-root-deck intertwiner.

## 1. The selected common and cofiber sheets

Put

```text
T = 297836897838480,
I = [142004992589460,142005019034340),
w = 27581135604,
c = 103478815440.
```

The two coefficient-effective right-cofiber sheets are

```text
J_- = I - T/13^5
    = [142004190428100,142004216872980),

J_+ = I + 96T/13^5
    = [142082000080020,142082026524900).
```

Their full `13^6` address shifts from `I` are `-13` and `+1248=13*96`.
The raw unweighted `R=B\A` carrier has `241` pieces in this cell, but after
the inherited physical coefficient cuts only `J_-` and `J_+` remain.
Both have weight `w` and length `26444880`.

The object typing is sharp:

```text
I   lies in M=A intersect B,
J_- and J_+ lie in R=B minus A.
```

Thus the two copies explain a cofiber coefficient.  They are not two
alternative representatives of the common object.

## 2. What the translations preserve and destroy

Use factor order

```text
(E3, clock, q1, q2, c2, c3).
```

The complete masks are

```text
I:
  source native   111111
  source pulled   111111
  target native   111111
  target adjacent 111111

J_-:
  source native   011111
  source pulled   111111
  target native   111111
  target adjacent 011111

J_+:
  source native   011101
  source pulled   111111
  target native   111111
  target adjacent 011101.
```

Therefore the first lost coordinate is native `E3` on both copies.
`J_+` additionally loses native `c2`.  In particular the two cofiber
copies are not exchangeable as fully typed atoms, even though their
coefficient content agrees.

All three intervals do preserve the complete literal THM-2584 ancestry
sets:

```text
|U|=966606,
|V|= 28534,
|U||V|=27581135604=w,
combined digest =
15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd.
```

The sets agree by identity, not merely cardinality.  The supplied path

```text
(a,b,e')=(59162,26,56658)
```

is active on all three intervals.

The unit delayed profiles also agree:

```text
source: only carry 12 is nonzero, with (root0,root1)=(0,c),
target: only carry  6 is nonzero, with (root0,root1)=(0,c),
(source root,target root,deep digit)=(5,6,12).
```

Thus ancestry and semantic integration alone cannot distinguish the
common sheet from the cofiber copies.  The native/pulled factor sidecar can.

## 3. Carrier and endpoint obstruction

The thirteen carrier-twist masks are

```text
I:   source=delta_0, target=delta_0,
J_-: source=0,       target=delta_0,
J_+: source=0,       target=delta_0.
```

Hence no cofiber copy has any source allocation state, at any twist.

Across the full `169` inherited endpoint-present addresses, the
zero/one support counts are

```text
              source          target
I             88/81           88/81
J_-           169/0           88/81
J_+           79/90           70/99.
```

The `J_-` target mask is literally the `I` target mask.  This tempting
one-sided equality does not extend across the source: `J_-` has no source
endpoint support.  The `J_+` source and target masks differ from `I` in
`27` and `18` addresses respectively.  Exhausting all `169`
`F_13^2` address translations finds no translation from either common mask
to either nonidentical cofiber mask.

Consequently the interval translations

```text
I -> J_-,       I -> J_+
```

are maps only after applying a forgetful functor that discards native
factor membership, carrier allocation, and endpoint masks.  That quotient
is too coarse for the THM-2806/2771 physical bridge.

## 4. Exact factor-two and Bockstein split

The arithmetic identities are

```text
57068*k1 = w,                 k1=483303,
57068*c  = g0,                g0=5905329039529920,
g0*k1    = w*c
          =2854063240791928925760.
```

Each of `J_-` and `J_+` has individual delayed clock vector

```text
(0,wc,wc,wc,wc,wc,wc).
```

Their disjoint sum is

```text
(0,2wc,2wc,2wc,2wc,2wc,2wc),
```

and division by `g0` gives THM-2771's primitive coefficient

```text
2k1=966606.
```

Since `v_13(wc)=1`, the intrinsic copywise coefficient Bockstein is

```text
beta(J_-)=beta(J_+)=wc/13=9 mod13,
beta(R)=9+9=5 mod13.
```

This exactly explains the local THM-2771 raw Bockstein entry.  It does not
produce the target convolution `K_beta`: that operation mixes target
labels, while the present copy split is internal to one target cell.

After removing the common factor `13` from the full-address shifts, their
depth-one residues are `12` and `5`.  The phase-refined local germ is

```text
9(z^12+z^5)
 =5+10 epsilon+...       in F_13[C_13], epsilon=z-1.
```

Its augmentation is `5`, so it is a local-algebra unit, not a uniformizer.
THM-2771's pure-target uniformizer therefore arises only after the global
cross-column/clock assembly; it is absent in this local two-copy germ.

## 5. All 28 nonzero cells

Every coefficient-effective interval has the same length `26444880`.
There are two ancestry/weight chambers:

```text
physical clock 1:
  |U|=966606, |V|=28534,
  w1=27581135604,
  copy value g0*k1=2854063240791928925760,
  copy Bockstein=9.

physical clocks 2 and 3:
  |U|=966574, |V|=28534,
  w2=27580222516,
  copy value g0*k2=2853968755527296447040,
  k2=483287,
  copy Bockstein=2.
```

Within every one of the `28` nonzero cells, all raw intervals lie in one
wall-free `U` chamber and one wall-free `V` chamber.  Hence the live copies
have identical ancestry sets throughout that cell.  The supplied path
`(59162,26,56658)` remains active in both chambers.

For every nonexceptional target column the raw/effective partition is

```text
2 raw = 2 live + 0 zero.
```

For `t=12` it is

```text
(e,t)=(1,12): 241 raw =121 live +120 zero,
(e,t)=(2,12): 528 raw =265 live +263 zero,
(e,t)=(3,12): 506 raw =254 live +252 zero.
```

Thus the primitive multipliers are literal positive live-copy counts:

```text
2, 121, 265, 254.
```

No negative piece and no cancellation occurs.

## 6. Why the exceptional counts are 121, 265, and 254

The `t=12` pieces form linearly ordered blocks.  Consecutive left endpoints
inside a block differ by

```text
T/(2*13^5)=401080680.
```

The delayed terminal selector is exactly the alternating word

```text
1,0,1,0,...
```

starting with a live copy in every block.  The exact block profiles
`(length,live,zero)` are

```text
clock 1: (145,73,72), ( 96,48,48),
clock 2: (143,72,71), (289,145,144), (96,48,48),
clock 3: (143,72,71), (289,145,144), (74,37,37).
```

Therefore

```text
121 = ceil(145/2)+ceil( 96/2),
265 = ceil(143/2)+ceil(289/2)+ceil(96/2),
254 = ceil(143/2)+ceil(289/2)+ceil(74/2).
```

This is the structural origin of the previously opaque large primitive
coefficients: they are parity counts on a small number of interval chains.

The first raw failure of uniform semantic content is `(1,12)`, where
`121` copies carry `c` and `120` carry zero.  After deleting the zero
copies, no cellwise uniformity failure remains.

The first change from the clock-one global ancestry prototype is `(2,2)`:
the `U` set loses exactly `32` labels and the weight changes from `w1` to
`w2`; `V` and the supplied path remain unchanged.  This is a chamber change,
not a failure of the copy-count mechanism.

## 7. Sharp boundary and next object

What is proved in scratch:

```text
coefficient-effective R
  = disjoint positive live interval copies,

raw THM-2771 coefficient
  = sum of their identical clock vectors,

intrinsic Bockstein
  = sum of their copywise Bocksteins.
```

What is not constructed:

```text
M common sheet
  -> R copy with native factor membership,
  -> source carrier allocation,
  -> transported endpoint mask,
  -> target convolution,
  -> target-role/root-deck character.
```

The cheapest next object is not another scalar identity.  It is a typed
boundary path that records crossings of the `E3` and `c2` walls together
with endpoint origin.  Without that wall-incidence sidecar, the exact
interval translations live only in a nonfaithful quotient.

## 8. Reproduction

From the repository root:

```bash
python3 .scratch/lrc_cofiber_copy_bridge/audit.py
python3 -O .scratch/lrc_cofiber_copy_bridge/audit.py

python3 .scratch/lrc_cofiber_copy_bridge/all_cells.py
python3 -O .scratch/lrc_cofiber_copy_bridge/all_cells.py
```

The ordinary/optimized pairs byte-match.  Current LF-normalized hashes are

```text
a154a4fba6e4901264e2228b29b980eadff3702d69b04d0061c77c695c2462cf
  audit.py
4d2458e2a5f71ba32580484aed9076b42280144e5d196bd52c07d32e006aef33
  audit.out

d7d6b86ffb8ac61249c5309f0369c0cde630a80afbb954081b79c8b17540d2eb
  all_cells.py
43d73c8878aab073094d79ba4a1f675d5e98273786a190df44ab681a8dc28bc9
  all_cells.out
```
