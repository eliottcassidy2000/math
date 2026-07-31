# Full endpoint rectangles and the local augmentation obstruction

**Status: FINITE-EXACT SCRATCH AUDIT.  Not canon and not an LRC(14)
closure.**

This audit answers a deliberately narrow question left by the THM-2825
half-step collar: after decorating all `587` cofiber-rooted paths by their
full `169`-bit **source** endpoint masks, can the masks be described by a
lawful endpoint-address `0`-potential or `1`-coboundary?

The answer depends on which of three inequivalent meanings of
“coboundary” is intended:

1. as a cochain with values in the ambient module
   \(\mathbb Z^{\mathbb F_{13}^2}\), it is tautologically exact;
2. as a physical translation gauge
   \(E_y=T_{\delta_e}E_x\), it fails on every nonconstant adjacent edge;
3. as an address derivative
   \(1_{E_y}-1_{E_x}=(T_\delta-I)\phi\), it likewise fails over
   \(\mathbb Z\) or \(\mathbb F_{13}\) on every nonconstant adjacent edge.

The failure is local and stronger than a holonomy failure.  The path
complex is a forest, so any lawful translation labels would integrate.
Most edges never acquire a lawful label because their endpoint-mask
cardinality changes.

## 1. Independent reconstruction

The companion
[`lrc_endpoint_coboundary_audit.py`](lrc_endpoint_coboundary_audit.py)
pins the THM-2825 primary companion at

```text
bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a
```

and imports only its physical interval geometry.  It does **not**
materialize the large union of present intervals for each of the `169`
endpoint addresses.  Instead it evaluates the defining periodic
inequalities exactly on every selected interval, subdividing at all
relevant comb boundaries.

This evaluator exactly reproduces the stored THM-2825 root/M1/M2 source
censuses and both stored translation-pair censuses:

```text
R -> M1 Hamming classes: 0:74, 9:187, 10:245, 81:81
R -> M2 Hamming classes: 0:74, 9:187, 10:245, 81:81
```

It then evaluates all `55,341` vertices of the `587` paths.  These use
`1,192` distinct physical intervals and only `73` distinct endpoint masks.

## 2. The hidden object is a moving rectangle

Every full source endpoint mask is exactly a Cartesian product

\[
E_x=A_x\times B_x\subseteq\mathbb F_{13}^2.
\]

This is not an empirical coincidence.  The inherited address
representative is

```text
REPS(a,b)=(0,-a,-b,0,0,0,a,b,0).
```

Consequently the address-dependent inequalities split into an
`a`-family and a `b`-family; all remaining inequalities are independent of
the address.  The complete vertex census by
\((|A_x|,|B_x|,|E_x|)\) is

```text
(0,0,0):       81
(9,9,81):   17,286
(9,10,90):   4,798
(9,11,99):   2,306
(10,9,90):  20,335
(10,10,100): 7,798
(10,11,110): 2,737
```

The full rectangle census has SHA-256

```text
f35de3c4d0fe5e38d612cc32ee8fc574662225ba4f83d038487b4e861d3b80e3
```

Thus a more faithful state object than “one endpoint bit” is the ordered
pair of one-dimensional complements \((A_x,B_x)\), together with the
physical interval carrying it.

## 3. Adjacent edges: the augmentation is complete

For the `54,754` `+h` edges:

```text
Hamming:      0:50,725, 9:1,806, 10:2,142, 81:81
augmentation: -10:973, -9:848, 0:50,725, +9:958, +10:1,169, +81:81
relation:     equal:50,725, birth:2,208, death:1,821
changed axis: first:63, second:3,885, empty boundary:81
```

For the `54,167` `+2h` edges:

```text
Hamming:      0:46,622, 9:3,425, 10:4,039, 81:81
augmentation: -10:1,946, -9:1,696, 0:46,622,
              +9:1,729, +10:2,093, +81:81
relation:     equal:46,622, birth:3,903, death:3,642
changed axis: first:112, second:7,352, empty boundary:81
```

Every nonempty nonconstant defect is a single row or column slice.
Crucially,

> an adjacent edge has zero augmentation if and only if its two full masks
> are literally equal.

This one scalar is therefore a complete local classifier on the audited
edge universe.

Any permutation of the `169` addresses preserves cardinality.  Hence no
nonconstant adjacent edge can be a permutation conjugacy, much less a
translation in \(\mathbb F_{13}^2\) or \(C_{169}\).  The maximal
permutation/translation-gauge-flat `+h` subforest has exactly `50,725`
edges and `4,616` components.  Only `14/587` complete paths remain in one
translation orbit, and those paths have a constant mask.

## 4. What “coboundary” can and cannot mean

Let \(\Phi(x)=1_{E_x}\).

### Ambient additive cochain

\[
D(x\to y)=\Phi(y)-\Phi(x)
\]

is, by definition, the coboundary \(d\Phi\).  This is useful bookkeeping,
but it does not produce a physical carrier action.

### Translation gauge

A split-plane gauge would require

\[
E_y=T_{\delta_e}E_x,\qquad \delta_e\in\mathbb F_{13}^2.
\]

Only the equal adjacent masks pass, with \(\delta_e=(0,0)\).  The same
verdict holds for cyclic translations in both digit orders
`n=13a+q` and `n=a+13q`.

Because the underlying graph is a forest, this is not an integrability or
loop-curvature obstruction.  Edgewise orbit equivalence fails before
there is a cochain to integrate.

### Additive translation derivative

For a nonzero split direction \(d\), membership in
\(\operatorname{im}(T_d-I)\) over the integers is equivalent to zero sum
on every order-`13` `d`-line.  In particular it has zero total
augmentation.  The exact audit finds:

```text
+h:  4,029 nonconstant defects admit 0 directions;
     50,725 zero defects admit all 14 projective directions.

+2h: 7,545 nonconstant defects admit 0 directions;
     46,622 zero defects admit all 14 projective directions.
```

The result remains true over \(\mathbb F_{13}\), since the nonzero
augmentations are `±9`, `±10`, or `81`, all nonzero modulo `13`.

### The characteristic-two survivor

Over \(\mathbb F_2\), and only after forgetting signs and positivity, the
even ten-point slices become exact one-axis derivatives:

```text
+h direction counts:  0:1,887, 1:2,142, 14:50,725
+2h direction counts: 0:3,506, 1:4,039, 14:46,622
```

The odd nine- and `81`-point defects still fail.  This survivor is a
coefficient-module fact; it does not yield mask conjugacy, preserve
endpoint origin, or supply a physical carrier.

## 5. Split endpoint plane versus the nonsplit carry

The endpoint key set is the split plane \(\mathbb F_{13}^2\), whereas the
THM-2851 carry sidecar is the nonsplit torsor \(C_{169}\).  The companion
checks both cyclic digit orders and exhaustively enumerates
\(\operatorname{AGL}(2,13)\):

```text
linear parts with f^13=id:   169
affine maps with f^13=id: 28,561
linear parts with f^169=id:  169
affine maps with f^169=id: 28,561
```

Thus no affine-plane element has order `169`.  Conceptually, a
characteristic-`13` affine map embeds in `GL(3,13)`; its `13`-primary
unipotent blocks have size at most `3`, so their exponent is `13`, not
`169`.

The THM-2851 horn has lawful local lifts and a nonzero carry `2`-cocycle.
The endpoint collar usually fails one level earlier: it changes the
support that a group action would have to transport.  Recharting the
address set cannot repair support birth or death.

## 6. The exact horn hinge

At the horn cell `(clock,s,t)=(1,0,5)` and root index `0`,

```text
R  = (142004190428100, 142004216872980)
M2 = (142004992589460, 142005019034340).
```

The source endpoint mask on `R` is empty.  The mask on `M2` is the
`81`-point rectangle

```text
A={0,1,2,3,4,5,6,7,12}
B={0,1,2,3,4,5,8,9,10}.
```

No permutation of addresses can transport the first mask to the second.
This is the minimal conceptual obstruction at the THM-2851 hinge:
physical passage from `R` to `M2` creates endpoint support; it does not
merely select a basepoint in pre-existing support.

A smallest nonempty `+h` witness occurs in cell `(1,0,3)`, root `0`,
source index `12`: the target adds exactly the nine-point slice

```text
{(a,7): a in {0,1,2,3,4,5,6,7,12}}.
```

Its augmentation is `+9`, so it already obstructs every address
permutation and every integral or characteristic-`13` translation
derivative.

## 7. A real nonlocal \(Z^8\) reference, but no carry triangle

The adjacent obstruction does not imply that all longer-range translated
pairs are absent.  Exhausting all `4,208,957` ordered pairs `i<j` on one
rooted path gives `798,322` translated pairs, of which `248,470` have
nonzero split-plane delta.  The complete nonzero delta census is

```text
(0,1): 50,712
(0,2): 14,812
(0,3): 11,109
(0,7):119,560
(0,8): 42,547
(0,9):  9,730
```

No horizontal delta occurs, even after the complete coordinate-swapped
search.  In particular `(0,4)`, `(4,0)`, `(8,0)`, and `(9,0)` never occur.

For the distinguished atom

```text
I=M2=(142004992589460,142005019034340),
```

the twelve horn paths with `s in {0,3,12}` and
`t in {5,6,9,10}` reach index `68`.  On every one,

```text
E_(index 68)=T_(0,8) E_I.
```

The gap from `I` at index `2` is `66h`.  Moreover indices `68,...,90`
carry the same `(0,8)` translate whenever the path is long enough.
No forward state on any of the twenty distinguished `I` paths is an
`(a,1)` translate of `E_I`.

Under the chosen split endpoint notation this is a genuine
\(Z^8=(0,8)\) **source-mask** reference.  A concurrent ancestry audit
reports that the selected source and target atoms also have identical
ancestry sets, which is the additional typing needed before calling this
a faithful \(T^8\) reference.  That extra ancestry equality is not
rederived by this companion.

This positive does not realize the THM-2851 carry triangle.  Although
\(8+9=4\pmod {13}\), the exhaustive within-path search finds:

```text
(0,8) pairs: 42,547 on 187 paths
(0,9) pairs:  9,730 on 49 paths
(0,4) pairs:      0
composable 8,9 versus 4 triangles: 0
```

Nor does either coordinate-swapped orientation occur.  The numerical
match between horn degree difference `q3 -> q11 = 8` and this endpoint
delta is therefore an untyped analogy until a map preserving the relevant
objects is exhibited.  The strongest present statement is: a typed
long-scale \(T^8\) reference exists, while the normalized \(T\) reference
and the `8+9=4` composition are absent from this collar path object.

## 8. Allocation-shifted macro boundary

For

\[
I_q=I+q(T/13),\qquad q\in\{0,3,7,11\},
\]

the independently evaluated full endpoint masks have the same source and
target mask at each allocation, with masses

```text
q0:  81
q3:  90
q7:   0
q11: 81
```

The `q0` and `q11` rectangles have equal cardinality but are **not**
translates in \(\mathbb F_{13}^2\), on either source or target.  Their
factors are

```text
q0:
 A={0,1,2,3,4,5,6,7,12}
 B={0,1,2,3,4,5,8,9,10}

q11:
 A={0,1,2,3,4,5,8,9,12}
 B={0,1,2,3,4,5,8,11,12}.
```

This is a sharper macro-complement boundary than a single occupancy bit:
`q7` is literal endpoint extinction, `q3` changes mass, and `q11`
preserves mass while changing translation orbit.  Any macro clutch must
therefore carry at least the rectangle factors or an equivalent sidecar;
allocation degree alone loses the obstruction.

## 9. Connection contract and next test

```text
source:
  the THM-2825 rooted physical collar decorated by full source endpoint
  masks E_x=A_x*B_x;

target:
  a lawful THM-2851/THM-2857 carry or semilinear endpoint clutch;

candidate map:
  physical interval step -> endpoint-address translation delta;

preserved on the positive long-range lane:
  full 169-bit source mask up to Z^8, and (from the independent sidecar)
  ancestry support on the selected atom;

destroyed on adjacent and allocation lanes:
  mask cardinality, translation orbit, endpoint support, and sometimes an
  entire rectangle factor;

needed sidecar:
  physical source/target typing, ancestry action, semantic survival,
  E3/complement truth, and the coefficient semilinearity required by
  THM-2857;

cheapest decisive next test:
  decorate the twelve typed I -> index-68 Z^8 pairs by the exact endpoint
  scalar and ancestry operator, then test whether the map is semilinear
  under one carried step; do not infer a normalized T step or an 8+9=4
  triangle from the Z^8 numerical label.
```

## 10. Reproduction boundary

Run:

```bash
python3 .scratch/lrc_endpoint_coboundary_20260728/lrc_endpoint_coboundary_audit.py
python3 -O .scratch/lrc_endpoint_coboundary_20260728/lrc_endpoint_coboundary_audit.py
```

The script contains no executable Python `assert`; it uses explicit
runtime checks.  The normal and optimized transcripts are required to be
byte-identical to
[`lrc_endpoint_coboundary_audit.out`](lrc_endpoint_coboundary_audit.out).
That three-way byte comparison passed.  LF-normalized SHA-256:

```text
script  418fa31e87093a03f63fc62b961da375b4d4c256be0ffeaac96c918c4f4d6171
output  c7c10df10e21bacebb9934f37bf2e48023bdc5a93246fbce59aeaafbc8b7540d
```
